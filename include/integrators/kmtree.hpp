#ifndef KMTREE_HPP
#define KMTREE_HPP

#include <vector>
#include <queue>
#include <algorithm>
#include <functional>

#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>
#include "kmeans.hpp"

template <class descr_t, class index_t>
class kmeans_tree
{

    public:

    class node_t
    {
        public:

        descr_t              _center;
        std::vector<node_t>  _childs;
        std::vector<index_t> _indices;
        double               _radius;

        bool is_leaf() const { return _childs.size() == 0; }

        friend class boost::serialization::access;
        template <class archive_t> void serialize(archive_t& ar, unsigned int /* version */)
        {
            ar & _center & _childs & _indices & _radius;
        }
    };

    kmeans_tree() {}

    kmeans_tree(std::size_t branching) : _branching(branching) {}

    const node_t& root() const { return _root; }

    int depth() const { return depth(root()); }

    int depth(const node_t& node) const
    {
        if (node.is_leaf()) return 0;
        else
        {
            int max_child_depth = 0;

            for (int i = 0; i < node._childs.size(); i++)
            {
                max_child_depth = std::max(max_child_depth, depth(node._childs[i]) + 1);
            }
            return max_child_depth;
        }
    }

    template <class collection_t, class dist_fn>
    kmeans_tree(const collection_t& collection, std::size_t branching, KmeansInitAlgorithm kminit, const dist_fn& distfn)
     : _branching(branching), _maxchecks(1000)
    {
        build(collection, distfn, kminit);
    }

    template <class collection_t, class dist_fn>
    void build(const collection_t& collection, const dist_fn& distfn, KmeansInitAlgorithm kminit, size_t limit)
    {
        _root = node_t();
        std::vector<index_t> indices(collection.size());
        for (std::size_t i = 0; i < collection.size(); i++) indices[i] = i;
        build_recursive(_root, collection, indices, _branching, distfn, kminit, limit, 0);
    }

    template <class collection_t, class dist_fn>
    void build(const collection_t& collection, const dist_fn& distfn, KmeansInitAlgorithm kminit)
    {
        build(collection, distfn, kminit, _branching);
    }

    template <class collection_t, class result_t, class dist_fn>
    void query(const descr_t& descr, const collection_t& collection, result_t& result, std::size_t rsize, const dist_fn& distfn) const
    {
        typedef std::pair<double, const node_t*> distnode_t;
        std::priority_queue<distnode_t, std::vector<distnode_t>, std::greater<distnode_t> > queue;

        queue.push(std::make_pair(0.0, &_root));

        while (!queue.empty() && result.size() < rsize)
        {
            const node_t* current = queue.top().second;
            queue.pop();

            if (!current->_indices.empty())
            {
                // leaf node: add all indices of leaf to result
                for (std::size_t i = 0; i < current->_indices.size(); i++)
                {
                    index_t index = current->_indices[i];
                    double dist = distfn(descr, collection[index]);
                    result.push_back(std::make_pair(dist, index));
                }
            }
            else
            {
                // inner node: add child nodes to priority queue
                for (std::size_t i = 0; i < current->_childs.size(); i++)
                {
                    const node_t& child = current->_childs[i];
                    double dist = distfn(descr, child._center);

                    queue.push(std::make_pair(dist, &child));
                }
            }
        }

        std::sort(result.begin(), result.end());
    }

    template <class result_t, class dist_fn>
    void query_no_pp(const descr_t& descr, result_t& result, std::size_t rsize, dist_fn& distfn) const
    {
        typedef std::pair<double, const node_t*> distnode_t;
        std::priority_queue<distnode_t, std::vector<distnode_t>, std::greater<distnode_t> > queue;

        queue.push(std::make_pair(0.0, &_root));

        while (!queue.empty() && result.size() < rsize)
        {
            const node_t* current = queue.top().second;
            queue.pop();

            if (!current->_indices.empty())
            {
                // leaf node: add all indices of leaf to result
                std::copy(current->_indices.begin(), current->_indices.end(), std::back_inserter(result));
            }
            else
            {
                // inner node: add child nodes to priority queue
                for (std::size_t i = 0; i < current->_childs.size(); i++)
                {
                    const node_t& child = current->_childs[i];
                    double dist = distfn(descr, child._center);

                    queue.push(std::make_pair(dist, &child));
                }
            }
        }
    }

    void set_maxchecks(std::size_t maxchecks)
    {
        _maxchecks = maxchecks;
    }

    void get_leaf_cluster(std::vector<index_t>& clusterindex) const
    {
        std::queue<const node_t*> queue;

        size_t cn = 0;

        queue.push(&_root);
        while (!queue.empty())
        {
            const node_t* n = queue.front();
            queue.pop();

            if (n->_childs.empty())
            {
                for (size_t i = 0; i < n->_indices.size(); i++)
                {
                    size_t index = n->_indices[i];
                    if (clusterindex.size() <= index) clusterindex.resize(index+1);
                    clusterindex[index] = cn;
                }

                ++cn;
            }
            else
            {
                for (size_t i = 0; i < n->_childs.size(); i++)
                {
                    queue.push(&n->_childs[i]);
                }
            }
        }
    }

    private:

    template <class collection_t, class dist_fn>
    void build_recursive(node_t& node,
                         const collection_t& collection,
                         const std::vector<index_t>& indices,
                         std::size_t branching,
                         const dist_fn& distfn,
                         KmeansInitAlgorithm kminit,
                         size_t limit,
                         std::size_t level)
    {
        typedef kmeans<collection_t, dist_fn> cluster_fn;

//        std::cout << "build recursive: level=" << level << " indices.size=" << indices.size() << std::endl;

        if (indices.size() <= limit)
        {
            node._indices = indices;
            std::cout << "level " << level << ": " << indices.size() << std::endl;
            return;
        }

        collection_t localcoll(indices.size());
        for (std::size_t i = 0; i < localcoll.size(); i++) localcoll[i] = collection[indices[i]];

        cluster_fn clusterfn(localcoll, branching, kminit, distfn);

        std::vector<std::vector<index_t> > clusters;

        clusterfn.run(30, 0.01);
        clusterfn.make_cluster_table(clusters);

        const std::vector<std::size_t>& ci(clusterfn.clusters());

        // intercept degenerated clusterings, i.e. all points belong to a single cluster
        // this happens only when all points are equal
        if (std::find_if(ci.begin(), ci.end(), std::bind1st(std::not_equal_to<std::size_t>(), ci.front())) == ci.end())
        {
            node._indices = indices;
            return;
        }

        for (std::size_t i = 0; i < clusters.size(); i++)
        {
            // intercept empty clusters
            if (clusters[i].size() == 0) continue;

            // make a child node
            node_t current;

            current._center = clusterfn.centers()[i];

            // recursive call
            std::vector<index_t> currentindices(clusters[i].size());
            for (std::size_t k = 0; k < clusters[i].size(); k++) currentindices[k] = indices[clusters[i][k]];
            build_recursive(current, collection, currentindices, branching, distfn, kminit, limit, level+1);

            // add node
            node._childs.push_back(current);
        }
    }

    friend class boost::serialization::access;
    template <class archive_t> void serialize(archive_t& ar, unsigned int /* version */)
    {
        ar & _root & _branching & _maxchecks;
    }

    node_t      _root;
    std::size_t _branching;
    std::size_t _maxchecks;
};

#endif // KMTREE_HPP
