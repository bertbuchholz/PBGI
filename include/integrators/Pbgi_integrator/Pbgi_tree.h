#ifndef PBGI_TREE_H
#define PBGI_TREE_H

#include <queue>

#include <core_api/bound.h>

__BEGIN_YAFRAY


struct Tree_subdivision_decider
{
    // returns true if a subdivision should be made
    bool operator() (bound_t const& bound, int const num_points) const
    {
        if (num_points == 1) return false;

        if (max_num_points == 0) // check for bounds
        {
            int const longest_axis = bound.largestAxis();

            return bound.get_length(longest_axis) > max_size;
        }
        else
        {
            return (num_points > max_num_points);
        }

        return false;
    }

    float max_size;
    int   max_num_points;
};







struct Compare_vector_by_axis
{
    Compare_vector_by_axis(int const a) : axis(a) {}

    bool operator() (vector3d_t const& a, vector3d_t const& b)
    {
        return a[axis] < b[axis];
    }

    int axis;
};



template <class Data>
class Node
{
public:
    Node() : _initialized(false) {}

    virtual ~Node() {}

    virtual void average() {}

    virtual void add_points(std::vector<vector3d_t> const& points, bound_t const& bound, int const depth, Tree_subdivision_decider const& subdiv_decider) = 0;

    virtual void get_post_order_queue(std::vector<Node*> & list, int const depth, int const max_depth = -1) = 0;

    virtual bool is_leaf() const = 0;

    virtual Node *      get_child(int const i) = 0;
    virtual Node const* get_child(int const i) const = 0;

    // virtual void average_node(Averager const& averager) = 0;

    // virtual bound_t get_bound() const = 0;

    Data &      get_data()       { assert(_initialized); return _data; }
    Data const& get_data() const { assert(_initialized); return _data; }

    void set_data(Data const& data) { _data = data; _initialized = true; }

    virtual int get_shortest_distance_to_leaf(int const height) const = 0;

    template <class T >
    // T const* get_derived() const { return static_cast<T const*>(this); }
    T const* get_derived() const { T const* ptr = dynamic_cast<T const*>(this); assert(ptr); return ptr; }
    template <class T >
    // T      * get_derived()       { return static_cast<T*>(this); }
    T      * get_derived()       { T* ptr = dynamic_cast<T*>(this); assert(ptr); return ptr; }

protected:
    Data _data;
    bool _initialized;
};

template <class Data, class Leaf_data>
class Leaf_node : public Node<Data>
{
public:
    typedef Node<Data> Tree_node;

    void add_points(std::vector<vector3d_t> const& points, bound_t const& /* bound */, int const depth, Tree_subdivision_decider const& /* subdiv_decider */)
    {
        // assert(points.size() == 1);
        for (size_t i = 0; i < points.size(); ++i)
        {
            Leaf_data d;
            d.pos = points[i];
            _surfels.push_back(d);
        }

        // this->_data = points[0].data;
        // this->_depth = depth; // FIXME: re-add
    }

    void get_post_order_queue(std::vector<Tree_node*> & list, int const /* depth */, int const /* max_depth = -1 */)
    {
        list.push_back(this);
    }

    bool is_leaf() const { return true; }

//    void average_node(Averager const& averager)
//    {
//        // this->_data = averager.average_leaf(_surfels);
//    }

    std::vector<Leaf_data> const& get_surfels() const { return _surfels; }
    std::vector<Leaf_data>      & get_surfels()       { return _surfels; }

    int get_shortest_distance_to_leaf(int const height) const { return height; }

    // shouldn't be called
    bound_t get_bound() const { assert(false); return bound_t(); }

    Tree_node *      get_child(int const i)       { return NULL; }
    Tree_node const* get_child(int const i) const { return NULL; }

private:
    std::vector<Leaf_data> _surfels;
};


template <typename T>
inline
T log_base(T const value, T const base)
{
    return std::log(value) / std::log(base);
}

inline
int log2_int(int const value)
{
    int result = -1;
    int tmp_val = value;
    while (tmp_val != 0) {
        tmp_val >>= 1;
        result++;
    }
    return result;
}

template <class Data, class Leaf_data, int Arity = 2>
class Inner_node : public Node<Data>
{
public:
    typedef Node<Data> Tree_node;

    Inner_node()
    {
        /*
        for (int i = 0; i < Arity; ++i)
        {
            _children[i] = NULL;
        }
        */
    }


    void add_points_any_arity(std::vector<vector3d_t> const& points, bound_t const& bound, int const depth, Tree_subdivision_decider const& subdiv_decider)
    {
        // this->_depth = depth; // FIXME: re-add

        // _bound = bound;

        // vector3d_t split_positions = bound.center();

        int const num_splits = log2_int(Arity);

        std::vector<bound_t> new_bounds(Arity, bound);
        // new_bounds[0] = bound;

        for (int split = 1; split <= num_splits; ++split)
        {
            int const spacing = Arity / std::pow(2, split);

            int split_axis = split - 1;

            if (num_splits < 3)
            {
                split_axis += depth % (3 - (num_splits - 1));
                split_axis %= 3;
            }

            for (unsigned int j = 0; j < std::pow(2, split - 1); ++j)
            {
                int parentIndex = j * spacing * 2;
                int childIndex = parentIndex + spacing;
                bound_t parent_bound = new_bounds[parentIndex];

                new_bounds[childIndex] = new_bounds[parentIndex];

                new_bounds[parentIndex].g[split_axis] -= 0.5f * (parent_bound.g[split_axis] - parent_bound.a[split_axis]);
                new_bounds[childIndex]. a[split_axis] += 0.5f * (parent_bound.g[split_axis] - parent_bound.a[split_axis]);
            }
        }

        // volume assert
        float volume = 0.0f;

        for (size_t i = 0; i < new_bounds.size(); ++i)
        {
            volume += new_bounds[i].vol();
        }

        assert(std::abs(1 - volume / bound.vol()) < 0.001f);

        std::vector<std::vector<vector3d_t> > new_points(Arity);

        for (size_t j = 0; j < points.size(); ++j)
        {
            for (int i = 0; i < Arity; ++i)
            {
                if (new_bounds[i].includes(points[j]))
                {
                    new_points[i].push_back(points[j]);
                    break;
                }
            }
        }

        bool all_empty = true;
        for (int i = 0; i < Arity; ++i)
        {
            all_empty = new_points[i].size() == 0;

            if (!all_empty) break;
        }

        assert(!all_empty);

        for (int i = 0; i < Arity; ++i)
        {
            if (new_points[i].size() == 0) // empty child
            {
                continue;
            }
            else
            {
                Tree_node* child = NULL;

                if (!subdiv_decider(new_bounds[i], new_points[i].size()))
                {
                    child = new Leaf_node<Data, Leaf_data>;
                }
                else
                {
                    child = new Inner_node;
                }

                child->add_points(new_points[i], new_bounds[i], depth + 1, subdiv_decider);

                _children.push_back(child);

                std::vector<vector3d_t>().swap(new_points[i]);
            }
        }
    }

//    bound_t get_points_bound(std::vector<vector3d_t> const& points)
//    {

//    }

    void add_points_binary_median(std::vector<vector3d_t> const& points, bound_t const& bound, int const depth, Tree_subdivision_decider const& subdiv_decider)
    {
        assert(Arity == 2);

//        calc_points_bounding_box();
//        Point cell_extent = get_points_extent();

        this->_depth = depth;

        // _bound = bound;

        int const split_axis = bound.largestAxis();
        // int const split_axis = get_longest_axis_3d(cell_extent);

        std::vector<vector3d_t> points_sorted_by_axis = points;

        std::sort(points_sorted_by_axis.begin(), points_sorted_by_axis.end(), Compare_vector_by_axis(split_axis));

        int const median_index = int(points_sorted_by_axis.size() / 2.0f);

        std::vector<vector3d_t> left_points(points_sorted_by_axis.begin(), points_sorted_by_axis.begin() + median_index);
        std::vector<vector3d_t> right_points(points_sorted_by_axis.begin() + median_index, points_sorted_by_axis.end());

        std::vector<vector3d_t>().swap(points_sorted_by_axis);

        float const split_position =
                0.5f * (left_points.back())[split_axis] +
                0.5f * (right_points.front())[split_axis];

        bound_t left_bound = bound;
        // left_bound.setMax(split_position, split_axis);
        left_bound.setMax(left_points.back()[split_axis], split_axis);

        bound_t right_bound = bound;
        // right_bound.setMin(split_position, split_axis);
        right_bound.setMin(right_points.front()[split_axis], split_axis);

        _children.resize(2);

        if (left_points.size() == 1)
        {
            _children[0] = new Leaf_node<Data, Leaf_data>;
        }
        else
        {
            _children[0] = new Inner_node;
        }

        if (right_points.size() == 1)
        {
            _children[1] = new Leaf_node<Data, Leaf_data>;
        }
        else
        {
            _children[1] = new Inner_node;
        }

        _children[0]->set_parent(this);
        _children[1]->set_parent(this);

        _children[0]->add_points(left_points, left_bound, depth + 1, subdiv_decider);
        std::vector<vector3d_t>().swap(left_points);

        _children[1]->add_points(right_points, right_bound, depth + 1, subdiv_decider);
        std::vector<vector3d_t>().swap(right_points);
    }

    void add_points(std::vector<vector3d_t> const& points, bound_t const& bound, int const depth, Tree_subdivision_decider const& subdiv_decider)
    {
        add_points_any_arity(points, bound, depth, subdiv_decider);

        /*
        if (Arity == 8)
        {
            add_points_any_arity(points, bound, depth, subdiv_decider);
        }
        else
        {
            add_points_binary_median(points, bound, depth, subdiv_decider);
        }
        */
    }

    void get_post_order_queue(std::vector<Tree_node*> & list, int const depth, int const max_depth = -1)
    {
        if (max_depth == -1 || (max_depth != -1 && depth < max_depth))
        {
            for (size_t i = 0; i < _children.size(); ++i)
            {
                _children[i]->get_post_order_queue(list, depth + 1, max_depth);
            }
        }

        list.push_back(this);
    }

    bool is_leaf() const { return false; }

    Tree_node *      get_child(int const i)       { return _children[i]; }
    Tree_node const* get_child(int const i) const { return _children[i]; }


//    void average_node(Averager const& averager)
//    {
//        // TODO: change bounds to sum of childrens' BB -> bounds in GiPoint, so no bounds needed (?) in the tree

//        std::vector<Data*> children_data;

//        for (size_t i = 0; i < _children.size(); ++i)
//        {
//            children_data.push_back(&_children[i]->get_data());
//        }

//        this->_data = averager.average_node(children_data);
//    }

    int get_shortest_distance_to_leaf(int const height) const
    {
        std::vector<int> heights;

        for (size_t i = 0; i < _children.size(); ++i)
        {
            heights.push_back(_children[i]->get_shortest_distance_to_leaf(height + 1));
        }

        return *std::min_element(heights.begin(), heights.end());
    }

    std::vector< Tree_node* >      & get_children()       { return _children; }
    std::vector< Tree_node* > const& get_children() const { return _children; }

    std::vector< Data const* > get_children_data() const
    {
        std::vector< Data const* > result;

        for (size_t i = 0; i < _children.size(); ++i)
        {
            result.push_back(&_children[i]->get_data());
        }

        return result;
    }

private:
    // Node<Data, Averager> * _children[Arity];
    std::vector< Tree_node* > _children;
};

template <class Inner_node_data, class Leaf_data, int n_Arity>
class Pbgi_tree
{
public:
    static int const Arity = n_Arity;

    typedef Node<Inner_node_data> Tree_node;
    typedef Inner_node<Inner_node_data, Leaf_data, Arity> Tree_inner_node;
    typedef Leaf_node<Inner_node_data, Leaf_data> Tree_leaf_node;

    Pbgi_tree() : _root(NULL) {}

    ~Pbgi_tree()
    {
        if (_root)
        {
            clear();
        }
    }

    void build_tree(std::vector<vector3d_t> const& points, bound_t const& bound, Tree_subdivision_decider const& subdiv_decider)
    {
        _root = new Tree_inner_node;
        _root->add_points(points, bound, 0, subdiv_decider);
    }

    void clear()
    {
        std::vector<Tree_node*> nodes = get_post_order_queue();

        for (size_t i = 0; i < nodes.size(); ++i)
        {
            Tree_node * node = nodes[i];
            delete node;
        }

        _root = NULL;
    }

    std::vector<Tree_node*> get_post_order_queue(int const max_depth = -1)
    {
        std::vector<Tree_node*> result;
        _root->get_post_order_queue(result, 0, max_depth);
        return result;
    }

    Tree_node const* get_root() const
    {
        return _root;
    }

    std::vector<Tree_node const*> get_leafs() const
    {
        std::vector<Tree_node const*> result;

        std::queue<Tree_node const*> queue;
        queue.push(_root);

        while (!queue.empty())
        {
            Tree_node const* node = queue.front();
            queue.pop();

            if (node->is_leaf())
            {
                result.push_back(node);
            }
            else
            {
                Tree_inner_node const* inner_node = node->template get_derived<Tree_inner_node>();
                std::vector< Tree_node* > const& children = inner_node->get_children();

                for (int i = 0; i < children.size(); ++i)
                {
                    queue.push(children[i]);
                }
            }
        }

        return result;
    }

    std::vector<Tree_node const*> get_nodes() const
    {
        std::vector<Tree_node const*> result;

        std::queue<Tree_node*> queue;
        queue.push(_root);

        while (!queue.empty())
        {
            Tree_node * node = queue.front();
            queue.pop();

            result.push_back(node);

            if (!node->is_leaf())
            {
                Tree_inner_node const* inner_node = node->template get_derived<Tree_inner_node>();
                std::vector< Tree_node* > const& children = inner_node->get_children();

                for (size_t i = 0; i < children.size(); ++i)
                {
                    queue.push(children[i]);
                }
            }
        }

        return result;
    }

    /*
    std::vector<Tree_node*> get_nodes_at_depth(int const depth) const
    {
        std::vector<Tree_node*> result;

        std::queue<Tree_node*> queue;
        queue.push(_root);

        while (!queue.empty())
        {
            Tree_node * node = queue.front();
            queue.pop();

            if (node->get_depth() == depth)
            {
                result.push_back(node);
            }
            else if (!node->is_leaf() && node->get_depth() < depth)
            {
                Tree_inner_node const* inner_node = node->template get_derived<Tree_inner_node>();
                std::vector< Tree_node* > const& children = inner_node->get_children();

                for (size_t i = 0; i < children.size(); ++i)
                {
                    queue.push(children[i]);
                }
            }

            assert(node->get_depth() <= depth);
        }

        return result;
    }
    */

    void count_nodes(int & num_inner_nodes, int & num_leafs) const
    {
        num_inner_nodes = 0;
        num_leafs = 0;

        std::queue<Tree_node*> queue;
        queue.push(_root);

        while (!queue.empty())
        {
            Tree_node * node = queue.front();
            queue.pop();

            if (node->is_leaf())
            {
                ++num_leafs;
            }
            else
            {
                ++num_inner_nodes;
            }

            if (!node->is_leaf())
            {
                Tree_inner_node const* inner_node = node->template get_derived<Tree_inner_node>();
                std::vector< Tree_node* > const& children = inner_node->get_children();

                for (int i = 0; i < children.size(); ++i)
                {
                    queue.push(children[i]);
                }
            }
        }
    }

    struct Node_gatherer
    {
        void operator() (Tree_node const* node)
        {
            result.push_back(node);
        }

        std::vector<Tree_node const*> result;
    };

    std::vector<Tree_node const*> get_nodes_using_handler() const
    {
        Node_gatherer gatherer;

        traverse_tree(gatherer, gatherer);

        return gatherer.result;
    }


    template <typename Inner_node_handler, typename Leaf_handler>
    void traverse_tree(Inner_node_handler & inner_node_handler, Leaf_handler & leaf_handler) const
    {
        std::queue<Tree_node*> queue;
        queue.push(_root);

        while (!queue.empty())
        {
            Tree_node * node = queue.front();
            queue.pop();


            if (node->is_leaf())
            {
                leaf_handler(node);
            }
            else
            {
                inner_node_handler(node);

                Tree_inner_node const* inner_node = node->template get_derived<Tree_inner_node>();
                std::vector< Tree_node* > const& children = inner_node->get_children();

                for (size_t i = 0; i < children.size(); ++i)
                {
                    queue.push(children[i]);
                }
            }
        }
    }

private:
    Tree_node* _root;
};


__END_YAFRAY

#endif
