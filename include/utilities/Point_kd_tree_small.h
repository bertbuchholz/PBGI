#ifndef POINT_KD_TREE_SMALL_H
#define POINT_KD_TREE_SMALL_H

#include <core_api/bound.h>
#include <utilities/PointKdTree.h>

__BEGIN_YAFRAY

template <class Data>
struct Point_data
{
    // Point_data() {}
    Point_data(vector3d_t const& p, Data * d) : point(p), data(d) {}

    vector3d_t point;
    Data *     data;
};

template <class Data>
struct Compare_by_axis
{
    Compare_by_axis(int const a) : axis(a) {}

    bool operator() (Point_data<Data> const& a, Point_data<Data> const& b)
    {
        return a.point[axis] < b.point[axis];
    }

    int axis;
};



template <class Data, class Averager>
class Node
{
public:
    Node() : _data(NULL), _depth(-1), _parent(NULL)
    { }

    virtual void average() {}

    virtual void add_points(std::vector<vector3d_t> const& points, bound_t const& bound, int const depth) = 0;

    virtual void get_post_order_queue(std::vector<Node*> & list, int const depth, int const max_depth = -1) = 0;

    virtual bool is_leaf() const = 0;

    virtual Node *      get_child(int const i) = 0;
    virtual Node const* get_child(int const i) const = 0;

    virtual void average_node(Averager const& averager) = 0;

    // virtual bound_t get_bound() const = 0;

    Data * get_data() { return _data; }
    Data const* get_data() const { return _data; }

    void set_data(Data * data) { _data = data; }

    virtual std::vector<Data*> const* get_surfels() const { assert(false); return NULL; }
    virtual std::vector<Data*>      * get_surfels()       { assert(false); return NULL; }

    int get_depth() const { return _depth; }

    Node const* get_parent() const { assert(_parent); return _parent; }

    void set_parent(Node * parent) { _parent = parent; }

    virtual int get_shortest_distance_to_leaf(int const height) const = 0;

protected:
    Data * _data;
    int _depth;
    Node * _parent;
};

template <class Data, class Averager>
class Leaf_node : public Node<Data, Averager>
{
public:
    void add_points(std::vector<vector3d_t> const& points, bound_t const& /* bound */, int const depth)
    {
        // assert(points.size() == 1);
        for (size_t i = 0; i < points.size(); ++i)
        {
            Data * d = new Data;
            d->pos = points[i];
            _surfels.push_back(d);
        }

        // this->_data = points[0].data;
        this->_depth = depth;
    }

    void get_post_order_queue(std::vector<Node<Data, Averager>*> & list, int const /* depth */, int const /* max_depth = -1 */)
    {
        list.push_back(this);
    }

    bool is_leaf() const { return true; }

    void average_node(Averager const& averager)
    {
        // this->_data = averager.average_leaf(_surfels);
    }

    std::vector<Data*> const* get_surfels() const { return &_surfels; }
    std::vector<Data*>      * get_surfels()       { return &_surfels; }

    int get_shortest_distance_to_leaf(int const height) const { return height; }

    // shouldn't be called
    bound_t get_bound() const { assert(false); return bound_t(); }


    Node<Data, Averager> *      get_child(int const i)       { return NULL; }
    Node<Data, Averager> const* get_child(int const i) const { return NULL; }

private:
    std::vector<Data*> _surfels;
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

template <class Data, class Averager, int Arity = 2>
class Inner_node : public Node<Data, Averager>
{
public:
    Inner_node()
    {
        for (int i = 0; i < Arity; ++i)
        {
            _children[i] = NULL;
        }
    }


    void add_points_any_arity(std::vector<vector3d_t> const& points, bound_t const& bound, int const depth)
    {
        this->_depth = depth;

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
                _children[i] = NULL;
            }
            else
            {
                // if (new_points[i].size() == 1)
                int const longest_axis = new_bounds[i].largestAxis();

                if (new_bounds[i].get_length(longest_axis) < 0.1f || new_points[i].size() == 1)
                {
                    _children[i] = new Leaf_node<Data, Averager>;
                }
                else
                {
                    _children[i] = new Inner_node;
                }

                _children[i]->set_parent(this);

                _children[i]->add_points(new_points[i], new_bounds[i], depth + 1);
                std::vector<vector3d_t>().swap(new_points[i]);
            }
        }
    }

    void add_points_binary(std::vector<vector3d_t> const& points, bound_t const& bound, int const depth)
    // void add_points(std::vector<vector3d_t> const& points, std::vector<Data*> const& data_per_point, bound_t const& bound)
    {
//        calc_points_bounding_box();
//        Point cell_extent = get_points_extent();

        this->_depth = depth;

        // _bound = bound;

        int const split_axis = bound.largestAxis();
        // int const split_axis = get_longest_axis_3d(cell_extent);

        std::vector<Point_data<Data> > points_sorted_by_axis = points;

        std::sort(points_sorted_by_axis.begin(), points_sorted_by_axis.end(), Compare_by_axis<Data>(split_axis));

        int const median_index = int(points_sorted_by_axis.size() / 2.0f);

        std::vector<vector3d_t> left_points(points_sorted_by_axis.begin(), points_sorted_by_axis.begin() + median_index);
        std::vector<vector3d_t> right_points(points_sorted_by_axis.begin() + median_index, points_sorted_by_axis.end());

        std::vector<vector3d_t>().swap(points_sorted_by_axis);

        float const split_position =
                0.5f * (left_points.back())[split_axis] +
                0.5f * (right_points.front())[split_axis];

        bound_t left_bound = bound;
        left_bound.setMax(split_position, split_axis);

        bound_t right_bound = bound;
        right_bound.setMin(split_position, split_axis);


        if (left_points.size() == 1)
        {
            _children[0] = new Leaf_node<Data, Averager>;
        }
        else
        {
            _children[0] = new Inner_node<Data, Averager>;
        }

        if (right_points.size() == 1)
        {
            _children[1] = new Leaf_node<Data, Averager>;
        }
        else
        {
            _children[1] = new Inner_node<Data, Averager>;
        }

        _children[0]->set_parent(this);
        _children[1]->set_parent(this);

        _children[0]->add_points(left_points, left_bound, depth + 1);
        std::vector<vector3d_t>().swap(left_points);

        _children[1]->add_points(right_points, right_bound, depth + 1);
        std::vector<vector3d_t>().swap(right_points);
    }

    void add_points(std::vector<vector3d_t> const& points, bound_t const& bound, int const depth)
    {
        add_points_any_arity(points, bound, depth);
    }

    void get_post_order_queue(std::vector<Node<Data, Averager>*> & list, int const depth, int const max_depth = -1)
    {
        if (max_depth == -1 || (max_depth != -1 && depth < max_depth))
        {
            for (int i = 0; i < Arity; ++i)
            {
                if (_children[i])
                {
                    _children[i]->get_post_order_queue(list, depth + 1, max_depth);
                }
            }
        }

        list.push_back(this);
    }

    bool is_leaf() const { return false; }

    Node<Data, Averager> *      get_child(int const i)       { return _children[i]; }
    Node<Data, Averager> const* get_child(int const i) const { return _children[i]; }


    void average_node(Averager const& averager)
    {
        // TODO: change bounds to sum of childrens' BB -> bounds in GiPoint, so no bounds needed (?) in the tree

        std::vector<Data*> leaf_data;

        for (int i = 0; i < Arity; ++i)
        {
            if (_children[i])
            {
                assert(_children[i]->get_data() != NULL);
                leaf_data.push_back(_children[i]->get_data());
            }
        }

        this->_data = averager.average_node(leaf_data);
    }

    int get_shortest_distance_to_leaf(int const height) const
    {
        std::vector<int> heights;

        for (int i = 0; i < Arity; ++i)
        {
            if (_children[i])
            {
                heights.push_back(_children[i]->get_shortest_distance_to_leaf(height + 1));
            }
        }

        return *std::min_element(heights.begin(), heights.end());
    }

    // bound_t get_bound() const { return _bound; }

private:
    Node<Data, Averager> * _children[Arity];

    // bound_t _bound;
};

template <class Data, class Averager, int n_Arity>
class Point_kd_tree_small
{
public:
    static int const Arity = n_Arity;

    typedef Node<Data, Averager> Tree_node;
    typedef Inner_node<Data, Averager, Arity> Tree_inner_node;
    typedef Leaf_node<Data, Averager> Tree_leaf_node;

    void build_tree(std::vector<vector3d_t> const& points, bound_t const& bound)
    {
        _root = new Tree_inner_node;
        _root->add_points(points, bound, 0);
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
                for (int i = 0; i < Arity; ++i)
                {
                    if (node->get_child(i))
                    {
                        queue.push(node->get_child(i));
                    }
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
                for (int i = 0; i < Arity; ++i)
                {
                    if (node->get_child(i))
                    {
                        queue.push(node->get_child(i));
                    }
                }
            }
        }

        return result;
    }

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
                for (int i = 0; i < Arity; ++i)
                {
                    if (node->get_child(i))
                    {
                        queue.push(node->get_child(i));
                    }
                }
            }
        }
    }

private:
    Tree_node* _root;
};


template <class Data, class Averager>
bool is_node_data_behind_plane(Node<Data, Averager> const* node, vector3d_t const& pos, vector3d_t const& normal, float const bias)
{
    // XXX: possible problem when using threads
    static int corners[][3] = {
        {0, 0, 0},
        {0, 1, 0},
        {0, 0, 1},
        {0, 1, 1},
        {1, 0, 0},
        {1, 1, 0},
        {1, 0, 1},
        {1, 1, 1}
    };

    vector3d_t const minmax[] = { vector3d_t(node->get_bound().a), vector3d_t(node->get_bound().g) };

    for (int i = 0; i < 8; ++i)
    {
        int const x = corners[i][0];
        int const y = corners[i][1];
        int const z = corners[i][2];

        Point const corner((minmax[x])[0],
                           (minmax[y])[1],
                           (minmax[z])[2]);

        float const front = (corner - pos).normalize() * normal;

        if (front > bias) return false;
    }

    return true;
}



__END_YAFRAY

#endif // POINT_KD_TREE_SMALL_H
