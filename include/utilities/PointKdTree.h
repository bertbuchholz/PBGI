#ifndef POINTKDTREE_H
#define POINTKDTREE_H

#include <iostream>
#include <vector>
#include <cmath>
#include <queue>
#include <list>

#include <yafray_constants.h>

//#include <core_api/color.h>
#include <core_api/vector3d.h>
#include <utilities/mcqmc.h>

__BEGIN_YAFRAY

// - Dimensions, size of each dimension
// - Max depth or max amount of points per cell
// - Enter a point into the tree
// - When adding a point, check first if the max amount of points is
//   excessed:
//   - If not, add the point to the cell.
//   - If yes, split the cell and redistribute the points into the new
//     cells.

template <class Point>
inline int get_longest_axis_2d(Point const& cell)
{
    int longest_axis = 0;
    if      (std::abs(cell[1]) > std::abs(cell[0]))
    {
        longest_axis = 1;
    }
    return longest_axis;
}

template <class Point>
inline int get_longest_axis_3d(Point const& cell)
{
    int longest_axis = 0;
    if      (std::abs(cell[1]) > std::abs(cell[0]) &&
             std::abs(cell[1]) > std::abs(cell[2]))
    {
        longest_axis = 1;
    }
    else if (std::abs(cell[2]) > std::abs(cell[0]))
    {
        longest_axis = 2;
    }
    return longest_axis;
}


template < class Point, int dim, class Data = int >
class Point_kd_tree
{
public:
    typedef Point_kd_tree Tree_node;

        Point_kd_tree(Point const& min, Point const& max, int maxDepth, int maxPoints) :
        _cellMin(min),
        _cellMax(max),
        _parent(NULL),
        _isLeaf(true),
        _depth(0),
        _maxDepth(maxDepth),
        _maxPoints(maxPoints),
        _childrenCount(2)
        { }

    Point_kd_tree() :
        _parent(NULL),
        _isLeaf(true),
        _depth(0),
        _maxDepth(10),
        _maxPoints(1),
        _childrenCount(2)
    { }

    void getLeafs(std::vector<Point_kd_tree const*>& leafs) const
    {
        if (!_isLeaf)
        {
            for (int i = 0; i < _childrenCount; ++i)
            {
                _children[i].getLeafs(leafs);
            }
        }
        else
        {
            leafs.push_back(this);
        }
    }

    bool is_empty() const
    {
        return (_points.size() == 0);
    }

    bool has_children() const
    {
        return (_children.size() > 0);
    }

    void setDepth(int d)
    {
        _depth = d;
    }

    int get_depth() const
    {
        return _depth;
    }

    bool is_leaf() const
    {
        return _isLeaf;
    }

    void setMin(Point const& min)
    {
        _cellMin = min;
    }

    void setMax(Point const& max)
    {
        _cellMax = max;
    }

    Point const& getMin() const
    {
        return _cellMin;
    }

    Point const& getMax() const
    {
        return _cellMax;
    }

    Point_kd_tree * get_parent() const
    {
        return _parent;
    }

    std::vector<Point_kd_tree> const& getChildren() const
    {
        return _children;
    }

    std::vector<Point_kd_tree> & getChildren()
    {
        return _children;
    }

    Point_kd_tree* get_child(int const i)
    {
        return &_children[i];
    }

    Point_kd_tree const* get_child(int const i) const
    {
        return &_children[i];
    }

    std::vector<Point> const& getPoints() const
    {
        return _points;
    }

    std::vector<Point> & getPoints()
    {
        return _points;
    }

    std::vector<Data> & get_points_data()
    {
        return _data;
    }

    std::vector<Data> const& get_points_data() const
    {
        return _data;
    }

    float getRadius() const
    {
        return (_cellMax - _cellMin).length() * 0.5f;
    }

    float get_points_radius() const
    {
        return (_points_max - _points_min).length() * 0.5f;
    }

    Point get_points_center() const
    {
        return (_points_max + _points_min) * 0.5f;
    }

    Point getCenter() const
    {
        return (_cellMax + _cellMin) * 0.5f;
    }

    Point get_points_extent() const
    {
        return (_points_max - _points_min);
    }

    void get_points_bounding_box(Point & min, Point & max) const
    {
        min = Point( 1e10f);
        max = Point(-1e10f);

        for (unsigned int i = 0; i < _points.size(); ++i)
        {
            Point const& p = _points[i];

            for (int axis = 0; axis < dim; ++axis)
            {
                if (p[axis] < min[axis]) min[axis] = p[axis];
                if (p[axis] > max[axis]) max[axis] = p[axis];
            }
        }

        assert(min[0] <= max[0] && min[1] <= max[1] && min[2] <= max[2]);
    }

    void calc_points_bounding_box()
    {
        get_points_bounding_box(_points_min, _points_max);
    }

    bool isNodeBehindPlane(Point const& pos, Point const& normal) const
    {
        int corners[][3] = {
            {0, 0, 0},
            {0, 1, 0},
            {0, 0, 1},
            {0, 1, 1},
            {1, 0, 0},
            {1, 1, 0},
            {1, 0, 1},
            {1, 1, 1}
        };

        Point const* minmax[] = { &_cellMin, &_cellMax };
        // Point minmax[2];

        for (int i = 0; i < 8; ++i)
        {
            int const x = corners[i][0];
            int const y = corners[i][1];
            int const z = corners[i][2];

            Point corner((*minmax[x])[0],
                         (*minmax[y])[1],
                         (*minmax[z])[2]);


            float front = (corner - pos) * normal;

            if (front > 0.0f) return false;
        }

        return true;
    }

    bool isNodeDataBehindPlane(Point const& pos, Point const& normal, float const bias) const
    {
        int corners[][3] = {
            {0, 0, 0},
            {0, 1, 0},
            {0, 0, 1},
            {0, 1, 1},
            {1, 0, 0},
            {1, 1, 0},
            {1, 0, 1},
            {1, 1, 1}
        };

        Point const* minmax[] = { &_points_min, &_points_max };

        for (int i = 0; i < 8; ++i)
        {
            int const x = corners[i][0];
            int const y = corners[i][1];
            int const z = corners[i][2];

            Point corner((*minmax[x])[0],
                         (*minmax[y])[1],
                         (*minmax[z])[2]);


            float front = (corner - pos).normalize() * normal;

            if (front > bias) return false;
        }

        return true;
    }

    bool is_node_data_behind_plane_approx(Point const& pos, Point const& normal, float const radius, float const bias) const
    {
        Point p = (get_points_center() + radius * normal) - pos;
        p.normalize();

        float const cos_normal_point_dir = p * normal;

        if (cos_normal_point_dir > bias) return false;

        return true;
    }

    void addPoint(Point const& point, Data const& data = Data())
    {
        _points.push_back(point);
        _data.push_back(data);
    }

    void prepare_containers(std::size_t const size)
    {
        _points.resize(size);
        _data.resize(size);
    }

    void add_point_prepared(int const index, Point const& point, Data const& data = Data())
    {
        _points[index] = point;
        _data[index] = data;
    }

    struct Point_data
    {
        // Point_data() {}
        Point_data(Point const* p, Data const* d) : point(p), data(d) {}

        Point const* point;
        Data  const* data;
    };


    struct Compare_by_axis
    {
        Compare_by_axis(int const a) : axis(a) {}

        bool operator() (Point_data const& a, Point_data const& b)
        {
            return (*a.point)[axis] < (*b.point)[axis];
        }

        int axis;
    };

    // actually build the tree from the previously added points
    void finalize_points()
    {
        // std::cout << "Depth: " << _depth << " adding point: " << point << std::endl;

        assert(_points.size() == _data.size());

        if (_isLeaf)
        {
            assert(_children.size() == 0);

            // add point into cell directly or if full, split the cell
            // and add then
            if (int(_points.size()) > _maxPoints && _depth < _maxDepth)
            {
                // split();
                _children.resize(2, Point_kd_tree(_cellMin, _cellMax, _maxDepth, _maxPoints));

                get_points_bounding_box(_points_min, _points_max);

                // Point cell_extent = _cellMax - _cellMin;
                Point cell_extent = _points_max - _points_min;

                assert(dim == 3 || dim == 2);

                int split_axis;

                if (dim == 3)
                {
                    split_axis = get_longest_axis_3d(cell_extent);
                }
                else if (dim == 2)
                {
                    split_axis = get_longest_axis_2d(cell_extent);
                }

                //std::vector<Point_data> points_sorted_by_axis(_points.size());
                std::vector<Point_data> points_sorted_by_axis;

                for (unsigned int i = 0; i < _points.size(); ++i)
                {
                    points_sorted_by_axis.push_back(Point_data(&_points[i], &_data[i]));
                }

                std::sort(points_sorted_by_axis.begin(), points_sorted_by_axis.end(), Compare_by_axis(split_axis));

                int const median_index = int(points_sorted_by_axis.size() / 2.0f);

                for (int i = 0; i < median_index; ++i)
                {
                    _children[0].addPoint(*points_sorted_by_axis[i].point, *points_sorted_by_axis[i].data);
                }

                for (unsigned int i = median_index; i < points_sorted_by_axis.size(); ++i)
                {
                    _children[1].addPoint(*points_sorted_by_axis[i].point, *points_sorted_by_axis[i].data);
                }


                float split_position =
                        0.5f * (*points_sorted_by_axis[median_index - 1].point)[split_axis] +
                        0.5f * (*points_sorted_by_axis[median_index    ].point)[split_axis];


/*

                float split_position = cell_extent[split_axis] / 2.0f + _cellMin[split_axis];

                for (unsigned int i = 0; i < _points.size(); ++i)
                {
                    if (_points[i][split_axis] < split_position)
                    {
                        _children[0].addPoint(_points[i], _data[i]);
                    }
                    else
                    {
                        _children[1].addPoint(_points[i], _data[i]);
                    }
                }
*/

                _points.clear();
                _data.clear();
                _isLeaf = false;

                _children[0]._cellMax[split_axis] = split_position;
                _children[1]._cellMin[split_axis] = split_position;

                _children[0]._parent = this;
                _children[1]._parent = this;

                _children[0].setDepth(_depth + 1);
                _children[1].setDepth(_depth + 1);

                _children[0].finalize_points();
                _children[1].finalize_points();
            }

            assert(_points.size() == _data.size());
        }
        else
        {
            // don't handle points that are inserted after the first initialization
            assert(false);
        }
    }



    void split()
    {
        assert(_isLeaf && _points.size() == _data.size());

        _children.resize(2, Point_kd_tree(_cellMin, _cellMax, _maxDepth, _maxPoints));

        calc_points_bounding_box();
        Point cell_extent = get_points_extent();

        assert(dim == 3 || dim == 2);

        int split_axis;

        if (dim == 3)
        {
            split_axis = get_longest_axis_3d(cell_extent);
        }
        else if (dim == 2)
        {
            split_axis = get_longest_axis_2d(cell_extent);
        }

        //std::vector<Point_data> points_sorted_by_axis(_points.size());
        std::vector<Point_data> points_sorted_by_axis;

        for (unsigned int i = 0; i < _points.size(); ++i)
        {
            points_sorted_by_axis.push_back(Point_data(&_points[i], &_data[i]));
        }

        std::sort(points_sorted_by_axis.begin(), points_sorted_by_axis.end(), Compare_by_axis(split_axis));

        int const median_index = int(points_sorted_by_axis.size() / 2.0f);

        for (int i = 0; i < median_index; ++i)
        {
            _children[0].addPoint(*points_sorted_by_axis[i].point, *points_sorted_by_axis[i].data);
        }

        for (unsigned int i = median_index; i < points_sorted_by_axis.size(); ++i)
        {
            _children[1].addPoint(*points_sorted_by_axis[i].point, *points_sorted_by_axis[i].data);
        }


        float split_position =
                0.5f * (*points_sorted_by_axis[median_index - 1].point)[split_axis] +
                0.5f * (*points_sorted_by_axis[median_index    ].point)[split_axis];

        // _points.clear();
        std::vector<Point>().swap(_points);

        // _data.clear();
        std::vector<Data>().swap(_data);

        _isLeaf = false;

        _children[0]._cellMax[split_axis] = split_position;
        _children[1]._cellMin[split_axis] = split_position;

        _children[0]._parent = this;
        _children[1]._parent = this;

        _children[0].setDepth(_depth + 1);
        _children[1].setDepth(_depth + 1);
    }


    // actually build the tree from the previously added points
    static void finalize_points_iterative(Point_kd_tree * root)
    {
        // std::cout << "Depth: " << _depth << " adding point: " << point << std::endl;

        std::queue<Point_kd_tree*> queue;
        queue.push(root);

        while (!queue.empty())
        {
            Point_kd_tree * node = queue.front();
            queue.pop();

            if (node->is_leaf())
            {
                assert(node->getChildren().size() == 0);

                // add point into cell directly or if full, split the cell
                // and add then
                if (int(node->getPoints().size()) > node->get_tree_max_points_per_node() && node->get_depth() < node->get_tree_max_depth())
                {
                    node->split();

                    std::vector<Point_kd_tree> & children = node->getChildren();
                    queue.push(&children[0]);
                    queue.push(&children[1]);
                }
            }
            else
            {
                // don't handle points that are inserted after the first initialization
                assert(false);
            }
        }
    }




    bool isPointIn(Point const& point) const
    {
        // return (point >= _cellMin && point <= _cellMax);

        for (unsigned int i = 0; i < dim; ++i)
        {
            if (point[i] < _cellMin[i] || point[i] > _cellMax[i]) return false;
        }

        return true;
    }

    friend std::ostream& operator<< (std::ostream& stream, Point_kd_tree const& q)
    {
        // stream << q._depth << " " << q._cellMin << " " << q._cellMax << std::endl;
        stream << q._depth << " " << q._children.size() << std::endl;

        /*
  for (unsigned int i = 0; i < q._points.size(); ++i)
  {
            std::cout << "[" << q._points[i] << "] Data: " << q._data[i] << std::endl;
  }
        */

        if (!q._isLeaf)
        {
            std::cout << "Children:" << std::endl;

            for (int i = 0; i < q._childrenCount; ++i)
            {
                stream << q._children[i] << std::endl;
            }
        }

        return stream;
    }

    Point_kd_tree const& query(Point const& point) const
    {
        if (!_isLeaf)
        {
            for (int j = 0; j < _childrenCount; ++j)
            {
                if (_children[j].isPointIn(point))
                {
                    return _children[j].query(point);
                }
            }
        }

        return *this;
    }

    void get_post_order_queue(std::vector<Point_kd_tree*> & list, int const max_depth = -1)
    {
        if (max_depth == -1 || (max_depth != -1 && _depth < max_depth))
        {
            for (unsigned int i = 0; i < _children.size(); ++i)
            {
                _children[i].get_post_order_queue(list, max_depth);
            }
        }

        list.push_back(this);
    }

    void get_pre_order_queue(std::vector<Point_kd_tree*> & list, int const max_depth = -1)
    {
        list.push_back(this);

        if (max_depth == -1 || (max_depth != -1 && _depth < max_depth))
        {
            for (unsigned int i = 0; i < _children.size(); ++i)
            {
                _children[i].get_pre_order_queue(list, max_depth);
            }
        }
    }

    void printAveraged() const
    {
        std::cout << "p: d: " << _depth << ", " << _isLeaf << ", " << _data.size() << ")" << std::endl;

        for (unsigned int i = 0; i < _children.size(); ++i)
        {
            _children[i].printAveraged();
        }
    }

    template <class Averager>
    void averageData(Averager averager)
    {
        std::vector<Point_kd_tree*> nodeQueue;
        getPostOrderQueue(nodeQueue, this);

        std::cout << "postorder.size: " << nodeQueue.size() << std::endl;

        for (unsigned int i = 0; i < nodeQueue.size(); ++i)
        {
            Point_kd_tree & node = *nodeQueue[i];

            // std::cout << "node depth: " << node._depth << " node: " << &node << std::endl;

            if (node.is_leaf())
//            if (node.getIsLeaf() && node._data.size() > 0)
            {
                if (node._data.size() > 0)
                {
                    node._clusteredData = averager.average(node._data);
                }
                else
                {
                    node._clusteredData = NULL;
                }

                // assert(node._averagedData.radius > 0.0f);
            }
            else if (!node.is_leaf())
            {
                std::vector<Data> leafData;

                /*
                for (unsigned int j = 0; j < node._children.size(); ++j)
                {
                    if (node._children[j].getIsLeaf() && node._children[j]._data.size() > 0)
                    {
                        leafData.push_back(node._children[j]._averagedData);
                    }
                }
                */

                for (unsigned int j = 0; j < node._children.size(); ++j)
                {
                    if (node._children[j]._clusteredData == NULL) continue;
                    // if (node._children[j].getIsLeaf() && node._children[j]._data.size() == 0) continue;

                    leafData.push_back(node._children[j]._clusteredData);
                }

                node._clusteredData = averager.average(leafData);

                // node._clusteredData->pos = (node._cellMin + node._cellMax) * 0.5f;
            }
        }
    }


    template <class Averager>
    void average_node(Averager averager)
    {
        if (is_leaf())
        {
            assert(_data.size() == 1);

            _clusteredData = averager.average(_data);
        }
        else
        {
            assert(_children.size() == 2);
            assert(_children[0].get_data() != NULL);
            assert(_children[1].get_data() != NULL);

            std::vector<Data> leafData;
            leafData.push_back(_children[0]._clusteredData);
            leafData.push_back(_children[1]._clusteredData);

            _clusteredData = averager.average(leafData);
        }
    }


    Data const& get_data() const
    {
        return _clusteredData;
    }

    Data & get_data()
    {
        return _clusteredData;
    }


    void set_clustered_data(Data const& data)
    {
        _clusteredData = data;
    }


    // get nodes in specified depth or all nodes for depth == -1
    std::vector<Point_kd_tree const*> get_nodes(int const depth = -1) const
    {
        std::vector<Point_kd_tree const*> result;

        std::queue<Point_kd_tree const*> queue;
        queue.push(this);

        while (!queue.empty())
        {
            Point_kd_tree const* node = queue.front();
            queue.pop();

            if (depth == -1 || node->_depth < depth)
            {
                for (unsigned int i = 0; i < node->_children.size(); ++i)
                {
                    queue.push(&(node->_children[i]));
                }
            }

            if (depth == -1 || node->_depth == depth)
            {
                result.push_back(node);
            }

            assert(depth == -1 || node->_depth <= depth);
        }

        return result;
    }

    std::vector<Point_kd_tree*> get_nodes(int const depth = -1)
    {
        std::vector<Point_kd_tree*> result;

        std::queue<Point_kd_tree*> queue;
        queue.push(this);

        while (!queue.empty())
        {
            Point_kd_tree * node = queue.front();
            queue.pop();

            if (depth == -1 || node->_depth < depth)
            {
                for (unsigned int i = 0; i < node->_children.size(); ++i)
                {
                    queue.push(&(node->_children[i]));
                }
            }

            if (depth == -1 || node->_depth == depth)
            {
                result.push_back(node);
            }

            assert(depth == -1 || node->_depth <= depth);
        }

        return result;
    }

    int getMaxDepth() const
    {
        int maxDepth = -1;

        std::queue<Point_kd_tree const*> queue;
        queue.push(this);

        while (!queue.empty())
        {
            Point_kd_tree const* node = queue.front();
            queue.pop();

            if (node->_depth > maxDepth)
            {
                maxDepth = node->_depth;
            }

            for (unsigned int i = 0; i < node->_children.size(); ++i)
            {
                queue.push(&(node->_children[i]));
            }
        }

        return maxDepth;
    }

    // the leaf closest to the root
    int getMinLeafDepth() const
    {
        int minDepth = 10000;

        std::queue<Point_kd_tree const*> queue;
        queue.push(this);

        while (!queue.empty())
        {
            Point_kd_tree const* node = queue.front();
            queue.pop();

            if (node->is_leaf() && node->_depth < minDepth)
            {
                minDepth = node->_depth;
            }

            for (unsigned int i = 0; i < node->_children.size(); ++i)
            {
                queue.push(&(node->_children[i]));
            }
        }

        return minDepth;
    }

    int getNodeCount() const
    {
        int count = 0;

        std::queue<Point_kd_tree const*> queue;
        queue.push(this);

        while (!queue.empty())
        {
            Point_kd_tree const* node = queue.front();
            queue.pop();

            ++count;

            for (unsigned int i = 0; i < node->_children.size(); ++i)
            {
                queue.push(&(node->_children[i]));
            }
        }

        return count;
    }

    // longest path to a leaf
    int get_node_height(int const height) const
    {
        if (is_leaf())
        {
            return height;
        }

        int const height_0 = _children[0].get_node_height(height + 1);
        int const height_1 = _children[1].get_node_height(height + 1);

        return std::max(height_0, height_1);
    }

    int get_shortest_distance_to_leaf(int const height = 0) const
    {
        if (is_leaf())
        {
            return height;
        }

        int const height_0 = _children[0].get_shortest_distance_to_leaf(height + 1);
        int const height_1 = _children[1].get_shortest_distance_to_leaf(height + 1);

        return std::min(height_0, height_1);
    }

    int get_tree_max_depth() const
    {
        return _maxDepth;
    }

    int get_tree_max_points_per_node() const
    {
        return _maxPoints;
    }


    static void test()
    {
        typedef Point_kd_tree<vector3d_t, 3> MyTree;
        MyTree tree(vector3d_t(0.0f, 0.0f, 0.0f), vector3d_t(1.0f, 0.0f, 0.0f), 10, 1);


        random_t my_random;

        for (int i = 0; i < 20; ++i)
        {
            tree.addPoint(vector3d_t(my_random(), my_random(), my_random()));
        }

        /*
        tree.addPoint(vector3d_t(0.0f, 0.0f, 0.0f));
        tree.addPoint(vector3d_t(0.5f, 0.0f, 0.0f));
        tree.addPoint(vector3d_t(1.0f, 0.0f, 0.0f));
        */

        tree.finalize_points();

        std::vector<MyTree*> nodes_in_depth_2 = tree.get_nodes(2);

        std::vector<MyTree*> nodes_with_depth_2_as_leafs;
        tree.get_post_order_queue(nodes_with_depth_2_as_leafs, 2);

        // std::vector<MyTree*> nodes = tree.getNodes();


        for (unsigned int i = 0; i < nodes_in_depth_2.size(); ++i)
        {
            MyTree * node = nodes_in_depth_2[i];
            std::cout << "i: " << i << " " << node << " " << node->get_depth() << std::endl;
        }

        for (unsigned int i = 0; i < nodes_with_depth_2_as_leafs.size(); ++i)
        {
            MyTree * node = nodes_with_depth_2_as_leafs[i];
            std::cout << "i: " << i << " " << node << " " << node->get_depth() << std::endl;
        }
    }

private:
    Point _cellMin;
    Point _cellMax;
    std::vector<Point_kd_tree> _children;
    Point_kd_tree* _parent;

    std::vector<Point> _points;
    std::vector<Data> _data;

    Point _points_min; // extent of the actual points data
    Point _points_max;

    Data _clusteredData;

    bool _isLeaf;
    int _depth;
    int _maxDepth;
    int _maxPoints;

    int _childrenCount;
};

/*
template <class Point>
void draw(QPainter & p, RegularBspTree<Point, 2> const& quadTree, Point const& min, Point const& max)
{
	// std::cout << "Drawing" << std::endl;

	p.setPen(Qt::black);

	float xRange = max[0] - min[0];
	float yRange = max[1] - min[1];

	QRect rect(
                std::floor((quadTree.getMin()[0] - min[0]) / xRange * p.window().width() + 0.5f),
                std::floor((quadTree.getMin()[1] - min[1]) / yRange * p.window().height() + 0.5f),
                std::floor((quadTree.getMax()[0] - quadTree.getMin()[0]) / xRange * p.window().width() + 0.5f),
                std::floor((quadTree.getMax()[1] - quadTree.getMin()[1]) / yRange * p.window().width() + 0.5f)
                );

	p.drawRect(rect);

	std::vector<Point> const& points = quadTree.getPoints();

	p.setPen(Qt::red);

	for (unsigned int i = 0; i < points.size(); ++i)
	{
		p.drawEllipse(
                    (points[i][0] - min[0]) / xRange * p.window().width(),
                    (points[i][1] - min[1]) / xRange * p.window().width(),
                    2, 2);
	}

	if (!quadTree.getIsLeaf())
	{
		std::vector<RegularBspTree<Point, 2> > const& children = quadTree.getChildren();

		for (unsigned int i = 0; i < children.size(); ++i)
		{
			draw(p, children[i], min, max);
		}
	}
}
*/


template <class Node>
bool is_node_data_behind_plane(Node const* node, vector3d_t const& pos, vector3d_t const& normal, float const bias)
{
    return node->isNodeDataBehindPlane(pos, normal, bias);
}

__END_YAFRAY

#endif // POINTKDTREE_H
