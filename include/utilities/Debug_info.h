#ifndef DEBUG_INFO_H
#define DEBUG_INFO_H

#include <tr1/unordered_set>

#include <yafray_config.h>

#include <yafraycore/timer.h>

__BEGIN_YAFRAY

class Cube_raster_buffer;
class GiPoint;



struct Node_weight_pair
{
    Node_weight_pair(GiPoint const* n, float w, int group_i = -1, float group_w = 0.0f) :
        node(n),
        weight(w),
        group_index(group_i),
        group_weight(group_w)
    {}

    GiPoint const* node;
    float weight;
    int group_index;
    float group_weight;
};

struct Comp_Node_weight_pair_by_weight
{
    bool operator() (Node_weight_pair const& i, Node_weight_pair const& j) { return i.weight < j.weight; }
};

struct Comp_Node_weight_pair_by_group_weight
{
    bool operator() (Node_weight_pair const& i, Node_weight_pair const& j) { return i.group_weight < j.group_weight; }
};


struct Debug_info
{
    Debug_info() :
        used_nodes(0),
        used_near_surfels(0),
        used_far_surfels(0),
        node_depth(-1),
        node_height(-1),
        color_by_depth(false),
        result_fb(NULL),
        checked_cells(0),
        time_add_point(0.0f),
        show_octree_node_area(false),
        num_node_checks(0),
        cube_plane(-1),
        cube_x(-1),
        cube_y(-1)
    {
        my_timer = timer_t();
        my_timer.addEvent("Surfel");
        my_timer.addEvent("Node");
        my_timer.addEvent("Traversal");
        my_timer.addEvent("Accumulating");
        my_timer.addEvent("NodeCheck");
    }

    // returned debug info, needs to be reset before use
    int used_nodes;
    int used_near_surfels;
    int used_far_surfels;
    int node_depth;
    int node_height;
    bool color_by_depth;
    std::tr1::unordered_set<GiPoint const*> gi_points;
    std::vector<Node_weight_pair> single_pixel_contributors;
    Cube_raster_buffer * result_fb;
    int checked_cells;
    float time_add_point;
    bool show_octree_node_area; // contrary to visible area
    int num_node_checks;
    timer_t my_timer;

    // debug settings, not to be reset
    int cube_plane, cube_x, cube_y;
};


__END_YAFRAY

#endif // DEBUG_INFO_H
