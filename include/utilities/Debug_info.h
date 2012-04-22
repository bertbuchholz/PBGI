#ifndef DEBUG_INFO_H
#define DEBUG_INFO_H

#include <tr1/unordered_set>

#include <vector>

#include <yafray_config.h>

#include <yafraycore/timer.h>
#include <core_api/color.h>

#include <integrators/Gi_point_info.h>

__BEGIN_YAFRAY

template <class Data>
class Cube_raster_buffer;
class Splat_cube_raster_buffer;
class GiPoint;
class Gi_point_base;



struct Node_weight_pair
{
    Node_weight_pair(Gi_point_info const& p_info, float w, color_t const& c, int group_i = -1, float group_w = 0.0f) :
        point_info(p_info),
        weight(w),
        color(c),
        group_index(group_i),
        group_weight(group_w)
    {}

    Gi_point_info point_info;
    float weight;
    color_t color;
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
        num_node_checks(0),
        debug_return_type(""),
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
    std::vector<Gi_point_info> gi_point_infos;
    std::vector<Node_weight_pair> single_pixel_contributors;
    Splat_cube_raster_buffer * result_fb;
    int checked_cells;
    float time_add_point;
    int num_node_checks;
    timer_t my_timer;

    std::string debug_return_type;

    // debug settings, not to be reset
    int cube_plane, cube_x, cube_y;
};


__END_YAFRAY

#endif // DEBUG_INFO_H
