#ifndef DEBUG_INFO_H
#define DEBUG_INFO_H

__BEGIN_YAFRAY

class Cube_raster_buffer;
class GiPoint;

struct Node_weight_pair
{
    Node_weight_pair(GiPoint const* n, float w) : node(n), weight(w) {}

    GiPoint const* node;
    float weight;
};

struct Debug_info
{
    Debug_info() :
        used_nodes(0),
        used_leafs_rays(0),
        used_leafs_splat(0),
        node_depth(-1),
        color_by_depth(false),
        cube_plane(-1),
        cube_x(-1),
        cube_y(-1),
        result_fb(NULL),
        checked_cells(0),
        time_add_point(0.0f)
    {}

    // returned debug info, needs to be reset before use
    int used_nodes;
    int used_leafs_rays;
    int used_leafs_splat;
    int node_depth;
    bool color_by_depth;
    std::vector<GiPoint const*> gi_points;
    std::vector<Node_weight_pair> single_pixel_contributors;
    Cube_raster_buffer * result_fb;
    int checked_cells;
    float time_add_point;

    // debug settings, not to be reset
    int cube_plane, cube_x, cube_y;
};


__END_YAFRAY

#endif // DEBUG_INFO_H
