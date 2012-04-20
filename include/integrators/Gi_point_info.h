#ifndef GI_POINT_INFO_H
#define GI_POINT_INFO_H

#include <yafray_config.h>

//#include <utilities/spherical_harmonics.h>

__BEGIN_YAFRAY

class Spherical_node_representation;
class Gi_point_base;

struct Gi_point_info
{
#ifdef DEBUG
    Gi_point_info() : splatted(false)
    {}
#endif

    enum Type { Node = 0, Far_surfel, Near_surfel };

    Type type;
    color_t color;
    vector3d_t position;
    vector3d_t receiver_position;
    float depth;
    vector3d_t direction; // from the receiving point to the GI point
    float solid_angle;
    vector3d_t disc_normal;
    float radius;
    float weight; // used for mixing when doing variational splatting
    float distance_from_center;
    Gi_point_base const* gi_point;
    Spherical_node_representation const* spherical_function;

#ifdef DEBUG
    mutable bool splatted;
#endif
};

__END_YAFRAY

#endif // GI_POINT_INFO_H
