#ifndef GI_POINT_INFO_H
#define GI_POINT_INFO_H

#include <yafray_config.h>

//#include <utilities/spherical_harmonics.h>

__BEGIN_YAFRAY

class Spherical_node_representation;

struct Gi_point_info
{
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
    Spherical_node_representation const* spherical_function;
};

__END_YAFRAY

#endif // GI_POINT_INFO_H
