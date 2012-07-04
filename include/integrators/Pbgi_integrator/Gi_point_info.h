#ifndef GI_POINT_INFO_H
#define GI_POINT_INFO_H

#include <yafray_config.h>
#include <core_api/color.h>
#include <core_api/vector3d.h>


__BEGIN_YAFRAY

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
};

__END_YAFRAY

#endif // GI_POINT_INFO_H
