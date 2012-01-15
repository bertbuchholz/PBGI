#ifndef SERIALIZATIONHELPER_H
#define SERIALIZATIONHELPER_H

#include <boost/serialization/vector.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/export.hpp>

#include <core_api/vector3d.h>
#include <core_api/color.h>

namespace boost {
namespace serialization {

template<class Archive>
void serialize(Archive & ar, yafaray::vector3d_t & o, const unsigned int version)
{
    ar & o.x;
    ar & o.y;
    ar & o.z;
}

template<class Archive>
void serialize(Archive & ar, yafaray::color_t & o, const unsigned int version)
{
    ar & o.R;
    ar & o.G;
    ar & o.B;
}

}//namespace serialization
}//namespace boost

#endif // SERIALIZATIONHELPER_H
