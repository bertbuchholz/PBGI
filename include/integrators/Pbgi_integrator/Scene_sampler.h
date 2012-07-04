#ifndef SCENE_SAMPLER_H
#define SCENE_SAMPLER_H

#include <yafray_config.h>

#include <core_api/scene.h>
#include <core_api/surface.h>
#include <utilities/mcqmc.h>

__BEGIN_YAFRAY




struct pbgi_sample_t
{
    triangle_t const* tri_pointer;
    point3d_t position;
    float area;
    intersectData_t intersect_data;
};




float get_total_scene_area(std::vector<triangle_t const*> const& triangles);
std::vector<triangle_t const*> get_scene_triangles(std::map<objID_t, objData_t> const& meshes);
float generate_histogram(std::vector<pbgi_sample_t> const& samples, float const min_radius);



class Scene_sampler
{
    public:
    virtual std::vector<pbgi_sample_t> generate_samples(std::vector<triangle_t const*> const& triangles) = 0;
    virtual float get_widest_gap() { assert(false); return -1.0f; }
};


class Scene_sampler_cdf : public Scene_sampler
{
public:
    Scene_sampler_cdf(int const number_of_samples) :
        _number_of_samples(number_of_samples)
    { }

    std::vector<pbgi_sample_t> generate_samples(std::vector<triangle_t const*> const& triangles);

private:
    int   _number_of_samples;
};



__END_YAFRAY

#endif // SCENE_SAMPLER_H
