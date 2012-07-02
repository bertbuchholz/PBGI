#ifndef SCENE_SAMPLER_H
#define SCENE_SAMPLER_H

#include <yafray_config.h>

#include <core_api/scene.h>
#include <core_api/surface.h>
#include <utilities/mcqmc.h>
#include <integrators/pointbased_gi.h>

__BEGIN_YAFRAY


float get_total_scene_area(std::vector<triangle_t const*> const& triangles);
std::vector<triangle_t const*> get_scene_triangles(std::map<objID_t, objData_t> const& meshes);
float generate_histogram(std::vector<pbgi_sample_t> const& samples, float const min_radius);





class Scene_sampler
{
    public:
    virtual ~Scene_sampler() {}
    virtual std::vector<pbgi_sample_t> generate_samples(std::vector<triangle_t const*> const& triangles) = 0;
    virtual float get_widest_gap() { assert(false); return -1.0f; }
};


class Scene_sampler_darts_hash : public Scene_sampler
{
public:
    Scene_sampler_darts_hash(float const min_radius, int const number_of_samples, float const cell_size_factor = 2.0f, float const bin_count_factor = 0.5f) :
        _min_radius(min_radius),
        _number_of_samples(number_of_samples),
        _cell_size_factor(cell_size_factor),
        _bin_count_factor(bin_count_factor)
    { }

    std::vector<pbgi_sample_t> generate_samples(std::vector<triangle_t const*> const& triangles);
    virtual float get_widest_gap() { return _widest_gap; }

private:
    float _min_radius;
    int   _number_of_samples;
    float _cell_size_factor;
    float _bin_count_factor;
    float _widest_gap;
};

class Scene_sampler_darts : public Scene_sampler
{
public:
    Scene_sampler_darts(float const min_radius, int const number_of_samples) :
        _min_radius(min_radius),
        _number_of_samples(number_of_samples)
    { }

    std::vector<pbgi_sample_t> generate_samples(std::vector<triangle_t const*> const& triangles);

private:
    float _min_radius;
    int   _number_of_samples;
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

class Scene_sampler_suicide : public Scene_sampler
{
public:
    Scene_sampler_suicide(int const number_of_samples, float const desired_radius) :
        _number_of_samples(number_of_samples),
        _desired_radius(desired_radius)
    { }

    std::vector<pbgi_sample_t> generate_samples(std::vector<triangle_t const*> const& triangles);

private:
    int _number_of_samples;
    float _desired_radius;
};

class Scene_sampler_reyes : public Scene_sampler
{
public:
    Scene_sampler_reyes(float const max_solid_angle, vector3d_t const& camera_pos, std::map<objID_t, objData_t> & meshes, std::vector<vector3d_t> * debug_new_triangles) :
        _meshes(meshes),
        _camera_pos(camera_pos),
        _max_solid_angle(max_solid_angle),
        _debug_new_triangles(debug_new_triangles)
    { }

    std::vector<pbgi_sample_t> generate_samples(std::vector<triangle_t const*> const& triangles);

private:
    std::map<objID_t, objData_t> const& _meshes;
    vector3d_t _camera_pos;
    float _max_solid_angle;
    std::vector<vector3d_t> * _debug_new_triangles;
};






__END_YAFRAY

#endif // SCENE_SAMPLER_H
