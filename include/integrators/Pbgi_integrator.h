#ifndef PBGI_INTEGRATOR_H
#define PBGI_INTEGRATOR_H


#include <tr1/unordered_map>

#include <core_api/environment.h>
#include <core_api/material.h>
#include <core_api/mcintegrator.h>
#include <core_api/background.h>
#include <core_api/light.h>
#include <yafraycore/meshtypes.h>
#include <utilities/mcqmc.h>
#include <integrators/Pbgi_integrator/spherical_harmonics.h>

#include <integrators/Pbgi_integrator/Pbgi_tree_node_data.h>
#include <integrators/Pbgi_integrator/Pbgi_tree.h>
#include <integrators/Pbgi_integrator/spherical_harmonics.h>

__BEGIN_YAFRAY

class pbgi_sample_t;


class Vector_hash {
public:
    size_t operator()(vector3d_t const& v) const
    {
        size_t h1 = std::tr1::hash<float>()(v.x);
        size_t h2 = std::tr1::hash<float>()(v.y);
        size_t h3 = std::tr1::hash<float>()(v.z);
        return h1 ^ ( (h2 ^ (h3 << 1)) << 2 );
    }
};


class YAFRAYPLUGIN_EXPORT Pbgi_integrator_t : public mcIntegrator_t
{
public:
    typedef Pbgi_tree<Gi_point_inner_node, Gi_point_surfel, 8> MyTree;

    Pbgi_integrator_t(bool transpShad=false, int shadowDepth=4, int rayDepth=6);
    virtual bool preprocess();
    virtual colorA_t integrate(renderState_t &state, diffRay_t &ray) const;
    static integrator_t* factory(paraMap_t &params, renderEnvironment_t &render);

    void generate_gi_points(MyTree & tree, int const number_of_samples);

    void generate_gi_points_data_2(std::vector<MyTree::Tree_node*> const& nodes,
                                   std::tr1::unordered_map<vector3d_t, pbgi_sample_t*, Vector_hash> const& point_to_sampling_point_map,
                                   float const area,
                                   float const radius,
                                   int const treat_depth_as_leaf = -1);


    virtual void cleanup();

    // debug hack
    MyTree _point_tree;

private:
    int surfel_samples;
    bool indirectOnly;

    float _solid_angle_factor;
    int raster_buffer_resolution;

    float surfel_near_threshold;

    Spherical_harmonics_factory<color_t> _spherical_function_color_factory;
    Spherical_harmonics_factory<float>   _spherical_function_area_factory;

    bool _variational;


    Cube_raster_buffer< Spherical_harmonics<float> > _precalculated_sf;

    float _disc_scale_factor;
};

__END_YAFRAY




#endif // PBGI_INTEGRATOR_H
