#include <sstream>
#include <cassert>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <map>


#include <omp.h>

#include <bert/shared/Parameter.h>

#include <yafray_config.h>
#include <core_api/environment.h>
#include <core_api/material.h>
#include <core_api/mcintegrator.h>
#include <core_api/background.h>
#include <core_api/light.h>
#include <core_api/camera.h>
#include <yafraycore/meshtypes.h>
#include <utilities/mcqmc.h>
#include <utilities/spherical_harmonics.h>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/export.hpp>

#include <utilities/RegularBspTree.h>
#include <utilities/CubeRasterBuffer.h>
#include <utilities/sample_utils.h>
#include <utilities/Dictionary_generator.h>
#include <integrators/Dictionary_converter.h>
#include <integrators/pointbased_gi.h>
//#include <integrators/KMeans.h>
#include <integrators/kmeans.hpp>

#include <yafraycore/timer.h>

#include <integrators/Scene_sampler.h>

__BEGIN_YAFRAY


const color_t colors[] = {
    color_t(1.0f, 0.0f, 0.0f), // red
    color_t(0.0f, 1.0f, 0.0f), // blue
    color_t(1.0f, 1.0f, 0.0f), // yellow
    color_t(0.0f, 1.0f, 1.0f), // turquois
    color_t(1.0f, 0.0f, 1.0f), // violet
    color_t(0.5f, 0.5f, 0.5f), // gray
    color_t(1.0f, 1.0f, 1.0f), // white
    color_t(1.0f, 0.0f, 0.0f)  // red
};

std::vector<color_t> pbLighting_t::debug_colors(colors, colors + sizeof(colors)/sizeof(*colors));


pbLighting_t::pbLighting_t(bool transpShad, int shadowDepth, int rayDepth) : solid_angle_factor(0.5f), _bspTree(NULL)
{
    type = SURFACE;
    causRadius = 0.25;
    causDepth = 10;
    nCausPhotons = 100000;
    nCausSearch = 100;
    trShad = transpShad;
    usePhotonCaustics = false;
    sDepth = shadowDepth;
    rDepth = rayDepth;
    intpb = 0;
    integratorName = "PB";
    integratorShortName = "PBGI";
    do_load_gi_points = false;
    _color_dictionary_generator = NULL;
}


color_t pbLighting_t::estimateIncomingLight(renderState_t & state, light_t *light, const surfacePoint_t &sp, const unsigned int &loffs) const
{
    color_t col(0.f);
    bool shadowed;
    // unsigned int l_offs = loffs * loffsDelta;
    // const material_t *material = sp.material;
    ray_t lightRay;
    lightRay.from = sp.P;
    color_t lcol(0.f), scol;
    //float lightPdf;

    // handle lights with delta distribution, e.g. point and directional lights
    if( light->diracLight() )
    {
        if( light->illuminate(sp, lcol, lightRay) )
        {
            // ...shadowed...
            lightRay.tmin = YAF_SHADOW_BIAS; // < better add some _smart_ self-bias value...this is bad.
            shadowed = scene->isShadowed(state, lightRay);
            if (!shadowed)
            {
                // if(trShad) lcol *= scol;
                // color_t surfCol = material->eval(state, sp, wo, lightRay.dir, BSDF_ALL);
                // color_t transmitCol = scene->volIntegrator->transmittance(state, lightRay);
                // col += lcol * std::fabs(sp.N * lightRay.dir);
                col += lcol * std::max(0.0f, sp.N * lightRay.dir);
            }
        }
    }
    /*
    else // area light and suchlike
    {
        Halton hal2(2);
        Halton hal3(3);
        int n = light->nSamples();
        if(state.rayDivision > 1) n = std::max(1, n/state.rayDivision);
        float invNS = 1.f / (float)n;
        unsigned int offs = n * state.pixelSample + state.samplingOffs + l_offs;
        bool canIntersect=light->canIntersect();
        color_t ccol(0.0);
        lSample_t ls;

        hal2.setStart(offs-1);
        hal3.setStart(offs-1);

        for(int i=0; i<n; ++i)
        {
            // ...get sample val...
            ls.s1 = hal2.getNext();
            ls.s2 = hal3.getNext();

            if( light->illumSample (sp, ls, lightRay) )
            {
                // ...shadowed...
                lightRay.tmin = YAF_SHADOW_BIAS; // < better add some _smart_ self-bias value...this is bad.
                shadowed = (trShad) ? scene->isShadowed(state, lightRay, sDepth, scol) : scene->isShadowed(state, lightRay);

                if(!shadowed && ls.pdf > 1e-6f)
                {
                    if(trShad) ls.col *= scol;
                    color_t transmitCol = scene->volIntegrator->transmittance(state, lightRay);
                    ls.col *= transmitCol;
                    color_t surfCol = material->eval(state, sp, wo, lightRay.dir, BSDF_ALL);
                    if( canIntersect)
                    {
                        float mPdf = material->pdf(state, sp, wo, lightRay.dir, BSDF_GLOSSY | BSDF_DIFFUSE | BSDF_DISPERSIVE | BSDF_REFLECT | BSDF_TRANSMIT);
                        if(mPdf > 1e-6f)
                        {
                            float l2 = ls.pdf * ls.pdf;
                            float m2 = mPdf * mPdf;
                            float w = l2 / (l2 + m2);
                            ccol += surfCol * ls.col * std::fabs(sp.N*lightRay.dir) * w / ls.pdf;
                        }
                        else
                        {
                            ccol += surfCol * ls.col * std::fabs(sp.N*lightRay.dir) / ls.pdf;
                        }
                    }
                    else ccol += surfCol * ls.col * std::fabs(sp.N*lightRay.dir) / ls.pdf;
                }
            }
        }

        col += ccol * invNS;

        if(canIntersect) // sample from BSDF to complete MIS
        {
            color_t ccol2(0.f);

            hal2.setStart(offs-1);
            hal3.setStart(offs-1);

            for(int i=0; i<n; ++i)
            {
                ray_t bRay;
                bRay.tmin = MIN_RAYDIST; bRay.from = sp.P;

                float s1 = hal2.getNext();
                float s2 = hal3.getNext();
                float W = 0.f;

                sample_t s(s1, s2, BSDF_GLOSSY | BSDF_DIFFUSE | BSDF_DISPERSIVE | BSDF_REFLECT | BSDF_TRANSMIT);
                color_t surfCol = material->sample(state, sp, wo, bRay.dir, s, W);
                if( s.pdf>1e-6f && light->intersect(bRay, bRay.tmax, lcol, lightPdf) )
                {
                    shadowed = (trShad) ? scene->isShadowed(state, bRay, sDepth, scol) : scene->isShadowed(state, bRay);
                    if(!shadowed && lightPdf > 1e-6f)
                    {
                        if(trShad) lcol *= scol;
                        color_t transmitCol = scene->volIntegrator->transmittance(state, lightRay);
                        lcol *= transmitCol;
                        float lPdf = 1.f/lightPdf;
                        float l2 = lPdf * lPdf;
                        float m2 = s.pdf * s.pdf;
                        float w = m2 / (l2 + m2);
                        ccol2 += surfCol * lcol * w * W;
                    }
                }
            }
            col += ccol2 * invNS;
        }
    }
    */

    return col;
}


/*
  classical_sampling()

  INPUT: number_of_samples, triangle_areas

  for each t in triangles:
      triangles_areas[t.index] = t.area

  cdf_triangles
  cdf_triangles[0] = 0

  for each t in triangles, t > 0:
      cdf_triangles[t.index] += cdf_triangles[t.index - 1]

  area_sum = cdf_triangles.last

  for each s in number_of_samples:
      ksi = hammersley(3d)
      triangle_index = lower_bound(cdf_triangles.begin, cdf_triangles.end, ksi[0] * area_sum)

      sample_pos = sample_triangle(triangle_index, ksi[1], ksi[2])
  */





void pbLighting_t::generate_spherical_function(renderState_t & state, GiPoint * gi_point, surfacePoint_t const& sp, std::vector<light_t*> const& lights)
{
    ray_t lightRay;
    lightRay.from = gi_point->pos;

    color_t lcol(0.f);

    for (std::vector<light_t *>::const_iterator light = lights.begin(); light != lights.end(); ++light)
    {
        if ((*light)->diracLight())
        {
            bool shadowed = true;

            if ((*light)->illuminate(sp, lcol, lightRay))
            {
                // ...shadowed...
                lightRay.tmin = YAF_SHADOW_BIAS; // < better add some _smart_ self-bias value...this is bad.
                shadowed = scene->isShadowed(state, lightRay);
                float const cos_sp_N_lightray_dir = sp.N * lightRay.dir;
                if (!shadowed && cos_sp_N_lightray_dir > 0.0f)
                {
                    color_t const incoming_light = lcol;
                    // color_t const reflected_light = lcol * std::max(0.0f, cos_sp_N_lightray_dir);
                    gi_point->color += incoming_light;
                }
                else
                {
                    lcol = color_t(0.0f);
                }
            }

            if (!shadowed)
            {
                Abstract_spherical_function_estimator<color_t> * color_estimator = new Spherical_function_light_color_estimator<color_t>(state, sp, lightRay, lcol);
                gi_point->sf_representation.color->calc_coefficients_random(color_estimator);
            }
        }
    }

    if (sp.material->getFlags() & BSDF_EMIT)
    {
        Abstract_spherical_function_estimator<color_t> * emit_estimator = new Spherical_function_emit_color_estimator<color_t>(state, sp);
        gi_point->sf_representation.color->calc_coefficients_random(emit_estimator);
    }

    Abstract_spherical_function_estimator<float> * area_estimator = new Spherical_function_area_estimator<float>(sp.N, gi_point->area);
    gi_point->sf_representation.area->calc_coefficients_random(area_estimator);
}


bound_t generate_disc_bounding_box(vector3d_t const& p, vector3d_t const& nu, vector3d_t const& nv, float const radius)
{
    bound_t bounding_box;

    vector3d_t disc_bound[4];
    disc_bound[0] = p + nu * radius + nv * radius;
    disc_bound[1] = p + nu * radius - nv * radius;
    disc_bound[2] = p - nu * radius - nv * radius;
    disc_bound[3] = p - nu * radius + nv * radius;

    bounding_box = bound_t(disc_bound[0], disc_bound[1]);
    bounding_box.include(disc_bound[2]);
    bounding_box.include(disc_bound[3]);

    return bounding_box;
}




void pbLighting_t::generate_gi_points_data_2(std::vector<MyTree*> const& nodes,
                                             std::tr1::unordered_map<GiPoint*, pbgi_sample_t*> const& point_to_sampling_point_map,
                                             float const area,
                                             float const radius,
                                             int const treat_depth_as_leaf)
{
    std::cout << "pbLighting_t::generate_gi_points_data(), nodes: " << nodes.size() << std::endl;

    Gi_point_averager averager(_spherical_function_color_factory, _spherical_function_area_factory);

    renderState_t state;
    unsigned char userdata[USER_DATA_SIZE+7];
    state.userdata = (void *)( &userdata[7] - ( ((size_t)&userdata[7])&7 ) ); // pad userdata to 8 bytes
    state.cam = scene->getCamera();

    for (unsigned int i = 0; i < nodes.size(); ++i)
    {
        MyTree* node = nodes[i];

        if ((i % 1000) == 0)
        {
            std::cout << "." << std::flush;
        }

        if (node->getIsLeaf())
        {
//            std::cout << "pbLighting_t::generate_gi_points_data: leaf" << std::endl;

            std::vector<GiPoint*> node_data = node->getData();
            assert(node_data.size() == 1);

            for (unsigned int j = 0; j < node_data.size(); ++j)
            {
                GiPoint * gi_point = node_data[j];

                std::tr1::unordered_map<GiPoint*, pbgi_sample_t*>::const_iterator iter = point_to_sampling_point_map.find(gi_point);
                assert(iter != point_to_sampling_point_map.end());

                pbgi_sample_t & sampling_point = *(iter->second);

                point3d_t const& hitPoint = sampling_point.position;
                triangle_t const& tri = *sampling_point.tri_pointer;
                intersectData_t & iData = sampling_point.intersect_data;

                surfacePoint_t sp;

                tri.getSurface(sp, hitPoint, iData);

                assert(!std::isnan(sp.P.x));

                gi_point->pos          = vector3d_t(sp.P);
                gi_point->normal       = sp.N;
                //giPoint->area = area_per_sample;
                gi_point->area         = area;
                // gi_point->area         = area;
                //giPoint->debug_radius = needed_radius;
                gi_point->debug_radius = 0.0f;

                gi_point->bounding_box = generate_disc_bounding_box(vector3d_t(sp.P), sp.NU, sp.NV, radius);

                unsigned char userdata[USER_DATA_SIZE];
                state.userdata = (void *) userdata;

                material_t const* material = sp.material;
                BSDF_t bsdfs;
                material->initBSDF(state, sp, bsdfs);

                gi_point->sf_representation.color = _spherical_function_color_factory->create();
                gi_point->sf_representation.area  = _spherical_function_area_factory->create();
                generate_spherical_function(state, gi_point, sp, lights);

                gi_point->color *= material->getDiffuseAtPoint(state, sp);
                gi_point->color += material->emission(state, sp, vector3d_t());

                node->average_node(averager);
            }
        }
        else
        {
//            std::cout << "pbLighting_t::generate_gi_points_data: node" << std::endl;

            bool const treat_as_leaf = (node->getDepth() == treat_depth_as_leaf);

            if (!treat_as_leaf)
            {
                node->average_node(averager);
            }

            std::vector<MyTree> & children = node->getChildren();

            for (unsigned int j = 0; j < children.size(); ++j)
            {
                // convert the children, not the node itself (yet)
                if (_area_converter && !treat_as_leaf)
                {
                    Spherical_function<float> * tmp;

                    if (children[j].getIsLeaf())
                    {
                        assert(children[j].getData().size() == 1);
                        //tmp = children[j].getData()[0]->sf_representation.area;
                        //children[j].getData()[0]->sf_representation.area = _area_converter->convert(children[j].getData()[0]->sf_representation.area);
                        //delete tmp;
                        delete children[j].getData()[0]->sf_representation.area;
                        children[j].getData()[0]->sf_representation.area = NULL;
                    }

                    tmp = children[j].getClusteredData()->sf_representation.area;
                    children[j].getClusteredData()->sf_representation.area = _area_converter->convert(children[j].getClusteredData()->sf_representation.area);
                    delete tmp;
                }
                else
                {
                    if (children[j].getIsLeaf())
                    {
                        assert(children[j].getData().size() == 1);
                        delete children[j].getData()[0]->sf_representation.area;
                        children[j].getData()[0]->sf_representation.area = NULL;
                    }
                }

                if (_color_converter && !treat_as_leaf)
                {
                    //                std::cout << "pbLighting_t::generate_gi_points_data: convert children" << std::endl;

                    Spherical_function<color_t> * tmp;

                    if (children[j].getIsLeaf())
                    {
                        //                        std::cout << "pbLighting_t::generate_gi_points_data: convert leaf" << std::endl;

                        //                        tmp = children[j].getData()[0]->sf_representation.color;
                        //                        children[j].getData()[0]->sf_representation.color = _color_converter->convert(children[j].getData()[0]->sf_representation.color);
                        //                        delete tmp;

                        assert(children[j].getData().size() == 1);
                        delete children[j].getData()[0]->sf_representation.color;
                        children[j].getData()[0]->sf_representation.color = NULL;
                    }

                    tmp = children[j].getClusteredData()->sf_representation.color;
                    children[j].getClusteredData()->sf_representation.color = _color_converter->convert(children[j].getClusteredData()->sf_representation.color);
                    delete tmp;
                }
                else
                {
                    if (children[j].getIsLeaf())
                    {
                        assert(children[j].getData().size() == 1);
                        delete children[j].getData()[0]->sf_representation.color;
                        children[j].getData()[0]->sf_representation.color = NULL;
                    }
                }
            }
        }
    }

    std::cout << "generate_gi_points_data: used cells/all cells: " << Cube_spherical_function<color_t>::_num_used_cells << "/" << Cube_spherical_function<color_t>::_num_all_cells << std::endl;
}




void pbLighting_t::generate_gi_points_data(renderState_t & state,
                                           std::vector<MyTree*> const& nodes,
                                           std::tr1::unordered_map<GiPoint*, surfacePoint_t> const& point_to_surface_point_map,
                                           int const treat_depth_as_leaf)
{
    std::cout << "pbLighting_t::generate_gi_points_data(), nodes: " << nodes.size() << std::endl;

    Gi_point_averager averager(_spherical_function_color_factory, _spherical_function_area_factory);

    for (unsigned int i = 0; i < nodes.size(); ++i)
    {
        MyTree* node = nodes[i];

        if ((i % 1000) == 0)
        {
            std::cout << "." << std::flush;
        }

        if (node->getIsLeaf())
        {
//            std::cout << "pbLighting_t::generate_gi_points_data: leaf" << std::endl;

            std::vector<GiPoint*> node_data = node->getData();
            assert(node_data.size() == 1);

            for (unsigned int j = 0; j < node_data.size(); ++j)
            {
                GiPoint * gi_point = node_data[j];

                std::tr1::unordered_map<GiPoint*, surfacePoint_t>::const_iterator iter = point_to_surface_point_map.find(gi_point);
                assert(iter != point_to_surface_point_map.end());

                surfacePoint_t const& sp = iter->second;

                unsigned char userdata[USER_DATA_SIZE];
                state.userdata = (void *) userdata;

                material_t const* material = sp.material;
                BSDF_t bsdfs;
                material->initBSDF(state, sp, bsdfs);

                gi_point->sf_representation.color = _spherical_function_color_factory->create();
                gi_point->sf_representation.area  = _spherical_function_area_factory->create();
                generate_spherical_function(state, gi_point, sp, lights);

                gi_point->color *= material->getDiffuseAtPoint(state, sp);
                gi_point->color += material->emission(state, sp, vector3d_t());

                node->average_node(averager);
            }
        }
        else
        {
//            std::cout << "pbLighting_t::generate_gi_points_data: node" << std::endl;

            bool const treat_as_leaf = (node->getDepth() == treat_depth_as_leaf);

            if (!treat_as_leaf)
            {
                node->average_node(averager);
            }

            // convert the children, not the node itself (yet)
            if (_area_converter && !treat_as_leaf)
            {
                std::vector<MyTree> & children = node->getChildren();

                for (unsigned int j = 0; j < children.size(); ++j)
                {
                    Spherical_function<float> * tmp;

                    if (children[j].getIsLeaf())
                    {
                        assert(children[j].getData().size() == 1);
                        tmp = children[j].getData()[0]->sf_representation.area;
                        children[j].getData()[0]->sf_representation.area = _area_converter->convert(children[j].getData()[0]->sf_representation.area);
                        delete tmp;
                    }

                    tmp = children[j].getClusteredData()->sf_representation.area;
                    children[j].getClusteredData()->sf_representation.area = _area_converter->convert(children[j].getClusteredData()->sf_representation.area);
                    delete tmp;
                }
            }

            if (_color_converter && !treat_as_leaf)
            {
//                std::cout << "pbLighting_t::generate_gi_points_data: convert children" << std::endl;

                std::vector<MyTree> & children = node->getChildren();

                for (unsigned int j = 0; j < children.size(); ++j)
                {
                    Spherical_function<color_t> * tmp;

                    if (children[j].getIsLeaf())
                    {
//                        std::cout << "pbLighting_t::generate_gi_points_data: convert leaf" << std::endl;

                        assert(children[j].getData().size() == 1);
                        tmp = children[j].getData()[0]->sf_representation.color;
                        children[j].getData()[0]->sf_representation.color = _color_converter->convert(children[j].getData()[0]->sf_representation.color);
                        delete tmp;
                    }

                    tmp = children[j].getClusteredData()->sf_representation.color;
                    children[j].getClusteredData()->sf_representation.color = _color_converter->convert(children[j].getClusteredData()->sf_representation.color);
                    delete tmp;
                }
            }
        }
    }

    std::cout << "generate_gi_points_data: used cells/all cells: " << Cube_spherical_function<color_t>::_num_used_cells << "/" << Cube_spherical_function<color_t>::_num_all_cells << std::endl;
}

void test_dart_params(float const radius, int const number_of_samples, std::vector<triangle_t const*> const& triangles)
{
    for (float cell_size_factor = 1.0f; cell_size_factor < 2.01f; cell_size_factor += 0.5f)
    {
        for (float bin_count_factor = 0.5f; bin_count_factor < 2.01f; bin_count_factor += 0.5f)
        {
            timer_t my_timer;
            my_timer.addEvent("t1");

            std::cout << "test: cell_size_factor: " << cell_size_factor << " bin_count_factor: " << bin_count_factor << std::endl;

            Scene_sampler_darts_hash sampler(radius, number_of_samples, cell_size_factor, bin_count_factor);

            my_timer.start("t1");
            std::vector<pbgi_sample_t> sampling_points = sampler.generate_samples(triangles);
            my_timer.stop("t1");

            std::cout << "samples/desired_samples: " << (sampling_points.size() / float(number_of_samples)) << " time: " << my_timer.getTime("t1") << std::endl;
        }
    }
}

void pbLighting_t::generate_gi_points(renderState_t & state, MyTree * tree, int const number_of_samples)
{
    std::vector<triangle_t const*> triangles = get_scene_triangles(scene->meshes);

    float const scene_area = get_total_scene_area(triangles);

    float const area_per_sample = scene_area / float(number_of_samples);

    std::cout << "generate_gi_points(): scene_area: " << scene_area << " area_per_sample: " << area_per_sample << std::endl;

    float const radius = std::sqrt(area_per_sample / M_PI) / std::sqrt(2.0f);
    std::cout << "generate_gi_points(): best radius: " << radius << std::endl;

//    test_dart_params(radius, number_of_samples, get_scene_triangles(scene->meshes));

    timer_t timer;
    timer.addEvent("t1");
    timer.start("t1");

    bool use_cdf_sampler = true;
    bool do_find_widest_gap = false;

    Scene_sampler * sampler = NULL;

    if (use_cdf_sampler)
    {
        sampler = new Scene_sampler_cdf(number_of_samples);
    }
    else
    {
        sampler = new Scene_sampler_darts_hash(radius, number_of_samples);
    }

    std::vector<pbgi_sample_t> sampling_points = sampler->generate_samples(triangles);


    // Scene_sampler_darts_hash sampler(radius, number_of_samples, 1.0f, 2.0f);


    // Scene_sampler_cdf sampler(number_of_samples);

    // float const sampling_solid_angle = 1.0f / float(number_of_samples);
    // Scene_sampler_reyes sampler(sampling_solid_angle, vector3d_t(scene->getCamera()->getPosition()), scene->meshes, &_debug_new_triangles);

    // std::vector<pbgi_sample_t> sampling_points = generate_samples_cdf(number_of_samples, get_scene_triangles(scene->meshes));
    // std::vector<pbgi_sample_t> sampling_points = generate_samples_darts(0.05f, number_of_samples, get_scene_triangles(scene->meshes));
    // std::vector<pbgi_sample_t> sampling_points = generate_samples_darts_hash(radius * 0.5f, number_of_samples, get_scene_triangles(scene->meshes));
    // std::vector<pbgi_sample_t> sampling_points = generate_samples_suicide(number_of_samples, radius, get_scene_triangles(scene->meshes));



    timer.stop("t1");

    std::cout << "generate_gi_points(): samples: " << sampling_points.size() << std::endl;

    std::cout << "samples/desired_samples: " << (sampling_points.size() / float(number_of_samples)) << " time: " << timer.getTime("t1") << std::endl;

    // float const needed_radius = radius;


    float needed_radius;

    if (do_find_widest_gap)
    {
        float const largest_distance = generate_histogram(sampling_points, radius);
        // float const largest_distance = sampler.get_widest_gap();
        needed_radius = largest_distance * std::sqrt(2.0f) / 2.0f; // correct factor
        // float const needed_radius = largest_distance * std::sqrt(2.0f) / 4.0f; // FIXME: correct factor: std::sqrt(2.0f) / 2.0f but leads to too fast fill degree
    }
    else
    {
        float const real_covered_area = sampling_points.size() * radius * radius * M_PI;
        float const missing_area_factor = real_covered_area / scene_area; // <= 1
        float const needed_area_per_sample = area_per_sample / missing_area_factor;
        needed_radius = std::sqrt(needed_area_per_sample / M_PI);

        if (use_cdf_sampler)
        {
            needed_radius *= 2.0f; // additional factor of 2.0 when using the simple CDF sampling
        }
    }

    float const area_estimation = M_PI * needed_radius * needed_radius;

    std::cout << "generate_gi_points(): needed radius: " << needed_radius << std::endl;

    // std::vector<GiPoint*> points(sampling_points.size());

    //std::tr1::unordered_map<GiPoint*, surfacePoint_t> point_to_surface_point_map;
    std::tr1::unordered_map<GiPoint*, pbgi_sample_t*> gi_point_to_sampling_point_map;

    std::cout << "generate_gi_points(): creating points" << std::endl;

    tree->prepare_containers(sampling_points.size());

    timer.addEvent("point_creation");
    timer.start("point_creation");

    // -------------------------------------------
    // Put the GI points into the tree and create a correspondency map between
    // point and its surface point
    // #pragma omp parallel for
    for (int i = 0; i < int(sampling_points.size()); ++i)
    {
        GiPoint * gi_point = new GiPoint(NULL, NULL);
        gi_point_to_sampling_point_map[gi_point] = &sampling_points[i];
        tree->add_point_prepared(i, vector3d_t(sampling_points[i].position), gi_point);


        /*
        point3d_t const& hitPoint = sampling_points[i].position;
        triangle_t const& tri = *sampling_points[i].tri_pointer;
        // float const area = sampling_points[i].area * 2.0f;
        intersectData_t & iData = sampling_points[i].intersect_data;

        surfacePoint_t sp;

        tri.getSurface(sp, hitPoint, iData);

        assert(!std::isnan(sp.P.x));

        // GiPoint * gi_point = new GiPoint(_spherical_function_color_factory, _spherical_function_area_factory);
        GiPoint * gi_point = new GiPoint(NULL, NULL);

        gi_point->pos          = vector3d_t(sp.P);
        gi_point->normal       = sp.N;
        //giPoint->area = area_per_sample;
        gi_point->area         = area_estimation;
        // gi_point->area         = area;
        //giPoint->debug_radius = needed_radius;
        gi_point->debug_radius = 0.0f;

        gi_point->bounding_box = generate_disc_bounding_box(vector3d_t(sp.P), sp.NU, sp.NV, needed_radius);

        point_to_surface_point_map[gi_point] = sp;

        tree->add_point_prepared(i, vector3d_t(sp.P), gi_point);
*/
    }

    std::cout << "generate_gi_points(): finish tree" << std::endl;
    // tree->finalize_points();
    MyTree::finalize_points_iterative(tree);
    timer.stop("point_creation");
    // ------------------------------------------


    // ------------------------------------------
    // go through all nodes in post order, generate the SF representation in the children,
    // then go up to the parent, use the children to generate its SF and convert it and delete the children's unconverted etc.

    timer.addEvent("node_data_creation");
    timer.start("node_data_creation");

    bool use_parallel_averaging = true;

    if (use_parallel_averaging)
    {
        int const num_procs = omp_get_num_procs();
        int const splitting_depth = std::ceil(std::log(num_procs) / std::log(2.0f));
        int const num_threads = std::pow(2, splitting_depth);

        std::cout << "generate_gi_points(): parallel averaging, num_threads: " << num_threads << " depth: " << splitting_depth << std::endl;

        std::vector<MyTree*> nodes_in_splitting_depth = tree->getNodes(splitting_depth);

        assert(int(nodes_in_splitting_depth.size()) == num_threads);


        #pragma omp parallel for
        for (int i = 0; i < num_threads; ++i)
        {
            renderState_t parallel_state;
            unsigned char userdata[USER_DATA_SIZE+7];
            parallel_state.userdata = (void *)( &userdata[7] - ( ((size_t)&userdata[7])&7 ) ); // pad userdata to 8 bytes
            parallel_state.cam = scene->getCamera();

            MyTree * split_root = nodes_in_splitting_depth[i];

            std::vector<MyTree*> split_nodes;
            split_root->get_post_order_queue(split_nodes);

            //generate_gi_points_data(parallel_state, split_nodes, point_to_surface_point_map);
            generate_gi_points_data_2(split_nodes, gi_point_to_sampling_point_map, area_estimation, needed_radius);
        }

        std::cout << "generate_gi_points(): finished lower tree part" << std::endl;

        std::vector<MyTree*> nodes;
        tree->get_post_order_queue(nodes, splitting_depth);

        //generate_gi_points_data(state, nodes, point_to_surface_point_map, splitting_depth);
        generate_gi_points_data_2(nodes, gi_point_to_sampling_point_map, area_estimation, needed_radius, splitting_depth);
    }
    else
    {
        std::vector<MyTree*> nodes;
        tree->get_post_order_queue(nodes);

        //generate_gi_points_data(state, nodes, point_to_surface_point_map);
        generate_gi_points_data_2(nodes, gi_point_to_sampling_point_map, area_estimation, needed_radius);
    }

    timer.stop("node_data_creation");

    // ------------------------------------------

    /*
    std::string fileName = "/tmp/pbgi_points_store";
    std::ofstream out_file(fileName.c_str(), std::ios_base::binary);
    //boost::archive::binary_oarchive ar(out_file);
    boost::archive::text_oarchive ar(out_file);
    ar.register_type<GiSphericalHarmonics<vector3d_t, color_t> >();
    ar.register_type<Cube_spherical_function<vector3d_t, color_t> >();

    ar << points;

    out_file.close();
    */

    std::cout << "surfel count: " << sampling_points.size() <<
                 " timing, points: " << timer.getTime("point_creation") <<
                 " node data: "      << timer.getTime("node_data_creation") <<
                 std::endl;

    _log_file << "Num surfels: " << sampling_points.size() << std::endl;
    _log_file << "Sample point generation time: " << timer.getTime("point_creation") << std::endl;
    _log_file << "Node data generation time: " << timer.getTime("node_data_creation") << std::endl;

    delete sampler;
}

// FIXME: broken due to fucked up serialisation
yafaray::pbLighting_t::MyTree * load_gi_points()
{
    std::cout << "Loading points" << std::endl;

    std::vector<GiPoint*> points;

    std::string fileName = "/tmp/pbgi_points_store";
    std::ifstream in_file(fileName.c_str(), std::ios_base::binary);

    /*
    boost::archive::text_iarchive ar(in_file);
    ar.register_type<GiSphericalHarmonics>();
    ar.register_type<Cube_spherical_function>();
    */

    // ar >> points;

    in_file.close();

    std::cout << "load_gi_points(): surfel count: " << points.size() << std::endl;

    yafaray::point3d_t min(1e10f, 1e10f, 1e10f), max(-1e10f, -1e10f, -1e10f);

    for (unsigned int i = 0; i < points.size(); ++i)
    {
        GiPoint * p = points[i];

        min.x = std::min(min.x, p->pos.x);
        min.y = std::min(min.y, p->pos.y);
        min.z = std::min(min.z, p->pos.z);

        max.x = std::max(max.x, p->pos.x);
        max.y = std::max(max.y, p->pos.y);
        max.z = std::max(max.z, p->pos.z);
    }

    pbLighting_t::MyTree* tree = new pbLighting_t::MyTree(vector3d_t(min), vector3d_t(max), 30, 1);

    for (unsigned int i = 0; i < points.size(); ++i)
    {
        tree->addPoint(points[i]->pos, points[i]);
    }

    return tree;
}



std::vector< Spherical_function<color_t> *> generate_dictionary_from_tree(pbLighting_t::MyTree const* tree,
                                                                          Dictionary_generator const* dictionary_generator,
                                                                          Spherical_function_factory<color_t> const* spherical_function_factory,
                                                                          int const dict_num_centers)
{
    std::vector<pbLighting_t::MyTree const*> nodes = tree->getNodes();

    std::vector<Word> word_list;

    for (unsigned int i = 0; i < nodes.size(); ++i)
    {
        GiPoint const* gi_point = nodes[i]->getClusteredData();

        assert(gi_point);

        if (gi_point)
        {
            Word word = gi_point->sf_representation.color->to_vector();

            word_list.push_back(word);
        }
    }

    std::vector<Word> unpacked_dictionary = dictionary_generator->generate(word_list, dict_num_centers);

    std::vector<Spherical_function<color_t> *> dictionary;

    for (unsigned int i = 0; i < unpacked_dictionary.size(); ++i)
    {
        // Spherical_function<Data> * sf = spherical_function_factory->create();
        Spherical_function<color_t> * sf = spherical_function_factory->create();
        sf->from_vector(unpacked_dictionary[i]);

        dictionary.push_back(sf);
    }

    return dictionary;
}



void push_light(renderState_t & state, pbLighting_t::MyTree * node, Splat_cube_raster_buffer const& fb_in)
{
    // add light to the node
    node->getClusteredData()->sf_representation.color->add_light(fb_in);

    if (!node->has_children()) return;

    int const local_raster_buffer_resolution = fb_in.get_resolution();

    Cube_raster_buffer<float> fb_area_0;
    fb_area_0.setup_surfel_buffer(local_raster_buffer_resolution);

    Cube_raster_buffer<float> fb_area_1;
    fb_area_1.setup_surfel_buffer(local_raster_buffer_resolution);

    Splat_cube_raster_buffer fb_in_0;
    fb_in_0.setup_surfel_buffer(local_raster_buffer_resolution);

    Splat_cube_raster_buffer fb_in_1;
    fb_in_1.setup_surfel_buffer(local_raster_buffer_resolution);

    pbLighting_t::MyTree* child_0 = &(node->getChildren()[0]);
    pbLighting_t::MyTree* child_1 = &(node->getChildren()[1]);

    GiPoint* child_0_data = child_0->getClusteredData();
    GiPoint* child_1_data = child_1->getClusteredData();

    std::vector<Cube_cell> const& cells = fb_area_0.get_cube_cells();

    for (std::size_t i = 0; i < cells.size(); ++i)
    {
        float const area_0 = child_0_data->sf_representation.area->get_value(fb_area_0.get_cell_direction(cells[i]));
        float const area_1 = child_1_data->sf_representation.area->get_value(fb_area_1.get_cell_direction(cells[i]));

        color_t const color_in_0 = fb_in.get_data(cells[i]) * area_0 / (area_0 + area_1);
        color_t const color_in_1 = fb_in.get_data(cells[i]) * area_1 / (area_0 + area_1);

        fb_in_0.set_data(cells[i], color_in_0);
        fb_in_1.set_data(cells[i], color_in_1);
    }

    // push down
    push_light(state, child_0, fb_in_0);
    push_light(state, child_1, fb_in_1);
}


void pbLighting_t::add_bounce(renderState_t & state, pbLighting_t::MyTree * tree)
{
//    std::cout << "pbLighting_t::add_bounce();" << std::endl;

//    // get a certain level of nodes
//    // for each node: gather indirect light incoming to this node's position
//    // add to node and push the (new) light down the tree

//    int const node_level = 12;
//    int const local_raster_buffer_resolution = 16;
//    float const local_solid_angle_factor = 16.0f;

//    std::vector<pbLighting_t::MyTree*> nodes = tree->getNodes(node_level);

//    for (std::size_t i = 0; i < nodes.size(); ++i)
//    {
//        if (i % 100 == 0)
//        {
//            std::cout << "." << std::flush;
//        }

//        pbLighting_t::MyTree* node = nodes[i];

//        Splat_cube_raster_buffer fb_in =
//                gather_light_at_point(tree,
//                                      node->getClusteredData()->pos,
//                                      local_solid_angle_factor,
//                                      _disc_scale_factor,
//                                      local_raster_buffer_resolution,
//                                      raster_buffer_type,
//                                      node_splat_type,
//                                      surfel_far_splat_type,
//                                      surfel_near_splat_type,
//                                      surfel_near_threshold);

//        push_light(state, node, fb_in);

//    }
}



std::string extract_initials(std::string const& name)
{
    std::string res;
    res += name[0];

    size_t pos = 0;

    while (true)
    {
        pos = name.find("_", pos + 1);

        if (pos == std::string::npos)
        {
            break;
        }

        res += name[pos + 1];
    }

    return res;
}

template <class T>
std::string get_string_from_type_map(T const& type, std::map<std::string, T> const& type_map)
{
    typename std::map<std::string, T>::const_iterator iter;

    for (iter = type_map.begin(); iter != type_map.end(); ++iter)
    {
        if (iter->second == type)
        {
            return iter->first;
        }
    }

    return "X";
}


template <class Data>
class Spherical_function_print_visitor : public Spherical_function_visitor<Data>
{
public:
    virtual void visit(GiSphericalHarmonics<Data> * sf) {}
    virtual void visit(Cube_spherical_function<Data> * sf) {}
    virtual void visit(Mises_Fisher_spherical_function<Data> * sf)
    {
        if (sf->get_lobes().size() == 0)
        {
            std::cout << "no lobes" << std::endl;
            return;
        }

        std::cout << sf->get_lobes()[0] << std::endl;
    }

    virtual void identify() { std::cout << "Spherical_function_print_visitor" << std::endl; }
};


float dictionary_mean_absolute_error(pbLighting_t::MyTree const* tree)
{
    std::vector<pbLighting_t::MyTree const*> nodes = tree->getNodes();

    float error = 0.0f;

    for (std::size_t i = 0; i < nodes.size(); ++i)
    {
        GiPoint * gi_point = nodes[i]->getClusteredData();

        Indexed_spherical_function<color_t> * indexed_sf = dynamic_cast<Indexed_spherical_function<color_t> *>(gi_point->sf_representation.color);

        if (!indexed_sf) continue;

        std::vector<float> original_values = indexed_sf->get_original_sf()->to_vector();
        std::vector<float> indexed_values  = indexed_sf->get_indexed_sf()->to_vector();

        Distance_function distance_fnc;

        float const error_i = distance_fnc(original_values, indexed_values);

        error += error_i;
    }

    error /= float(nodes.size());

    return error;
}




bool pbLighting_t::preprocess()
{
    Y_INFO << "PBGI Preprocess" << std::endl;

    _log_file.open("/tmp/yafaray_pbgi_debug_log.txt");

    bool success = true;
    settings = "";

    std::stringstream ss;
    ss << "#s: "    << surfel_samples <<
          " | rr: " << raster_buffer_resolution;

    if (variational) ss << " (V)";

    ss << ", "      << extract_initials(get_string_from_type_map(raster_buffer_type, Splat_cube_raster_buffer::enum_type_map)) <<
          ", "      << extract_initials(get_string_from_type_map(node_splat_type, Splat_cube_raster_buffer::enum_splat_type_map)) <<
          ", "      << extract_initials(get_string_from_type_map(surfel_far_splat_type, Splat_cube_raster_buffer::enum_splat_type_map)) <<
          ", "      << extract_initials(get_string_from_type_map(surfel_near_splat_type, Splat_cube_raster_buffer::enum_splat_type_map)) <<
          " | ds: " << _disc_scale_factor <<
          " | sa: " << solid_angle_factor <<
          " | sfc: " << _spherical_function_color_factory->get_name();

    if (_color_converter) ss << ", " << _color_converter->get_name();
    if (_color_dictionary_generator) ss << ", " << extract_initials(_color_dictionary_generator->identify());

    ss << ", "      << _sf_resolution <<
          ", sfa: " << _spherical_function_area_factory->get_name();
          //" | dict: " << _dictionary_type <<
          //" | samp frac: " << _dictionary_sample_fraction;

    settings = ss.str();

    _log_file << "Receiving point buffer: " <<
                 "resolution: " << raster_buffer_resolution << ", " <<
                 "type: " << get_string_from_type_map(raster_buffer_type, Splat_cube_raster_buffer::enum_type_map) << ", " <<
                 "node: " << get_string_from_type_map(node_splat_type, Splat_cube_raster_buffer::enum_splat_type_map) << ", " <<
                 "surfel far: " << get_string_from_type_map(surfel_far_splat_type, Splat_cube_raster_buffer::enum_splat_type_map) << ", " <<
                 "surfel near: " << get_string_from_type_map(surfel_near_splat_type, Splat_cube_raster_buffer::enum_splat_type_map) << ", " <<
                 "surfel threshold: " << surfel_near_threshold << std::endl;

    _log_file << "Disc scale factor: " << _disc_scale_factor << std::endl;
    _log_file << "Solid angle factor: " << solid_angle_factor << std::endl;
    _log_file << "SF color: " << _spherical_function_color_factory->get_name() << ", " << _spherical_function_color_factory->get_settings() << std::endl;
    _log_file << "SF area: " << _spherical_function_area_factory->get_name() << ", " << _spherical_function_area_factory->get_settings() << std::endl;

    if (_color_converter) _log_file << "Color converter: " << _color_converter->get_name() << std::endl;
    if (_color_dictionary_generator) _log_file << "Color Dict. converter: " << _color_dictionary_generator->identify() << ", " <<
                                                  "num centers: " << _dict_num_centers << ", " <<
                                                  "sample fraction: " << _dictionary_sample_fraction << std::endl;


    background = scene->getBackground();
    lights = scene->lights;

    renderState_t state;
    unsigned char userdata[USER_DATA_SIZE+7];
    state.userdata = (void *)( &userdata[7] - ( ((size_t)&userdata[7])&7 ) ); // pad userdata to 8 bytes
    state.cam = scene->getCamera();

    bound_t const& sceneBound = scene->getSceneBound();


    if (_color_dictionary_generator)
    {
        MyTree * dict_tree = new MyTree(vector3d_t(sceneBound.a), vector3d_t(sceneBound.g), 100, 1);
        generate_gi_points(state, dict_tree, surfel_samples * _dictionary_sample_fraction);


        timer_t my_timer;
        my_timer.addEvent("Dictionary_creation");
        my_timer.start("Dictionary_creation");

        _color_dictionary = generate_dictionary_from_tree(dict_tree, _color_dictionary_generator, _spherical_function_color_factory, _dict_num_centers);

        my_timer.stop("Dictionary_creation");

        _log_file << "Dictionary creation time: " << my_timer.getTime("Dictionary_creation");

        _color_converter = new Spherical_function_indexed_converter<color_t>(&_color_dictionary, _do_dictionary_stats);

        _bspTree = new MyTree(vector3d_t(sceneBound.a), vector3d_t(sceneBound.g), 100, 1);
        generate_gi_points(state, _bspTree, surfel_samples);

        if (_do_dictionary_stats)
        {
            float const ma_error = dictionary_mean_absolute_error(_bspTree);

            _log_file << "Dictionary mean absolute error: " << ma_error << std::endl;

            _log_file << "Dictionary stats: " << _color_dictionary_generator->get_stat_results() << std::endl;

            std::cout << "_dict_num_centers: " << _dict_num_centers << " " <<
                         "_dictionary_sample_fraction: " << _dictionary_sample_fraction << " " <<
                         "ma_error: " << ma_error << std::endl;
        }

        delete dict_tree;
    }
    else
    {
        if (do_load_gi_points)
        {
            _bspTree = load_gi_points();
        }
        else
        {
            _bspTree = new MyTree(vector3d_t(sceneBound.a), vector3d_t(sceneBound.g), 100, 1);
            generate_gi_points(state, _bspTree, surfel_samples);

            if (!_use_sf_files)
            {
                // framebuffer_nodes_to_text(_bspTree, _sf_resolution);
            }
        }
    }

    // add_bounce(state, _bspTree);


//    { // debug
//        Spherical_function_visitor<color_t> * visitor = new Spherical_function_print_visitor<color_t>;

//        std::vector<MyTree *> nodes = _bspTree->getNodes();
//        for (std::size_t i = 0; i < nodes.size(); ++i)
//        {
//            MyTree const* n = nodes[i];

//            if (n->has_children())
//            {
//                std::cout << "i: "  << i << std::endl;

//                std::cout << "p: ";
//                n->getClusteredData()->sf_representation.color->accept(visitor);
//                std::cout << "c0: ";
//                n->getChildren()[0].getClusteredData()->sf_representation.color->accept(visitor);
//                std::cout << "c1: ";
//                n->getChildren()[1].getClusteredData()->sf_representation.color->accept(visitor);
//            }
//        }
//    }


    // Y_INFO << "PBGI: BSP: " << *_bspTree << std::endl;

    Y_INFO << "PBGI: solid angle: "
           << solid_angle_factor << " "
           << debugTreeDepth << " "
           << "max depth: " << _bspTree->getMaxDepth() << " "
           << "min leaf depth: " << _bspTree->getMinLeafDepth() << " "
           << settings << " "
           << std::endl;

    _log_file.close();

    return success;
}



void process_surfel(
        GiPoint const& gi_point,
        surfacePoint_t const& rp,
        Cube_raster_buffer<color_t> & frame_buffer,
        float const surfel_near_threshold,
        float const mix_amount,
        float const disc_scale_factor,
        bool const handle_rp_as_surface,
        std::vector<Gi_point_info> & point_infos,
        std::vector<GiPoint*> & debug_points,
        Debug_info * debug_info)
{
    yafaray::vector3d_t surfel_to_rp = (vector3d_t(rp.P) - gi_point.pos);

    float const distance = surfel_to_rp.length();

    surfel_to_rp.normalize();

    if (handle_rp_as_surface)
    {
        float const cos_rp_n_gip = rp.N * (-surfel_to_rp);
        if (cos_rp_n_gip < 0.05f) return;
    }

    // float const cos_sp_gip = std::max(gi_point.normal * giToSp, 0.0f);
    float cos_sp_gip = gi_point.normal * surfel_to_rp;
    bool const back_facing = cos_sp_gip < 0.0f;

    if (back_facing) return;

    cos_sp_gip = std::abs(cos_sp_gip);


    // float const cos_sp_gip = std::abs(giP.normal * giToSp);

    float const visible_area = cos_sp_gip * gi_point.area;
    float const max_area = gi_point.area;

    float const solid_angle_real = visible_area / std::max(0.001f, distance * distance);

    float const disc_radius = std::sqrt(max_area / M_PI);
    float const visible_radius = std::sqrt(visible_area / M_PI);

    yafaray::color_t contribution;

    if (debug_info && debug_info->color_by_depth)
    {
        contribution = pbLighting_t::debug_colors.front();
    }
    else
    {
        // evaluate actual BRDF? need to use the spherical function here as well to capture anything else than the diffuse light
        // if (!back_facing)
        {
            contribution = gi_point.sf_representation.color->get_value(surfel_to_rp);
        }
    }



    Gi_point_info point_info;
    point_info.type = Gi_point_info::Far_surfel;

    if (distance < disc_radius * surfel_near_threshold)
    {
        point_info.type = Gi_point_info::Near_surfel;
    }

    point_info.color = contribution;
    point_info.disc_normal = gi_point.normal;
    point_info.direction = -surfel_to_rp;
    point_info.radius = disc_radius * disc_scale_factor;
    point_info.depth = distance;
    point_info.position = gi_point.pos;
    point_info.receiver_position = vector3d_t(rp.P);
    point_info.solid_angle = solid_angle_real * disc_scale_factor;
    point_info.weight             = mix_amount;
    // point_info.spherical_function = gi_point.sf_representation;
    point_info.spherical_function = &gi_point.sf_representation;

    if (debug_info)
    {
        if (point_info.type == Gi_point_info::Near_surfel)
        {
            ++debug_info->used_near_surfels;
        }
        else
        {
            ++debug_info->used_far_surfels;
        }

        yafaray::GiPoint * debug_point = new GiPoint(gi_point);

        //debug_point->area = max_area;
        debug_point->area = visible_area * disc_scale_factor;
        // debug_point->area = gi_point.area;
        debug_point->color = contribution;
        debug_point->depth = debug_info->node_depth;
        debug_point->debug_radius = visible_radius; // debug only!

        debug_point->debug_mix_amount = mix_amount;

        debug_points.push_back(debug_point);
    }

    // frame_buffer.add_point(point_info, debug_point);
    point_infos.push_back(point_info);
}



void process_node(
        pbLighting_t::MyTree const* node,
        GiPoint const& gi_point,
        surfacePoint_t const& receiving_point,
        float const distance,
        vector3d_t const& node_to_rp,
        float const node_radius,
        float const mix_amount,
        float const disc_scale_factor,
        std::vector<Gi_point_info> & point_infos,
        std::vector<GiPoint*> & debug_points,
        Debug_info * debug_info)
{
    vector3d_t const& position = gi_point.pos;

    float const visible_area = std::max(0.0f, gi_point.sf_representation.area->get_value(node_to_rp));

    float const real_solid_angle = visible_area / (distance * distance);

    yafaray::color_t cluster_contribution;

    if (debug_info && debug_info->color_by_depth)
    {
        cluster_contribution = pbLighting_t::debug_colors[std::min(node->get_shortest_distance_to_leaf(), int(pbLighting_t::debug_colors.size()) - 1)];
    }
    else
    {
        cluster_contribution = gi_point.sf_representation.color->get_value(node_to_rp);
    }

    if (debug_info)
    {
        yafaray::GiPoint * debug_point = new yafaray::GiPoint(gi_point);

        debug_point->pos = position;

        debug_point->area = visible_area;

        debug_point->color        = cluster_contribution;
        debug_point->depth        = node->getDepth();
        debug_point->debug_radius = std::sqrt(visible_area / M_PI); // debug only!

        ++debug_info->used_nodes;

        debug_point->debug_mix_amount = mix_amount;

        debug_points.push_back(debug_point);
    }

    Gi_point_info point_info;
    point_info.type               = Gi_point_info::Node;
    point_info.color              = cluster_contribution;
    point_info.direction          = -node_to_rp;
    point_info.depth              = distance;
    point_info.position           = position;
    point_info.receiver_position  = vector3d_t(receiving_point.P);
    // point_info.solid_angle        = real_solid_angle * 2.0f;
    point_info.solid_angle        = real_solid_angle * disc_scale_factor;
    point_info.radius             = node_radius * disc_scale_factor;
    point_info.weight             = mix_amount;
    // point_info.radius             = std::sqrt(visible_area / M_PI);
    point_info.spherical_function = &gi_point.sf_representation;

    // frame_buffer.add_point(point_info, debug_point);
    point_infos.push_back(point_info);
}



struct Compare_point_info_by_depth
{
    bool operator() (Gi_point_info const& p1, Gi_point_info const& p2)
    {
        return (p1.depth < p2.depth);
    }
};



struct Point_and_debug_info {
    GiPoint*       debug_point;
    Gi_point_info* point_info;

    bool operator() (Point_and_debug_info const& p1, Point_and_debug_info const& p2)
    {
        return (p1.point_info->depth < p2.point_info->depth);
    }
};


color_t doPointBasedGiTree_sh_fb(
        pbLighting_t::MyTree const* tree,
        renderState_t & state,
        surfacePoint_t const& receiving_point,
        float const solid_angle_factor,
        background_t * background,
        vector3d_t const& wo,
        int const raster_buffer_resolution,
        Splat_cube_raster_buffer::Buffer_type const raster_buffer_type,
        Splat_cube_raster_buffer::Splat_type const node_splat_type,
        Splat_cube_raster_buffer::Splat_type const surfel_far_splat_type,
        Splat_cube_raster_buffer::Splat_type const surfel_near_splat_type,
        float const surfel_near_threshold,
        bool const variational,
        float const disc_scale_factor,
        Parameter_list const* parameters,
        Debug_info * debug_info
        )
{
    Splat_cube_raster_buffer frame_buffer;

    if (parameters)
    {
        frame_buffer.setup(*parameters);
    }
    else
    {
        frame_buffer.setup(raster_buffer_type, raster_buffer_resolution, node_splat_type, surfel_far_splat_type, surfel_near_splat_type);
    }

    std::vector<Gi_point_info> point_infos;
    std::vector<GiPoint*> debug_points;

    point_infos.reserve(10000);

    if (debug_info)
    {
        debug_points.reserve(10000);
    }

    float const fixed_max_solid_angle = frame_buffer.get_solid_angle(vector3d_t(1, 0, 0)) * solid_angle_factor;

    if (debug_info)
    {
        debug_info->my_timer.start("Traversal");
    }

    // traverse tree, if solid angle of node > max, traverse into the children
    std::queue<pbLighting_t::MyTree const*> queue;
    queue.push(tree);

    while (!queue.empty())
    {
        pbLighting_t::MyTree const* node = queue.front();
        queue.pop();

        if (debug_info)
        {
            debug_info->node_depth  = node->getDepth();
            // debug_info->node_height = node->get_node_height(0); // this  is really expensive for every node!
        }

        // if we are here and we have a leaf, we sample all points in the leaf
        if (node->getIsLeaf())
        {
            if (debug_info)
            {
                debug_info->my_timer.start("Surfel");
            }

            std::vector<GiPoint*> const& points = node->getData();

            assert(points.size() == 1); // true as long as we force 1 surfel per leaf in the tree

            // only use the clustered data, the SF stuff on the leaf's vector<Data> itself was deleted
            process_surfel(*node->getClusteredData(), receiving_point, frame_buffer, surfel_near_threshold, 1.0f, disc_scale_factor, true, point_infos, debug_points, debug_info);

//            for (unsigned int i = 0; i < points.size(); ++i)
//            {
//                yafaray::GiPoint const& giP = *points[i];
//                process_surfel(giP, receiving_point, frame_buffer, surfel_near_threshold, 1.0f, disc_scale_factor, true, point_infos, debug_points, debug_info);
//            }

            if (debug_info)
            {
                debug_info->my_timer.stop("Surfel");
            }
        }
        else if (node->has_children())
        {
            if (debug_info)
            {
                debug_info->my_timer.start("NodeCheck");
                ++debug_info->num_node_checks;
            }

            if (node->isNodeDataBehindPlane(vector3d_t(receiving_point.P), receiving_point.N, 0.1f))
            // if (node->is_node_data_behind_plane_approx(vector3d_t(receiving_point.P), receiving_point.N, bb_radius, 0.2f))
            {
                if (debug_info)
                {
                    debug_info->my_timer.stop("NodeCheck");
                }

                continue;
            }

            GiPoint const& gi_point = *node->getClusteredData();
            float const bb_radius = gi_point.bounding_box.get_enclosing_radius();

            vector3d_t node_to_rp = (vector3d_t(receiving_point.P) - gi_point.pos);

            float const distance = node_to_rp.length();

            node_to_rp.normalize();

            // float const node_radius = node->getRadius();
            float const node_radius = bb_radius;
            float const max_visible_area = node_radius * node_radius * M_PI;
            float const max_node_solid_angle = max_visible_area / (distance * distance);

            // using correction term to remove overshooting values at the "back"
            // float const visible_area = std::max(0.0f, (gi_point.sh_representation->get_area(giToSp) - 0.3f * max_visible_area) * 1.3f);

            // float const direction_max_solid_angle = frame_buffer.get_solid_angle(node_to_rp) * solid_angle_factor;
            float const direction_max_solid_angle = fixed_max_solid_angle;

            if (debug_info)
            {
                debug_info->my_timer.stop("NodeCheck");
            }

            if (max_node_solid_angle > direction_max_solid_angle /* solid_angle_threshold */ || distance < node_radius * 1.5f)
            // if (real_solid_angle > solid_angle_threshold || distance < node->getRadius() * 1.5f) // || distance < node->getRadius() * 1.5f)
            {
                std::vector<pbLighting_t::MyTree> const& children = node->getChildren();
                for (unsigned int i = 0; i < children.size(); ++i)
                {
                    queue.push(&children[i]);
                }
            }
            else
            {
                if (debug_info)
                {
                    debug_info->my_timer.start("Node");
                }


                if (variational)
                {
                    float const parent_radius = node->get_parent()->getClusteredData()->bounding_box.get_enclosing_radius();
                    float const parent_max_visible_area = parent_radius * parent_radius * M_PI;

                    float const dist_next_split = std::sqrt(max_visible_area / direction_max_solid_angle);

                    //                float const dist_previous_split = dist_next_split * 4.0f; // FIXME: improve by taking actual distance using the node's parent
                    float const dist_previous_split = std::sqrt(parent_max_visible_area / direction_max_solid_angle);

                    // FIXME: the clamping should not be necessary, somehow, the split distances are not correct yet
                    float const mix_amount = std::max(0.0f, std::min(1.0f, (dist_previous_split - distance) / (dist_previous_split - dist_next_split)));

                    if (debug_info && !(dist_previous_split > distance && distance > dist_next_split))
                    {
                        std::cout << "doPointBasedGiTree_sh_fb_variational(): " <<
                                     dist_previous_split << " " <<
                                     dist_next_split << " " <<
                                     distance << " " <<
                                     std::endl;
                    }

                    process_node(node, gi_point, receiving_point, distance, node_to_rp, node_radius, 1.0f - mix_amount, disc_scale_factor, point_infos, debug_points, debug_info);

                    std::vector<pbLighting_t::MyTree> const& children = node->getChildren();
                    for (unsigned int i = 0; i < children.size(); ++i)
                    {
                        pbLighting_t::MyTree const* child_node = &children[i];

                        if (child_node->getIsLeaf())
                        {
                            std::vector<GiPoint*> const& points = child_node->getData();

                            assert(points.size() == 1); // true as long as we force 1 surfel per leaf in the tree

                            process_surfel(*node->getClusteredData(), receiving_point, frame_buffer, surfel_near_threshold, mix_amount, disc_scale_factor, true, point_infos, debug_points, debug_info);

//                            for (unsigned int i = 0; i < points.size(); ++i)
//                            {
//                                yafaray::GiPoint const& giP = *points[i];

//                                process_surfel(giP, receiving_point, frame_buffer, surfel_near_threshold, mix_amount, disc_scale_factor, true, point_infos, debug_points, debug_info);
//                            }
                        }
                        else
                        {
                            GiPoint const& child_gi_point = *child_node->getClusteredData();
                            float const child_node_radius = child_gi_point.bounding_box.get_enclosing_radius();

                            vector3d_t const& child_position = child_gi_point.pos;

                            vector3d_t child_node_to_rp = (vector3d_t(receiving_point.P) - child_position);

                            float const child_distance = child_node_to_rp.length();
                            child_node_to_rp.normalize();

                            process_node(child_node, child_gi_point, receiving_point, child_distance, child_node_to_rp, child_node_radius, mix_amount, disc_scale_factor, point_infos, debug_points, debug_info);
                        }
                    }
                }
                else // non-variational
                {
                    process_node(node, gi_point, receiving_point, distance, node_to_rp, node_radius, 1.0f, disc_scale_factor, point_infos, debug_points, debug_info);
                }


                if (debug_info)
                {
                    debug_info->my_timer.stop("Node");
                }
            } // end else (bad solid angle)
        }
        else
        {
            assert(false);
        }
    }

    if (debug_info)
    {
        debug_info->my_timer.stop("Traversal");
        debug_info->my_timer.start("Accumulating");
    }

    if (debug_info)
    {
        assert(point_infos.size() == debug_points.size());

        std::vector<Point_and_debug_info> pd_infos(point_infos.size());

        for (unsigned int i = 0; i < point_infos.size(); ++i)
        {
            Point_and_debug_info & pdi = pd_infos[i];
            pdi.debug_point = debug_points[i];
            pdi.point_info  = &point_infos[i];
        }

        std::sort(pd_infos.begin(), pd_infos.end(), Point_and_debug_info());

        for (unsigned int i = 0; i < pd_infos.size(); ++i)
        {
            frame_buffer.add_point(*pd_infos[i].point_info, pd_infos[i].debug_point);
        }
    }
    else
    {
        // sort and splat, then integrate
        std::sort(point_infos.begin(), point_infos.end(), Compare_point_info_by_depth());

        for (unsigned int i = 0; i < point_infos.size(); ++i)
        {
            frame_buffer.add_point(point_infos[i]);
            //Gi_point_info const& pi = point_infos[i];
            //col += pi.color * std::max(0.0f, pi.direction * receiving_point.N) * pi.solid_angle;
        }

        //col *= 1.0f / (2.0f * M_PI);
    }

    // frame_buffer.add_background(background);


    color_t col(0.0f);

    if (receiving_point.material)
    {
//        col = frame_buffer.get_diffuse(receiving_point.N);
//        color_t surfCol = material->getDiffuseAtPoint(state, receiving_point); //material->eval(state, sp, wo, vector3d_t(0.0f), BSDF_ALL);
//        col *= surfCol;

        col = frame_buffer.get_brdf_response(state, receiving_point, wo);
    }
    else
    {
        col = frame_buffer.get_diffuse(receiving_point.N);
    }

    if (debug_info)
    {
        debug_info->my_timer.stop("Accumulating");

        if (debug_info->result_fb)
        {
            *debug_info->result_fb = frame_buffer;
        }

        std::cout <<
                     " NodeCheck: " << debug_info->my_timer.getTime("NodeCheck") <<
                     " Node: " << debug_info->my_timer.getTime("Node") <<
                     " Surfel: " << debug_info->my_timer.getTime("Surfel") <<
                     " Traversal: " << debug_info->my_timer.getTime("Traversal") <<
                     " Accumulating: " << debug_info->my_timer.getTime("Accumulating") <<
                     " num_node_checks: " << debug_info->num_node_checks <<
                     std::endl ;
    }

    return col;
}


#if 0

Splat_cube_raster_buffer gather_light_at_point(
        pbLighting_t::MyTree const* tree,
        vector3d_t const& receiving_point,
        float const solid_angle_factor,
        float const disc_scale_factor,
        int const raster_buffer_resolution,
        Splat_cube_raster_buffer::Buffer_type const raster_buffer_type,
        Splat_cube_raster_buffer::Splat_type const node_splat_type,
        Splat_cube_raster_buffer::Splat_type const surfel_far_splat_type,
        Splat_cube_raster_buffer::Splat_type const surfel_near_splat_type,
        float const surfel_near_threshold)
{
    Splat_cube_raster_buffer frame_buffer;
    frame_buffer.setup(raster_buffer_type, raster_buffer_resolution, node_splat_type, surfel_far_splat_type, surfel_near_splat_type);

    std::vector<Gi_point_info> point_infos;
    std::vector<GiPoint*> debug_points;

    surfacePoint_t receiving_point_sp; // dummy containing the position only
    receiving_point_sp.P = point3d_t(receiving_point);

    point_infos.reserve(10000);

    float const fixed_max_solid_angle = frame_buffer.get_solid_angle(vector3d_t(1, 0, 0)) * solid_angle_factor;


    // traverse tree, if solid angle of node > max, traverse into the children
    std::queue<pbLighting_t::MyTree const*> queue;
    queue.push(tree);

    while (!queue.empty())
    {
        pbLighting_t::MyTree const* node = queue.front();
        queue.pop();

        // if we are here and we have a leaf, we sample all points in the leaf
        if (node->getIsLeaf())
        {
            std::vector<GiPoint*> const& points = node->getData();

            assert(points.size() == 1); // true as long as we force 1 surfel per leaf in the tree

            for (unsigned int i = 0; i < points.size(); ++i)
            {
                yafaray::GiPoint const& giP = *points[i];

                process_surfel(giP, receiving_point_sp, frame_buffer, surfel_near_threshold, 1.0f, disc_scale_factor, false, point_infos, debug_points);
            }
        }
        else if (node->has_children())
        {
            GiPoint const& gi_point = *node->getClusteredData();
            float const bb_radius = gi_point.bounding_box.get_enclosing_radius();

            // vector3d_t const position = node->getCenter();
            vector3d_t const& position = gi_point.pos;

            vector3d_t node_to_rp = (vector3d_t(receiving_point) - position);

            float const distance = node_to_rp.length();

            node_to_rp.normalize();

            // float const node_radius = node->getRadius();
            float const node_radius = bb_radius;
            float const max_visible_area = node_radius * node_radius * M_PI;
            float const max_node_solid_angle = max_visible_area / (distance * distance);

            // using correction term to remove overshooting values at the "back"
            // float const visible_area = std::max(0.0f, (gi_point.sh_representation->get_area(giToSp) - 0.3f * max_visible_area) * 1.3f);

            // float const direction_max_solid_angle = frame_buffer.get_solid_angle(node_to_rp) * solid_angle_threshold;
            float const direction_max_solid_angle = fixed_max_solid_angle;

            if (max_node_solid_angle > direction_max_solid_angle /* solid_angle_threshold */ || distance < node_radius * 1.5f)
            // if (real_solid_angle > solid_angle_threshold || distance < node->getRadius() * 1.5f) // || distance < node->getRadius() * 1.5f)
            {
                std::vector<pbLighting_t::MyTree> const& children = node->getChildren();
                for (unsigned int i = 0; i < children.size(); ++i)
                {
                    queue.push(&children[i]);
                }
            }
            else
            {
                process_node(node, gi_point, receiving_point_sp, distance, node_to_rp, node_radius, 1.0f, disc_scale_factor, point_infos, debug_points, NULL);
            }
        }
        else
        {
            assert(false);
        }
    }

    // sort and splat, then integrate
    std::sort(point_infos.begin(), point_infos.end(), Compare_point_info_by_depth());

    for (unsigned int i = 0; i < point_infos.size(); ++i)
    {
        frame_buffer.add_point(point_infos[i]);
    }

    // frame_buffer.add_background(background);

    return frame_buffer;
}

#endif



colorA_t pbLighting_t::integrate(renderState_t &state, diffRay_t &ray) const
{
    if (render_single_pixel && (pixel_x != state.pixel_x || pixel_y != state.pixel_y))
    {
        return color_t(0.0);
    }

    color_t col(0.0);
    float alpha = 0.0;
    surfacePoint_t sp;
    void *o_udat = state.userdata;
    bool oldIncludeLights = state.includeLights;

    if(scene->intersect(ray, sp)) // If it hits
    {
        unsigned char userdata[USER_DATA_SIZE];
        const material_t *material = sp.material;
        BSDF_t bsdfs;

        state.userdata = (void *) userdata;
        vector3d_t wo = -ray.dir;
        if(state.raylevel == 0) state.includeLights = true;

        material->initBSDF(state, sp, bsdfs);

        if(bsdfs & BSDF_EMIT)
        {
            col += material->emission(state, sp, wo);
        }

        // if(bsdfs & BSDF_DIFFUSE)
        {
            if (!indirectOnly)
            {
                col += estimateAllDirectLight(state, sp, wo);
            }

            if (debug_type == Tree_sh_fb)
            {

                col += doPointBasedGiTree_sh_fb(_bspTree,
                                                state,
                                                sp,
                                                solid_angle_factor,
                                                background,
                                                wo,
                                                raster_buffer_resolution,
                                                raster_buffer_type,
                                                node_splat_type,
                                                surfel_far_splat_type,
                                                surfel_near_splat_type,
                                                surfel_near_threshold,
                                                variational,
                                                _disc_scale_factor
                                                );
            }
            else
            {
                assert(false);
            }
        }

        if (!indirectOnly && !(bsdfs & BSDF_GLOSSY))
        {
            recursiveRaytrace(state, ray, bsdfs, sp, wo, col, alpha);
        }

        float m_alpha = material->getAlpha(state, sp, wo);
        alpha = m_alpha + (1.f - m_alpha) * alpha;
    }
    else // Nothing hit, return background if any
    {
        if (background) col += (*background)(ray, state, false);
    }

    state.userdata = o_udat;
    state.includeLights = oldIncludeLights;
    return colorA_t(col, alpha);
}

integrator_t* pbLighting_t::factory(paraMap_t &params, renderEnvironment_t &render)
{
    bool transpShad=false;
    bool caustics=false;
    bool do_AO=false;
    int shadowDepth=5;
    int raydepth=5, cDepth=10;
    int search=100, photons=500000;
    int AO_samples = 32;
    double cRad = 0.25;
    double AO_dist = 1.0;
    color_t AO_col(1.f);
    int samples = 10;
    bool debug = false;
    bool indirectOnly = false;
    float maxSolidAngle = 1.0f;
    int debugTreeDepth = 2;
    bool debugOutputPointsToFile = false;
    std::string debug_type = "NoTree";
    std::map<std::string, Debug_type> debug_type_map;
    debug_type_map["NoTree"] = NoTree;
    debug_type_map["Tree"] = Tree;
    debug_type_map["Tree_sh"] = Tree_sh;
    debug_type_map["Tree_sh_fb"] = Tree_sh_fb;
    debug_type_map["Tree_sh_leafs"] = Tree_sh_leafs;
    bool render_single_pixel = false;
    int pixel_x = -1;
    int pixel_y = -1;
    bool do_load_gi_points = false;
    int fb_resolution = 8;
    std::string fb_type = Splat_cube_raster_buffer::enum_type_map.begin()->first;
    std::string node_splat_type = Splat_cube_raster_buffer::enum_splat_type_map.begin()->first;
    std::string surfel_far_splat_type = Splat_cube_raster_buffer::enum_splat_type_map.begin()->first;
    std::string surfel_near_splat_type = Splat_cube_raster_buffer::enum_splat_type_map.begin()->first;
    float surfel_near_threshold = 2.0f;
    std::string dictionary_type_str = "No_dict";
    int sf_resolution = 4;
    int dict_num_centers = 100;
    float dictionary_sample_fraction = 0.1f;
    bool do_dictionary_stats = true;
    bool use_sf_files = false;
    bool variational = false;
    bool enable_conversion = false;
    float disc_scale_factor = 1.0f;

    params.getParam("raydepth", raydepth);
    params.getParam("transpShad", transpShad);
    params.getParam("shadowDepth", shadowDepth);
    params.getParam("caustics", caustics);
    params.getParam("photons", photons);
    params.getParam("caustic_mix", search);
    params.getParam("caustic_depth", cDepth);
    params.getParam("caustic_radius", cRad);
    params.getParam("do_AO", do_AO);
    params.getParam("AO_samples", AO_samples);
    params.getParam("AO_distance", AO_dist);
    params.getParam("AO_color", AO_col);

    params.getParam("samplesPerArea", samples);
    params.getParam("debug", debug);
    params.getParam("indirectOnly", indirectOnly);
    params.getParam("maxSolidAngle", maxSolidAngle);
    params.getParam("debugTreeDepth", debugTreeDepth);
    params.getParam("debugOutputPointsToFile", debugOutputPointsToFile);
    params.getParam("debug_type", debug_type);
    params.getParam("render_single_pixel", render_single_pixel);
    params.getParam("pixel_x", pixel_x);
    params.getParam("pixel_y", pixel_y);
    params.getParam("do_load_gi_points", do_load_gi_points);
    params.getParam("fb_type", fb_type);
    params.getParam("fb_resolution", fb_resolution);
    params.getParam("node_splat_type", node_splat_type);
    params.getParam("surfel_far_splat_type", surfel_far_splat_type);
    params.getParam("surfel_near_splat_type", surfel_near_splat_type);
    params.getParam("surfel_near_threshold", surfel_near_threshold);
    params.getParam("dictionary_type", dictionary_type_str);
    params.getParam("dict_num_centers", dict_num_centers);
    params.getParam("dict_sample_fraction", dictionary_sample_fraction);
    params.getParam("do_dict_stats", do_dictionary_stats);
    params.getParam("sf_resolution", sf_resolution);
    params.getParam("use_sf_files", use_sf_files);
    params.getParam("variational", variational);
    params.getParam("enable_conversion", enable_conversion);
    params.getParam("disc_scale_factor", disc_scale_factor);



    pbLighting_t *inte = new pbLighting_t(transpShad, shadowDepth, raydepth);
    // caustic settings
    inte->usePhotonCaustics = caustics;
    inte->nCausPhotons = photons;
    inte->nCausSearch = search;
    inte->causDepth = cDepth;
    inte->causRadius = cRad;
    // AO settings
    inte->useAmbientOcclusion = do_AO;
    inte->aoSamples = AO_samples;
    inte->aoDist = AO_dist;
    inte->aoCol = AO_col;

    inte->surfel_samples = samples;
    inte->debug = debug;
    inte->indirectOnly = indirectOnly;
    inte->debugTreeDepth = debugTreeDepth;
    inte->debugOutputPointsToFile = debugOutputPointsToFile;
    if (debug_type_map.find(debug_type) != debug_type_map.end())
    {
        inte->debug_type = debug_type_map[debug_type];
    }

    std::cout << "Debug type: " << inte->debug_type << " " << debug_type << std::endl;
    inte->debug_type_str = debug_type;

    inte->render_single_pixel = render_single_pixel;
    inte->pixel_x             = pixel_x;
    inte->pixel_y             = pixel_y;
    inte->variational         = variational;

    inte->do_load_gi_points = do_load_gi_points;

    inte->raster_buffer_resolution = fb_resolution;
    inte->raster_buffer_type       = Splat_cube_raster_buffer::enum_type_map[fb_type];
    inte->node_splat_type          = Splat_cube_raster_buffer::enum_splat_type_map[node_splat_type];
    inte->surfel_far_splat_type    = Splat_cube_raster_buffer::enum_splat_type_map[surfel_far_splat_type];
    inte->surfel_near_splat_type   = Splat_cube_raster_buffer::enum_splat_type_map[surfel_near_splat_type];
    inte->surfel_near_threshold    = surfel_near_threshold;

    inte->solid_angle_factor = maxSolidAngle;

    std::string spherical_function_factory_type = "SH";
    params.getParam("spherical_function_type", spherical_function_factory_type);
    if (spherical_function_factory_type == "Cube")
    {
        inte->_spherical_function_color_factory = new Cube_spherical_function_factory<color_t>(sf_resolution, false);
    }
    else if (spherical_function_factory_type == "SH")
    {
        inte->_spherical_function_color_factory = new Spherical_harmonics_factory<color_t>(3, true);
    }

    // inte->_spherical_function_area_factory = new Cube_spherical_function_factory<float>(4, false);
    inte->_spherical_function_area_factory = new Spherical_harmonics_factory<float>(3, true);

    inte->_area_converter  = NULL;
    inte->_color_converter = NULL;

    if (enable_conversion)
    {
        inte->_color_converter = new Cube_to_mises_fisher_converter<color_t>(1);
    }

    inte->_sf_resolution              = sf_resolution;
    inte->_disc_scale_factor          = disc_scale_factor;
    inte->_use_sf_files               = use_sf_files;
    inte->_dict_num_centers           = dict_num_centers;
    inte->_dictionary_sample_fraction = dictionary_sample_fraction;
    inte->_do_dictionary_stats        = do_dictionary_stats;

    if (dictionary_type_str == "No_dict")
    {
        inte->_color_dictionary_generator = NULL;
    }
    else if (dictionary_type_str == "Random_dict")
    {
        inte->_color_dictionary_generator = new Random_dictionary_generator;
    }
    else if (dictionary_type_str == "Kmeans_dict")
    {
        inte->_color_dictionary_generator = new Kmeans_dictionary_generator(true);
    }


    Y_INFO <<
              "maxSolidAngle: " << inte->solid_angle_factor << " " <<
              "raster_buffer_type: " << inte->raster_buffer_type << " " <<
              "fb_resolution: " << inte->raster_buffer_resolution << " " <<
              "node_splat_type: " << inte->node_splat_type << " " <<
              "surfel_far_splat_type: " << inte->surfel_far_splat_type << " " <<
              "surfel_near_splat_type: " << inte->surfel_near_splat_type << " " <<
              "surfel_near_threshold: " << inte->surfel_near_threshold << " " <<
              "dict_num_centers: " << inte->_dict_num_centers << " " <<
              "dict_sample_fraction: " << inte->_dictionary_sample_fraction << " " <<
              "sf_resolution: " << inte->_sf_resolution << " " <<
              std::endl;


    return inte;
}

extern "C"
{

YAFRAYPLUGIN_EXPORT void registerPlugin(renderEnvironment_t &render)
{
    render.registerFactory("pbgi", pbLighting_t::factory);
}

}

__END_YAFRAY
