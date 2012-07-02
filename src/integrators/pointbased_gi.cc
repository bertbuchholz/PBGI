#include <sstream>
#include <cassert>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <map>


#include <omp.h>
#include <valgrind/callgrind.h>

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


/*
vector3d_t const& Gi_point_base::get_normal()
{
    assert(_is_surfel);

    Gi_point_surfel * n = static_cast<Gi_point_surfel*>(this);
    return n->get_normal();
}

float Gi_point_base::get_area(vector3d_t const& dir)
{
    if (_is_surfel)
    {
        Gi_point_surfel * n = static_cast<Gi_point_surfel*>(this);
        return n->get_area();
    }
    else
    {
        Gi_point_inner_node * n = static_cast<Gi_point_inner_node*>(this);
        return n->get_area(dir);
    }
}

color_t Gi_point_base::get_color(vector3d_t const& dir)
{
    if (_is_surfel)
    {
        Gi_point_surfel * n = static_cast<Gi_point_surfel*>(this);
        return n->get_color();
    }
    else
    {
        Gi_point_inner_node * n = static_cast<Gi_point_inner_node*>(this);
        return n->get_color(dir);
    }
}
*/

pbLighting_t::pbLighting_t(bool transpShad, int shadowDepth, int rayDepth) : _solid_angle_factor(0.5f)
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


Cube_raster_buffer<Spherical_function<float> *> precalculate_spherical_harmonics(Spherical_function_factory<float> const* _area_factory)
{
    std::cout << "precalculate_spherical_harmonics()" << std::endl;

    int const resolution = 32;

    Cube_raster_buffer<Spherical_function<float> *> normal_cube;
    normal_cube.setup_surfel_buffer(resolution);

    std::vector<Cube_cell> const& cells = normal_cube.get_cube_cells();

    for (size_t i = 0; i < cells.size(); ++i)
    {
        vector3d_t normal = normal_cube.get_cell_direction(cells[i]);
        Abstract_spherical_function_estimator<float> * area_estimator = new Spherical_function_area_estimator<float>(normal, 1.0f);

        // Spherical_function<float> * sf = new Spherical_harmonics<float>(exact_sh, 3);
        Spherical_function<float> * sf = _area_factory->create();
        sf->calc_coefficients_random(area_estimator);

        normal_cube.set_data(cells[i], sf);

        delete area_estimator;
    }

    return normal_cube;
}


void pbLighting_t::generate_spherical_function(renderState_t & state, GiPoint * gi_point, surfacePoint_t const& sp, std::vector<light_t*> const& lights)
{
    ray_t lightRay;
    lightRay.from = gi_point->pos;

    color_t light_color(0.f);

    for (std::vector<light_t *>::const_iterator light = lights.begin(); light != lights.end(); ++light)
    {
        if ((*light)->diracLight())
        {
            bool shadowed = true;

            if ((*light)->illuminate(sp, light_color, lightRay))
            {
                // ...shadowed...
                lightRay.tmin = YAF_SHADOW_BIAS; // < better add some _smart_ self-bias value...this is bad.
                shadowed = scene->isShadowed(state, lightRay);
                float const cos_sp_N_lightray_dir = sp.N * lightRay.dir;
                if (!shadowed && cos_sp_N_lightray_dir > 0.0f)
                {
                    color_t const incoming_light = light_color;
                    // color_t const reflected_light = lcol * std::max(0.0f, cos_sp_N_lightray_dir);
                    gi_point->color += incoming_light;
                }
                else
                {
                    light_color = color_t(0.0f);
                }
            }

            if (!shadowed)
            {
                if (_use_precalculated_sf)
                {
                    float const lambert_scale = lightRay.dir * sp.N;
                    color_t scale = sp.material->getDiffuseAtPoint(state, sp) * light_color * lambert_scale;
                    gi_point->sf_representation.color->get_precalculated_coefficients(_precalculated_sf, scale, sp.N);
                }
                else
                {
                    Abstract_spherical_function_estimator<color_t> * color_estimator = new Spherical_function_light_color_estimator<color_t>(state, sp, lightRay, light_color, gi_point->area);
                    gi_point->sf_representation.color->calc_coefficients_random(color_estimator);
                    delete color_estimator;
                }
            }
        }
    }

    if (sp.material->getFlags() & BSDF_EMIT)
    {
        Abstract_spherical_function_estimator<color_t> * emit_estimator = new Spherical_function_emit_color_estimator<color_t>(state, sp, gi_point->area);
        gi_point->sf_representation.color->calc_coefficients_random(emit_estimator);
        delete emit_estimator;
    }

    if (_use_precalculated_sf)
    {
        gi_point->sf_representation.area->get_precalculated_coefficients(_precalculated_sf, gi_point->area, sp.N);
    }
    else
    {
        Abstract_spherical_function_estimator<float> * area_estimator = new Spherical_function_area_estimator<float>(sp.N, gi_point->area);
        gi_point->sf_representation.area->calc_coefficients_random(area_estimator);
        delete area_estimator;
    }
}



Spherical_node_representation pbLighting_t::generate_spherical_function_2(renderState_t & state, Gi_point_surfel & surfel, surfacePoint_t const& sp, std::vector<light_t*> const& lights)
{
    ray_t lightRay;
    lightRay.from = surfel.pos;



    Spherical_node_representation sf_representation;

    sf_representation.color = _spherical_function_color_factory->create();
    sf_representation.area  = _spherical_function_area_factory->create();

    color_t accumulated_reflected_light(0.0f);

    for (std::vector<light_t *>::const_iterator light = lights.begin(); light != lights.end(); ++light)
    {
        if ((*light)->diracLight())
        {
            bool shadowed = true;
            bool back_facing = true;

            color_t light_color(0.f);

            if ((*light)->illuminate(sp, light_color, lightRay))
            {
                // ...shadowed...
                lightRay.tmin = YAF_SHADOW_BIAS; // < better add some _smart_ self-bias value...this is bad.
                shadowed = scene->isShadowed(state, lightRay);
                float const cos_sp_N_lightray_dir = sp.N * lightRay.dir;

                back_facing = cos_sp_N_lightray_dir < 0.0f;

                if (!shadowed && !back_facing)
                {
#ifndef SURFEL_HAS_COLOR_SF
                    color_t const incoming_light = light_color;
                    // color_t const reflected_light = lcol * std::max(0.0f, cos_sp_N_lightray_dir);
                    float const lambert_scale = lightRay.dir * sp.N;
                    surfel.color += incoming_light * lambert_scale;
#endif
                }
                else
                {
                    light_color = color_t(0.0f);
                }
            }

            if (!shadowed && !back_facing)
            {
                if (_use_precalculated_sf)
                {
                    float const lambert_scale = lightRay.dir * sp.N;
                    accumulated_reflected_light += lambert_scale * light_color;
                    //color_t scale = sp.material->getDiffuseAtPoint(state, sp) * light_color * lambert_scale * surfel.area;
                    //sf_representation.color->get_precalculated_coefficients(_precalculated_sf, scale, sp.N);
                }
                else
                {
                    Abstract_spherical_function_estimator<color_t> * color_estimator = new Spherical_function_light_color_estimator<color_t>(state, sp, lightRay, light_color, surfel.area);
                    sf_representation.color->calc_coefficients_random(color_estimator);
                    delete color_estimator;
                }
            }
        }
    }

    // FIXME: emit should be added, right now it replaces the color if a material has diffuse and emit
    if (!_use_precalculated_sf && (sp.material->getFlags() & BSDF_EMIT))
    {
        Abstract_spherical_function_estimator<color_t> * emit_estimator = new Spherical_function_emit_color_estimator<color_t>(state, sp, surfel.area);
        sf_representation.color->calc_coefficients_random(emit_estimator);
        delete emit_estimator;
    }

    if (_use_precalculated_sf)
    {
        sf_representation.area->get_precalculated_coefficients(_precalculated_sf, surfel.area, sp.N);

        color_t const scale = (sp.material->getDiffuseAtPoint(state, sp) * accumulated_reflected_light + sp.material->emission(state, sp, vector3d_t())) * surfel.area;
        sf_representation.color->get_precalculated_coefficients(_precalculated_sf, scale, sp.N);
    }
    else
    {
        Abstract_spherical_function_estimator<float> * area_estimator = new Spherical_function_area_estimator<float>(sp.N, surfel.area);
        sf_representation.area->calc_coefficients_random(area_estimator);
        delete area_estimator;
    }

#ifndef SURFEL_HAS_COLOR_SF
    surfel.color *= sp.material->getDiffuseAtPoint(state, sp);
    surfel.color += sp.material->emission(state, sp, vector3d_t());
#endif

    return sf_representation;
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




#ifdef USE_FAT_TREE
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

        if (node->is_leaf())
        {
//            std::cout << "pbLighting_t::generate_gi_points_data: leaf" << std::endl;

            std::vector<GiPoint*> node_data = node->get_points_data();
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

            bool const treat_as_leaf = (node->get_depth() == treat_depth_as_leaf);

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

                    if (children[j].is_leaf())
                    {
                        assert(children[j].get_points_data().size() == 1);
                        //tmp = children[j].getData()[0]->sf_representation.area;
                        //children[j].getData()[0]->sf_representation.area = _area_converter->convert(children[j].getData()[0]->sf_representation.area);
                        //delete tmp;
                        delete children[j].get_points_data()[0]->sf_representation.area;
                        children[j].get_points_data()[0]->sf_representation.area = NULL;
                    }

                    tmp = children[j].get_data()->sf_representation.area;
                    children[j].get_data()->sf_representation.area = _area_converter->convert(children[j].get_data()->sf_representation.area);
                    delete tmp;
                }
                else
                {
                    if (children[j].is_leaf())
                    {
                        assert(children[j].get_points_data().size() == 1);
                        delete children[j].get_points_data()[0]->sf_representation.area;
                        children[j].get_points_data()[0]->sf_representation.area = NULL;
                    }
                }

                if (_color_converter && !treat_as_leaf)
                {
                    //                std::cout << "pbLighting_t::generate_gi_points_data: convert children" << std::endl;

                    Spherical_function<color_t> * tmp;

                    if (children[j].is_leaf())
                    {
                        //                        std::cout << "pbLighting_t::generate_gi_points_data: convert leaf" << std::endl;

                        //                        tmp = children[j].getData()[0]->sf_representation.color;
                        //                        children[j].getData()[0]->sf_representation.color = _color_converter->convert(children[j].getData()[0]->sf_representation.color);
                        //                        delete tmp;

                        assert(children[j].get_points_data().size() == 1);
                        delete children[j].get_points_data()[0]->sf_representation.color;
                        children[j].get_points_data()[0]->sf_representation.color = NULL;
                    }

                    tmp = children[j].get_data()->sf_representation.color;
                    children[j].get_data()->sf_representation.color = _color_converter->convert(children[j].get_data()->sf_representation.color);
                    delete tmp;
                }
                else
                {
                    if (children[j].is_leaf())
                    {
                        assert(children[j].get_points_data().size() == 1);
                        delete children[j].get_points_data()[0]->sf_representation.color;
                        children[j].get_points_data()[0]->sf_representation.color = NULL;
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

        if (node->is_leaf())
        {
//            std::cout << "pbLighting_t::generate_gi_points_data: leaf" << std::endl;

            std::vector<GiPoint*> node_data = node->get_points_data();
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

            bool const treat_as_leaf = (node->get_depth() == treat_depth_as_leaf);

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

                    if (children[j].is_leaf())
                    {
                        assert(children[j].get_points_data().size() == 1);
                        tmp = children[j].get_points_data()[0]->sf_representation.area;
                        children[j].get_points_data()[0]->sf_representation.area = _area_converter->convert(children[j].get_points_data()[0]->sf_representation.area);
                        delete tmp;
                    }

                    tmp = children[j].get_data()->sf_representation.area;
                    children[j].get_data()->sf_representation.area = _area_converter->convert(children[j].get_data()->sf_representation.area);
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

                    if (children[j].is_leaf())
                    {
//                        std::cout << "pbLighting_t::generate_gi_points_data: convert leaf" << std::endl;

                        assert(children[j].get_points_data().size() == 1);
                        tmp = children[j].get_points_data()[0]->sf_representation.color;
                        children[j].get_points_data()[0]->sf_representation.color = _color_converter->convert(children[j].get_points_data()[0]->sf_representation.color);
                        delete tmp;
                    }

                    tmp = children[j].get_data()->sf_representation.color;
                    children[j].get_data()->sf_representation.color = _color_converter->convert(children[j].get_data()->sf_representation.color);
                    delete tmp;
                }
            }
        }
    }

    std::cout << "generate_gi_points_data: used cells/all cells: " << Cube_spherical_function<color_t>::_num_used_cells << "/" << Cube_spherical_function<color_t>::_num_all_cells << std::endl;
}
#else

void pbLighting_t::generate_gi_points_data_2(std::vector<MyTree::Tree_node*> const& nodes,
                                             std::tr1::unordered_map<vector3d_t, pbgi_sample_t*, My_hash> const& point_to_sampling_point_map,
                                             float const area,
                                             float const radius,
                                             int const treat_depth_as_leaf)
{
    std::cout << "pbLighting_t::generate_gi_points_data_2(), small tree, nodes: " << nodes.size() << std::endl;

    Gi_point_averager averager(_spherical_function_color_factory, _spherical_function_area_factory);

    renderState_t state;
    unsigned char userdata[USER_DATA_SIZE+7];
    state.userdata = (void *)( &userdata[7] - ( ((size_t)&userdata[7])&7 ) ); // pad userdata to 8 bytes
    state.cam = scene->getCamera();

    for (unsigned int node_index = 0; node_index < nodes.size(); ++node_index)
    {
        MyTree::Tree_node* node = nodes[node_index];

        if ((node_index % 1000) == 0)
        {
            std::cout << "." << std::flush;
        }

        if (node->is_leaf())
        {
//            std::cout << "pbLighting_t::generate_gi_points_data: leaf" << std::endl;

            MyTree::Tree_leaf_node * leaf_node = node->get_derived<MyTree::Tree_leaf_node>();

//            Gi_point_inner_node gi_point(_spherical_function_color_factory, _spherical_function_area_factory);
//            gi_point.bounding_box = bound_t(point3d_t(1e10f, 1e10f, 1e10f), point3d_t(-1e10f, -1e10f, -1e10f));

            std::vector<Gi_point_surfel> & surfels = leaf_node->get_surfels();
            assert(surfels.size() > 0);

            for (size_t surfel_index = 0; surfel_index < surfels.size(); ++surfel_index)
            {

                Gi_point_surfel & surfel = surfels[surfel_index];

                std::tr1::unordered_map<vector3d_t, pbgi_sample_t*>::const_iterator iter = point_to_sampling_point_map.find(surfel.pos);
                assert(iter != point_to_sampling_point_map.end());

                pbgi_sample_t & sampling_point = *(iter->second);

                point3d_t const& hitPoint = sampling_point.position;
                triangle_t const& tri = *sampling_point.tri_pointer;
                intersectData_t & iData = sampling_point.intersect_data;

                surfacePoint_t sp;

                tri.getSurface(sp, hitPoint, iData);

                assert(!std::isnan(sp.P.x));

                surfel.normal       = sp.N;
                surfel.area         = area;
#ifndef SURFEL_HAS_COLOR_SF
                surfel.color        = color_t(0.0f);
#endif



                unsigned char userdata[USER_DATA_SIZE];
                state.userdata = (void *) userdata;

                material_t const* material = sp.material;
                BSDF_t bsdfs;
                material->initBSDF(state, sp, bsdfs);

                surfel.color = estimateAllDirectLight(state, sp, surfel.normal);

                if(bsdfs & BSDF_EMIT)
                {
                    surfel.color += material->emission(state, sp, surfel.normal);
                }

#ifdef PBGI_DEBUG
                surfel.debug_radius = 0.0f;
                surfel.diffuse_color = material->getDiffuseAtPoint(state, sp);
#endif

                /*
                Spherical_node_representation sf_representation = generate_spherical_function_2(state, surfel, sp, lights);

                gi_point.sf_representation.area->add(sf_representation.area);
                gi_point.sf_representation.color->add(sf_representation.color);

                bound_t bounding_box = generate_disc_bounding_box(vector3d_t(sp.P), sp.NU, sp.NV, radius);

                for (int j = 0; j < 3; ++j)
                {
                    gi_point.bounding_box.a[j] = std::min(gi_point.bounding_box.a[j], bounding_box.a[j]);
                    gi_point.bounding_box.g[j] = std::max(gi_point.bounding_box.g[j], bounding_box.g[j]);
                }

                gi_point.pos += surfel.pos;

#ifdef SURFEL_HAS_COLOR_SF
                if (_color_converter)
                {
                    Spherical_function<color_t> * tmp = sf_representation.color;
                    surfel.color = _color_converter->convert(sf_representation.color, sf_representation.area);
                    delete tmp;
                }
                else
                {
                    surfel.color = sf_representation.color;
                }
#endif

                delete sf_representation.area;

#ifndef SURFEL_HAS_COLOR_SF
                delete sf_representation.color;
#endif
                */
            }

//            gi_point.pos *= 1.0f / float(surfels.size());
//            node->set_data(gi_point);

            Gi_point_inner_node const leaf_node_data = averager.average_leaf(surfels, _precalculated_sf);
            node->set_data(leaf_node_data);

#ifdef SURFEL_HAS_COLOR_SF
            // color sh is computed in
            if (_color_converter)
            {
                for (size_t surfel_index = 0; surfel_index < surfels.size(); ++surfel_index)
                {
                    Gi_point_surfel & surfel = surfels[surfel_index];

                    Spherical_function<color_t> * tmp = sf_representation.color;
                    surfel.color = _color_converter->convert(sf_representation.color, sf_representation.area);
                    delete tmp;
                }
            }
#endif
        }
        else
        {
//            std::cout << "pbLighting_t::generate_gi_points_data: node" << std::endl;

            bool const treat_as_leaf = (node->get_depth() == treat_depth_as_leaf);

            if (!treat_as_leaf)
            {
                node->average_node(averager);
            }

            MyTree::Tree_inner_node const* inner_node = node->get_derived<MyTree::Tree_inner_node>();
            std::vector<MyTree::Tree_node*> const& children = inner_node->get_children();

            for (size_t i = 0; i < children.size(); ++i)
            {
                MyTree::Tree_node* child = children[i];

                Gi_point_inner_node & gi_point = child->get_data();

                // convert the children, not the node itself (yet)
                if (_area_converter && !treat_as_leaf)
                {
                    Spherical_function<float> * tmp;

                    tmp = gi_point.sf_representation.area;
                    gi_point.sf_representation.area = _area_converter->convert(gi_point.sf_representation.area);
                    delete tmp;
                }

                if (_color_converter && !treat_as_leaf)
                {
                    //                std::cout << "pbLighting_t::generate_gi_points_data: convert children" << std::endl;

                    Spherical_function<color_t> * tmp;

                    tmp = gi_point.sf_representation.color;
                    gi_point.sf_representation.color = _color_converter->convert(gi_point.sf_representation.color, gi_point.sf_representation.area);
                    delete tmp;
                }

                if (child->is_leaf())
                {
#ifndef SURFEL_HAS_COLOR_SF
                    delete gi_point.sf_representation.color;
                    gi_point.sf_representation.color = NULL;
#endif
                    delete gi_point.sf_representation.area;
                    gi_point.sf_representation.area = NULL;
                }
            }
        }
    }

    std::cout << "generate_gi_points_data: used cells/all cells: " << Cube_spherical_function<color_t>::_num_used_cells << "/" << Cube_spherical_function<color_t>::_num_all_cells << std::endl;
}
#endif

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

// create an enclosing cube shaped bound
bound_t create_fitting_cube(bound_t const& bound)
{
    int const longest_axis = bound.longestAxis();
    float const longest_extent_2 = bound.get_length(longest_axis) / 2.0f;

    vector3d_t center = vector3d_t(bound.center());

    bound_t result;

    result.a = center - vector3d_t(longest_extent_2 * 1.001f);
    result.g = center + vector3d_t(longest_extent_2 * 1.001f);

    return result;
}

void pbLighting_t::generate_gi_points(renderState_t & state, MyTree & tree, int const number_of_samples)
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

    bool const use_cdf_sampler = true;
    bool const do_find_widest_gap = false;

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


    std::cout << "generate_gi_points(): creating points" << std::endl;

    // std::tr1::unordered_map<Gi_point_base*, pbgi_sample_t*> gi_point_to_sampling_point_map;
    std::tr1::unordered_map<vector3d_t, pbgi_sample_t*, My_hash> point_to_sampling_point_map;

#ifdef USE_FAT_TREE
    tree.prepare_containers(sampling_points.size());
#else
    std::vector<vector3d_t> point_data;
    point_data.reserve(sampling_points.size());
#endif

    timer.addEvent("point_creation");
    timer.start("point_creation");

    // -------------------------------------------
    // Put the GI points into the tree and create a correspondency map between
    // point and its surface point
    for (int i = 0; i < int(sampling_points.size()); ++i)
    {
        // Gi_point_base * gi_point = new GiPoint(NULL, NULL);
        // gi_point_to_sampling_point_map[gi_point] = &sampling_points[i];
        point_to_sampling_point_map[vector3d_t(sampling_points[i].position)] = &sampling_points[i];
#ifdef USE_FAT_TREE
        tree.add_point_prepared(i, vector3d_t(sampling_points[i].position), gi_point);
#else
        point_data.push_back(vector3d_t(sampling_points[i].position));
        // point_data.push_back(Point_data<GiPoint>(vector3d_t(sampling_points[i].position), gi_point));
#endif
    }

    std::cout << "generate_gi_points(): finish tree" << std::endl;

#ifdef USE_FAT_TREE
    MyTree::finalize_points_iterative(&tree);
#else
    My_tree_subdivision_decider subdivision_decider;

    subdivision_decider.max_num_points = 0;
    subdivision_decider.max_size = needed_radius * 1.0f;

    if (MyTree::Arity == 2)
    {
        tree.build_tree(point_data, scene->getSceneBound(), subdivision_decider);
    }
    else
    {
        bound_t cube_bound = create_fitting_cube(scene->getSceneBound());
        tree.build_tree(point_data, cube_bound, subdivision_decider);
    }
#endif

    timer.stop("point_creation");
    // ------------------------------------------


    // ------------------------------------------
    // go through all nodes in post order, generate the SF representation in the children,
    // then go up to the parent, use the children to generate its SF and convert it and delete the children's unconverted etc.

    timer.addEvent("node_data_creation");
    timer.start("node_data_creation");

    bool const use_parallel_averaging = false;

#ifdef USE_FAT_TREE
    if (use_parallel_averaging)
    {
        int const num_procs = omp_get_num_procs();
        int const splitting_depth = std::ceil(std::log(num_procs) / std::log(2.0f));
        int const num_threads = std::pow(2, splitting_depth);

        std::cout << "generate_gi_points(): parallel averaging, num_threads: " << num_threads << " depth: " << splitting_depth << std::endl;

        std::vector<MyTree*> nodes_in_splitting_depth = tree.get_nodes(splitting_depth);

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
        tree.get_post_order_queue(nodes, splitting_depth);

        //generate_gi_points_data(state, nodes, point_to_surface_point_map, splitting_depth);
        generate_gi_points_data_2(nodes, gi_point_to_sampling_point_map, area_estimation, needed_radius, splitting_depth);
    }
    else
    {
        std::vector<MyTree*> nodes;
        tree.get_post_order_queue(nodes);

        generate_gi_points_data_2(nodes, gi_point_to_sampling_point_map, area_estimation, needed_radius);
    }
#else
    if (use_parallel_averaging)
    {
        int const num_procs = omp_get_num_procs();
        int const splitting_depth = std::ceil(std::log(num_procs) / std::log(2.0f));
        int const num_threads = std::pow(2, splitting_depth);

        std::cout << "generate_gi_points(): parallel averaging, num_threads: " << num_threads << " depth: " << splitting_depth << std::endl;

        std::vector<MyTree::Tree_node*> nodes_in_splitting_depth = tree.get_nodes_at_depth(splitting_depth);

        assert(int(nodes_in_splitting_depth.size()) == num_threads);

        #pragma omp parallel for
        for (int i = 0; i < num_threads; ++i)
        {
            renderState_t parallel_state;
            unsigned char userdata[USER_DATA_SIZE+7];
            parallel_state.userdata = (void *)( &userdata[7] - ( ((size_t)&userdata[7])&7 ) ); // pad userdata to 8 bytes
            parallel_state.cam = scene->getCamera();

            MyTree::Tree_node* split_root = nodes_in_splitting_depth[i];

            std::vector<MyTree::Tree_node*> split_nodes;
            split_root->get_post_order_queue(split_nodes, -1);

            //generate_gi_points_data(parallel_state, split_nodes, point_to_surface_point_map);
            generate_gi_points_data_2(split_nodes, point_to_sampling_point_map, area_estimation, needed_radius);
        }

        std::cout << "generate_gi_points(): finished lower tree part" << std::endl;

        std::vector<MyTree::Tree_node*> nodes = tree.get_post_order_queue(splitting_depth);

        generate_gi_points_data_2(nodes, point_to_sampling_point_map, area_estimation, needed_radius, splitting_depth);
    }
    else
    {
        std::vector<MyTree::Tree_node*> nodes = tree.get_post_order_queue();

        // generate_gi_points_data_2(nodes, gi_point_to_sampling_point_map, area_estimation, needed_radius);
        generate_gi_points_data_2(nodes, point_to_sampling_point_map, area_estimation, needed_radius);

    }
#endif

    timer.stop("node_data_creation");

    // ------------------------------------------

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
pbLighting_t::MyTree load_gi_points()
{
#ifdef USE_FAT_TREE
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

    pbLighting_t::MyTree tree(vector3d_t(min), vector3d_t(max), 30, 1);

    for (unsigned int i = 0; i < points.size(); ++i)
    {
        tree.addPoint(points[i]->pos, points[i]);
    }

    return tree;
#else
    return pbLighting_t::MyTree();
#endif
}


std::vector< Spherical_function<color_t> *> generate_dictionary_from_tree(pbLighting_t::MyTree const& tree,
                                                                          Dictionary_generator const* dictionary_generator,
                                                                          Spherical_function_factory<color_t> const* spherical_function_factory,
                                                                          int const dict_num_centers)
{
    std::vector<pbLighting_t::MyTree::Tree_node const*> nodes = tree.get_nodes();

    std::vector<Word> word_list;

    for (unsigned int i = 0; i < nodes.size(); ++i)
    {
#ifndef SURFEL_HAS_COLOR_SF
        if (nodes[i]->is_leaf()) continue;
#endif
        Gi_point_inner_node const& gi_point = nodes[i]->get_data();

        Word word = gi_point.sf_representation.color->to_vector();
        word_list.push_back(word);
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



void push_light(renderState_t & state, pbLighting_t::MyTree::Tree_node * node, Splat_cube_raster_buffer const& fb_in)
{
    // add light to the node
    node->get_data().sf_representation.color->add_light(fb_in);

    if (node->is_leaf()) return;

    int const local_raster_buffer_resolution = fb_in.get_resolution();

    Cube_raster_buffer<float> fb_area_0;
    fb_area_0.setup_surfel_buffer(local_raster_buffer_resolution);

    Cube_raster_buffer<float> fb_area_1;
    fb_area_1.setup_surfel_buffer(local_raster_buffer_resolution);

    Splat_cube_raster_buffer fb_in_0;
    fb_in_0.setup_surfel_buffer(local_raster_buffer_resolution);

    Splat_cube_raster_buffer fb_in_1;
    fb_in_1.setup_surfel_buffer(local_raster_buffer_resolution);

    pbLighting_t::MyTree::Tree_node* child_0 = node->get_child(0);
    pbLighting_t::MyTree::Tree_node* child_1 = node->get_child(1);

    Gi_point_inner_node const&  child_0_data = child_0->get_data();
    Gi_point_inner_node const&  child_1_data = child_1->get_data();

    std::vector<Cube_cell> const& cells = fb_area_0.get_cube_cells();

    for (std::size_t i = 0; i < cells.size(); ++i)
    {
        float const area_0 = child_0_data.sf_representation.area->get_value(fb_area_0.get_cell_direction(cells[i]));
        float const area_1 = child_1_data.sf_representation.area->get_value(fb_area_1.get_cell_direction(cells[i]));

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
    virtual void visit(Spherical_harmonics<Data> * sf) {}
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


float dictionary_mean_absolute_error(pbLighting_t::MyTree const& tree)
{
    std::vector<pbLighting_t::MyTree::Tree_node const*> nodes = tree.get_nodes();

    float error = 0.0f;

    for (std::size_t i = 0; i < nodes.size(); ++i)
    {
        Gi_point_inner_node const& gi_point = nodes[i]->get_data();

        Indexed_spherical_function<color_t> * indexed_sf = dynamic_cast<Indexed_spherical_function<color_t> *>(gi_point.sf_representation.color);

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

    if (_variational) ss << " (V)";

    ss << ", "      << extract_initials(get_string_from_type_map(raster_buffer_type, Splat_cube_raster_buffer::enum_type_map)) <<
          ", "      << extract_initials(get_string_from_type_map(node_splat_type, Splat_cube_raster_buffer::enum_splat_type_map)) <<
          ", "      << extract_initials(get_string_from_type_map(surfel_far_splat_type, Splat_cube_raster_buffer::enum_splat_type_map)) <<
          ", "      << extract_initials(get_string_from_type_map(surfel_near_splat_type, Splat_cube_raster_buffer::enum_splat_type_map)) <<
          " | ds: " << _disc_scale_factor <<
          " | sa: " << _solid_angle_factor <<
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
    _log_file << "Solid angle factor: " << _solid_angle_factor << std::endl;
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

    if (_use_precalculated_sf)
    {
        _precalculated_sf = precalculate_spherical_harmonics(_spherical_function_area_factory);
    }

    if (_color_dictionary_generator)
    {
        MyTree dict_tree;
#ifdef USE_FAT_TREE
        bound_t const& sceneBound = scene->getSceneBound();
        dict_tree = MyTree(vector3d_t(sceneBound.a), vector3d_t(sceneBound.g), 100, 1);
#endif
        generate_gi_points(state, dict_tree, surfel_samples * _dictionary_sample_fraction);

        timer_t my_timer;
        my_timer.addEvent("Dictionary_creation");
        my_timer.start("Dictionary_creation");

        _color_dictionary = generate_dictionary_from_tree(dict_tree, _color_dictionary_generator, _spherical_function_color_factory, _dict_num_centers);

        my_timer.stop("Dictionary_creation");

        _log_file << "Dictionary creation time: " << my_timer.getTime("Dictionary_creation");

        _color_converter = new Spherical_function_indexed_converter<color_t>(&_color_dictionary, _use_ann, _do_dictionary_stats);

#ifdef USE_FAT_TREE
        _point_tree = MyTree(vector3d_t(sceneBound.a), vector3d_t(sceneBound.g), 100, 1);
#endif

        generate_gi_points(state, _point_tree, surfel_samples);

        if (_do_dictionary_stats)
        {
            float const ma_error = dictionary_mean_absolute_error(_point_tree);

            _log_file << "Dictionary mean absolute error: " << ma_error << std::endl;

            _log_file << "Dictionary stats: " << _color_dictionary_generator->get_stat_results() << std::endl;

            std::cout << "_dict_num_centers: " << _dict_num_centers << " " <<
                         "_dictionary_sample_fraction: " << _dictionary_sample_fraction << " " <<
                         "ma_error: " << ma_error << std::endl;
        }
    }
    else
    {
        if (do_load_gi_points)
        {
            _point_tree = load_gi_points();
        }
        else
        {
#ifdef USE_FAT_TREE
            bound_t const& sceneBound = scene->getSceneBound();
            _point_tree = MyTree(vector3d_t(sceneBound.a), vector3d_t(sceneBound.g), 100, 1);
#endif
            generate_gi_points(state, _point_tree, surfel_samples);

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
           << _solid_angle_factor << " "
           << debugTreeDepth << " "
           //<< "max depth: " << _bspTree->getMaxDepth() << " "
           //<< "min leaf depth: " << _bspTree->getMinLeafDepth() << " "
           << settings << " "
           << std::endl;

    _log_file.close();

    return success;
}


void pbLighting_t::cleanup()
{
    _point_tree.clear();

    delete _spherical_function_color_factory;
    delete _spherical_function_area_factory;

    std::vector<Cube_cell> const& cells = _precalculated_sf.get_cube_cells();

    for (size_t i = 0; i < cells.size(); ++i)
    {
        Cube_cell const& c = cells[i];
        Spherical_function<float> * sf = _precalculated_sf.get_data(c);
        delete sf;
    }

    delete _color_dictionary_generator;
    delete _color_converter;
    delete _area_converter;

    for (size_t i = 0; i < _color_dictionary.size(); ++i)
    {
        delete _color_dictionary[i];
    }

    _color_dictionary.clear();

    for (size_t i = 0; i < _area_dictionary.size(); ++i)
    {
        delete _area_dictionary[i];
    }

    _area_dictionary.clear();
}


void process_surfel(
        Gi_point_surfel const& surfel,
        surfacePoint_t const& rp,
        Cube_raster_buffer<color_t> & frame_buffer,
        float const surfel_near_threshold,
        float const mix_amount,
        float const disc_scale_factor,
        bool const handle_rp_as_surface,
        std::vector<Gi_point_info> & point_infos,
        Debug_info * debug_info)
{
    yafaray::vector3d_t surfel_to_rp = (vector3d_t(rp.P) - surfel.pos);

    float const distance = surfel_to_rp.length();

    surfel_to_rp.normalize();

    if (handle_rp_as_surface)
    {
        float const cos_rp_n_gip = rp.N * (-surfel_to_rp);
        if (cos_rp_n_gip < -0.1f) return;
    }

    // float const cos_sp_gip = std::max(gi_point.normal * giToSp, 0.0f);
    float cos_rp_surfel = surfel.normal * surfel_to_rp;
    bool const back_facing = cos_rp_surfel < 0.0f;

    if (back_facing) return;

    cos_rp_surfel = std::abs(cos_rp_surfel);


    // float const cos_sp_gip = std::abs(giP.normal * giToSp);

    float const visible_area = cos_rp_surfel * surfel.area;
    float const max_area = surfel.area;

    float const solid_angle_real = visible_area / (distance * distance); // std::max(0.001f, distance * distance);

    float const disc_radius = std::sqrt(max_area / M_PI);
    // float const visible_radius = std::sqrt(visible_area / M_PI);

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
#ifndef SURFEL_HAS_COLOR_SF
            contribution = surfel.color;
#else
            contribution = gi_point_surfel.color->get_value(surfel_to_rp) / visible_area;
#endif
        }
    }

    Gi_point_info point_info;

    point_info.type = Gi_point_info::Far_surfel;
    if (distance < disc_radius * surfel_near_threshold)
    {
        point_info.type = Gi_point_info::Near_surfel;
    }

    point_info.color             = contribution;
    point_info.disc_normal       = surfel.normal;
    point_info.direction         = -surfel_to_rp;
    point_info.radius            = disc_radius * disc_scale_factor;
    point_info.depth             = distance;
    point_info.position          = surfel.pos;
    point_info.receiver_position = vector3d_t(rp.P);
    point_info.solid_angle       = solid_angle_real * disc_scale_factor * disc_scale_factor;
    point_info.weight            = mix_amount;
    point_info.gi_point          = &surfel;
#ifndef SURFEL_HAS_COLOR_SF
    point_info.spherical_function = NULL;
#else
    Spherical_node_representation * snr = new Spherical_node_representation; // FIXME: Memory leak!
    snr->color = gi_point_surfel.color;
    snr->area  = NULL;
    point_info.spherical_function = snr;
#endif
    // point_info.spherical_function = &gi_point.sf_representation; // FIXME: maybe need to readd?

    if (debug_info)
    {
        if (point_info.type == Gi_point_info::Far_surfel)
        {
            ++debug_info->used_far_surfels;
        }
        else
        {
            ++debug_info->used_near_surfels;
        }
    }

    point_infos.push_back(point_info);
}



void process_node(
        Gi_point_inner_node const* gi_point,
        surfacePoint_t const& receiving_point,
        float const distance,
        vector3d_t const& node_to_rp,
        float const node_radius,
        float const mix_amount,
        float const disc_scale_factor,
        std::vector<Gi_point_info> & point_infos,
        Debug_info * debug_info)
{
    vector3d_t const& position = gi_point->pos;

    float const visible_area = std::max(0.0f, gi_point->sf_representation.area->get_value(node_to_rp));

    // if (visible_area < 0.001f) return;

    float const real_solid_angle = visible_area / (distance * distance);

    yafaray::color_t cluster_contribution;

    if (debug_info && debug_info->color_by_depth)
    {
        // cluster_contribution = pbLighting_t::debug_colors[std::min(node->get_shortest_distance_to_leaf(), int(pbLighting_t::debug_colors.size()) - 1)];
        cluster_contribution = pbLighting_t::debug_colors[std::min(debug_info->node_depth, int(pbLighting_t::debug_colors.size()) - 1)];
    }
    else
    {
        cluster_contribution = gi_point->sf_representation.color->get_value(node_to_rp);
        cluster_contribution *= 1.0f / visible_area;
    }

    Gi_point_info point_info;
    point_info.type               = Gi_point_info::Node;
    point_info.color              = cluster_contribution;
    point_info.direction          = -node_to_rp;
    point_info.depth              = distance;
    point_info.position           = position;
    point_info.receiver_position  = vector3d_t(receiving_point.P);
    // point_info.solid_angle        = real_solid_angle * 2.0f;
    point_info.solid_angle        = real_solid_angle * disc_scale_factor * disc_scale_factor;
    point_info.radius             = node_radius * disc_scale_factor;
    point_info.weight             = mix_amount;
    // point_info.radius             = std::sqrt(visible_area / M_PI);
    point_info.spherical_function = &gi_point->sf_representation;
    point_info.gi_point           = gi_point;

#ifdef PBGI_DEBUG
    point_info.visible_area       = visible_area;
#endif


    if (debug_info)
    {
        ++debug_info->used_nodes;
    }

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
    Gi_point_base* debug_point;
    Gi_point_info* point_info;

    bool operator() (Point_and_debug_info const& p1, Point_and_debug_info const& p2)
    {
        return (p1.point_info->depth < p2.point_info->depth);
    }
};




// #ifdef USE_FAT_TREE
#if 0
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
    typedef pbLighting_t::My_small_tree::Tree_node Node;
    //typedef pbLighting_t::MyTree Node;

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
//    std::queue<pbLighting_t::MyTree const*> queue;
//    queue.push(tree);

    std::queue<Node const*> queue;
    queue.push(tree.get_root());

    while (!queue.empty())
    {
        // pbLighting_t::MyTree const* node = queue.front();
        Node const* node = queue.front();

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
#else
color_t doPointBasedGiTree_sh_fb(
        pbLighting_t::MyTree const& tree,
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
    typedef pbLighting_t::MyTree::Tree_node Node;

    Splat_cube_raster_buffer frame_buffer;

    if (parameters)
    {
        frame_buffer.setup(*parameters->get_child("receiving_fb"));
    }
    else
    {
        frame_buffer.setup(raster_buffer_type, raster_buffer_resolution, node_splat_type, surfel_far_splat_type, surfel_near_splat_type);
    }

    std::vector<Gi_point_info> point_infos;

    point_infos.reserve(10000);

    float const fixed_max_solid_angle = frame_buffer.get_solid_angle(vector3d_t(1, 0, 0)) * solid_angle_factor;

    if (debug_info)
    {
        debug_info->my_timer.start("Traversal");
    }

    float const culling_bias = parameters->get_parameter("culling_bias")->get_value<float>();
    bool const always_use_nodes = parameters->get_parameter("always_use_nodes")->get_value<bool>();
    bool const use_background = parameters->get_parameter("use_background")->get_value<bool>();

    // traverse tree, if solid angle of node > max, traverse into the children
//    std::queue<pbLighting_t::MyTree const*> queue;
//    queue.push(tree);

    std::queue<Node const*> queue;
#ifdef USE_FAT_TREE
    queue.push(&tree);
#else
    queue.push(tree.get_root());
#endif


    while (!queue.empty())
    {
        // pbLighting_t::MyTree const* node = queue.front();
        Node const* node = queue.front();

        queue.pop();

        if (debug_info)
        {
            debug_info->node_depth  = node->get_depth();
            // debug_info->node_height = node->get_node_height(0); // this  is really expensive for every node!
        }

        // if we are here and we have a leaf, we sample all points in the leaf
        if (node->is_leaf())
        {
            if (debug_info)
            {
                debug_info->my_timer.start("Surfel");
            }

            pbLighting_t::MyTree::Tree_leaf_node const* leaf_node = node->get_derived<pbLighting_t::MyTree::Tree_leaf_node>();

            if (always_use_nodes)
            {
                Gi_point_inner_node const& gi_point = leaf_node->get_data();

                float const node_radius = gi_point.bounding_box.get_enclosing_radius();

                vector3d_t node_to_rp = (vector3d_t(receiving_point.P) - gi_point.pos);

                float const distance = node_to_rp.length();

                node_to_rp.normalize();

                process_node(&gi_point, receiving_point, distance, node_to_rp, node_radius, 1.0f, disc_scale_factor, point_infos, debug_info);
            }
            else
            {
                std::vector<Gi_point_surfel> const& surfels = leaf_node->get_surfels();

                for (size_t i = 0; i < surfels.size(); ++i)
                {
                    Gi_point_surfel const& surfel = surfels[i];
                    process_surfel(surfel, receiving_point, frame_buffer, surfel_near_threshold, 1.0f, disc_scale_factor, true, point_infos, debug_info);
                }
            }

            if (debug_info)
            {
                debug_info->my_timer.stop("Surfel");
            }
        }
        else
        {
            if (debug_info)
            {
                debug_info->my_timer.start("NodeCheck");
                ++debug_info->num_node_checks;
            }

            Gi_point_inner_node const& gi_point = node->get_data();

            // if (is_node_data_behind_plane(node, vector3d_t(receiving_point.P), receiving_point.N, 0.1f))
            if (is_bound_behind_plane(gi_point.bounding_box, vector3d_t(receiving_point.P), receiving_point.N, culling_bias))
            // if (node->is_node_data_behind_plane_approx(vector3d_t(receiving_point.P), receiving_point.N, bb_radius, 0.2f))
            {
                if (debug_info)
                {
                    debug_info->my_timer.stop("NodeCheck");
                }

                continue;
            }

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

                pbLighting_t::MyTree::Tree_inner_node const* inner_node = node->get_derived<pbLighting_t::MyTree::Tree_inner_node>();
                std::vector<pbLighting_t::MyTree::Tree_node*> const& children = inner_node->get_children();

                for (size_t i = 0; i < children.size(); ++i)
                {
                    queue.push(children[i]);
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
                    float const parent_radius = node->get_parent()->get_data().bounding_box.get_enclosing_radius();
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

                    process_node(&gi_point, receiving_point, distance, node_to_rp, node_radius, 1.0f - mix_amount, disc_scale_factor, point_infos, debug_info);

                    pbLighting_t::MyTree::Tree_inner_node const* inner_node = node->get_derived<pbLighting_t::MyTree::Tree_inner_node>();
                    std::vector<Node*> const& children = inner_node->get_children();

                    for (size_t i = 0; i < children.size(); ++i)
                    {
                        Node const* child_node = children[i];

                        if (!child_node) continue;

                        Gi_point_inner_node const& child_gi_point = child_node->get_data();
                        float const child_node_radius = child_gi_point.bounding_box.get_enclosing_radius();

                        if (child_node->is_leaf()) // && distance < child_node_radius * surfel_near_threshold)
                        {
                            pbLighting_t::MyTree::Tree_leaf_node const* child_leaf_node = child_node->get_derived<pbLighting_t::MyTree::Tree_leaf_node>();
                            std::vector<Gi_point_surfel> const& child_surfels = child_leaf_node->get_surfels();

                            for (size_t i = 0; i < child_surfels.size(); ++i)
                            {
                                Gi_point_surfel const& child_surfel = child_surfels[i];
                                process_surfel(child_surfel, receiving_point, frame_buffer, surfel_near_threshold, mix_amount, disc_scale_factor, true, point_infos, debug_info);
                            }



//                            std::vector<Gi_point_base*> const* child_surfels = child_node->get_surfels();

//                            for (size_t i = 0; i < child_surfels->size(); ++i)
//                            {
//                                Gi_point_surfel const* child_surfel = (*child_surfels)[i]->get_derived<Gi_point_surfel>();
//                                process_surfel(child_surfel, receiving_point, frame_buffer, surfel_near_threshold, mix_amount, disc_scale_factor, true, point_infos, debug_info);
//                            }
                        }
                        else
                        {
//                            Gi_point_inner_node const* child_gi_point = child_node->get_data()->get_derived<Gi_point_inner_node>();
//                            float const child_node_radius = child_gi_point->bounding_box.get_enclosing_radius();

                            vector3d_t const& child_position = child_gi_point.pos;

                            vector3d_t child_node_to_rp = (vector3d_t(receiving_point.P) - child_position);

                            float const child_distance = child_node_to_rp.length();
                            child_node_to_rp.normalize();

                            process_node(&child_gi_point, receiving_point, child_distance, child_node_to_rp, child_node_radius, mix_amount, disc_scale_factor, point_infos, debug_info);
                        }
                    }
                }
                else // non-variational
                {
                    process_node(&gi_point, receiving_point, distance, node_to_rp, node_radius, 1.0f, disc_scale_factor, point_infos, debug_info);
                }


                if (debug_info)
                {
                    debug_info->my_timer.stop("Node");
                }
            } // end else (bad solid angle)
        } // inner node handling
    } // queue empty

    if (debug_info)
    {
        debug_info->my_timer.stop("Traversal");
        debug_info->my_timer.start("Accumulating");
    }


    // sort and splat, then integrate
    std::sort(point_infos.begin(), point_infos.end(), Compare_point_info_by_depth());

    for (unsigned int i = 0; i < point_infos.size(); ++i)
    {
        frame_buffer.add_point(point_infos[i]);
    }

    if (use_background)
    {
        frame_buffer.add_background(background);
    }


    color_t col(0.0f);

    if (receiving_point.material)
    {
        col = frame_buffer.get_brdf_response(state, receiving_point, wo);
    }
    else
    {
        col = frame_buffer.get_diffuse(receiving_point.N);
    }

    if (debug_info)
    {
        debug_info->gi_point_infos = point_infos;

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
#endif


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

        if (!indirectOnly)
        {
            if(bsdfs & BSDF_EMIT)
            {
                col += material->emission(state, sp, wo);
            }

            col += estimateAllDirectLight(state, sp, wo);
        }

        CALLGRIND_START_INSTRUMENTATION;

        col += doPointBasedGiTree_sh_fb(
                    _point_tree,
                    state,
                    sp,
                    _solid_angle_factor,
                    background,
                    wo,
                    raster_buffer_resolution,
                    raster_buffer_type,
                    node_splat_type,
                    surfel_far_splat_type,
                    surfel_near_splat_type,
                    surfel_near_threshold,
                    _variational,
                    _disc_scale_factor,
                    &_parameters
                    );

        CALLGRIND_STOP_INSTRUMENTATION;

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

    inte->do_load_gi_points = do_load_gi_points;

    inte->_variational         = variational;
    inte->_parameters.add_parameter(new Parameter("variational", variational));

    inte->_solid_angle_factor = maxSolidAngle;
    inte->_parameters.add_parameter(new Parameter("solid_angle_factor", maxSolidAngle, 0.0f, 256.0f));

    inte->_disc_scale_factor = disc_scale_factor;
    inte->_parameters.add_parameter(new Parameter("disc_scale_factor", disc_scale_factor, 0.01f, 20.0f));

    // Parameter_list parameters;

    Parameter_list * l;

    inte->_parameters.add_parameter(new Parameter("culling_bias", 0.1f, 0.0f, 0.5f));
    inte->_parameters.add_parameter(new Parameter("always_use_nodes", false));
    inte->_parameters.add_parameter(new Parameter("use_background", false));

    Parameter_list * receiving_fb_list = inte->_parameters.add_child("receiving_fb");
    Parameter_registry< Abstract_frame_buffer<color_t> >::create_single_select_instance(receiving_fb_list, "fb_type");
    receiving_fb_list->get_child("fb_type")->get_parameter("type")->set_value(fb_type);

    l = receiving_fb_list->get_child("fb_type")->get_child("Accumulating_frame_buffer_without_queue");
    (*l)["resolution"]->set_value(fb_resolution);
    (*l)["use_visibility"]->set_value(true);
    (*l)["use_depth_modulation"]->set_value(false);
    (*l)["normalize"]->set_value(false);

    Parameter_registry<Splat_strategy>::create_single_select_instance(receiving_fb_list, "node_splat_type");
    Parameter_registry<Splat_strategy>::create_single_select_instance(receiving_fb_list, "surfel_far_splat_type");
    Parameter_registry<Splat_strategy>::create_single_select_instance(receiving_fb_list, "surfel_near_splat_type");

    receiving_fb_list->get_child("node_splat_type")->get_parameter("type")->set_value(node_splat_type);
    receiving_fb_list->get_child("surfel_far_splat_type")->get_parameter("type")->set_value(surfel_far_splat_type);
    receiving_fb_list->get_child("surfel_near_splat_type")->get_parameter("type")->set_value(surfel_near_splat_type);

    l = receiving_fb_list->get_child("node_splat_type")->get_child("Gaussian_splat_strategy");
    (*l)["wendland_integral"]->set_value(false);

    l = receiving_fb_list->get_child("surfel_far_splat_type")->get_child("Gaussian_splat_strategy");
    (*l)["wendland_integral"]->set_value(false);

    inte->raster_buffer_resolution = fb_resolution;
    inte->raster_buffer_type       = Splat_cube_raster_buffer::enum_type_map[fb_type];
    inte->node_splat_type          = Splat_cube_raster_buffer::enum_splat_type_map[node_splat_type];
    inte->surfel_far_splat_type    = Splat_cube_raster_buffer::enum_splat_type_map[surfel_far_splat_type];
    inte->surfel_near_splat_type   = Splat_cube_raster_buffer::enum_splat_type_map[surfel_near_splat_type];
    inte->surfel_near_threshold    = surfel_near_threshold;

    bool const use_exact_sh = true;

    std::string spherical_function_factory_type = "SH";
    params.getParam("spherical_function_type", spherical_function_factory_type);
    if (spherical_function_factory_type == "Cube")
    {
        inte->_spherical_function_color_factory = new Cube_spherical_function_factory<color_t>(sf_resolution, false);
        inte->_use_precalculated_sf       = false;
    }
    else if (spherical_function_factory_type == "SH")
    {
        inte->_spherical_function_color_factory = new Spherical_harmonics_factory<color_t>(3, use_exact_sh);
        inte->_use_precalculated_sf       = true;
    }

    // inte->_spherical_function_area_factory = new Cube_spherical_function_factory<float>(4, false);
    inte->_spherical_function_area_factory = new Spherical_harmonics_factory<float>(3, use_exact_sh);

    inte->_area_converter  = NULL;
    inte->_color_converter = NULL;

    if (enable_conversion)
    {
        inte->_color_converter = new Cube_to_mises_fisher_converter<color_t>(1);
        inte->_parameters["always_use_nodes"]->set_value(true);
    }

    inte->_sf_resolution              = sf_resolution;
    inte->_use_sf_files               = use_sf_files;

    // dictionary settings
    inte->_dict_num_centers           = dict_num_centers;
    inte->_dictionary_sample_fraction = dictionary_sample_fraction;
    inte->_do_dictionary_stats        = do_dictionary_stats;
    inte->_use_ann                    = true;

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
        inte->_color_dictionary_generator = new Kmeans_dictionary_generator(inte->_use_ann, true);
    }


    Y_INFO <<
              "maxSolidAngle: " << inte->_solid_angle_factor << " " <<
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
