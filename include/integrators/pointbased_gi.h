#ifndef POINTBASED_GI_H
#define POINTBASED_GI_H

#include <yafray_config.h>
#include <core_api/environment.h>
#include <core_api/material.h>
#include <core_api/mcintegrator.h>
#include <core_api/background.h>
#include <core_api/light.h>
#include <yafraycore/meshtypes.h>
#include <utilities/mcqmc.h>
#include <utilities/spherical_harmonics.h>

#include <utilities/RegularBspTree.h>
#include <utilities/PointKdTree.h>
#include <utilities/CubeRasterBuffer.h>
#include <utilities/Debug_info.h>
#include <utilities/SerializationHelper.h>

#include <integrators/Dictionary_converter.h>

#include "distance.hpp"

__BEGIN_YAFRAY

struct GiPoint
{
    GiPoint(Spherical_function_factory const* sf_factory) :
        pos(0.0f),
        normal(0.0f),
        color(0.0f),
        area(0.0f),
        // energy(0.0f),
        depth(-1),
        bounding_box(point3d_t(0.0f), point3d_t(0.0f)),
        is_surfel(true)
    {
        if (sf_factory)
        {
            sh_representation = sf_factory->create();
        }
    }

    vector3d_t pos;
    vector3d_t normal;
    color_t    color;
    float      area;
    int        depth;
    bound_t    bounding_box;

    bool is_surfel;
    mutable float debug_radius;

    Spherical_function * sh_representation;

    /*
    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & pos;
        ar & normal;
        ar & color;
        ar & area;
        ar & energy;
        ar & depth;
        ar & is_surfel;
        ar & sh_representation;
    }
    */
};


class Gi_point_averager
{
    public:
    Gi_point_averager(Spherical_function_factory const* spherical_function_factory) :
        _spherical_function_factory(spherical_function_factory)
    { }

    GiPoint * average(std::vector<GiPoint*> const& points)
    {
        assert(points.size() == 2 || points.size() == 1);

        GiPoint * result = new GiPoint(_spherical_function_factory);

        if (points.size() == 0) return result;

        // SH summation
        for (unsigned int i = 0; i < points.size(); ++i)
        {
            GiPoint const& p = *points[i];

            result->sh_representation->add(p.sh_representation);
        }

        result->sh_representation->normalize(1.0f / float(points.size()));

        // -------------------- test

        vector3d_t dir(0.4f, 0.3f, 0.5f);
        dir.normalize();

        float area_sum = 0;

        for (unsigned int i = 0; i < points.size(); ++i)
        {
            GiPoint const& p = *points[i];

            area_sum += p.sh_representation->get_area(dir);
        }

        if (!is_in_range(0.9f, 1.1f, result->sh_representation->get_area(dir) / area_sum))
        {
            std::cout << "surfel?: " << points[0]->is_surfel << ", not in range: " << area_sum << " " << result->sh_representation->get_area(dir) << std::endl;
        }
        /*
        else
        {
            std::cout << "in range: " << area_sum << " " << result->sh_representation->get_area(dir) << std::endl;
        }
        */

        // -------------------------

        for (unsigned int i = 0; i < points.size(); ++i)
        {
            GiPoint const* p = points[i];

            result->pos    += p->pos;
            result->normal += p->normal;
            result->color  += p->color;
            // result->energy += p->energy;
        }

        if (points.size() == 2)
        {
            result->bounding_box = bound_t(points[0]->bounding_box, points[1]->bounding_box);

            assert(
                        result->bounding_box.a.x <= points[0]->bounding_box.a.x &&
                        result->bounding_box.a.x <= points[1]->bounding_box.a.x &&
                        result->bounding_box.a.y <= points[0]->bounding_box.a.y &&
                        result->bounding_box.a.y <= points[1]->bounding_box.a.y &&
                        result->bounding_box.a.z <= points[0]->bounding_box.a.z &&
                        result->bounding_box.a.z <= points[1]->bounding_box.a.z &&
                        result->bounding_box.g.x >= points[0]->bounding_box.g.x &&
                        result->bounding_box.g.x >= points[1]->bounding_box.g.x &&
                        result->bounding_box.g.y >= points[0]->bounding_box.g.y &&
                        result->bounding_box.g.y >= points[1]->bounding_box.g.y &&
                        result->bounding_box.g.z >= points[0]->bounding_box.g.z &&
                        result->bounding_box.g.z >= points[1]->bounding_box.g.z
                        );
        }
        else
        {
            result->bounding_box = points[0]->bounding_box;
        }


        result->pos    *= 1.0f / float(points.size());
        result->color  *= 1.0f / float(points.size());
        // result->energy *= 1.0f / float(points.size());


        result->normal.normalize();

        result->is_surfel = false;

        return result;
    }

    private:
    Spherical_function_factory const* _spherical_function_factory;
};



class Indexed_gi_point_averager
{
    public:
    Indexed_gi_point_averager(Spherical_function_factory const* spherical_function_factory, Dictionary_converter const* converter) :
        _spherical_function_factory(spherical_function_factory),
        _converter(converter)
    { }

    GiPoint * average(std::vector<GiPoint*> const& points)
    {
        //    assert(points.size() > 0);

        GiPoint * result = new GiPoint(_spherical_function_factory);

        if (points.size() == 0) return result;

        // convert the point sfs from indexed to cube sf
        std::vector<Spherical_function*> converted_functions;

        for (unsigned int i = 0; i < points.size(); ++i)
        {
            GiPoint const& p = *points[i];

            Spherical_function * csf = _converter->reconvert(p.sh_representation);
            converted_functions.push_back(csf);
        }

        // SH summation
        for (unsigned int i = 0; i < converted_functions.size(); ++i)
        {
            result->sh_representation->add(converted_functions[i]);
        }

        result->sh_representation->normalize(1.0f / float(points.size()));

        Spherical_function * tmp = result->sh_representation;
        result->sh_representation = _converter->convert(result->sh_representation);
        delete tmp;

        for (unsigned int i = 0; i < converted_functions.size(); ++i)
        {
            delete converted_functions[i];
        }

        // -------------------- test

        vector3d_t dir(0.4f, 0.3f, 0.5f);
        dir.normalize();

        float area_sum = 0;

        for (unsigned int i = 0; i < points.size(); ++i)
        {
            GiPoint const& p = *points[i];

            area_sum += p.sh_representation->get_area(dir);
        }

        if (!is_in_range(0.9f, 1.1f, result->sh_representation->get_area(dir) / area_sum))
        {
            std::cout << "Indexed_gi_point_averager::average(): is surfel: " << points[0]->is_surfel << ", not in range: " << area_sum << " " << result->sh_representation->get_area(dir) << std::endl;
        }
        else
        {
            // std::cout << "in range: " << area_sum << " " << result->sh_representation->get_area(dir) << std::endl;
        }

        // -------------------------

        for (unsigned int i = 0; i < points.size(); ++i)
        {
            GiPoint const* p = points[i];

            result->pos    += p->pos;
            result->normal += p->normal;
            result->color  += p->color;
            // result->energy += p->energy;
        }

        result->pos    *= 1.0f / float(points.size());
        result->color  *= 1.0f / float(points.size());
        // result->energy *= 1.0f / float(points.size());

        result->normal.normalize();

        result->is_surfel = false;

        return result;
    }

private:
    Spherical_function_factory const* _spherical_function_factory;
    Dictionary_converter const* _converter;
};


class YAFRAYPLUGIN_EXPORT pbLighting_t: public mcIntegrator_t
{
public:
    enum Dictionary_type { No_dict, Random_dict, Kmeans_dict };


    static std::vector<color_t> debug_colors;

    // typedef RegularBspTree<vector3d_t, 3, GiPoint*> MyTree;
    typedef Point_kd_tree<vector3d_t, 3, GiPoint*> MyTree;

    pbLighting_t(bool transpShad=false, int shadowDepth=4, int rayDepth=6);
    virtual bool preprocess();
    virtual colorA_t integrate(renderState_t &state, diffRay_t &ray) const;
    static integrator_t* factory(paraMap_t &params, renderEnvironment_t &render);

    color_t estimateIncomingLight(renderState_t & state, light_t *light, const surfacePoint_t &sp, const unsigned int &loffs) const;
    color_t doPointBasedGiTree(renderState_t & state, surfacePoint_t const& sp, vector3d_t const& wo) const;
    color_t doPointBasedGiTreeSH(renderState_t & state, surfacePoint_t const& sp, vector3d_t const& wo) const;
    color_t doPointBasedGiTreeSH_leafs_only(renderState_t & state, surfacePoint_t const& sp, vector3d_t const& wo) const;

    void generate_spherical_function(renderState_t & state, GiPoint * gi_point, surfacePoint_t const& sp, std::vector<light_t*> const& lights);
    void generate_gi_points         (renderState_t & state, MyTree * tree, int const number_of_samples, Spherical_function_converter const* converter = NULL);
    void generate_gi_points_data    (renderState_t & state,
                                     std::vector<MyTree*> const& nodes,
                                     std::map<GiPoint*, surfacePoint_t> const& point_to_surface_point_map,
                                     Spherical_function_converter const* converter = NULL,
                                     int const treat_depth_as_leaf = -1);

    float get_max_solid_angle() { return maxSolidAngle; }
    int get_raster_buffer_resolution() const { return raster_buffer_resolution; }
    Cube_raster_buffer::Type get_raster_buffer_type() const { return raster_buffer_type; }

    float get_surfel_near_threshold() const { return surfel_near_threshold; }
    Cube_raster_buffer::Splat_type get_node_splat_type() const { return node_splat_type; }
    Cube_raster_buffer::Splat_type get_surfel_far_splat_type() const { return surfel_far_splat_type; }
    Cube_raster_buffer::Splat_type get_surfel_near_splat_type() const { return surfel_near_splat_type; }

    MyTree* get_tree() { return _bspTree; }

    Dictionary_type get_dictionary_type() { return _dictionary_type; }

    std::vector<Spherical_function*> const& get_dictionary() { return _dictionary; }

    void set_load_gi_points(bool const b) { do_load_gi_points = b; }

    int get_sf_resolution() const { return _sf_resolution; }

private:
    enum Debug_type { NoTree, Tree, Tree_sh, Tree_sh_fb, Tree_sh_leafs };

    // std::vector<GiPoint> giPoints;
    int surfel_samples;
    bool debug;
    bool indirectOnly;
    int debugTreeDepth;
    bool debugOutputPointsToFile;
    Debug_type debug_type;
    std::string debug_type_str;

    bool render_single_pixel;
    int pixel_x, pixel_y;

    bool do_load_gi_points;

    float maxSolidAngle;
    int raster_buffer_resolution;
    Cube_raster_buffer::Type raster_buffer_type;

    float surfel_near_threshold;
    Cube_raster_buffer::Splat_type node_splat_type;
    Cube_raster_buffer::Splat_type surfel_far_splat_type;
    Cube_raster_buffer::Splat_type surfel_near_splat_type;

    Spherical_function_factory const* _spherical_function_factory;

    bool variational;

    MyTree* _bspTree;

    std::vector<Spherical_function*> _dictionary;

    Dictionary_type _dictionary_type;
    int _dict_num_centers;
    float _dictionary_sample_fraction; // fraction of samples to generate the dictionary from the scene samples

    int _sf_resolution; // spherical function resolution, for cube buffer, used for debug

    bool _use_sf_files;
};


pbLighting_t::MyTree* load_gi_points();
GiPoint * averageGiPoints(std::vector<GiPoint*> const& points);

color_t doPointBasedGiTree_sh_fb(
        pbLighting_t::MyTree const* tree,
        renderState_t & state,
        surfacePoint_t const& sp,
        float const maxSolidAngle,
        vector3d_t const& wo,
        int const raster_buffer_resolution,
        Cube_raster_buffer::Type const raster_buffer_type,
        Cube_raster_buffer::Splat_type const node_splat_type,
        Cube_raster_buffer::Splat_type const surfel_far_splat_type,
        Cube_raster_buffer::Splat_type const surfel_near_splat_type,
        float const surfel_near_threshold,
        Debug_info * debug_info = NULL);

color_t doPointBasedGiTree_sh_fb_variational(
        pbLighting_t::MyTree const* tree,
        renderState_t & state,
        surfacePoint_t const& receiving_point,
        float const solid_angle_threshold,
        vector3d_t const& wo,
        int const raster_buffer_resolution,
        Cube_raster_buffer::Type const raster_buffer_type,
        Cube_raster_buffer::Splat_type const node_splat_type,
        Cube_raster_buffer::Splat_type const surfel_far_splat_type,
        Cube_raster_buffer::Splat_type const surfel_near_splat_type,
        float const surfel_near_threshold,
        Debug_info * debug_info = NULL
        );

void process_surfel(
        GiPoint const& gi_point,
        surfacePoint_t const& sp,
        Cube_raster_buffer & frame_buffer,
        float const surfel_near_threshold,
        float const mix_amount_node_points,
        std::vector<Gi_point_info> & point_infos,
        std::vector<GiPoint*> & debug_points,
        Debug_info * debug_info = NULL);


std::vector<Spherical_function*> generate_dictionary_from_gi_points_kmeans(std::vector<GiPoint const*> const& gi_points, Spherical_function_factory const* spherical_function_factory, const int dict_num_centers);

__END_YAFRAY

#endif // POINTBASED_GI_H
