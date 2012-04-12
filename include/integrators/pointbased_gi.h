#ifndef POINTBASED_GI_H
#define POINTBASED_GI_H

#include <tr1/unordered_map>
#include <fstream>

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
#include <utilities/Dictionary_generator.h>


#include "distance.hpp"

__BEGIN_YAFRAY


struct GiPoint
{
    GiPoint(Spherical_function_factory<color_t> const* sf_color_factory, Spherical_function_factory<float> const* sf_area_factory) :
        pos(0.0f),
        normal(0.0f),
        color(0.0f),
        area(0.0f),
        // energy(0.0f),
        depth(-1),
        bounding_box(point3d_t(0.0f), point3d_t(0.0f)),
        is_surfel(true)
    {
        if (sf_color_factory)
        {
            sf_representation.color = sf_color_factory->create();
        }

        if (sf_area_factory)
        {
            sf_representation.area = sf_area_factory->create();
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
    mutable float debug_mix_amount;

    Spherical_node_representation sf_representation;

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
    Gi_point_averager(Spherical_function_factory<color_t> const* sf_color_factory, Spherical_function_factory<float> const* sf_area_factory) :
        _sf_color_factory(sf_color_factory),
        _sf_area_factory(sf_area_factory)
    { }

    GiPoint * average(std::vector<GiPoint*> const& points)
    {
        assert(points.size() == 2 || points.size() == 1);

        GiPoint * result = new GiPoint(_sf_color_factory, _sf_area_factory);

        if (points.size() == 0) return result;

        // SH summation
        for (unsigned int i = 0; i < points.size(); ++i)
        {
            GiPoint const& p = *points[i];

            result->sf_representation.area->add(p.sf_representation.area);
            result->sf_representation.color->add(p.sf_representation.color);
        }

        result->sf_representation.color->normalize(1.0f / float(points.size()));

        // -------------------- test

        vector3d_t dir(0.4f, 0.3f, 0.5f);
        dir.normalize();

        float area_sum = 0;

        for (unsigned int i = 0; i < points.size(); ++i)
        {
            GiPoint const& p = *points[i];

            area_sum += p.sf_representation.area->get_value(dir);
        }

        if (result->sf_representation.area->get_value(dir) > 0.0f && !is_in_range(0.9f, 1.1f, result->sf_representation.area->get_value(dir) / area_sum))
        {
            std::cout << "Gi_point_averager::average(): surfel?: " << points[0]->is_surfel << ", not in range: " << area_sum << " " << result->sf_representation.area->get_value(dir) << std::endl;
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
            result->area   += p->area;
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

        result->normal.normalize();

        result->is_surfel = (points.size() == 1);

        return result;
    }

private:
    Spherical_function_factory<color_t> const* _sf_color_factory;
    Spherical_function_factory<float> const* _sf_area_factory;
};


struct pbgi_sample_t
{
    triangle_t const* tri_pointer;
    point3d_t position;
    float area;
    intersectData_t intersect_data;
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
    void generate_gi_points         (renderState_t & state, MyTree * tree, int const number_of_samples);
    void generate_gi_points_data    (renderState_t & state,
                                     std::vector<MyTree*> const& nodes,
                                     std::tr1::unordered_map<GiPoint*, surfacePoint_t> const& point_to_surface_point_map,
                                     int const treat_depth_as_leaf = -1);
    void generate_gi_points_data_2(std::vector<MyTree*> const& nodes,
                                   std::tr1::unordered_map<GiPoint*, pbgi_sample_t*> const& point_to_sampling_point_map,
                                   float const area,
                                   float const radius,
                                   int const treat_depth_as_leaf = -1);

    void add_bounce(renderState_t & state, pbLighting_t::MyTree * tree);

    float get_max_solid_angle() { return _solid_angle_factor; }
    int get_raster_buffer_resolution() const { return raster_buffer_resolution; }
    Splat_cube_raster_buffer::Buffer_type get_raster_buffer_type() const { return raster_buffer_type; }

    float get_surfel_near_threshold() const { return surfel_near_threshold; }
    Splat_cube_raster_buffer::Splat_type get_node_splat_type() const { return node_splat_type; }
    Splat_cube_raster_buffer::Splat_type get_surfel_far_splat_type() const { return surfel_far_splat_type; }
    Splat_cube_raster_buffer::Splat_type get_surfel_near_splat_type() const { return surfel_near_splat_type; }

    MyTree* get_tree() { return _bspTree; }

    std::vector<Spherical_function<color_t> *> const& get_dictionary_color() { return _color_dictionary; }

    Dictionary_generator const* get_color_dictionary_generator() { return _color_dictionary_generator; }

    void set_load_gi_points(bool const b) { do_load_gi_points = b; }

    int get_sf_resolution() const { return _sf_resolution; }

    bool get_variational() const { return _variational; }

    Parameter_list const& get_parameters() const { return _parameters; }

    std::vector<vector3d_t> _debug_new_triangles;

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

    float _solid_angle_factor;
    int raster_buffer_resolution;
    Splat_cube_raster_buffer::Buffer_type raster_buffer_type;

    float surfel_near_threshold;
    Splat_cube_raster_buffer::Splat_type node_splat_type;
    Splat_cube_raster_buffer::Splat_type surfel_far_splat_type;
    Splat_cube_raster_buffer::Splat_type surfel_near_splat_type;

    Spherical_function_factory<color_t> const* _spherical_function_color_factory;
    Spherical_function_factory<float>   const* _spherical_function_area_factory;

    Spherical_function_converter<color_t> const* _color_converter;
    Spherical_function_converter<float>   const* _area_converter;

    bool _variational;

    MyTree* _bspTree;

    Cube_raster_buffer<Spherical_function<float> *> _precalculated_sf;

    Dictionary_generator * _color_dictionary_generator;

    std::vector< Spherical_function<color_t> *> _color_dictionary;
    std::vector< Spherical_function<float>   *> _area_dictionary;

    bool _do_dictionary_stats;
    int _dict_num_centers;
    float _dictionary_sample_fraction; // fraction of samples to generate the dictionary from the scene samples
    bool _use_ann; // use ann to accelerate the closest cluster search

    int _sf_resolution; // spherical function resolution, for cube buffer, used for debug

    bool _use_sf_files;

    float _disc_scale_factor;

    std::ofstream _log_file;

    Parameter_list _parameters;
};


pbLighting_t::MyTree* load_gi_points();
GiPoint * averageGiPoints(std::vector<GiPoint*> const& points);

color_t doPointBasedGiTree_sh_fb(
        pbLighting_t::MyTree const* tree,
        renderState_t & state,
        surfacePoint_t const& sp,
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
        Parameter_list const* parameters = NULL,
        Debug_info * debug_info = NULL);

void process_surfel(
        GiPoint const& gi_point,
        surfacePoint_t const& sp,
        Cube_raster_buffer<color_t> & frame_buffer,
        float const surfel_near_threshold,
        float const mix_amount_node_points,
        float const disc_scale_factor,
        bool const handle_rp_as_surface,
        std::vector<Gi_point_info> & point_infos,
        std::vector<GiPoint*> & debug_points,
        Debug_info * debug_info = NULL);

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
        float const surfel_near_threshold);

std::vector<Spherical_function<color_t>*> generate_dictionary_from_gi_points_kmeans(std::vector<GiPoint const*> const& gi_points, Spherical_function_factory<color_t> const* spherical_function_factory, const int dict_num_centers);

__END_YAFRAY

#endif // POINTBASED_GI_H
