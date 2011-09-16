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
#include <utilities/CubeRasterBuffer.h>
#include <utilities/Debug_info.h>
#include <utilities/SerializationHelper.h>

__BEGIN_YAFRAY

struct GiPoint
{
    GiPoint() :
        pos(0.0f),
        normal(0.0f),
        color(0.0f),
        area(0.0f),
        energy(0.0f),
        depth(-1),
        is_surfel(true)
    {
        // sh_representation = new GiSphericalHarmonics<vector3d_t, color_t>(true, 3);
        sh_representation = new Cube_spherical_function();
    }

    vector3d_t pos;
    vector3d_t normal;
    color_t    color;
    float      area;
    color_t    energy;
    int        depth;

    bool is_surfel;
    mutable float debug_radius;

    //GiSphericalHarmonics<vector3d_t, color_t> sh_representation;
    Spherical_function * sh_representation;

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
};



class YAFRAYPLUGIN_EXPORT pbLighting_t: public mcIntegrator_t
{
public:
    static std::vector<color_t> debug_colors;

    typedef RegularBspTree<vector3d_t, 3, GiPoint*> MyTree;

    pbLighting_t(bool transpShad=false, int shadowDepth=4, int rayDepth=6);
    virtual bool preprocess();
    virtual colorA_t integrate(renderState_t &state, diffRay_t &ray) const;
    static integrator_t* factory(paraMap_t &params, renderEnvironment_t &render);

    color_t estimateIncomingLight(renderState_t & state, light_t *light, const surfacePoint_t &sp, const unsigned int &loffs) const;
    color_t doPointBasedGiTree(renderState_t & state, surfacePoint_t const& sp, vector3d_t const& wo) const;
    color_t doPointBasedGiTreeSH(renderState_t & state, surfacePoint_t const& sp, vector3d_t const& wo) const;
    // color_t doPointBasedGiTree_sh_fb(renderState_t & state, surfacePoint_t const& sp, vector3d_t const& wo) const;
    color_t doPointBasedGiTreeSH_leafs_only(renderState_t & state, surfacePoint_t const& sp, vector3d_t const& wo) const;

    void generate_gi_points(renderState_t & state);

    float get_max_solid_angle() { return maxSolidAngle; }
    int get_raster_buffer_resolution() const { return raster_buffer_resolution; }
    Cube_raster_buffer::Type get_raster_buffer_type() const { return raster_buffer_type; }

    float get_surfel_near_threshold() const { return surfel_near_threshold; }
    Cube_raster_buffer::Splat_type get_node_splat_type() const { return node_splat_type; }
    Cube_raster_buffer::Splat_type get_surfel_far_splat_type() const { return surfel_far_splat_type; }
    Cube_raster_buffer::Splat_type get_surfel_near_splat_type() const { return surfel_near_splat_type; }

    MyTree* get_tree() { return _bspTree; }

    void set_load_gi_points(bool const b) { do_load_gi_points = b; }

private:
    enum Debug_type { NoTree, Tree, Tree_sh, Tree_sh_fb, Tree_sh_leafs };

    // std::vector<GiPoint> giPoints;
    int samplesPerArea;
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

    MyTree* _bspTree;
};


pbLighting_t::MyTree* load_gi_points();
GiPoint * averageGiPoints(std::vector<GiPoint*> const& points);
color_t doPointBasedGiTree_sh_fb(
    pbLighting_t::MyTree const* tree, renderState_t & state,
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

void process_surfel(
    GiPoint const& gi_point,
    surfacePoint_t const& sp,
    Cube_raster_buffer & frame_buffer,
    float const surfel_near_threshold,
    Debug_info * debug_info = NULL);

__END_YAFRAY

#endif // POINTBASED_GI_H
