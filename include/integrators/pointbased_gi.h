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

__BEGIN_YAFRAY

struct GiPoint
{
    GiPoint() :
        pos(0.0f),
        normal(0.0f),
        color(0.0f),
        area(0.0f),
        energy(0.0f),
        sh_representation(GiSphericalHarmonics<vector3d_t, color_t>(true, 3))
    { }

    vector3d_t pos;
    vector3d_t normal;
    color_t    color;
    float      area;
    color_t    energy;

    GiSphericalHarmonics<vector3d_t, color_t> sh_representation;

    friend std::ostream & operator<<(std::ostream & s, GiPoint const& p)
    {
        s <<
             p.pos.x << " " << p.pos.y << " " << p.pos.z << " " <<
             p.normal.x << " " << p.normal.y << " " << p.normal.z << " " <<
             p.color.R << " " << p.color.G << " " << p.color.B << " " <<
             p.area << " " <<
             p.energy.R << " " << p.energy.G << " " << p.energy.B << " " <<
             p.sh_representation
             ;

        return s;
    }

    friend std::istream & operator>>(std::istream & s, GiPoint & p)
    {
        s >>
             p.pos.x >> p.pos.y >> p.pos.z >>
             p.normal.x >> p.normal.y >> p.normal.z >>
             p.color.R >> p.color.G >> p.color.B >>
             p.area >>
             p.energy.R >> p.energy.G >> p.energy.B >>
             p.sh_representation
             ;

        return s;
    }
};


class YAFRAYPLUGIN_EXPORT pbLighting_t: public mcIntegrator_t
{
public:
    typedef RegularBspTree<vector3d_t, 3, GiPoint> MyTree;

    pbLighting_t(bool transpShad=false, int shadowDepth=4, int rayDepth=6);
    virtual bool preprocess();
    virtual colorA_t integrate(renderState_t &state, diffRay_t &ray) const;
    static integrator_t* factory(paraMap_t &params, renderEnvironment_t &render);

    color_t estimateIncomingLight(renderState_t & state, light_t *light, const surfacePoint_t &sp, const unsigned int &loffs) const;
    color_t doPointBasedGi(renderState_t & state, surfacePoint_t const& sp, vector3d_t const& wo) const;
    color_t doPointBasedGiTree(renderState_t & state, surfacePoint_t const& sp, vector3d_t const& wo) const;
    color_t doPointBasedGiTreeSH(renderState_t & state, surfacePoint_t const& sp, vector3d_t const& wo) const;
    color_t doPointBasedGiTreeSH_leafs_only(renderState_t & state, surfacePoint_t const& sp, vector3d_t const& wo) const;

private:
    enum Debug_type { NoTree, Tree, Tree_sh, Tree_sh_leafs };

    std::vector<GiPoint> giPoints;
    int samplesPerArea;
    bool debug;
    bool indirectOnly;
    float maxSolidAngle;
    int debugTreeDepth;
    bool debugOutputPointsToFile;
    Debug_type debug_type;

    bool render_single_pixel;
    int pixel_x, pixel_y;

    MyTree* _bspTree;
    std::vector<MyTree *> leafNodes;
};

__END_YAFRAY

#endif // POINTBASED_GI_H
