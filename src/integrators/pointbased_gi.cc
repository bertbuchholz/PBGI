#include <sstream>
#include <cassert>
#include <algorithm>
#include <fstream>
#include <iomanip>

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
#include <integrators/pointbased_gi.h>

__BEGIN_YAFRAY


GiPoint * averageGiPoints(std::vector<GiPoint*> const& points)
{
//    assert(points.size() > 0);

    GiPoint * result = new GiPoint();

    if (points.size() == 0) return result;

    // SH summation
    for (unsigned int i = 0; i < points.size(); ++i)
    {
        GiPoint const& p = *points[i];

        result->sh_representation = result->sh_representation + p.sh_representation;

        result->pos += p.pos;
    }

    result->sh_representation = result->sh_representation / float(points.size());

    result->pos *= 1.0f / float(points.size());

    // sanity check
    /*
    yafaray::vector3d_t mainDir(0.5f, 0.5f, 0.5f);
    mainDir.normalize();

    float accAreaChildren = 0;
    float accAreaParent = result.sh_representation.get_sh_area(mainDir);

    for (int i = 0; i < points.size(); ++i)
    {
        accAreaChildren += points[i].sh_representation.get_sh_area(mainDir);
    }

    // if (std::abs(accAreaParent / accAreaChildren) < 0.5f || std::abs(accAreaParent / accAreaChildren) > 2.0f)
    {
        std::cout << "parent: " << accAreaParent << " children: " << accAreaChildren << std::endl;
    }
    */


    return result;

#if 0
    if (points.size() == 0) return result;

    bool invalidChild = false;

    for (unsigned int i = 0; i < points.size(); ++i)
    {
        GiPoint const& p = points[i];

        if (p.area < 0.0f)
        {
            invalidChild = true;
            break;
        }
    }

    for (unsigned int i = 0; i < points.size(); ++i)
    {
        GiPoint const& p = points[i];

        result.pos    += p.pos;
        result.normal += p.normal;
        result.color  += p.color;
        result.area   += p.area;
        result.energy += p.energy;
    }

    // if the normals vary too much, "invalidate" the node as being an average
    bool invalidNormals = result.normal.length() / float(points.size()) < 0.2f;

    result.pos    *= 1.0f / float(points.size());
    result.color  *= 1.0f / float(points.size());
    result.energy *= 1.0f / float(points.size());

    result.normal.normalize();

    // radius needs to be averaged by adding together the area of all discs and then taking the radius of that area
    // assert(area > 0.0f);
    if (invalidNormals || invalidChild)
    {
        result.area = -1.0f;
    }

    return result;
#endif
}


pbLighting_t::pbLighting_t(bool transpShad, int shadowDepth, int rayDepth) : maxSolidAngle(0.5f), _bspTree(NULL)
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
    integratorName = "PointBased";
    integratorShortName = "PBGI";
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
                col += lcol * std::fabs(sp.N * lightRay.dir);
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


void pbLighting_t::generate_gi_points(renderState_t & state)
{
    std::string fileName = "/tmp/pbgi_points_store";

    std::ofstream fileStream(fileName.c_str());

    for (std::map<objID_t, objData_t>::const_iterator iter = scene->meshes.begin(); iter != scene->meshes.end(); ++iter)
    {
        triangleObject_t const* obj = iter->second.obj;

        Y_INFO << "PBGI: Processing obj: " << obj << std::endl;

        std::vector<triangle_t> const& triangles = obj->getTriangles();

        for (unsigned int i = 0; i < triangles.size(); ++i)
        {
            triangle_t const& tri = triangles[i];

            float t;
            intersectData_t iData;

            ray_t ray;

            float triArea = tri.surfaceArea();

            int sampleCount = samplesPerArea * triArea;
            if (sampleCount < 1) sampleCount = 1;

            float radius = std::sqrt(triArea / (sampleCount * M_PI));

            std::vector<point3d_t> samplingPoints;

            /*
            int uCount = std::sqrt(sampleCount) + 1;
            int vCount = uCount;

            for (int u = 0; u < uCount; ++u)
            {
                for (int v = 0; v < vCount; ++v)
                {
                    point3d_t p;
                    vector3d_t n;

                    float s1 = (u + 1) / float(uCount + 1);
                    float s2 = (v + 1) / float(vCount + 1);

                    tri.sample(s1, s2, p, n);

                    samplingPoints.push_back(p);
                    samplingNormals.push_back(n);
                }
            }
            */

            vector3d_t edge0 = tri.getVertex(1) - tri.getVertex(0);

            // make a rectangular sampling pattern on the triangle
            vector3d_t edge1 = tri.getNormal() ^ edge0;

            float e0Length = edge0.length();
            float e1Length = edge1.length();

            // edge0.normalize();
            // edge1.normalize();

            float d = 2.0f * radius / std::sqrt(2.0f);

            int uCount = e0Length / d;
            int vCount = e1Length / d;

            for (int u = 0; u < uCount; ++u)
            {
                for (int v = 0; v < vCount; ++v)
                {
                    point3d_t p =
                            (u + 0.5f) / float(uCount) * edge0 +
                            (v + 0.5f) / float(vCount) * edge1 +
                            tri.getVertex(0);

                    vector3d_t const& normal = vector3d_t(tri.getNormal());
                    point3d_t const& center  = p;

                    ray.from = center + normal;
                    ray.dir  = center - ray.from;

                    if (tri.intersect(ray, &t, iData))
                    {
                        samplingPoints.push_back(p);
                    }
                }
            }

            float accDiscAreas = samplingPoints.size() * M_PI * radius * radius;

            float ratio = triArea / accDiscAreas;

            // radius *= std::sqrt(ratio * 2.0f);
            float singleDiscArea = M_PI * radius * radius;

            std::cout << "samples: " << samplingPoints.size() << " radius; " << radius << std::endl;
            std::cout << "ratio: " << ratio << " area: " << triArea << " accDiscArea: " << accDiscAreas << std::endl;

            for (unsigned int i = 0; i < samplingPoints.size(); ++i)
            {
                point3d_t const& hitPoint = samplingPoints[i];

                surfacePoint_t sp;

                tri.getSurface(sp, hitPoint, iData);

                GiPoint * giPoint = new GiPoint();
                giPoint->pos = vector3d_t(sp.P);
                giPoint->normal = sp.N;
                giPoint->area = singleDiscArea;

                const material_t *material = sp.material;
                BSDF_t bsdfs;
                material->initBSDF(state, sp, bsdfs);
                giPoint->color = material->getDiffuseAtPoint(state, sp);

                color_t incomingLight;
                unsigned int loffs = 0;
                for(std::vector<light_t *>::const_iterator l=lights.begin(); l!=lights.end(); ++l)
                {
                    incomingLight += estimateIncomingLight(state, *l, sp, loffs);
                    loffs++;
                }

                giPoint->energy = incomingLight + material->emission(state, sp, vector3d_t());

                giPoint->sh_representation.calc_coefficients_random(giPoint->normal, giPoint->color, giPoint->energy, giPoint->area);

                // giPoints.push_back(giPoint);

                _bspTree->addPoint(giPoint->pos, giPoint);

                fileStream << *giPoint << std::endl;
            }
        }
    }

    fileStream.close();
}


void load_gi_points(pbLighting_t::MyTree* tree)
{
    std::string fileName = "/tmp/pbgi_points_store";

    std::ifstream fileStream(fileName.c_str());

    while (fileStream.good())
    {
        GiPoint * p = new GiPoint();
        fileStream >> *p;
        tree->addPoint(p->pos, p);
    }

    fileStream.close();
}


bool pbLighting_t::preprocess()
{
    Y_INFO << "PBGI Preprocess" << std::endl;

	bool success = true;
	settings = "";

    std::stringstream ss;
    ss << "type: " << debug_type_str << " | samples: " << samplesPerArea << " | solid angle: " << maxSolidAngle;

    settings = ss.str();

    background = scene->getBackground();
	lights = scene->lights;

    renderState_t state;
    unsigned char userdata[USER_DATA_SIZE+7];
    state.userdata = (void *)( &userdata[7] - ( ((size_t)&userdata[7])&7 ) ); // pad userdata to 8 bytes
    state.cam = scene->getCamera();


    bound_t const& sceneBound = scene->getSceneBound();
    _bspTree = new MyTree(vector3d_t(sceneBound.a), vector3d_t(sceneBound.g), 1000, 16);

    if (do_load_gi_points)
    {
        load_gi_points(_bspTree);
    }
    else
    {
        generate_gi_points(state);
    }

    MyTree::averageData(_bspTree);

    if (debugOutputPointsToFile)
    {
        std::cout << "debugOutputPointsToFile, maxDepth " << _bspTree->getMaxDepth() << " nodes: " << _bspTree->getNodeCount() << std::endl;

        std::string fileName = "/tmp/pbgi_points";

        std::ofstream fileStream(fileName.c_str());

        for (int depth = 0; depth < 20; ++depth)
        {
            std::vector<MyTree const*> nodes = _bspTree->getNodes(depth);
            if (nodes.size() == 0) break;

            for (unsigned int i = 0; i < nodes.size(); ++i)
            {
                GiPoint p = *nodes[i]->getClusteredData();
                p.energy = 1.0f;
                p.area = p.sh_representation.get_area_coefficient(0);
                p.depth = nodes[i]->getDepth();

                fileStream << p << std::endl;
            }
        }

        std::vector<MyTree *> leafs;
        _bspTree->getLeafs(leafs);

        fileName = "/tmp/pbgi_points_leafs";

        std::ofstream fileStream_leafs(fileName.c_str());

        for (unsigned int i = 0; i < leafs.size(); ++i)
        {
            GiPoint p = *leafs[i]->getClusteredData();
            p.area = p.sh_representation.get_area_coefficient(0);
            p.color = p.sh_representation.get_color_coefficient(0);
            p.energy = 1.0f;
            p.depth = 0;

            fileStream_leafs << p << std::endl;

            // SH representation of surfels at depth 2
            for (unsigned int j = 0; j < leafs[i]->getData().size(); ++j)
            {
                GiPoint p = *leafs[i]->getData()[j];
                p.area = p.sh_representation.get_area_coefficient(0);
                p.color = p.sh_representation.get_color_coefficient(0);
                p.energy = 1.0f;
                p.depth = 1;
                fileStream_leafs << p << std::endl;
            }

            // surfels at depth 2
            for (unsigned int j = 0; j < leafs[i]->getData().size(); ++j)
            {
                GiPoint p = *leafs[i]->getData()[j];
                p.depth = 2;
                fileStream_leafs << p << std::endl;
            }
        }

        return false;
    }

    // empty out debug file
    std::ofstream file_stream_fb("/tmp/pbgi_frame_buffer");
    file_stream_fb.close();

//    _bspTree->printAveraged();

    // Y_INFO << "PBGI: BSP: " << *_bspTree << std::endl;

    Y_INFO << "PBGI: Created " << giPoints.size() << " gi points, solid angle: "
           << maxSolidAngle << " "
           << debugTreeDepth << " "
           << "max depth: " << _bspTree->getMaxDepth() << " "
           << "min leaf depth: " << _bspTree->getMinLeafDepth() << " "
           << std::endl;

	return success;
}


color_t pbLighting_t::doPointBasedGiTree(renderState_t & state, surfacePoint_t const& sp, vector3d_t const& wo) const
{
    color_t col(0.0f);
    const material_t *material = sp.material;

    // traverse tree, if solid angle of node > max, traverse into the children
    std::queue<MyTree const*> queue;
    queue.push(_bspTree);

    int shadingNodes = 0;
    int shadingDiscs = 0;

    while (!queue.empty())
    {
        MyTree const* node = queue.front();
        queue.pop();

        // if we are here and we have a leaf, we sample all points in the leaf
        if (node->getIsLeaf())
        {
            std::vector<GiPoint*> const& points = node->getData();

            for (unsigned int i = 0; i < points.size(); ++i)
            {
                GiPoint const& giP = *points[i];

                vector3d_t giToSp = (vector3d_t(sp.P) - giP.pos);

                float distance = giToSp.length();

                giToSp.normalize();

                float cos_normal_gip = -giToSp * sp.N;
                float cos_sp_gip = giP.normal * giToSp;

                float solidAngle = cos_sp_gip * giP.area / (distance * distance);

                ray_t raySpToGiP(giP.pos, giToSp, 0.0001f, distance);

                if (!scene->isShadowed(state, raySpToGiP) && !(cos_sp_gip <= 0.0f) && !(cos_normal_gip <= 0.0f))
                {
                    color_t surfCol = material->eval(state, sp, wo, -giToSp, BSDF_ALL);

                    col += solidAngle * giP.color * cos_normal_gip * giP.energy * surfCol;
                }

                ++shadingDiscs;
            }
        }
        else
        {
            GiPoint const& giP = *node->getClusteredData();

            if (giP.area > 0.0f)
            {
                vector3d_t giToSp = (vector3d_t(sp.P) - giP.pos);

                float distance = giToSp.length();

                giToSp.normalize();

                float cos_sp_gip = giP.normal * giToSp;

                float solidAngle = cos_sp_gip * giP.area / (distance * distance);
                // float solidAngle = calcSolidAngle(giP.radius, distance);

                if (std::abs(solidAngle) > maxSolidAngle)
                // if (true)
                {
                    std::vector<MyTree> const& children = node->getChildren();
                    for (unsigned int i = 0; i < children.size(); ++i)
                    {
                        queue.push(&children[i]);
                    }
                }
                else
                {
                    float cos_normal_gip = -giToSp * sp.N;

                    ray_t raySpToGiP(giP.pos, giToSp, 0.0001f, distance);

                    if (!scene->isShadowed(state, raySpToGiP) && !(cos_sp_gip <= 0.0f) && !(cos_normal_gip <= 0.0f))
                    {
                        color_t surfCol = material->eval(state, sp, wo, -giToSp, BSDF_ALL);

                        col += solidAngle * giP.color * cos_normal_gip * giP.energy * surfCol;
                    }

                    ++shadingNodes;
                }

            }
            else
            {
                std::vector<MyTree> const& children = node->getChildren();
                for (unsigned int i = 0; i < children.size(); ++i)
                {
                    queue.push(&children[i]);
                }
            }

        }
    }

//    std::cout << "shadingDiscs: " << shadingDiscs << " shadingNodes: " << shadingNodes << std::endl;

    return col;
}



color_t pbLighting_t::doPointBasedGiTreeSH(renderState_t & state, surfacePoint_t const& sp, vector3d_t const& wo) const
{
    color_t col(0.0f);
    const material_t *material = sp.material;

    Halton hal2(2);
    Halton hal3(3);

    std::vector<GiPoint const*> debugStorageShadingClusters;

    // traverse tree, if solid angle of node > max, traverse into the children
    std::queue<MyTree const*> queue;
    queue.push(_bspTree);

    while (!queue.empty())
    {
        MyTree const* node = queue.front();
        queue.pop();

        // if we are here and we have a leaf, we sample all points in the leaf
        if (node->getIsLeaf() && node->getData().size() > 0)
        {
            float leaf_contribution = 0.0f;
            float discs_contribution = 0.0f;

            if (pixel_x == state.pixel_x && pixel_y == state.pixel_y)
            {
                GiPoint const& giP = *node->getClusteredData();

                vector3d_t giToSp = (vector3d_t(sp.P) - giP.pos);

                float distance = giToSp.length();

                giToSp.normalize();

                // float cos_sp_gip = giP.normal * giToSp;

                float const area = giP.sh_representation.get_sh_area(giToSp);

                float solidAngle = std::max(0.0f, area / (distance * distance));

                float cos_normal_gip = (-giToSp) * sp.N;

                ray_t raySpToGiP(giP.pos, giToSp, 0.0001f, distance);

                if (cos_normal_gip > 0.0f && !scene->isShadowed(state, raySpToGiP))
                {
                    color_t cluster_contribution = giP.sh_representation.get_sh_color(giToSp) * solidAngle;

                    leaf_contribution = cluster_contribution.energy();
                }
            }



            std::vector<GiPoint*> const& points = node->getData();

            for (unsigned int i = 0; i < points.size(); ++i)
            {
                GiPoint const& giP = *points[i];

                vector3d_t giToSp = (vector3d_t(sp.P) - giP.pos);

                float distance = giToSp.length();

                giToSp.normalize();

                float cos_normal_gip = std::max((-giToSp) * sp.N, 0.0f);
                float cos_sp_gip = std::max(giP.normal * giToSp, 0.0f);

                // float const area = cos_sp_gip * giP.area;
                float const area = giP.sh_representation.get_sh_area(giToSp);

                float solidAngle = std::max(0.0f, area / (distance * distance));

                ray_t raySpToGiP(giP.pos, giToSp, 0.0001f, distance);

                if (cos_sp_gip > 0.0f && cos_normal_gip > 0.0f && !scene->isShadowed(state, raySpToGiP))
                {
                    color_t const surfCol = material->eval(state, sp, wo, -giToSp, BSDF_ALL);
                    // color_t const contribution = solidAngle * giP.color * giP.energy;
                    color_t const contribution = giP.sh_representation.get_sh_color(giToSp) * solidAngle;
                    col += contribution * cos_normal_gip * surfCol;

                    if (pixel_x == state.pixel_x && pixel_y == state.pixel_y)
                    {
                        discs_contribution += contribution.energy();

                        GiPoint * p = new GiPoint(giP);
                        p->area = area;;
                        p->color = contribution;
                        p->energy = 1.0f;
                        p->depth = node->getDepth();

                        debugStorageShadingClusters.push_back(p);
                    }
                }
            }

            if (pixel_x == state.pixel_x && pixel_y == state.pixel_y)
            {
                discs_contribution /= float(points.size());

                std::cout << "visible leaf/discs: " << leaf_contribution << " / " << discs_contribution << " ratio: " <<
                             (leaf_contribution / discs_contribution) << std::endl;
            }

        }
        else
        {
            if (node->isNodeBehindPlane(vector3d_t(sp.P), sp.N)) continue;

            GiPoint const& giP = *node->getClusteredData();

            vector3d_t giToSp = (vector3d_t(sp.P) - giP.pos);

            float distance = giToSp.length();

            giToSp.normalize();

            // float cos_sp_gip = giP.normal * giToSp;

            float const area = giP.sh_representation.get_sh_area(giToSp);

            float solidAngle = std::max(0.0f, area / (distance * distance));
            // float solidAngle = area / (distance * distance);

            if (solidAngle > maxSolidAngle || distance < node->getRadius() * 1.5f)
            {
                std::vector<MyTree> const& children = node->getChildren();
                for (unsigned int i = 0; i < children.size(); ++i)
                {
                    queue.push(&children[i]);
                }
            }
            else
            {
                float cos_normal_gip = (-giToSp) * sp.N;

                ray_t raySpToGiP(giP.pos, giToSp, 0.0001f, distance);

                if (cos_normal_gip > 0.0f && !scene->isShadowed(state, raySpToGiP))
                    // if (!scene->isShadowed(state, raySpToGiP))
                {
                    color_t surfCol = material->eval(state, sp, wo, -giToSp, BSDF_ALL);

                    // color_t cluster_contribution = giP.sh_representation.get_sh_color(giToSp) / (distance * distance);
                    color_t cluster_contribution = giP.sh_representation.get_sh_color(giToSp) * solidAngle;
                    col += surfCol * cluster_contribution * cos_normal_gip;

                    if (pixel_x == state.pixel_x && pixel_y == state.pixel_y)
                    {
                        GiPoint * p = new GiPoint(giP);
                        p->area = area;
                        p->color = cluster_contribution;
                        p->energy = 1.0f;
                        p->depth = node->getDepth();

                        debugStorageShadingClusters.push_back(p);
                    }
                }
            }
        }
    }

    if (pixel_x == state.pixel_x && pixel_y == state.pixel_y)
    {
        std::string fileName = "/tmp/pbgi_points_debug_pixel_shading_clusters";

        std::ofstream fileStream2(fileName.c_str());

        for (unsigned int i = 0; i < debugStorageShadingClusters.size(); ++i)
        {
            fileStream2 << *debugStorageShadingClusters[i] << std::endl;
        }
    }

    return col;
}




color_t pbLighting_t::doPointBasedGiTree_sh_fb(renderState_t & state, surfacePoint_t const& sp, vector3d_t const& wo) const
{
    color_t col(0.0f);
    const material_t *material = sp.material;

    // traverse tree, if solid angle of node > max, traverse into the children
    std::queue<MyTree const*> queue;
    queue.push(_bspTree);

    int shadingNodes = 0;
    int shadingDiscs = 0;

    std::vector<GiPoint const*> debugStorageAllClusters;
    std::vector<GiPoint const*> debugStorageShadingClusters;

    cube_raster_buffer_type frame_buffer;

    float total_energy = 0.0f;

    std::vector<color_t> colors;
    colors.push_back(color_t(1.0f, 0.0f, 0.0f)); // red
    colors.push_back(color_t(0.0f, 1.0f, 0.0f)); // green
    colors.push_back(color_t(0.0f, 0.0f, 1.0f)); // blue
    colors.push_back(color_t(1.0f, 1.0f, 0.0f)); // yellow
    colors.push_back(color_t(0.0f, 1.0f, 1.0f)); // turquois
    colors.push_back(color_t(1.0f, 0.0f, 1.0f)); // violet
    colors.push_back(color_t(0.5f, 0.5f, 0.5f)); // gray
    colors.push_back(color_t(1.0f, 1.0f, 1.0f)); // white (7)


    while (!queue.empty())
    {
        MyTree const* node = queue.front();
        queue.pop();

        if (col.col2bri() < 0)
        {
//            std::cout << "col 0 < 0: " << col << " " << state.pixel_x << " " << state.pixel_y << std::endl;
        }

        // if we are here and we have a leaf, we sample all points in the leaf
        if (node->getIsLeaf() && node->getData().size() > 0)
        {
            std::vector<GiPoint*> const& points = node->getData();

            for (unsigned int i = 0; i < points.size(); ++i)
            {
                GiPoint const& giP = *points[i];

                vector3d_t giToSp = (vector3d_t(sp.P) - giP.pos);

                float const distance = giToSp.length();

                giToSp.normalize();

                float const cos_sp_gip = std::max(giP.normal * giToSp, 0.0f);

                float const area = cos_sp_gip * giP.area;
                // float const area = std::max(0.0f, giP.sh_representation.get_sh_area(giToSp));

                float const solidAngle = area / (distance * distance);

                color_t const contribution = giP.color * giP.energy;
                // color_t const contribution = colors[7];
                // color_t const contribution = giP.sh_representation.get_sh_color(giToSp);

                // bool use_rays = (distance < std::sqrt(area / M_PI) * 4.0f);
                bool const use_rays = false;

                frame_buffer.add_point(-giToSp, contribution, solidAngle, distance, use_rays);

                if (pixel_x == state.pixel_x && pixel_y == state.pixel_y)
                {
                    GiPoint * p = new GiPoint(giP);
                    p->area = area; //  * cos_normal_gip;
                    p->color = contribution;
                    p->energy = 1.0f;
                    p->depth = node->getDepth();

                    debugStorageShadingClusters.push_back(p);

                    debugStorageAllClusters.push_back(p);

                    total_energy += contribution.col2bri();
                }

                ++shadingDiscs;
            }
        }
        else
        {
            if (node->isNodeBehindPlane(vector3d_t(sp.P), sp.N)) continue;

            GiPoint const& giP = *node->getClusteredData();

            vector3d_t giToSp = (vector3d_t(sp.P) - giP.pos);

            float const distance = giToSp.length();

            giToSp.normalize();

            // float cos_sp_gip = giP.normal * giToSp;

            // float const area = std::max(0.0f, giP.sh_representation.get_sh_area(giToSp));
            float const area = node->getRadius() * node->getRadius() * M_PI;

            float const solidAngle = area / (distance * distance);

            if (solidAngle > maxSolidAngle || distance < node->getRadius() * 1.5f)
            // if (1.0f / (distance) >= maxSolidAngle || distance < node->getRadius() * 1.5f)
            // if (node->getDepth() < debugTreeDepth || distance < node->getRadius() * 2.0f)
            // if (false)
            // if (solidAngle >= maxSolidAngle || distance < node->getRadius() * 4.0f)
            {
                std::vector<MyTree> const& children = node->getChildren();
                for (unsigned int i = 0; i < children.size(); ++i)
                {
                    queue.push(&children[i]);
                }
            }
            else
            {
                // if (solidAngle < 0.0001f) continue;

                // float const cos_normal_gip = -giToSp * sp.N;

                color_t cluster_contribution = giP.sh_representation.get_sh_color(giToSp);
                // color_t const cluster_contribution = colors[node->getDepth()];

                bool const use_rays = false;
                frame_buffer.add_point(-giToSp, cluster_contribution, solidAngle, distance, use_rays);

                if (pixel_x == state.pixel_x && pixel_y == state.pixel_y)
                {
                    GiPoint * p = new GiPoint(giP);
                    p->area = area;
                    p->color = cluster_contribution;
                    p->energy = 1.0f;
                    p->depth = node->getDepth();

                    debugStorageShadingClusters.push_back(p);

                    total_energy += cluster_contribution.col2bri();
                }

                ++shadingNodes;
            } // end else (bad solid angle)


            if (pixel_x == state.pixel_x && pixel_y == state.pixel_y)
            {
                GiPoint * p = new GiPoint(giP);
                p->area = area;
                p->color = giP.sh_representation.get_sh_color(giToSp);
                p->energy = 1.0f;
                p->depth = node->getDepth();

                debugStorageAllClusters.push_back(p);
            }
        }
    }


    frame_buffer.accumulate();
    color_t surfCol = material->getDiffuseAtPoint(state, sp); //material->eval(state, sp, wo, vector3d_t(0.0f), BSDF_ALL);
    col = frame_buffer.get_diffuse(sp.N) * surfCol;

    if (pixel_x == state.pixel_x && pixel_y == state.pixel_y)
    {
        std::cout << "shadingDiscs: " << shadingDiscs << " shadingNodes: " << shadingNodes << " total energy: " << total_energy
                  << " total energy in fb: " << frame_buffer.get_total_energy() << std::endl;

        std::string fileName = "/tmp/pbgi_points_debug_pixel";

        std::ofstream fileStream(fileName.c_str());

        for (unsigned int i = 0; i < debugStorageAllClusters.size(); ++i)
        {
            fileStream << *debugStorageAllClusters[i] << std::endl;
        }

        fileName = "/tmp/pbgi_points_debug_pixel_shading_clusters";

        std::ofstream fileStream2(fileName.c_str());

        for (unsigned int i = 0; i < debugStorageShadingClusters.size(); ++i)
        {
            fileStream2 << *debugStorageShadingClusters[i] << std::endl;
        }

        /*
        std::cout << "points for " << pixel_x << ", " << pixel_y << ":" << std::endl;
        for (unsigned int i = 0; i < debugStorage.size(); ++i)
        {
            std::cout << *debugStorage[i] << std::endl;
        }
        */

    }

    if (state.pixel_x >= pixel_x && state.pixel_x <= pixel_x + 100 && pixel_y == state.pixel_y)
    {
        std::string fileName = "/tmp/pbgi_frame_buffer";

        std::ofstream file_stream_fb(fileName.c_str(), std::ios::app);

        file_stream_fb << state.pixel_x << " " << frame_buffer << std::endl;

        file_stream_fb.close();
    }

    return col;
}




color_t pbLighting_t::doPointBasedGiTreeSH_leafs_only(renderState_t & state, surfacePoint_t const& sp, vector3d_t const& wo) const
{
    color_t col(0.0f);
    const material_t *material = sp.material;

    // traverse tree, if solid angle of node > max, traverse into the children
    std::queue<MyTree const*> queue;
    queue.push(_bspTree);

    int shadingNodes = 0;
    int shadingDiscs = 0;

    while (!queue.empty())
    {
        MyTree const* node = queue.front();
        queue.pop();

        // if we are here and we have a leaf, we sample all points in the leaf
        if (node->getIsLeaf())
        {
            // if (node->isNodeBehindPlane(vector3d_t(sp.P), sp.N)) continue;

            GiPoint const& giP = *node->getClusteredData();

            vector3d_t giToSp = (vector3d_t(sp.P) - giP.pos);

            float distance = giToSp.length();

            giToSp.normalize();

            float cos_normal_gip = std::max(-giToSp * sp.N, 0.0f);
            float cos_sp_gip = std::max(giP.normal * giToSp, 0.0f);

            ray_t raySpToGiP(giP.pos, giToSp, 0.0001f, distance);

            // if (!scene->isShadowed(state, raySpToGiP) && !(cos_normal_gip <= 0.0f))
            if (!scene->isShadowed(state, raySpToGiP))
            {
                color_t surfCol = material->eval(state, sp, wo, -giToSp, BSDF_ALL);

                color_t cluster_contribution = giP.sh_representation.get_sh_color(giToSp);

                col += cos_normal_gip * surfCol * cluster_contribution / (distance * distance);
            }

            ++shadingNodes;
        }
        else
        {
//            if (node->isNodeBehindPlane(vector3d_t(sp.P), sp.N)) continue;

            std::vector<MyTree> const& children = node->getChildren();
            for (unsigned int i = 0; i < children.size(); ++i)
            {
                queue.push(&children[i]);
            }
        }
    }

//    std::cout << "shadingDiscs: " << shadingDiscs << " shadingNodes: " << shadingNodes << std::endl;

    return col;
}



color_t pbLighting_t::doPointBasedGi(renderState_t & state, surfacePoint_t const& sp, vector3d_t const& wo) const
{
    color_t col(0.0f);
    const material_t *material = sp.material;

    for (unsigned int i = 0; i < giPoints.size(); ++i)
    {
        GiPoint const& giP = giPoints[i];

        vector3d_t giToSp = (vector3d_t(sp.P) - giP.pos);

        float distance = giToSp.length();

        giToSp.normalize();

        ray_t raySpToGiP(giP.pos, giToSp, 0.0001f, distance);

        if (scene->isShadowed(state, raySpToGiP)) continue;

        float distanceThreshold = 0.0f;

        if (distance < distanceThreshold)
        {
            color_t colAcc(0.0f);

            Halton hal2(2);
            Halton hal3(3);

            int sampleCount = 100.0f * (distanceThreshold - distance) / distanceThreshold;
            if (sampleCount < 1) sampleCount = 1;

            // sampleCount = 10;

            vector3d_t tangent, bitangent;
            createCS(giP.normal, tangent, bitangent);

            for (int i = 0; i < sampleCount; ++i)
            {
                hal2.setStart(i);
                hal3.setStart(i);

                float radiusSamples = std::sqrt(hal2.getNext()) * std::sqrt(giP.area / (2.0f * M_PI));
                float angleSample = hal3.getNext() * 2.0f * M_PI;

                vector3d_t baseVector(radiusSamples * std::cos(angleSample), radiusSamples * std::sin(angleSample), 0.0f);

                vector3d_t transformedVector = tangent * baseVector.x + bitangent * baseVector.y;

                vector3d_t newGiPos = giP.pos + transformedVector;

                giToSp = (vector3d_t(sp.P) - newGiPos);
                distance = giToSp.length();
                giToSp.normalize();

                float cos_sp_gip = giP.normal * giToSp;

                if (cos_sp_gip <= 0.0f) continue;

                float cos_normal_gip = -giToSp * sp.N;

                if (cos_normal_gip <= 0.0f) continue;

                color_t surfCol = material->eval(state, sp, wo, -giToSp, BSDF_ALL);

                colAcc += giP.color * cos_normal_gip * giP.energy * surfCol / (distance * distance);
            }

            // if (colAcc.abscol2bri() > 5.0f) colAcc = colAcc * 5.0f / colAcc.abscol2bri();

            col += colAcc * giP.area / float(sampleCount);
        }
        else
        {
            float cos_sp_gip = giP.normal * giToSp;

            if (cos_sp_gip <= 0.0f) continue;

            float cos_normal_gip = -giToSp * sp.N;

            if (cos_normal_gip <= 0.0f) continue;

            float solidAngle = cos_sp_gip * giP.area / (distance * distance);

            color_t surfCol = material->eval(state, sp, wo, -giToSp, BSDF_ALL);

            col += solidAngle * giP.color * cos_normal_gip * giP.energy * surfCol;
        }
    }

    return col;
}


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

	// Shoot ray into scene
    bool useTree = true;
	
	if(scene->intersect(ray, sp)) // If it hits
	{
        if (debug && !useTree)
        {
            // find closest gi point
            float minDist = 1e10f;
            int index = -1;

            for (unsigned int i = 0; i < giPoints.size(); ++i)
            {
                float distSqr = (giPoints[i].pos - sp.P).lengthSqr();

                if (distSqr < minDist)
                {
                    minDist = distSqr;
                    index = i;
                }
            }

            col = giPoints[index].color * giPoints[index].energy;

            return colorA_t(col, 1.0f);
        }
        else if (debug && useTree)
        {
            // find closest gi point
            float minDist = 1e10f;
            int index = -1;

            std::vector<MyTree const*> points = _bspTree->getNodes(debugTreeDepth);

            for (unsigned int i = 0; i < points.size(); ++i)
            {
                GiPoint const& p = *points[i]->getClusteredData();

                float distSqr = (p.pos - sp.P).lengthSqr();

                if (distSqr < minDist)
                {
                    minDist = distSqr;
                    index = i;
                }
            }

            GiPoint const& p = *points[index]->getClusteredData();

            /*
            MyTree const& t = _bspTree->query(vector3d_t(sp.P));
            GiPoint const& p = t.getAveragedData();
            */

            col = p.color * p.energy;

            return colorA_t(col, 1.0f);
        }

		unsigned char userdata[USER_DATA_SIZE];
		const material_t *material = sp.material;
		BSDF_t bsdfs;

		state.userdata = (void *) userdata;
		vector3d_t wo = -ray.dir;
		if(state.raylevel == 0) state.includeLights = true;
		
		material->initBSDF(state, sp, bsdfs);
		
        if(bsdfs & BSDF_EMIT) col += material->emission(state, sp, wo);
		
		if(bsdfs & BSDF_DIFFUSE)
		{
            if (!indirectOnly)
            {
                col += estimateAllDirectLight(state, sp, wo);
            }

            if (debug_type == NoTree)
            {
                col += doPointBasedGi(state, sp, wo);
            }
            else if (debug_type == Tree)
            {
                col += doPointBasedGiTree(state, sp, wo);
            }
            else if (debug_type == Tree_sh)
            {
                col += doPointBasedGiTreeSH(state, sp, wo);
            }
            else if (debug_type == Tree_sh_fb)
            {
                col += doPointBasedGiTree_sh_fb(state, sp, wo);
            }
            else if (debug_type == Tree_sh_leafs)
            {
                col += doPointBasedGiTreeSH_leafs_only(state, sp, wo);
            }
        }
		
        // recursiveRaytrace(state, ray, bsdfs, sp, wo, col, alpha);
		
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
    float maxSolidAngle = 0.5f;
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

    inte->samplesPerArea = samples;
    inte->debug = debug;
    inte->indirectOnly = indirectOnly;
    inte->maxSolidAngle = maxSolidAngle;
    inte->debugTreeDepth = debugTreeDepth;
    inte->debugOutputPointsToFile = debugOutputPointsToFile;
    if (debug_type_map.find(debug_type) != debug_type_map.end())
    {
        inte->debug_type = debug_type_map[debug_type];
    }

    std::cout << "Debug type: " << inte->debug_type << " " << debug_type << std::endl;
    inte->debug_type_str = debug_type;

    inte->render_single_pixel = render_single_pixel;
    inte->pixel_x = pixel_x;
    inte->pixel_y = pixel_y;

    inte->do_load_gi_points = do_load_gi_points;

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
