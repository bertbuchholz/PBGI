#include <sstream>
#include <cassert>
#include <algorithm>
#include <fstream>

#include <yafray_config.h>
#include <core_api/environment.h>
#include <core_api/material.h>
#include <core_api/mcintegrator.h>
#include <core_api/background.h>
#include <core_api/light.h>
#include <yafraycore/meshtypes.h>
#include <utilities/mcqmc.h>

#include <utilities/RegularBspTree.h>

__BEGIN_YAFRAY

float calcSolidAngle(float radius, float distance)
{
    return 2.0f * M_PI * radius * radius / (distance * distance);
}

struct GiPoint
{
    GiPoint() :
        pos(0.0f),
        normal(0.0f),
        color(0.0f),
        radius(0.0f),
        energy(0.0f)
    {}

    vector3d_t pos;
    vector3d_t normal;
    color_t    color;
    float      radius;
    color_t    energy;

    friend std::ostream & operator<<(std::ostream & s, GiPoint const& p)
    {
        s <<
             p.pos.x << " " << p.pos.y << " " << p.pos.z << " " <<
             p.normal.x << " " << p.normal.y << " " << p.normal.z << " " <<
             p.color.R << " " << p.color.G << " " << p.color.B << " " <<
             p.radius << " " <<
             p.energy.R << " " << p.energy.G << " " << p.energy.B
             ;

        return s;
    }
};



GiPoint averageGiPoints(std::vector<GiPoint> const& points)
{
    GiPoint result;

    if (points.size() == 0) return result;

    bool invalidChild = false;

    for (unsigned int i = 0; i < points.size(); ++i)
    {
        GiPoint const& p = points[i];

        if (p.radius < 0.0f)
        {
            invalidChild = true;
            break;
        }
    }

    float area = 0.0f;

    for (unsigned int i = 0; i < points.size(); ++i)
    {
        GiPoint const& p = points[i];

        result.pos    += p.pos;
        result.normal += p.normal;
        result.color  += p.color;
        // result.radius += p.radius;
        result.energy += p.energy;

        area += 2.0f * M_PI * p.radius * p.radius;
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
        result.radius = -1.0f;
    }
    else
    {
        result.radius = std::sqrt(area / (2.0f * M_PI));
    }

    /*
    if (result.pos.x > -2.9f)
    {
        std::cout << "Wrong averaging: " << std::endl;

        for (unsigned int j = 0; j < points.size(); ++j)
        {
            std::cout << points[j] << std::endl;
        }
    }
    */

    return result;
}


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

private:
    std::vector<GiPoint> giPoints;
    int samplesPerArea;
    bool debug;
    bool indirectOnly;
    float maxSolidAngle;
    int debugTreeDepth;
    bool debugOutputPointsToFile;

    MyTree* _bspTree;
    std::vector<MyTree *> leafNodes;
};

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



bool pbLighting_t::preprocess()
{
    Y_INFO << "PBGI Preprocess" << std::endl;

	bool success = true;
	settings = "";
	
	background = scene->getBackground();
	lights = scene->lights;

    renderState_t state;
    unsigned char userdata[USER_DATA_SIZE+7];
    state.userdata = (void *)( &userdata[7] - ( ((size_t)&userdata[7])&7 ) ); // pad userdata to 8 bytes
    state.cam = scene->getCamera();

    Halton hal2(2);
    Halton hal3(3);

    bound_t const& sceneBound = scene->getSceneBound();
    _bspTree = new RegularBspTree<vector3d_t, 3, GiPoint>(vector3d_t(sceneBound.a), vector3d_t(sceneBound.g), 1000, 16);

    for (std::map<objID_t, objData_t>::const_iterator iter = scene->meshes.begin(); iter != scene->meshes.end(); ++iter)
    {
        triangleObject_t const* obj = iter->second.obj;

        Y_INFO << "PBGI: Processing obj: " << obj << std::endl;

        std::vector<triangle_t> const& triangles = obj->getTriangles();

        for (unsigned int i = 0; i < triangles.size(); ++i)
        {
            hal2.setStart(i);
            hal3.setStart(i);

            triangle_t const& tri = triangles[i];

            float t;
            intersectData_t iData;

            ray_t ray;

            float area = tri.surfaceArea();

            int sampleCount = samplesPerArea * area;
            if (sampleCount < 1) sampleCount = 1;

            float radius = std::sqrt(area / (sampleCount * M_PI * 2.0f));

            std::cout << "samples: " << sampleCount << " radius; " << radius << std::endl;

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

            float discArea = samplingPoints.size() * 2.0f * M_PI * radius * radius;

            float ratio = area / discArea;

            radius *= std::sqrt(ratio * 2.0f);

            std::cout << "ratio: " << ratio << " area: " << area << " discArea: " << discArea << std::endl;


            for (unsigned int i = 0; i < samplingPoints.size(); ++i)
            {
                point3d_t const& hitPoint  = samplingPoints[i];

                surfacePoint_t sp;

                tri.getSurface(sp, hitPoint, iData);

                GiPoint giPoint;
                giPoint.pos = vector3d_t(sp.P);
                giPoint.normal = sp.N;
                giPoint.radius = radius;

                const material_t *material = sp.material;
                BSDF_t bsdfs;
                material->initBSDF(state, sp, bsdfs);
                giPoint.color = material->getDiffuseAtPoint(state, sp);

                color_t incomingLight;
                unsigned int loffs = 0;
                for(std::vector<light_t *>::const_iterator l=lights.begin(); l!=lights.end(); ++l)
                {
                    incomingLight += estimateIncomingLight(state, *l, sp, loffs);
                    loffs++;
                }

                giPoint.energy = incomingLight;

                giPoints.push_back(giPoint);

                _bspTree->addPoint(giPoint.pos, giPoint);
            }
        }
    }

    MyTree::averageData(_bspTree);

    _bspTree->getLeafs(leafNodes);

    if (debugOutputPointsToFile)
    {
        std::cout << "debugOutputPointsToFile, maxDepth " << _bspTree->getMaxDepth() << " nodes: " << _bspTree->getNodeCount() << std::endl;

        std::string fileName = "/tmp/pbgi_points";

        std::ofstream fileStream(fileName.c_str());

        for (unsigned int i = 0; i < giPoints.size(); ++i)
        {
            fileStream << giPoints[i] << std::endl;
        }

        for (int depth = 0; depth < 20; ++depth)
        {
            std::vector<MyTree const*> nodes = _bspTree->getNodes(depth);
            if (nodes.size() == 0) break;

            std::stringstream ss;
            ss << "/tmp/pbgi_points_" << depth;

            fileName = ss.str();

            std::ofstream fileStream_depth(fileName.c_str());

            for (unsigned int i = 0; i < nodes.size(); ++i)
            {
                fileStream_depth << nodes[i]->getAveragedData() << std::endl;
            }
        }

        std::vector<MyTree *> leafs;
        _bspTree->getLeafs(leafs);

        fileName = "/tmp/pbgi_points_leafs";

        std::ofstream fileStream_leafs(fileName.c_str());

        for (unsigned int i = 0; i < leafs.size(); ++i)
        {
            fileStream_leafs << leafs[i]->getAveragedData() << std::endl;
        }

        return false;
    }

//    _bspTree->printAveraged();

    // Y_INFO << "PBGI: BSP: " << *_bspTree << std::endl;

    Y_INFO << "PBGI: Created " << giPoints.size() << " gi points, solid angle: "
           << maxSolidAngle << " "
           << debugTreeDepth << " "
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
            std::vector<GiPoint> const& points = node->getData();

            for (unsigned int i = 0; i < points.size(); ++i)
            {
                GiPoint const& giP = points[i];

                vector3d_t giToSp = (vector3d_t(sp.P) - giP.pos);

                float distance = giToSp.length();

                giToSp.normalize();

                float cos_normal_gip = -giToSp * sp.N;
                float cos_sp_gip = giP.normal * giToSp;

                float solidAngle = cos_sp_gip * calcSolidAngle(giP.radius, distance);

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
            GiPoint const& giP = node->getAveragedData();

            if (giP.radius > 0.0f)
            {
                vector3d_t giToSp = (vector3d_t(sp.P) - giP.pos);

                float distance = giToSp.length();

                giToSp.normalize();

                float cos_sp_gip = giP.normal * giToSp;

                float solidAngle = cos_sp_gip * calcSolidAngle(giP.radius, distance);
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

        float distanceThreshold = 1.0f;

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

                float radiusSamples = std::sqrt(hal2.getNext()) * giP.radius;
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

            col += colAcc * giP.radius / float(sampleCount);
        }
        else
        {
            float cos_sp_gip = giP.normal * giToSp;

            if (cos_sp_gip <= 0.0f) continue;

            float cos_normal_gip = -giToSp * sp.N;

            if (cos_normal_gip <= 0.0f) continue;

            float solidAngle = cos_sp_gip * calcSolidAngle(giP.radius, distance);

            color_t surfCol = material->eval(state, sp, wo, -giToSp, BSDF_ALL);

            col += solidAngle * giP.color * cos_normal_gip * giP.energy * surfCol;
        }
    }

    return col;
}


colorA_t pbLighting_t::integrate(renderState_t &state, diffRay_t &ray) const
{
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
                GiPoint const& p = points[i]->getAveragedData();

                float distSqr = (p.pos - sp.P).lengthSqr();

                if (distSqr < minDist)
                {
                    minDist = distSqr;
                    index = i;
                }
            }

            GiPoint const& p = points[index]->getAveragedData();

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
		
		if(bsdfs & BSDF_EMIT) col += material->emit(state, sp, wo);
		
		if(bsdfs & BSDF_DIFFUSE)
		{
            if (!indirectOnly)
            {
                col += estimateAllDirectLight(state, sp, wo);
            }

            if (useTree)
            {
                col += doPointBasedGiTree(state, sp, wo);
            }
            else
            {
                col += doPointBasedGi(state, sp, wo);
            }
		}
		
		recursiveRaytrace(state, ray, bsdfs, sp, wo, col, alpha);
		
		float m_alpha = material->getAlpha(state, sp, wo);
		alpha = m_alpha + (1.f - m_alpha) * alpha;
	}
	else // Nothing hit, return background if any
	{
		if(background) col += (*background)(ray, state, false);
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
