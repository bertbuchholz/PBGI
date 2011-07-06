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

#include <yafraycore/timer.h>

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
    }

    // result->sh_representation = result->sh_representation / float(points.size());

    result->sh_representation.normalize_color(1.0f / float(points.size()));

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


    for (unsigned int i = 0; i < points.size(); ++i)
    {
        GiPoint const* p = points[i];

        result->pos    += p->pos;
        result->normal += p->normal;
        result->color  += p->color;
        result->energy += p->energy;
    }

    result->pos    *= 1.0f / float(points.size());
    result->color  *= 1.0f / float(points.size());
    result->energy *= 1.0f / float(points.size());

    result->normal.normalize();

    return result;
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




float radicalInverse(int n, int base)
{
    float value = 0;
    float invBase = 1.0/(float)(base), invBi = invBase;

    while(n > 0)
    {
        int d_i = (n % base);
        value += d_i * invBi;
        n /= base;
        invBi *= invBase;

    }

    return value;
}


vector3d_t halton_3(int n)
{
    vector3d_t result;

    result[0] = radicalInverse(n, 2);
    result[1] = radicalInverse(n, 3);
    result[2] = radicalInverse(n, 5);

    return result;
}


vector3d_t hammersley_3(int const n, int const n_max)
{
    vector3d_t result;

    result[0] = n / float(n_max);
    result[1] = radicalInverse(n, 2);
    result[2] = radicalInverse(n, 3);

    return result;
}



struct pbgi_sample_t
{
    triangle_t const* tri_pointer;
    point3d_t position;
};


unsigned int hash_function(point3d_t const& sample, const int n, const float cell_size)
{
    /*
    hash(x,y,z) = ( x p1 xor y p2 xor z p3) mod n
    where p1, p2, p3 are large prime numbers, in
    our case 73856093, 19349663, 83492791
    */

    // float bb_size = 6.0f;

    unsigned int i_x = int((sample.x + 10000.0f) / cell_size);
    unsigned int i_y = int((sample.y + 10000.0f) / cell_size);
    unsigned int i_z = int((sample.z + 10000.0f) / cell_size);

    // std::cout << "i_x: " << (i_x % n) << " sample.x: " << sample.x << std::endl;

    unsigned int hash_value = ((i_x * 73856093) ^ (i_y * 19349663) ^ (i_z * 83492791)) % n;
    // unsigned int hash_value = (i_x + i_y + i_z) % n;

    return hash_value;

    // return i_x % n;

            /*

    unsigned int i_x = std::abs(int(sample.x / cell_size));
    unsigned int i_y = std::abs(int(sample.y / cell_size));
    unsigned int i_z = std::abs(int(sample.z / cell_size));

    int hash_value = ((i_x * 73856093) ^ (i_y * 19349663) ^ (i_z * 83492791)) % n;

    assert(hash_value >= 0 && hash_value < n);

    return hash_value;

    //return (i_x + i_y + i_z) % n;

*/

    // return rand() / RAND_MAX * (n - 1);
}

float generate_histogram(std::vector<pbgi_sample_t> const& samples, float const min_radius)
{
    std::cout << "generate_histogram() start" << std::endl;

    std::vector<float> distances(samples.size());

    int bin_count = samples.size() * 0.05f;
    float cell_size = min_radius * 20.0f;

    std::vector<std::vector<int> > hash_map(bin_count);

    for (unsigned int i = 0; i < samples.size(); ++i)
    {
        for (int x = -1; x <= 1; ++x)
        {
            for (int y = -1; y <= 1; ++y)
            {
                for (int z = -1; z <= 1; ++z)
                {
                    int const hash_value = hash_function(samples[i].position + point3d_t(x * cell_size, y * cell_size, z * cell_size), bin_count, cell_size);
                    hash_map[hash_value].push_back(i);
                }
            }
        }
    }

    for (unsigned int i = 0; i < samples.size(); ++i)
    {
        float smallest_distance_i = 1e10f;

        int const hash_value = hash_function(samples[i].position, bin_count, cell_size);
        std::vector<int> const& neighbors = hash_map[hash_value];

        for (unsigned int j = 0; j < neighbors.size(); ++j)
        {
            int sample_index = neighbors[j];

            if (sample_index == i) continue;

            float const distance = (samples[i].position - samples[sample_index].position).length();

            if (distance < smallest_distance_i)
            {
                smallest_distance_i = distance;
            }
        }

        distances[i] = smallest_distance_i;
    }

/*
    for (unsigned int i = 0; i < samples.size(); ++i)
    {
        float smallest_distance_i = 1e10f;

        for (unsigned int j = 0; j < samples.size(); ++j)
        {
            if (i == j) continue;

            float const distance = (samples[i].position - samples[j].position).length();

            if (distance < smallest_distance_i)
            {
                smallest_distance_i = distance;
            }
        }

        distances[i] = smallest_distance_i;
    }
*/

    float largest_distance = -1e10f;
    float smallest_distance = 1e10f;

    for (unsigned int i = 0; i < distances.size(); ++i)
    {
        if (distances[i] > largest_distance)
        {
            largest_distance = distances[i];
        }

        if (distances[i] < smallest_distance)
        {
            smallest_distance = distances[i];
        }
    }

    float const histo_lowest = smallest_distance;
    float const histo_highest = largest_distance;
    int bins = 100;

    std::vector<int> histogram(bins, 0);

    std::cout << "generate_histogram(): " << histo_highest << std::endl;

    for (unsigned int i = 0; i < distances.size(); ++i)
    {
        float dist = distances[i];

        int const bin = (dist - histo_lowest) / (histo_highest - histo_lowest) * bins;

        ++histogram[bin];
    }

    std::ofstream histo_file("/tmp/histogram");

    for (unsigned int i = 0; i < histogram.size(); ++i)
    {
        histo_file << ((histo_highest - histo_lowest) / float(bins) * i + histo_lowest) << " " << histogram[i] / float(samples.size()) << std::endl;
    }

    std::cout << "generate_histogram() finished" << std::endl;

    return largest_distance;
}

std::vector<triangle_t const*> get_scene_triangles(std::map<objID_t, objData_t> const& meshes)
{
    std::vector<triangle_t const*> triangles;

    for (std::map<objID_t, objData_t>::const_iterator iter = meshes.begin(); iter != meshes.end(); ++iter)
    {
        triangleObject_t const* obj = iter->second.obj;

        std::vector<triangle_t> const& obj_triangles = obj->getTriangles();

        for (unsigned int i = 0; i < obj_triangles.size(); ++i)
        {
            triangles.push_back(&obj_triangles[i]);
        }
    }

    return triangles;
}


float get_total_scene_area(std::map<objID_t, objData_t> const& meshes)
{
    float scene_area = 0.0f;

    for (std::map<objID_t, objData_t>::const_iterator iter = meshes.begin(); iter != meshes.end(); ++iter)
    {
        triangleObject_t const* obj = iter->second.obj;

        std::vector<triangle_t> const& triangles = obj->getTriangles();

        for (unsigned int i = 0; i < triangles.size(); ++i)
        {
            triangle_t const& tri = triangles[i];

            scene_area += tri.surfaceArea();
        }
    }

    return scene_area;
}


float get_total_scene_area(std::vector<triangle_t const*> const& triangles)
{
    float scene_area = 0.0f;

    for (unsigned int i = 0; i < triangles.size(); ++i)
    {
        triangle_t const& tri = *triangles[i];

        scene_area += tri.surfaceArea();
    }

    return scene_area;
}


std::vector<float> get_triangle_areas_cdf(std::vector<triangle_t const*> const& triangles)
{
    std::vector<float> triangle_areas(triangles.size());

    for (unsigned int i = 0; i < triangles.size(); ++i)
    {
        triangle_areas[i] = triangles[i]->surfaceArea();
    }

    std::vector<float> triangle_areas_cdf(triangle_areas.size());
    triangle_areas_cdf[0] = triangle_areas[0];

    for (unsigned int i = 1; i < triangle_areas_cdf.size(); ++i)
    {
        triangle_areas_cdf[i] = triangle_areas_cdf[i - 1] + triangle_areas[i];
    }

    return triangle_areas_cdf;
}



std::vector<pbgi_sample_t> generate_samples_darts_hash(float const min_radius, int const number_of_samples, std::vector<triangle_t const*> const& triangles)
{
    std::cout << "generate_samples_darts_hash() min_radius: " << min_radius << std::endl;

    random_t my_random;
    std::vector<pbgi_sample_t> sampling_points;

    /*
    for (float x = 0.0f; x < 5.0f; x += 0.1f)
    {
        hash_function(point3d_t(x, my_random(), my_random()), 400, 0.5f);
    }

    return sampling_points;
    */

    std::vector<float> triangle_areas_cdf = get_triangle_areas_cdf(triangles);

    float const scene_area = triangle_areas_cdf.back();


    int const max_tries = 10000;
    int bin_count = 320000;
    float cell_size = min_radius * 2.0f;

    std::vector<std::vector<int> > hash_map(bin_count);

    int rejected_count = 0;

    int debug_rejects = 0;

    while (rejected_count < max_tries)
    {
        float ksi[3] = { my_random(), my_random(), my_random() };

        int triangle_index = std::lower_bound(triangle_areas_cdf.begin(), triangle_areas_cdf.end(), ksi[0] * scene_area) - triangle_areas_cdf.begin();
        // int triangle_index = int(ksi[0] * triangles.size()) % triangles.size();

        vector3d_t sample_normal;
        point3d_t sample_point;

        triangles[triangle_index]->sample(ksi[1], ksi[2], sample_point, sample_normal);

        bool sample_rejected = false;

        int const hash_value = hash_function(sample_point, bin_count, cell_size);
        std::vector<int> const& neighbors = hash_map[hash_value];

        for (unsigned int i = 0; i < neighbors.size(); ++i)
        {
            int sample_index = neighbors[i];

            assert(sample_index < int(sampling_points.size()));

            if ((sampling_points[sample_index].position - sample_point).length() < min_radius * 2.0f)
            {
                sample_rejected = true;
                ++rejected_count;
                ++debug_rejects;
                break;
            }
        }


        /*
        bool sample_rejected_2 = false;

        for (unsigned int i = 0; i < sampling_points.size(); ++i)
        {
            if ((sampling_points[i].position - sample_point).length() < min_radius * 2.0f)
            {
                sample_rejected_2 = true;
                break;
            }
        }

        if (sample_rejected != sample_rejected_2)
        {
            if (!sample_rejected && sample_rejected_2)
            {
                std::cout << "close sample not found with hash" << std::endl;
                assert(false);
            }

            if (sample_rejected && !sample_rejected_2)
            {
                std::cout << "close sample found with hash IMPOSSIBLE ;)" << std::endl;
            }


        }
        */

        if (!sample_rejected)
        {
            int new_sample_index = sampling_points.size();

            for (int x = -1; x <= 1; ++x)
            {
                for (int y = -1; y <= 1; ++y)
                {
                    for (int z = -1; z <= 1; ++z)
                    {
                        int const hash_value = hash_function(sample_point + point3d_t(x * cell_size, y * cell_size, z * cell_size), bin_count, cell_size);
                        hash_map[hash_value].push_back(new_sample_index);
                    }
                }
            }

            pbgi_sample_t sample;
            sample.position = sample_point;
            sample.tri_pointer = triangles[triangle_index];

            sampling_points.push_back(sample);

            if (sampling_points.size() % 10000 == 0)
            {
                std::cout << "generate_samples_darts_hash(): accepted: " << sampling_points.size() << " " << rejected_count << std::endl;
            }

            rejected_count = 0;
        }
    }

    std::cout << "generate_samples_darts_hash() finished" << std::endl;

    std::ofstream hash_file("/tmp/hash_map");

    for (unsigned int i = 0; i < hash_map.size(); ++i)
    {
        hash_file << i << " " << hash_map[i].size() << std::endl;
    }

    return sampling_points;
}


std::vector<pbgi_sample_t> generate_samples_darts(float const min_radius, int const number_of_samples, std::vector<triangle_t const*> const& triangles)
{
    std::cout << "generate_samples_darts() start" << std::endl;

    random_t my_random;

    std::vector<pbgi_sample_t> sampling_points;

    std::vector<float> triangle_areas_cdf = get_triangle_areas_cdf(triangles);

    float const scene_area = triangle_areas_cdf.back();


    int const max_tries = 10000;

    for (int i = 0; i < number_of_samples; ++i)
    {
        bool sample_too_close;
        point3d_t sample_point;
        int triangle_index;

        for (int try_i = 0; try_i < max_tries; ++try_i)
        {
            sample_too_close = false;

            float ksi[3] = { my_random(), my_random(), my_random() };

            triangle_index = std::lower_bound(triangle_areas_cdf.begin(), triangle_areas_cdf.end(), ksi[0] * scene_area) - triangle_areas_cdf.begin();
            // triangle_index = int(ksi[0] * triangles.size()) % triangles.size();

            vector3d_t sample_normal;

            triangles[triangle_index]->sample(ksi[1], ksi[2], sample_point, sample_normal);

            for (unsigned int j = 0; j < sampling_points.size(); ++j)
            {
                if ((sampling_points[j].position - sample_point).length() < min_radius * 2.0f)
                {
                    sample_too_close = true;
                    break;
                }
            }

            if (!sample_too_close) break;
        }

        if (!sample_too_close)
        {
            pbgi_sample_t sample;
            sample.position = sample_point;
            sample.tri_pointer = triangles[triangle_index];

            sampling_points.push_back(sample);
        }
    }

    std::cout << "generate_samples_darts() finish" << std::endl;

    return sampling_points;
}


std::vector<pbgi_sample_t> generate_samples_cdf(int const number_of_samples, std::vector<triangle_t const*> const& triangles)
{
    std::cout << "generate_samples_cdf(): triangles: " << triangles.size() << std::endl;

    std::vector<float> triangle_areas_cdf = get_triangle_areas_cdf(triangles);

    float const area_sum = triangle_areas_cdf.back();

    random_t my_random;

    std::vector<pbgi_sample_t> sampling_points;

    for (int i = 0; i < number_of_samples; ++i)
    {
        // float ksi[3] = { my_random(), my_random(), my_random() };
        vector3d_t ksi = hammersley_3(i, number_of_samples);

        int triangle_index = std::lower_bound(triangle_areas_cdf.begin(), triangle_areas_cdf.end(), ksi[0] * area_sum) - triangle_areas_cdf.begin();
        // int triangle_index = std::lower_bound(triangle_areas_cdf.begin(), triangle_areas_cdf.end(), ksi[0] * (triangle_areas_cdf.size() - 1)) - triangle_areas_cdf.begin();
        // int triangle_index = int(my_random() * triangles.size()) % triangles.size();

        point3d_t sample_point;
        vector3d_t sample_normal;

        triangles[triangle_index]->sample(ksi[1], ksi[2], sample_point, sample_normal);

        pbgi_sample_t sample;
        sample.position = sample_point;
        sample.tri_pointer = triangles[triangle_index];

        sampling_points.push_back(sample);
    }

    return sampling_points;
}


/*
  "random point suicide"

  triangle_areas
  // number_of_samples
  desired_radius
  number_of_samples_2 (candidates)

  samples = classical_sampling(number_of_samples_2, triangle_areas)

  while samples.size > number_of_samples:
     sample = random from samples



  */

std::vector<pbgi_sample_t> generate_samples_suicide(int const number_of_samples, float const desired_radius, std::vector<triangle_t const*> const& triangles)
{
    std::cout << "generate_samples_suicide(): " << number_of_samples << " radius: " << desired_radius << std::endl;

    int number_of_candidates = number_of_samples * 20;

    std::vector<pbgi_sample_t> candidate_samples;
    candidate_samples.reserve(number_of_candidates);

    random_t my_random;



    int bin_count = number_of_candidates * 0.05f;
    std::vector<std::vector<int> > hash_map(bin_count);
    float cell_size = desired_radius * 2.0f;

    int const expected_bin_size = number_of_candidates / bin_count;

    for (unsigned int i = 0; i < hash_map.size(); ++i)
    {
        hash_map[i].reserve(expected_bin_size);
    }

    for (int i = 0; i < number_of_candidates; ++i)
    {
        float ksi[3] = { my_random(), my_random(), my_random() };
        // vector3d_t ksi = hammersley_3(i, number_of_candidates);

        // int triangle_index = std::lower_bound(triangle_areas_cdf.begin(), triangle_areas_cdf.end(), ksi[0] * area_sum) - triangle_areas_cdf.begin();
        // int triangle_index = std::lower_bound(triangle_areas_cdf.begin(), triangle_areas_cdf.end(), ksi[0] * (triangle_areas_cdf.size() - 1)) - triangle_areas_cdf.begin();
        int triangle_index = int(ksi[0] * triangles.size()) % triangles.size();

        point3d_t sample_point;
        vector3d_t sample_normal;

        triangles[triangle_index]->sample(ksi[1], ksi[2], sample_point, sample_normal);

        pbgi_sample_t sample;
        sample.position = sample_point;
        sample.tri_pointer = triangles[triangle_index];

        candidate_samples.push_back(sample);

        for (int x = -1; x <= 1; ++x)
        {
            for (int y = -1; y <= 1; ++y)
            {
                for (int z = -1; z <= 1; ++z)
                {
                    int const hash_value = hash_function(sample.position + point3d_t(x * cell_size, y * cell_size, z * cell_size), bin_count, cell_size);
                    hash_map[hash_value].push_back(i);
                }
            }
        }
    }




    std::vector<bool> dead_samples(candidate_samples.size(), false);

    /*
    for (unsigned int i = 0; i < candidate_samples.size(); ++i)
    {
        point3d_t const& suicide_pos = candidate_samples[i].position;

        if (i % 1000 == 0)
        {
            std::cout << "i: " << i << std::endl;
        }

        // if (dead_samples[i]) continue;



        int const hash_value = hash_function(suicide_pos, bin_count, cell_size);
        std::vector<int> const& neighbors = hash_map[hash_value];

        for (unsigned int j = 0; j < neighbors.size(); ++j)
        {
            int neighbor_index = neighbors[j];

            if (dead_samples[neighbor_index] || i == neighbor_index) continue;

            if ((candidate_samples[neighbor_index].position - suicide_pos).length() < desired_radius * 2.0f)
            {
                dead_samples[i] = true;
                break;
            }

        }
    }
    */


      // killer loop
    for (unsigned int i = 0; i < candidate_samples.size(); ++i)
    {
        point3d_t const& killer_pos = candidate_samples[i].position;

        if (i % 1000 == 0)
        {
            std::cout << "i: " << i << std::endl;
        }

        if (dead_samples[i]) continue;



        int const hash_value = hash_function(killer_pos, bin_count, cell_size);
        std::vector<int> const& neighbors = hash_map[hash_value];

        for (unsigned int j = 0; j < neighbors.size(); ++j)
        {
            int neighbor_index = neighbors[j];

            if (i == neighbor_index) continue;

            if ((candidate_samples[neighbor_index].position - killer_pos).length() < desired_radius * 2.0f)
            {
                dead_samples[neighbor_index] = true;
            }

        }



//        for (unsigned int j = 0; j < candidate_samples.size(); ++j)
//        {
//            // if (dead_samples[j] || i == j) continue;
//            if (i == j) continue;

//            if ((candidate_samples[j].position - candidate_samples[i].position).length() < desired_radius * 2.0f)
//            {
//                dead_samples[j] = true;
//            }
//        }

    }


    std::vector<pbgi_sample_t> survivors;

    for (unsigned int i = 0; i < candidate_samples.size(); ++i)
    {
        if (!dead_samples[i])
        {
            survivors.push_back(candidate_samples[i]);
        }
    }

    std::cout << "generate_samples_suicide(): generated samples: " << survivors.size() << std::endl;

    return survivors;



    /*
    while ()
    {
        std::vector<pbgi_sample_t> survivors;

        pbgi_sample_t const& sample = candidate_samples.front();

        for (unsigned int i = 1; i < candidate_samples.size(); ++i)
        {
            if ((candidate_samples[i].position - sample.position).length() >= desired_radius * 2.0f)
            {
                survivors.push_back(sample);
            }
        }

        candidate_samples = survivors;
        candidate_samples.push_back(sample);
    }
    */
}


void pbLighting_t::generate_gi_points(renderState_t & state)
{
    int node_count = 0;

    std::string fileName = "/tmp/pbgi_points_store";

    std::ofstream fileStream(fileName.c_str());

    float const scene_area = get_total_scene_area(scene->meshes);

    int const number_of_samples = samplesPerArea;
    float const area_per_sample = scene_area / float(number_of_samples);

    std::cout << "generate_gi_points(): scene_area: " << scene_area << " area_per_sample: " << area_per_sample << std::endl;

    float radius = std::sqrt(area_per_sample / M_PI);
    std::cout << "best radius: " << radius << std::endl;

    timer_t my_timer;
    my_timer.addEvent("t1");
    my_timer.start("t1");

    // std::vector<pbgi_sample_t> sampling_points = generate_samples_cdf(number_of_samples, get_scene_triangles(scene->meshes));
    // std::vector<pbgi_sample_t> sampling_points = generate_samples_darts(0.05f, number_of_samples, get_scene_triangles(scene->meshes));
    // std::vector<pbgi_sample_t> sampling_points = generate_samples_darts_hash(radius * 0.5f, number_of_samples, get_scene_triangles(scene->meshes));
    std::vector<pbgi_sample_t> sampling_points = generate_samples_darts_hash(radius, number_of_samples, get_scene_triangles(scene->meshes));
    //std::vector<pbgi_sample_t> sampling_points = generate_samples_suicide(number_of_samples, radius, get_scene_triangles(scene->meshes));

    my_timer.stop("t1");

    std::cout << "samples/desired_samples: " << (sampling_points.size() / float(number_of_samples)) << " time: " << my_timer.getTime("t1") << std::endl;

    float const largest_distance = generate_histogram(sampling_points, radius);
    // float const largest_distance = radius;
    float const needed_radius = largest_distance * std::sqrt(2.0f);
    float const area_estimation = M_PI * needed_radius * needed_radius;

    for (unsigned int i = 0; i < sampling_points.size(); ++i)
    {
        point3d_t const& hitPoint = sampling_points[i].position;
        triangle_t const& tri = *sampling_points[i].tri_pointer;

        surfacePoint_t sp;

        intersectData_t iData;
        tri.getSurface(sp, hitPoint, iData);

        GiPoint * giPoint = new GiPoint();
        giPoint->pos = vector3d_t(sp.P);
        giPoint->normal = sp.N;
        //giPoint->area = area_per_sample;
        giPoint->area = area_estimation;
        //giPoint->debug_radius = needed_radius;
        giPoint->debug_radius = radius;

        material_t const* material = sp.material;
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

        _bspTree->addPoint(giPoint->pos, giPoint);

        fileStream << *giPoint << std::endl;

        ++node_count;
    }

    std::cout << "surfel count: " << node_count << std::endl;

    fileStream.close();

}


yafaray::pbLighting_t::MyTree * load_gi_points()
{
    int node_count = 0;
    std::vector<GiPoint*> points;

    std::string fileName = "/tmp/pbgi_points_store";

    std::ifstream fileStream(fileName.c_str());

    yafaray::point3d_t min(1e10f, 1e10f, 1e10f), max(-1e10f, -1e10f, -1e10f);

    while (fileStream.good())
    {
        GiPoint * p = new GiPoint();
        fileStream >> *p;
        points.push_back(p);
        // tree->addPoint(p->pos, p);

        min.x = std::min(min.x, p->pos.x);
        min.y = std::min(min.y, p->pos.y);
        min.z = std::min(min.z, p->pos.z);

        max.x = std::max(max.x, p->pos.x);
        max.y = std::max(max.y, p->pos.y);
        max.z = std::max(max.z, p->pos.z);

        ++node_count;
    }

    fileStream.close();

    pbLighting_t::MyTree* tree = new pbLighting_t::MyTree(vector3d_t(min), vector3d_t(max), 30, 1);

    for (unsigned int i = 0; i < points.size(); ++i)
    {
        tree->addPoint(points[i]->pos, points[i]);
    }

    std::cout << "load_gi_points(): surfel count: " << node_count << std::endl;

    return tree;
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


    if (do_load_gi_points)
    {
        _bspTree = load_gi_points();
    }
    else
    {
        _bspTree = new MyTree(vector3d_t(sceneBound.a), vector3d_t(sceneBound.g), 30, 1);
        generate_gi_points(state);
    }

    MyTree::averageData(_bspTree);

    // empty out debug file
    std::ofstream file_stream_fb("/tmp/pbgi_frame_buffer");
    file_stream_fb.close();

//    _bspTree->printAveraged();

    // Y_INFO << "PBGI: BSP: " << *_bspTree << std::endl;

    Y_INFO << "PBGI: solid angle: "
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
            if (node->isNodeBehindPlane(vector3d_t(sp.P), sp.N) || !node->getClusteredData()) continue;

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



color_t doPointBasedGiTree_sh_fb(
    pbLighting_t::MyTree const* tree,
    renderState_t & state,
    surfacePoint_t const& sp,
    float const solid_angle_threshold,
    bool const color_by_depth,
    vector3d_t const& wo,
    pbLighting_t::cube_raster_buffer_type * result_fb,
    std::vector<yafaray::GiPoint const*> * gi_points)
{
    color_t col(0.0f);
    const material_t *material = sp.material;

    // traverse tree, if solid angle of node > max, traverse into the children
    std::queue<pbLighting_t::MyTree const*> queue;
    queue.push(tree);

    int shadingNodes = 0;
    int shading_discs_rays = 0;
    int shading_discs_square = 0;
    int shading_discs = 0;

    pbLighting_t::cube_raster_buffer_type frame_buffer;

    std::vector<color_t> colors;
    colors.push_back(color_t(1.0f, 0.0f, 0.0f)); // red
    colors.push_back(color_t(0.0f, 1.0f, 0.0f)); // green
    colors.push_back(color_t(0.0f, 0.0f, 1.0f)); // blue
    colors.push_back(color_t(1.0f, 1.0f, 0.0f)); // yellow
    colors.push_back(color_t(0.0f, 1.0f, 1.0f)); // turquois
    colors.push_back(color_t(1.0f, 0.0f, 1.0f)); // violet
    colors.push_back(color_t(0.5f, 0.5f, 0.5f)); // gray
    colors.push_back(color_t(1.0f, 1.0f, 1.0f)); // white (7)
    colors.push_back(color_t(1.0f, 0.0f, 0.0f)); // red, surfel (8)


    while (!queue.empty())
    {
        pbLighting_t::MyTree const* node = queue.front();
        queue.pop();

        // if we are here and we have a leaf, we sample all points in the leaf
        if (node->getIsLeaf() && node->getData().size() > 0)
        {
            std::vector<GiPoint*> const& points = node->getData();

            // assert(points.size() == 1); // true as long as we force 1 surfel per leaf in the tree

            for (unsigned int i = 0; i < points.size(); ++i)
            {
                yafaray::GiPoint const& giP = *points[i];

                yafaray::vector3d_t giToSp = (vector3d_t(sp.P) - giP.pos);

                float const distance = giToSp.length();

                giToSp.normalize();

                float const cos_sp_gip = std::max(giP.normal * giToSp, 0.0f);

                if (cos_sp_gip <= 0.001f) continue;

                float const visible_area = std::sqrt(cos_sp_gip) * giP.area;
                // float const visible_area = cos_sp_gip * giP.area;
                float const max_area = giP.area;
                // float const area = std::max(0.0f, giP.sh_representation.get_sh_area(giToSp));

                float const solid_angle_real = visible_area / (distance * distance);

                //frame_buffer.add_point(-giToSp, contribution, solidAngle, distance, use_rays);
                float const disc_radius = std::sqrt(max_area / M_PI);
                float const visible_radius = std::sqrt(visible_area / M_PI);

                yafaray::color_t contribution;

                if (color_by_depth)
                {
                    contribution = colors.back(); // * cos_sp_gip;

                    if (distance < disc_radius * 4.0f)
                    {
                        contribution = color_t(1.0f, 1.0f, 0.0f);
                    }
                }
                else
                {
                    // contribution = giP.sh_representation.get_sh_color(giToSp);
                    contribution = giP.color * giP.energy;
//                    contribution = giP.color * giP.energy * cos_sp_gip;
                }


                yafaray::GiPoint * debug_point = NULL;

                if (gi_points)
                {
                    debug_point = new yafaray::GiPoint(giP);

                    //debug_point->area = max_area;
                    debug_point->area = visible_area;
                    debug_point->color = contribution;
                    debug_point->energy = 1.0f;
                    debug_point->depth = node->getDepth();
                    debug_point->debug_radius = visible_radius; // debug only!
                }

                // float cos_normal_gip = -giToSp * sp.N; // lambert receiver

                // col += solidAngle * giP.color * giP.energy * cos_normal_gip;

                //frame_buffer.add_point_from_sphere_with_rays(-giToSp, contribution, radius * 1.2f, distance, debug_point);
                //frame_buffer.add_point_square_rasterization(-giToSp, contribution, solidAngle, distance, debug_point);
                //frame_buffer.add_point_exact(contribution, giP.normal, giP.pos - vector3d_t(sp.P), radius, distance, debug_point);
                // frame_buffer.add_point_exact(contribution, giP.normal, giP.pos, radius, distance, debug_point);
                //++shading_discs;


                //if (cos_sp_gip > 0.001f && distance > radius && distance < radius * 4.0f)
                if (distance < disc_radius * 4.0f)
                {
                    frame_buffer.add_point_exact(contribution, giP.normal, giP.pos - vector3d_t(sp.P), disc_radius, distance, debug_point);
                    ++shading_discs_rays;
                }
                else
                {
                    // frame_buffer.add_point_rays(-giToSp, contribution, solidAngle, distance);
                    frame_buffer.add_point_square_rasterization(-giToSp, contribution, solid_angle_real, distance, debug_point);
                    ++shading_discs_square;
                }
            }
        }
        else if (node->has_children())
        {
            if (node->isNodeBehindPlane(vector3d_t(sp.P), sp.N)) continue;

            GiPoint const& giP = *node->getClusteredData();
            // vector3d_t const position = node->getCenter();
            vector3d_t const position = giP.pos;

            vector3d_t giToSp = (vector3d_t(sp.P) - position);

            float const distance = giToSp.length();

            giToSp.normalize();

//            float const visible_area = std::max(0.0f, giP.sh_representation.get_sh_area(giToSp));
            float const max_visible_area = node->getRadius() * node->getRadius() * M_PI;
            float const visible_area = std::max(0.0f, giP.sh_representation.get_sh_area(giToSp));

            float const max_solid_angle = max_visible_area / (distance * distance);


            if (max_solid_angle > solid_angle_threshold || distance < node->getRadius() * 1.5f)
            {
                std::vector<pbLighting_t::MyTree> const& children = node->getChildren();
                for (unsigned int i = 0; i < children.size(); ++i)
                {
                    queue.push(&children[i]);
                }
            }
            else
            {
                // if (solidAngle < 0.0001f) continue;

                // float const cos_normal_gip = -giToSp * sp.N;

                yafaray::color_t cluster_contribution;

                if (color_by_depth)
                {
                    // float const cos_sp_gip = std::max(0.0f, giP.normal * giToSp);
                    cluster_contribution = colors[std::min(node->getDepth(), int(colors.size()) - 2)]; // * cos_sp_gip;
                }
                else
                {
                    cluster_contribution = giP.sh_representation.get_sh_color(giToSp);
                    // cluster_contribution = giP.color * giP.energy; // * cos_sp_gip;
                }

                // bool const use_rays = true;
                // frame_buffer.add_point(-giToSp, cluster_contribution, solidAngle, distance, use_rays);
                // frame_buffer.add_point_from_sphere_with_rays(-giToSp, cluster_contribution, node->getRadius(), distance, node->getClusteredData());

                yafaray::GiPoint * debug_point = NULL;

                if (gi_points)
                {
                    debug_point = new yafaray::GiPoint(giP);

                    debug_point->pos = position;
                    // debug_point->area = max_visible_area;
                    debug_point->area = visible_area;
                    debug_point->color = cluster_contribution;
                    debug_point->energy = 1.0f;
                    debug_point->depth = node->getDepth();
                    debug_point->debug_radius = std::sqrt(visible_area / M_PI); // debug only!
                }

                // float const solid_angle_real = cos_sp_gip * area / (distance * distance);
                float const real_solid_angle = std::max(0.0f, giP.sh_representation.get_sh_area(giToSp)) / (distance * distance);

                // float cos_normal_gip = -giToSp * sp.N; // lambert receiver
                // col += real_solid_angle * giP.color * giP.energy * cos_normal_gip;

                frame_buffer.add_point_square_rasterization(-giToSp, cluster_contribution, real_solid_angle, distance, debug_point);
                // frame_buffer.add_point_single_cell(-giToSp, cluster_contribution, distance, debug_point);

                ++shadingNodes;
            } // end else (bad solid angle)
        }
    }



    frame_buffer.accumulate(gi_points);

    color_t surfCol(1.0f);
    if (material)
    {
        surfCol = material->getDiffuseAtPoint(state, sp); //material->eval(state, sp, wo, vector3d_t(0.0f), BSDF_ALL);
        // surfCol = material->eval(state, sp, wo, vector3d_t(0.0f), BSDF_ALL);
        col = frame_buffer.get_diffuse(sp.N) * surfCol;
    }

    if (result_fb)
    {
        col = frame_buffer.get_diffuse(sp.N);
        *result_fb = frame_buffer;
    }

    // col = frame_buffer.get_diffuse(sp.N);

    if (gi_points)
    {
        std::cout << "gipoints: " << gi_points << std::endl;
        std::cout << "shading_discs_rays: " << shading_discs_rays
                  << " shading_discs_square: " << shading_discs_square
                  << " shadingNodes: " << shadingNodes
                  << " shading_discs: " << shading_discs
                  << std::endl;
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
		
        if(bsdfs & BSDF_EMIT) col += material->emission(state, sp, wo);
		
		if(bsdfs & BSDF_DIFFUSE)
		{
            if (!indirectOnly)
            {
                col += estimateAllDirectLight(state, sp, wo);
            }

            if (debug_type == Tree)
            {
                col += doPointBasedGiTree(state, sp, wo);
            }
            else if (debug_type == Tree_sh)
            {
                col += doPointBasedGiTreeSH(state, sp, wo);
            }
            else if (debug_type == Tree_sh_fb)
            {
                // col += doPointBasedGiTree_sh_fb(state, sp, wo);
                col += doPointBasedGiTree_sh_fb(_bspTree, state, sp, maxSolidAngle, false, wo, NULL, NULL);
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

    // inte->maxSolidAngle = maxSolidAngle;
    //float angle = std::atan(1.0f / float(pbLighting_t::cube_raster_buffer_type::resolution_2));
    //maxSolidAngle = 2.0f * M_PI * (1.0f - std::cos(angle));
    inte->maxSolidAngle = maxSolidAngle;

    Y_INFO << "maxSolidAngle: " << maxSolidAngle << std::endl;


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
