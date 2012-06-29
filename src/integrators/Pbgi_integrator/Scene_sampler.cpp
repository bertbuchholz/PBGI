#include <cassert>
#include <fstream>
#include <queue>
#include <vector>

#include <yafray_config.h>
#include <core_api/environment.h>
#include <yafraycore/meshtypes.h>
#include <utilities/sample_utils.h>

#include <yafraycore/timer.h>

#include <integrators/Pbgi_integrator/Scene_sampler.h>

__BEGIN_YAFRAY

unsigned int hash_function(point3d_t const& sample, const int n, const float cell_size)
{
    /*
    hash(x,y,z) = ( x p1 xor y p2 xor z p3) mod n
    where p1, p2, p3 are large prime numbers, in
    our case 73856093, 19349663, 83492791
    */

    // float bb_size = 6.0f;

    unsigned int i_x = std::floor(sample.x / cell_size);
    unsigned int i_y = std::floor(sample.y / cell_size);
    unsigned int i_z = std::floor(sample.z / cell_size);

    // std::cout << "i_x: " << (i_x % n) << " sample.x: " << sample.x << std::endl;

    unsigned int hash_value = ((i_x * 73856093u) ^ (i_y * 19349663u) ^ (i_z * 83492791u)) % unsigned(n);
    // unsigned int hash_value = (i_x + i_y + i_z) % n;

    return hash_value;
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

    #pragma omp parallel for
    for (int i = 0; i < int(samples.size()); ++i)
    {
        float smallest_distance_i = 1e10f;

        int const hash_value = hash_function(samples[i].position, bin_count, cell_size);
        std::vector<int> const& neighbors = hash_map[hash_value];

        for (unsigned int j = 0; j < neighbors.size(); ++j)
        {
            int sample_index = neighbors[j];

            if (sample_index == int(i)) continue;

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

    bool do_output_histogram = false;

    if (do_output_histogram)
    {
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

        // if (!obj->do_use_for_pbgi()) continue; // FIXME: add back in some way?

        std::vector<triangle_t> const& obj_triangles = obj->getTriangles();

        for (unsigned int i = 0; i < obj_triangles.size(); ++i)
        {
            triangles.push_back(&obj_triangles[i]);
        }
    }

    return triangles;
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




std::vector<pbgi_sample_t> Scene_sampler_cdf::generate_samples(std::vector<triangle_t const*> const& triangles)
{
    std::cout << "generate_samples_cdf(): triangles: " << triangles.size() << std::endl;

    std::vector<float> triangle_areas_cdf = get_triangle_areas_cdf(triangles);

    float const area_sum = triangle_areas_cdf.back();

    random_t my_random;

    std::vector<pbgi_sample_t> sampling_points;

    for (int i = 0; i < _number_of_samples; ++i)
    {
        // float ksi[3] = { my_random(), my_random(), my_random() };
        vector3d_t ksi = hammersley_3(i, _number_of_samples);

        int triangle_index = std::lower_bound(triangle_areas_cdf.begin(), triangle_areas_cdf.end(), ksi[0] * area_sum) - triangle_areas_cdf.begin();
        // int triangle_index = std::lower_bound(triangle_areas_cdf.begin(), triangle_areas_cdf.end(), ksi[0] * (triangle_areas_cdf.size() - 1)) - triangle_areas_cdf.begin();
        // int triangle_index = int(my_random() * triangles.size()) % triangles.size();

        point3d_t sample_point;
        vector3d_t sample_normal;
        intersectData_t data;

        triangles[triangle_index]->sample(ksi[1], ksi[2], sample_point, sample_normal, data);

        pbgi_sample_t sample;
        sample.position = sample_point;
        sample.tri_pointer = triangles[triangle_index];
        sample.intersect_data = data;

        sampling_points.push_back(sample);
    }

    return sampling_points;
}





__END_YAFRAY
