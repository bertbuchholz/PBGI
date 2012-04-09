#include <cassert>
#include <fstream>
#include <queue>
#include <vector>

#include <yafray_config.h>
#include <core_api/environment.h>
#include <yafraycore/meshtypes.h>
#include <utilities/sample_utils.h>

#include <yafraycore/timer.h>

#include <integrators/Scene_sampler.h>

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

        if (!obj->do_use_for_pbgi()) continue;

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


std::vector<pbgi_sample_t> Scene_sampler_darts_hash::generate_samples(std::vector<triangle_t const*> const& triangles)
{
    std::cout << "Scene_sampler_darts_hash::generate_samples(), min_radius: " << _min_radius << std::endl;

    random_t my_random;
    std::vector<pbgi_sample_t> sampling_points;

    /*
        // debug, hash value test
        for (float x = -3.0f; x < 3.0f; x += 0.3f)
        {
            std::cout << x << " hash: " << hash_function(point3d_t(x, my_random(), my_random()), 320000, 1.0f) << std::endl;
        }
        */

    std::vector<float> triangle_areas_cdf = get_triangle_areas_cdf(triangles);

    float const scene_area = triangle_areas_cdf.back();


    int const max_tries = 10000;
    int bin_count = _number_of_samples * _bin_count_factor;
    float cell_size = _min_radius * _cell_size_factor;

    std::vector<std::vector<int> > hash_map(bin_count);

    int rejected_count = 0;

    int debug_rejects = 0;

    while (rejected_count < max_tries)
    {
        float ksi[3] = { my_random(), my_random(), my_random() };
        // vector3d_t ksi = hammersley_3(i, number_of_samples);

        int triangle_index = std::lower_bound(triangle_areas_cdf.begin(), triangle_areas_cdf.end(), ksi[0] * scene_area) - triangle_areas_cdf.begin();
        // int triangle_index = int(ksi[0] * triangles.size()) % triangles.size();

        vector3d_t sample_normal;
        point3d_t sample_position;
        intersectData_t data;

        triangles[triangle_index]->sample(ksi[1], ksi[2], sample_position, sample_normal, data);

        bool sample_rejected = false;

        int const hash_value = hash_function(sample_position, bin_count, cell_size);
        std::vector<int> const& neighbors = hash_map[hash_value];

        for (unsigned int i = 0; i < neighbors.size(); ++i)
        {
            int neighbor_index = neighbors[i];

            assert(neighbor_index < int(sampling_points.size()));

            float const distance = (sampling_points[neighbor_index].position - sample_position).length();

            if (distance < _min_radius * 2.0f)
                // if ((sampling_points[sample_index].position - sample_point).length() < min_radius)
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
                    std::cout << "close sample not found with hash, accepted: " << sampling_points.size() << std::endl;
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
                        int const hash_value = hash_function(sample_position + point3d_t(x * cell_size, y * cell_size, z * cell_size), bin_count, cell_size);
                        hash_map[hash_value].push_back(new_sample_index);
                    }
                }
            }

            pbgi_sample_t sample;
            sample.position = sample_position;
            sample.tri_pointer = triangles[triangle_index];
            sample.intersect_data = data;

            sampling_points.push_back(sample);

            if (sampling_points.size() % 10000 == 0)
            {
                // std::cout << "generate_samples_darts_hash(): accepted: " << sampling_points.size() << " " << rejected_count << std::endl;
                std::cout << "." << std::flush;
            }

            rejected_count = 0;
        }
    }

    std::cout << "Scene_sampler_darts_hash::generate_samples(), sampling done" << std::endl;

    /*
    std::ofstream hash_file("/tmp/hash_map");

    for (unsigned int i = 0; i < hash_map.size(); ++i)
    {
        hash_file << i << " " << hash_map[i].size() << std::endl;
    }
    */

    /*
    std::vector<float> distances_squared(sampling_points.size());

    // #pragma omp parallel for
    for (int i = 0; i < int(sampling_points.size()); ++i)
    {
        float smallest_distance_i = 1e10f;

        point3d_t const& sampling_point_position = sampling_points[i].position;

        int const hash_value = hash_function(sampling_point_position, bin_count, cell_size);
        std::vector<int> const& neighbors = hash_map[hash_value];

        if (neighbors.size() <= 1) continue;

        for (unsigned int j = 0; j < neighbors.size(); ++j)
        {
            int neighbor_index = neighbors[j];

            if (neighbor_index == int(i)) continue;

            point3d_t const& neighbor_position = sampling_points[neighbor_index].position;

            float const distance = (sampling_point_position - neighbor_position).lengthSqr();

            if (distance < smallest_distance_i)
            {
                smallest_distance_i = distance;
            }
        }

        if (smallest_distance_i > 1e9f)
        {
            for (unsigned int j = 0; j < neighbors.size(); ++j)
            {
                int neighbor_index = neighbors[j];
                float const distance = (sampling_points[i].position - sampling_points[neighbor_index].position).lengthSqr();

                std::cout << "j: " << distance  << std::endl;
            }
        }

        distances_squared[i] = smallest_distance_i;
    }


    _widest_gap = 0.0f;

    for (std::size_t i = 0; i < distances_squared.size(); ++i)
    {
        if (distances_squared[i] > _widest_gap)
        {
            _widest_gap = distances_squared[i];
        }
    }

    _widest_gap = std::sqrt(_widest_gap);

    std::cout << "Scene_sampler_darts_hash::generate_samples(), widest_gap: " << _widest_gap << std::endl;
*/

    return sampling_points;
}





std::vector<pbgi_sample_t> Scene_sampler_darts::generate_samples(std::vector<triangle_t const*> const& triangles)
{
    std::cout << "generate_samples_darts() start" << std::endl;

    random_t my_random;

    std::vector<pbgi_sample_t> sampling_points;

    std::vector<float> triangle_areas_cdf = get_triangle_areas_cdf(triangles);

    float const scene_area = triangle_areas_cdf.back();


    int const max_tries = 10000;

    for (int i = 0; i < _number_of_samples; ++i)
    {
        bool sample_too_close;
        point3d_t sample_point;
        intersectData_t data;
        int triangle_index;

        for (int try_i = 0; try_i < max_tries; ++try_i)
        {
            sample_too_close = false;

            float ksi[3] = { my_random(), my_random(), my_random() };

            triangle_index = std::lower_bound(triangle_areas_cdf.begin(), triangle_areas_cdf.end(), ksi[0] * scene_area) - triangle_areas_cdf.begin();
            // triangle_index = int(ksi[0] * triangles.size()) % triangles.size();

            vector3d_t sample_normal;

            triangles[triangle_index]->sample(ksi[1], ksi[2], sample_point, sample_normal, data);

            for (unsigned int j = 0; j < sampling_points.size(); ++j)
            {
                if ((sampling_points[j].position - sample_point).length() < _min_radius * 2.0f)
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
            sample.intersect_data = data;
            sample.tri_pointer = triangles[triangle_index];

            sampling_points.push_back(sample);
        }
    }

    std::cout << "generate_samples_darts() finish" << std::endl;

    return sampling_points;
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

std::vector<pbgi_sample_t> Scene_sampler_suicide::generate_samples(std::vector<triangle_t const*> const& triangles)
{
    std::cout << "generate_samples_suicide(): " << _number_of_samples << " radius: " << _desired_radius << std::endl;

    int number_of_candidates = _number_of_samples * 20;

    std::vector<pbgi_sample_t> candidate_samples;
    candidate_samples.reserve(number_of_candidates);

    random_t my_random;



    int bin_count = number_of_candidates * 0.05f;
    std::vector<std::vector<int> > hash_map(bin_count);
    float cell_size = _desired_radius * 2.0f;

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
        intersectData_t data;

        triangles[triangle_index]->sample(ksi[1], ksi[2], sample_point, sample_normal, data);

        pbgi_sample_t sample;
        sample.position = sample_point;
        sample.tri_pointer = triangles[triangle_index];
        sample.intersect_data = data;

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

            if (int(i) == neighbor_index) continue;

            if ((candidate_samples[neighbor_index].position - killer_pos).length() < _desired_radius * 2.0f)
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


// http://en.wikipedia.org/wiki/Solid_angle#Tetrahedron
float triangle_solid_angle(triangle_t const& triangle, vector3d_t const& origin)
{
    float omega = 0.0f;

    vector3d_t a = vector3d_t(triangle.getVertex(0)) - origin;
    vector3d_t b = vector3d_t(triangle.getVertex(1)) - origin;
    vector3d_t c = vector3d_t(triangle.getVertex(2)) - origin;

    float const determ = std::abs(a * (b ^ c));

    float const al = a.length();
    float const bl = b.length();
    float const cl = c.length();

    float const div = al * bl * cl + (a * b) * cl + (a * c) * bl + (b * c) * al;
    float at = std::atan2(determ, div);

    if (at < 0) at += M_PI; // If det > 0 and div < 0 arctan2 returns < 0, so add pi.
    omega = 2.0f * at;

    return omega;
}


std::vector<float> triangle_angles(triangle_t const& triangle)
{
    std::vector<float> angles;

    {
        vector3d_t v0 = vector3d_t(triangle.getVertex(1) - triangle.getVertex(0));
        v0.normalize();
        vector3d_t v1 = vector3d_t(triangle.getVertex(2) - triangle.getVertex(0));
        v1.normalize();

        angles.push_back(std::acos(std::max(-1.0f, std::min(1.0f, v0 * v1))));
    }

    {
        vector3d_t v0 = vector3d_t(triangle.getVertex(0) - triangle.getVertex(1));
        v0.normalize();
        vector3d_t v1 = vector3d_t(triangle.getVertex(2) - triangle.getVertex(1));
        v1.normalize();

        angles.push_back(std::acos(std::max(-1.0f, std::min(1.0f, v0 * v1))));
    }

    {
        vector3d_t v0 = vector3d_t(triangle.getVertex(0) - triangle.getVertex(2));
        v0.normalize();
        vector3d_t v1 = vector3d_t(triangle.getVertex(1) - triangle.getVertex(2));
        v1.normalize();

        angles.push_back(std::acos(std::max(-1.0f, std::min(1.0f, v0 * v1))));
    }

    return angles;
}


int get_longest_index(std::vector<vector3d_t> const& edges)
{
    float max_length = 0;
    int longest_index = -1;

    for (std::size_t i = 0; i < edges.size(); ++i)
    {
        float const tmp_length = edges[i].length();

        if (tmp_length > max_length)
        {
            max_length = tmp_length;
            longest_index = i;
        }
    }

    return longest_index;
}


int get_shortest_index(std::vector<vector3d_t> const& edges)
{
    float min_length = 1e10f;
    int shortest_index = -1;

    for (std::size_t i = 0; i < edges.size(); ++i)
    {
        float const tmp_length = edges[i].length();

        if (tmp_length < min_length)
        {
            min_length = tmp_length;
            shortest_index = i;
        }
    }

    return shortest_index;
}



std::vector<pbgi_sample_t> Scene_sampler_reyes::generate_samples(std::vector<triangle_t const*> const& /* triangles */)
{
    std::vector<pbgi_sample_t> samples;

    for (std::map<objID_t, objData_t>::const_iterator iter = _meshes.begin(); iter != _meshes.end(); ++iter)
    {
        triangleObject_t * obj = new triangleObject_t(*iter->second.obj);


        std::vector<triangle_t> const& triangles = obj->getTriangles();

        std::vector<int> triangles_for_sampling;

        std::vector<triangle_t> queue;

        for (std::size_t i = 0; i < triangles.size(); ++i)
        {
            queue.push_back(triangles[i]);
        }

        while (!queue.empty())
        {
            triangle_t const triangle = queue.back();
            queue.pop_back();

            // float const solid_angle = tetrahedron_solid_angle(tri);
            vector3d_t const center = vector3d_t(triangle.getVertex(0) + triangle.getVertex(1) + triangle.getVertex(2)) / 3.0f;
            float const distance = (center - _camera_pos).length();
            float const solid_angle = triangle.surfaceArea() / (distance * distance);

            if (solid_angle > _max_solid_angle)
            {
                int index_a, index_b, index_c;
                triangle.getVertexIndices(index_a, index_b, index_c);

                int uv0, uv1, uv2;
                uv0 = triangle.get_uv_index(0);
                uv1 = triangle.get_uv_index(1);
                uv2 = triangle.get_uv_index(2);

                std::vector<vector3d_t> edges;
                edges.push_back(triangle.getVertex(1) - triangle.getVertex(0));
                edges.push_back(triangle.getVertex(2) - triangle.getVertex(1));
                edges.push_back(triangle.getVertex(0) - triangle.getVertex(2));

                int const longest_edge_index = get_longest_index(edges);
                int const shortest_edge_index = get_shortest_index(edges);
                int const middle_edge_index = ((longest_edge_index + 1) % 3) == shortest_edge_index ? (longest_edge_index + 2) % 3 : (longest_edge_index + 1) % 3;

                bool short_middle_long_order = ((longest_edge_index + 1) % 3) == shortest_edge_index ? true : false;

                assert(longest_edge_index != shortest_edge_index);
                assert(longest_edge_index != middle_edge_index);
                assert(middle_edge_index != shortest_edge_index);

                // int const length_ratio = edges[longest_edge_index].length() / edges[shortest_edge_index].length();
                // float const length_ratio = edges[longest_edge_index].length() / float(edges[shortest_edge_index].length());
                float const length_ratio = edges[middle_edge_index].length() / float(edges[shortest_edge_index].length());


                bound_t bound = triangle.getBound();

                if (length_ratio > 1.5f)
                {
                    float offset = 1.0f / float(length_ratio);

                    assert(!std::isnan(offset));
                    assert(!std::isinf(offset));

                    assert(offset < 1.0f && offset > 0.0f);

                    if (short_middle_long_order)
                    {
                        vector3d_t new_point_0 = vector3d_t(triangle.getVertex(middle_edge_index))  + edges[middle_edge_index]  * offset;
                        vector3d_t new_point_1 = vector3d_t(triangle.getVertex(longest_edge_index)) + edges[longest_edge_index] * (1.0f - offset);

                        assert(bound.includes(new_point_0));
                        assert(bound.includes(new_point_1));

                        int const new_index_0 = obj->add_point(new_point_0);
                        int const new_index_1 = obj->add_point(new_point_1);

                        triangle_t triangle_new_0(triangle.get_vertex_index(middle_edge_index), new_index_1, triangle.get_vertex_index(shortest_edge_index), obj);
                        triangle_new_0.setNormal(triangle.getNormal());
                        triangle_new_0.setMaterial(triangle.getMaterial());

                        //triangle_t triangle_new_1(new_index_1, new_index_0, triangle.get_vertex_index(shortest_edge_index), obj);
                        triangle_t triangle_new_1(triangle.get_vertex_index(middle_edge_index), new_index_0, new_index_1, obj);
                        triangle_new_1.setNormal(triangle.getNormal());
                        triangle_new_1.setMaterial(triangle.getMaterial());

                        obj->add_triangle_with_uv_indices(triangle_new_0, uv0, uv1, uv2);
                        triangle_new_0.setIndex(triangles.size() - 1);

                        obj->add_triangle_with_uv_indices(triangle_new_1, uv0, uv1, uv2);
                        triangle_new_1.setIndex(triangles.size() - 1);

                        queue.push_back(triangle_new_0);
                        queue.push_back(triangle_new_1);

                        triangle_t triangle_final(new_index_0, triangle.get_vertex_index(longest_edge_index), new_index_1, obj);
                        triangle_final.setNormal(triangle.getNormal());
                        triangle_final.setMaterial(triangle.getMaterial());

                        obj->add_triangle_with_uv_indices(triangle_final, uv0, uv1, uv2);
                        triangle_final.setIndex(triangles.size() - 1);

                        queue.push_back(triangle_final);
                    }
                    else
                    {
                        vector3d_t new_point_0 = vector3d_t(triangle.getVertex(longest_edge_index)) + edges[longest_edge_index] * offset;
                        vector3d_t new_point_1 = vector3d_t(triangle.getVertex(middle_edge_index))  + edges[middle_edge_index]  * (1.0f - offset);

                        assert(bound.includes(new_point_0));
                        assert(bound.includes(new_point_1));

                        int const new_index_0 = obj->add_point(new_point_0);
                        int const new_index_1 = obj->add_point(new_point_1);

                        triangle_t triangle_new_0(triangle.get_vertex_index(shortest_edge_index), triangle.get_vertex_index(longest_edge_index), new_index_0, obj);
                        triangle_new_0.setNormal(triangle.getNormal());
                        triangle_new_0.setMaterial(triangle.getMaterial());

                        triangle_t triangle_new_1(new_index_0, new_index_1, triangle.get_vertex_index(shortest_edge_index), obj);
                        triangle_new_1.setNormal(triangle.getNormal());
                        triangle_new_1.setMaterial(triangle.getMaterial());

                        obj->add_triangle_with_uv_indices(triangle_new_0, uv0, uv1, uv2);
                        triangle_new_0.setIndex(triangles.size() - 1);

                        obj->add_triangle_with_uv_indices(triangle_new_1, uv0, uv1, uv2);
                        triangle_new_1.setIndex(triangles.size() - 1);

                        queue.push_back(triangle_new_0);
                        queue.push_back(triangle_new_1);

                        triangle_t triangle_final(triangle.get_vertex_index(middle_edge_index), new_index_1, new_index_0, obj);
                        triangle_final.setNormal(triangle.getNormal());
                        triangle_final.setMaterial(triangle.getMaterial());

                        obj->add_triangle_with_uv_indices(triangle_final, uv0, uv1, uv2);
                        triangle_final.setIndex(triangles.size() - 1);

                        queue.push_back(triangle_final);
                    }
                }
                else
                {
                    // subdivide and push new triangles into queue
                    int index_a, index_b, index_c;
                    triangle.getVertexIndices(index_a, index_b, index_c);

                    int uv0, uv1, uv2;
                    uv0 = triangle.get_uv_index(0);
                    uv1 = triangle.get_uv_index(1);
                    uv2 = triangle.get_uv_index(2);

                    vector3d_t pos_x = vector3d_t(triangle.getVertex(0) + triangle.getVertex(1)) / 2.0f;
                    vector3d_t pos_y = vector3d_t(triangle.getVertex(1) + triangle.getVertex(2)) / 2.0f;
                    vector3d_t pos_z = vector3d_t(triangle.getVertex(2) + triangle.getVertex(0)) / 2.0f;

                    int index_x = obj->add_point(pos_x);
                    int index_y = obj->add_point(pos_y);
                    int index_z = obj->add_point(pos_z);

                    triangle_t triangle_axz(index_a, index_x, index_z, obj);
                    triangle_axz.setNormal(triangle.getNormal());
                    triangle_axz.setMaterial(triangle.getMaterial());

                    triangle_t triangle_xby(index_x, index_b, index_y, obj);
                    triangle_xby.setNormal(triangle.getNormal());
                    triangle_xby.setMaterial(triangle.getMaterial());

                    triangle_t triangle_xyz(index_x, index_y, index_z, obj);
                    triangle_xyz.setNormal(triangle.getNormal());
                    triangle_xyz.setMaterial(triangle.getMaterial());

                    triangle_t triangle_ycz(index_y, index_c, index_z, obj);
                    triangle_ycz.setNormal(triangle.getNormal());
                    triangle_ycz.setMaterial(triangle.getMaterial());

                    obj->add_triangle_with_uv_indices(triangle_axz, uv0, uv1, uv2);
                    triangle_axz.setIndex(triangles.size() - 1);

                    obj->add_triangle_with_uv_indices(triangle_xby, uv0, uv1, uv2);
                    triangle_xby.setIndex(triangles.size() - 1);

                    obj->add_triangle_with_uv_indices(triangle_xyz, uv0, uv1, uv2);
                    triangle_xyz.setIndex(triangles.size() - 1);

                    obj->add_triangle_with_uv_indices(triangle_ycz, uv0, uv1, uv2);
                    triangle_ycz.setIndex(triangles.size() - 1);

                    queue.push_back(triangle_axz);
                    queue.push_back(triangle_xby);
                    queue.push_back(triangle_xyz);
                    queue.push_back(triangle_ycz);
                }
            }
            else
            {
                triangles_for_sampling.push_back(triangle.getIndex());
                _debug_new_triangles->push_back(vector3d_t(triangle.getVertex(0)));
                _debug_new_triangles->push_back(vector3d_t(triangle.getVertex(1)));
                _debug_new_triangles->push_back(vector3d_t(triangle.getVertex(2)));
            }
        }

        for (std::size_t i = 0; i < triangles_for_sampling.size(); ++i)
        {
            // generate sample

            triangle_t const& triangle = triangles[triangles_for_sampling[i]];

            std::vector<float> angles = triangle_angles(triangle);

            intersectData_t data;

            data.b0 = angles[0] / M_PI;
            data.b1 = angles[1] / M_PI;
            data.b2 = angles[2] / M_PI;
//            data.b0 = 1.0f / 3.0f;
//            data.b1 = 1.0f / 3.0f;
//            data.b2 = 1.0f / 3.0f;


            vector3d_t const center = vector3d_t(triangle.getVertex(0) * data.b0 +
                                                 triangle.getVertex(1) * data.b1 +
                                                 triangle.getVertex(2) * data.b2);

            assert(!std::isnan(center.x));

            pbgi_sample_t sample;
            sample.position = center;
            sample.intersect_data = data;
            sample.area = triangle.surfaceArea();
            sample.tri_pointer = &triangle;

            samples.push_back(sample);
        }
    }

    return samples;
}



void tesselate_object(triangleObject_t * obj)
{
    std::vector<triangle_t> const& triangles = obj->getTriangles();

    std::vector<triangle_t> new_triangles;

    std::vector<triangle_t> queue;

    for (std::size_t i = 0; i < triangles.size(); ++i)
    {
        queue.push_back(triangles[i]);
    }

    while (!queue.empty())
    {
        triangle_t const triangle = queue.back();
        queue.pop_back();


        int uv0, uv1, uv2;
        uv0 = triangle.get_uv_index(0);
        uv1 = triangle.get_uv_index(1);
        uv2 = triangle.get_uv_index(2);

        std::vector<vector3d_t> edges;
        edges.push_back(triangle.getVertex(1) - triangle.getVertex(0));
        edges.push_back(triangle.getVertex(2) - triangle.getVertex(1));
        edges.push_back(triangle.getVertex(0) - triangle.getVertex(2));

        int const longest_edge_index = get_longest_index(edges);
        int const shortest_edge_index = get_shortest_index(edges);
        int const middle_edge_index = ((longest_edge_index + 1) % 3) == shortest_edge_index ? (longest_edge_index + 2) % 3 : (longest_edge_index + 1) % 3;

        bool short_middle_long_order = ((longest_edge_index + 1) % 3) == shortest_edge_index ? true : false;

        assert(longest_edge_index != shortest_edge_index);
        assert(longest_edge_index != middle_edge_index);
        assert(middle_edge_index != shortest_edge_index);

        // int const length_ratio = edges[longest_edge_index].length() / edges[shortest_edge_index].length();
        // float const length_ratio = edges[longest_edge_index].length() / float(edges[shortest_edge_index].length());
        float const length_ratio = edges[middle_edge_index].length() / float(edges[shortest_edge_index].length());


        bound_t bound = triangle.getBound();

        if (length_ratio > 1.5f && triangle.surfaceArea() > 0.001f)
        {
            float offset = 1.0f / float(length_ratio);

            assert(!std::isnan(offset));
            assert(!std::isinf(offset));

            assert(offset < 1.0f && offset > 0.0f);

            if (short_middle_long_order)
            {
                vector3d_t new_point_0 = vector3d_t(triangle.getVertex(middle_edge_index))  + edges[middle_edge_index]  * offset;
                vector3d_t new_point_1 = vector3d_t(triangle.getVertex(longest_edge_index)) + edges[longest_edge_index] * (1.0f - offset);

                assert(bound.includes(new_point_0));
                assert(bound.includes(new_point_1));

                int const new_index_0 = obj->add_point(new_point_0);
                int const new_index_1 = obj->add_point(new_point_1);

                triangle_t triangle_new_0(triangle.get_vertex_index(middle_edge_index), new_index_1, triangle.get_vertex_index(shortest_edge_index), obj);
                triangle_new_0.setNormal(triangle.getNormal());
                triangle_new_0.setMaterial(triangle.getMaterial());

                //triangle_t triangle_new_1(new_index_1, new_index_0, triangle.get_vertex_index(shortest_edge_index), obj);
                triangle_t triangle_new_1(triangle.get_vertex_index(middle_edge_index), new_index_0, new_index_1, obj);
                triangle_new_1.setNormal(triangle.getNormal());
                triangle_new_1.setMaterial(triangle.getMaterial());

                obj->add_triangle_with_uv_indices(triangle_new_0, uv0, uv1, uv2);
                triangle_new_0.setIndex(triangles.size() - 1);

                obj->add_triangle_with_uv_indices(triangle_new_1, uv0, uv1, uv2);
                triangle_new_1.setIndex(triangles.size() - 1);

                queue.push_back(triangle_new_0);
                queue.push_back(triangle_new_1);

                triangle_t triangle_final(new_index_0, triangle.get_vertex_index(longest_edge_index), new_index_1, obj);
                triangle_final.setNormal(triangle.getNormal());
                triangle_final.setMaterial(triangle.getMaterial());

                obj->add_triangle_with_uv_indices(triangle_final, uv0, uv1, uv2);
                triangle_final.setIndex(triangles.size() - 1);

                queue.push_back(triangle_final);
            }
            else
            {
                vector3d_t new_point_0 = vector3d_t(triangle.getVertex(longest_edge_index)) + edges[longest_edge_index] * offset;
                vector3d_t new_point_1 = vector3d_t(triangle.getVertex(middle_edge_index))  + edges[middle_edge_index]  * (1.0f - offset);

                assert(bound.includes(new_point_0));
                assert(bound.includes(new_point_1));

                int const new_index_0 = obj->add_point(new_point_0);
                int const new_index_1 = obj->add_point(new_point_1);

                triangle_t triangle_new_0(triangle.get_vertex_index(shortest_edge_index), triangle.get_vertex_index(longest_edge_index), new_index_0, obj);
                triangle_new_0.setNormal(triangle.getNormal());
                triangle_new_0.setMaterial(triangle.getMaterial());

                triangle_t triangle_new_1(new_index_0, new_index_1, triangle.get_vertex_index(shortest_edge_index), obj);
                triangle_new_1.setNormal(triangle.getNormal());
                triangle_new_1.setMaterial(triangle.getMaterial());

                obj->add_triangle_with_uv_indices(triangle_new_0, uv0, uv1, uv2);
                triangle_new_0.setIndex(triangles.size() - 1);

                obj->add_triangle_with_uv_indices(triangle_new_1, uv0, uv1, uv2);
                triangle_new_1.setIndex(triangles.size() - 1);

                queue.push_back(triangle_new_0);
                queue.push_back(triangle_new_1);

                triangle_t triangle_final(triangle.get_vertex_index(middle_edge_index), new_index_1, new_index_0, obj);
                triangle_final.setNormal(triangle.getNormal());
                triangle_final.setMaterial(triangle.getMaterial());

                obj->add_triangle_with_uv_indices(triangle_final, uv0, uv1, uv2);
                triangle_final.setIndex(triangles.size() - 1);

                queue.push_back(triangle_final);
            }
        }
        else
        {
            // put into good triangles list
            new_triangles.push_back(triangle);
        }
    }

    for (std::size_t i = 0; i < new_triangles.size(); ++i)
    {
        new_triangles[i].setIndex(i);
    }

    obj->setTriangles(triangles);
}



//std::vector<pbgi_sample_t> Scene_sampler_reyes::generate_samples(std::vector<triangle_t const*> const& /* triangles */)
//{
//    std::vector<pbgi_sample_t> samples;

//    for (std::map<objID_t, objData_t>::const_iterator iter = _meshes.begin(); iter != _meshes.end(); ++iter)
//    {
//        triangleObject_t * obj = new triangleObject_t(*iter->second.obj);

//        tesselate_object(obj);


//        std::vector<triangle_t> const& triangles = obj->getTriangles();

//        std::vector<int> triangles_for_sampling;

//        std::vector<triangle_t> queue;

//        for (std::size_t i = 0; i < triangles.size(); ++i)
//        {
//            queue.push_back(triangles[i]);
//        }

//        while (!queue.empty())
//        {
//            triangle_t const triangle = queue.back();
//            queue.pop_back();

//            // float const solid_angle = tetrahedron_solid_angle(tri);
//            vector3d_t const center = vector3d_t(triangle.getVertex(0) + triangle.getVertex(1) + triangle.getVertex(2)) / 3.0f;
//            float const distance = (center - _camera_pos).length();
//            float const solid_angle = triangle.surfaceArea() / (distance * distance);

//            if (solid_angle > _max_solid_angle)
//            {
//                // subdivide and push new triangles into queue
//                int index_a, index_b, index_c;
//                triangle.getVertexIndices(index_a, index_b, index_c);

//                int uv0, uv1, uv2;
//                uv0 = triangle.get_uv_index(0);
//                uv1 = triangle.get_uv_index(1);
//                uv2 = triangle.get_uv_index(2);

//                vector3d_t pos_x = vector3d_t(triangle.getVertex(0) + triangle.getVertex(1)) / 2.0f;
//                vector3d_t pos_y = vector3d_t(triangle.getVertex(1) + triangle.getVertex(2)) / 2.0f;
//                vector3d_t pos_z = vector3d_t(triangle.getVertex(2) + triangle.getVertex(0)) / 2.0f;

//                int index_x = obj->add_point(pos_x);
//                int index_y = obj->add_point(pos_y);
//                int index_z = obj->add_point(pos_z);

//                triangle_t triangle_axz(index_a, index_x, index_z, obj);
//                triangle_axz.setNormal(triangle.getNormal());
//                triangle_axz.setMaterial(triangle.getMaterial());

//                triangle_t triangle_xby(index_x, index_b, index_y, obj);
//                triangle_xby.setNormal(triangle.getNormal());
//                triangle_xby.setMaterial(triangle.getMaterial());

//                triangle_t triangle_xyz(index_x, index_y, index_z, obj);
//                triangle_xyz.setNormal(triangle.getNormal());
//                triangle_xyz.setMaterial(triangle.getMaterial());

//                triangle_t triangle_ycz(index_y, index_c, index_z, obj);
//                triangle_ycz.setNormal(triangle.getNormal());
//                triangle_ycz.setMaterial(triangle.getMaterial());

//                obj->add_triangle_with_uv_indices(triangle_axz, uv0, uv1, uv2);
//                triangle_axz.setIndex(triangles.size() - 1);

//                obj->add_triangle_with_uv_indices(triangle_xby, uv0, uv1, uv2);
//                triangle_xby.setIndex(triangles.size() - 1);

//                obj->add_triangle_with_uv_indices(triangle_xyz, uv0, uv1, uv2);
//                triangle_xyz.setIndex(triangles.size() - 1);

//                obj->add_triangle_with_uv_indices(triangle_ycz, uv0, uv1, uv2);
//                triangle_ycz.setIndex(triangles.size() - 1);

//                queue.push_back(triangle_axz);
//                queue.push_back(triangle_xby);
//                queue.push_back(triangle_xyz);
//                queue.push_back(triangle_ycz);
//            }
//            else
//            {
//                triangles_for_sampling.push_back(triangle.getIndex());
//                _debug_new_triangles->push_back(vector3d_t(triangle.getVertex(0)));
//                _debug_new_triangles->push_back(vector3d_t(triangle.getVertex(1)));
//                _debug_new_triangles->push_back(vector3d_t(triangle.getVertex(2)));
//            }
//        }

//        for (std::size_t i = 0; i < triangles_for_sampling.size(); ++i)
//        {
//            // generate sample

//            triangle_t const& triangle = triangles[triangles_for_sampling[i]];

//            std::vector<float> angles = triangle_angles(triangle);

//            intersectData_t data;

//            data.b0 = angles[0] / M_PI;
//            data.b1 = angles[1] / M_PI;
//            data.b2 = angles[2] / M_PI;
////            data.b0 = 1.0f / 3.0f;
////            data.b1 = 1.0f / 3.0f;
////            data.b2 = 1.0f / 3.0f;


//            vector3d_t const center = vector3d_t(triangle.getVertex(0) * data.b0 +
//                                                 triangle.getVertex(1) * data.b1 +
//                                                 triangle.getVertex(2) * data.b2);

//            assert(!std::isnan(center.x));

//            pbgi_sample_t sample;
//            sample.position = center;
//            sample.intersect_data = data;
//            sample.area = triangle.surfaceArea();
//            sample.tri_pointer = &triangle;

//            samples.push_back(sample);
//        }
//    }

//    return samples;
//}



__END_YAFRAY
