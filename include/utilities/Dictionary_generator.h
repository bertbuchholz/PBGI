#ifndef DICTIONARY_GENERATOR_H
#define DICTIONARY_GENERATOR_H

#include <ANN/ANN.h>

#include <utilities/spherical_harmonics.h>
#include <utilities/ANN_wrapper_functions.h>

__BEGIN_YAFRAY

class Dictionary_generator
{
    public:
    virtual std::vector<Word> generate(std::vector<Word> const& words,
                               int const dict_num_centers) const = 0;

    virtual std::string identify() const = 0;

    virtual std::string get_stat_results() const { return "No stats."; }
};



class Uniform_dictionary_generator : public Dictionary_generator
{
public:
    template <class T>
    bool contains(std::vector<T> v, T const& val) const
    {
        for (size_t i = 0; i < v.size(); ++i)
        {
            if (v[i] == val)
            {
                return true;
            }
        }

        return false;
    }

    std::vector<Word> generate(std::vector<Word> const& /* words */,
                               int const dict_num_centers) const
    {
        random_t my_random;

        std::vector<Word> dictionary;

        int const dimensions = 27;

        float base = 2.0f;
        float arg  = dict_num_centers;
        int const num_axis_with_2_samples = ceil(std::log(arg) / std::log(base));

        std::vector<int> num_samples_per_dimension(27, 1);

        for (int i = 0; i < num_axis_with_2_samples; )
        {
            int const r = my_random() * dimensions;

            if (num_samples_per_dimension[r] < 2)
            {
                num_samples_per_dimension[r] += 1;
                ++i;
            }
        }

        int const final_num_samples = std::pow(2, num_axis_with_2_samples);

        std::vector<int> used_num_samples_per_dimension = num_samples_per_dimension;

        for (std::size_t i = 0; i < used_num_samples_per_dimension.size(); ++i)
        {
            if (used_num_samples_per_dimension[i] == 1) used_num_samples_per_dimension[i] = 0;
        }


        for (int i = 0; i < final_num_samples; ++i)
        {
            Word w(dimensions);

            bool used_a_dimension = false;

            for (int i = 0; i < dimensions; ++i)
            {
                int const num_samples = num_samples_per_dimension[i];
                int const used_num_samples = used_num_samples_per_dimension[i];

                if (num_samples == 1)
                {
                    w[i] = 0.5f;
                }
                else
                {
                    if (used_num_samples > 0 && !used_a_dimension)
                    {
                        w[i] = 1.0f / (num_samples + 1) * used_num_samples;
                    }
                    else
                    {
                        w[i] = 1.0f / (num_samples + 1);
                    }
                }
            }

            dictionary.push_back(w);
        }



        return dictionary;
    }

    std::string identify() const
    {
        return "Uniform_dictionary";
    }
};




class Random_dictionary_generator : public Dictionary_generator
{
public:
    std::vector<Word> generate(std::vector<Word> const& words,
                               int const dict_num_centers) const
    {
        random_t my_random;

        std::vector<Word> dictionary;
        std::vector<int> available_indices;

        for (std::size_t i = 0; i < words.size(); ++i)
        {
            available_indices.push_back(i);
        }

        for (int i = 0; i < dict_num_centers; ++i)
        {
            int const index = available_indices[my_random() * available_indices.size()];
            std::swap(available_indices[index], available_indices.back());
            available_indices.pop_back();

            dictionary.push_back(words[index]);
        }

        return dictionary;
    }

    std::string identify() const
    {
        return "Random_dictionary";
    }
};



class Kmeans_dictionary_generator : public Dictionary_generator
{
public:
    Kmeans_dictionary_generator(bool const use_ann, bool const do_stats) : _use_ann(use_ann), _do_stats(do_stats)
    { }


    int find_closest_center_linear(Word const& word, std::vector<Word> const& centers, Distance_function const& dist_fnc) const
    {
        float closest_distance = 1e10f;
        int closest_center = -1;

        for (std::size_t i = 0; i < centers.size(); ++i)
        {
            float distance = dist_fnc(word, centers[i]);
            if (distance < closest_distance)
            {
                closest_center = i;
                closest_distance = distance;
            }
        }

        return closest_center;
    }


    Word calc_new_center(std::vector<Word> const& words, std::vector<int> const& cluster_assignment, random_t & my_random) const
    {
        Word new_center;

        if (cluster_assignment.empty())
        {
            int const random_word_index = my_random() * words.size();
            new_center = words[random_word_index];
        }
        else
        {
            new_center = words[cluster_assignment[0]];

            for (std::size_t j = 1; j < cluster_assignment.size(); ++j)
            {
                Word const& word = words[cluster_assignment[j]];

                for (std::size_t k = 0; k < new_center.size(); ++k)
                {
                    new_center[k] += word[k];
                }
            }

            for (std::size_t k = 0; k < new_center.size(); ++k)
            {
                new_center[k] /= cluster_assignment.size();
            }
        }

        return new_center;
    }




    std::vector<Word> brute_kmeans(std::vector<Word> const& words_all,
                                   int const dict_num_centers) const
    {
        Distance_function dist_fnc;
        random_t my_random;

        std::vector<Word> words;

        Word zero_word(27, 0.0f);

        for (std::size_t i = 0; i < words_all.size(); ++i)
        {
            if (dist_fnc(words_all[i], zero_word) > 0.0001f)
            {
                words.push_back(words_all[i]);
            }
        }

        int const iterations = 10;

        std::vector<Word> centers;

        std::vector<int> available_indices;

        for (std::size_t i = 0; i < words.size(); ++i)
        {
            available_indices.push_back(i);
        }

        for (int i = 0; i < dict_num_centers; ++i)
        {
            int const index = available_indices[my_random() * available_indices.size()];
            std::swap(available_indices[index], available_indices.back());
            available_indices.pop_back();

            centers.push_back(words[index]);
        }


        for (int iteration = 0; iteration < iterations; ++iteration)
        {
            std::cout << "Iteration: " << iteration << std::endl;

            ANN_wrapper ann_wrapper;

            if (_use_ann)
            {
                ann_wrapper.generate_tree_from_centers(centers);
            }

            std::vector< std::vector<int> > assignment = std::vector< std::vector<int> >(dict_num_centers);

            for (std::size_t i = 0; i < words.size(); ++i)
            {
                Word const& word = words[i];

                int closest_center;

                if (_use_ann)
                {
                    closest_center = ann_wrapper.find_closest_center_ann(word);
                }
                else
                {
                    closest_center = find_closest_center_linear(word, centers, dist_fnc);
                }

                assignment[closest_center].push_back(i);
            }

            assert(assignment.size() == centers.size());

            // find new centers as the assigned points' average
            for (std::size_t cluster_index = 0; cluster_index < assignment.size(); ++cluster_index)
            {
                std::vector<int> const& cluster_assignment = assignment[cluster_index];
                centers[cluster_index] = calc_new_center(words, cluster_assignment, my_random);
            }

        }

        centers.push_back(zero_word);

        return centers;
    }


    std::vector<Word> generate(std::vector<Word> const& words,
                               int const dict_num_centers) const
    {
        std::cout << "generate_dictionary_from_gi_points_kmeans()" << std::endl;

        // kmeans<std::vector<Word>, Distance_function> my_kmeans(words, dict_num_centers);
        // my_kmeans.run_default();
        // my_kmeans.run(10, 0.01);

        // std::vector<Word> const& centers = my_kmeans.centers();

        std::vector<Word> centers = brute_kmeans(words, dict_num_centers);

//        if (_do_stats)
//        {
//            calc_stats(my_kmeans, words);
//        }

        return centers;
    }

    std::string identify() const
    {
        return "Kmeans_dictionary";
    }

    std::string get_stat_results() const { return _stat_results; }

private:

//    float sphere_volume(int const dimension, float const radius) const
//    {
//        int const n = dimension;
//        float const R = radius;

//        // http://en.wikipedia.org/wiki/Deriving_the_volume_of_an_n-ball
//        // return std::pow(M_PI, n / 2.0f) * std::pow(R, n) / gamma_fnc((n / 2.0f) + 1.0f);
//    }

    Word find_min(Word const& w_0, Word const& w_1) const
    {
        Word result(w_0.size());

        for (std::size_t i = 0; i < result.size(); ++i)
        {
            result[i] = std::min(w_0[i], w_1[i]);
        }

        return result;
    }

    Word find_max(Word const& w_0, Word const& w_1) const
    {
        Word result(w_0.size());

        for (std::size_t i = 0; i < result.size(); ++i)
        {
            result[i] = std::max(w_0[i], w_1[i]);
        }

        return result;
    }

    double calc_volume(Word const& w_min, Word const& w_max) const
    {
        double result = 1.0f;
        for (std::size_t i = 0; i < w_min.size(); ++i)
        {
            result *= w_max[i] - w_min[i];
        }

        return result;
    }


    float cluster_variance(Word const& center, std::vector<int> const& assigned_node_indices, std::vector<Word> const& words) const
    {
        imdb::l2norm_squared<Word> distance_fnc;

        float variance = 0.0f;

        for (std::size_t i = 0; i < assigned_node_indices.size(); ++i)
        {
            Word const& node = words[assigned_node_indices[i]];

            float const squared_distance = distance_fnc(node, center);
            variance += squared_distance;
        }

        variance /= float(assigned_node_indices.size());

        return variance;
    }

    bool _use_ann;
    bool _do_stats;
    mutable std::string _stat_results;
};



__END_YAFRAY

#endif // DICTIONARY_GENERATOR_H
