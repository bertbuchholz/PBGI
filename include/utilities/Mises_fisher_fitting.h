#ifndef MISES_FISHER_FITTING_H
#define MISES_FISHER_FITTING_H

#include <utilities/CubeRasterBuffer.h>
#include <utilities/Mises_fisher.h>

__BEGIN_YAFRAY

template <class Data>
Mises_fisher_lobe<Data> fit(Cube_raster_buffer<Data> const& likelihoods)
{
    Mises_fisher_lobe<Data> lobe;

    std::vector<Cube_cell> const& cells = likelihoods.get_cube_cells();

    const int size = cells.size();

    Data likelihoodSum(0.0f);
    for(int k = 0; k < size; k++)
    {
        Cube_cell const& c = cells[k];
        likelihoodSum += likelihoods.get_data(c);
    }

    // Clamp to minimum
    Data const minLikelihoodSum(0.001f);
    if (likelihoodSum < minLikelihoodSum)
    {
        std::cout << "Likelihood sum loss" << std::endl;
        likelihoodSum = minLikelihoodSum;
        // likelihoods.print();
    }

    // Weight
    // const float newWeight = likelihoodSum / float(size);

    // Mean
    vector3d_t weightedDirectionSum(0.0f);

    for (int k = 0; k < size; k++)
    {
        Cube_cell const& c = cells[k];
        vector3d_t const direction = likelihoods.get_cell_direction(c);

        const Data likelihood = likelihoods.get_data(c);
        weightedDirectionSum += convert_to_float(likelihood) * direction;
    }
    weightedDirectionSum /= convert_to_float(likelihoodSum);


    vector3d_t newMean = weightedDirectionSum;
    newMean.normalize();

    // Concentration
    const float weightedDirectionSumLength = std::min(0.9999f, weightedDirectionSum.length());
    // assert(weightedDirectionSumLength < 0.9999f);

    const float newConcentration =
            (3.0f * weightedDirectionSumLength - weightedDirectionSumLength * weightedDirectionSumLength * weightedDirectionSumLength) /
            (1.0f - weightedDirectionSumLength * weightedDirectionSumLength);

    // std::cout << concentration << std::endl;

    // Copy
    //weight = newWeight;
    lobe.weight = Data(1.0f);
    lobe.mean_dir = newMean;
    lobe.concentration = newConcentration;

    return lobe;
}




template <class Data>
// std::vector<Mises_fisher_lobe<Data> > Mises_Fisher_spherical_function<Data>::generate_from_cube(Cube_raster_buffer<Data> const& distribution, int const num_lobes)
std::vector<Mises_fisher_lobe<Data> > generate_mises_fisher_from_cube(Cube_raster_buffer<Data> const& distribution,
                                                         int const num_lobes,
                                                         Spherical_function<float> const* area_sf)
{
    std::vector<Cube_cell> const& cells = distribution.get_cube_cells();

    const int size = cells.size();

    const int num_iterations = 5;

    std::vector<Mises_fisher_lobe<Data> > lobes(num_lobes);

    std::vector<Cube_raster_buffer<Data> > likelihoods(num_lobes);
    for (int i = 0; i < num_lobes; i++)
    {
        likelihoods[i].setup_surfel_buffer(distribution.get_resolution());
    }

    Cube_raster_buffer<Data> totalLikelihoods;
    totalLikelihoods.setup_surfel_buffer(distribution.get_resolution());

    // Random start guess
    random_t my_random;
    for (int i = 0; i < num_lobes; i++)
    {
        lobes[i].randomize(my_random);
        // std::cout << "Mises_Fisher_spherical_function::generate_from_cube(): " << i << " " << lobes[i] << std::endl;
    }

    // http://crsouza.blogspot.com/2010/10/gaussian-mixture-models-and-expectation.html
    for (int step = 0; step < num_iterations; ++step)
    {
        // E-Step

        // For every observation ..
        for (int j = 0; j < size; j++)
        {
            Cube_cell const& c = cells[j];

            const vector3d_t direction = distribution.get_cell_direction(c);

            // 1.) Evaluate each lobe
            for (int k = 0; k < num_lobes; k++)
            {
                likelihoods[k].set_data(c, lobes[k].evaluate(direction));
            }
        }

        for (int j = 0; j < size; j++)
        {
            // 2.) Sum over all lobes
            Cube_cell const& c = cells[j];

            totalLikelihoods.set_data(c, 0.0f);
            for(int k = 0; k < num_lobes; k++)
            {
                Data current_likelihood = totalLikelihoods.get_data(c);
                current_likelihood += likelihoods[k].get_data(c);
                totalLikelihoods.set_data(c, current_likelihood);
            }
        }

        for (int j = 0; j < size; j++)
        {
            Cube_cell const& c = cells[j];

            // 3.) Divide by total likelihood
            if (totalLikelihoods.get_data(c) > 0.0001f)
            {
                const Data pixel_value = distribution.get_data(c); // + 0.01f;
                const float solid_angle = distribution.get_solid_angle(c); // / (4.0f * M_PI);

                float area = 1.0f;

                if (area_sf)
                {
                    vector3d_t const cell_dir = distribution.get_cell_direction(c);
                    area = area_sf->get_value(cell_dir);
                }

                for (int k = 0; k < num_lobes; k++)
                {
                    likelihoods[k].set_data(c, pixel_value / area * solid_angle * likelihoods[k].get_data(c) / totalLikelihoods.get_data(c));
                }
            }
        }

        // M-step
        for (int lobe_index = 0; lobe_index < num_lobes; lobe_index++)
        {
            lobes[lobe_index] = fit(likelihoods[lobe_index]);

            // std::cout << "step: " << step << " " << lobes[lobe_index] << std::endl;
        }
    }

    // calculate lobe weight
    for (int cell_index = 0; cell_index < size; cell_index++)
    {
        // 2.) Sum over all lobes
        Cube_cell const& c = cells[cell_index];

        Data likelihood_sum(0.0f);

        for (int k = 0; k < num_lobes; k++)
        {
            likelihood_sum += likelihoods[k].get_data(c);
        }

        totalLikelihoods.set_data(c, likelihood_sum);

        for (int k = 0; k < num_lobes; k++)
        {
            if (totalLikelihoods.get_data(c) > Data(0.0001f))
            {
                Data likelihood_normalized = likelihoods[k].get_data(c) / totalLikelihoods.get_data(c);
                likelihoods[k].set_data(c, likelihood_normalized);
            }
            else
            {
                likelihoods[k].set_data(c, Data(0.0f));
            }
        }


//        for (int k = 0; k < num_lobes; k++)
//        {
//            Data current_likelihood = totalLikelihoods.get_data(c);
//            current_likelihood += likelihoods[k].get_data(c);
//            totalLikelihoods.set_data(c, current_likelihood);
//        }
    }

    /*
    for (int lobe_index = 0; lobe_index < num_lobes; lobe_index++)
    {
        vector3d_t const& mean = lobes[lobe_index].mean_dir;

        if (mean.length() < 0.001f)
        {
            lobes[lobe_index].weight = Data(0.0f);
            continue;
        }

        Cube_cell c = likelihoods[lobe_index].get_cell(mean);
        Data const likelihood = likelihoods[lobe_index].get_data_interpolated(c);
        Data const total_likelihood = totalLikelihoods.get_data_interpolated(c);
        Data const amplitude = distribution.get_data_interpolated(c);
        Data const lobe_max = lobes[lobe_index].evaluate(mean);

        if (total_likelihood > Data(0.00001f))
        {
            // lobes[lobe_index].weight = (likelihood / total_likelihood) * amplitude;
            // lobes[lobe_index].weight = (convert_to_float(likelihood) / convert_to_float(total_likelihood)) * amplitude;
            lobes[lobe_index].weight = (convert_to_float(likelihood) / convert_to_float(total_likelihood)) * amplitude / convert_to_float(lobe_max);
        }

        std::cout << "weight: " << lobes[lobe_index].weight.energy() << std::endl;
    }
    */


    if (area_sf)
    {
        for (int cell_index = 0; cell_index < size; cell_index++)
        {
            Cube_cell const& c = cells[cell_index];

            vector3d_t const cell_dir = distribution.get_cell_direction(c);
            float const area = area_sf->get_value(cell_dir);

            for (int k = 0; k < num_lobes; k++)
            {
                if (area > 0.0001f)
                {
                    Data likelihood_normalized = likelihoods[k].get_data(c) / area;
                    likelihoods[k].set_data(c, likelihood_normalized);
                }
                else
                {
                    likelihoods[k].set_data(c, Data(0.0f));
                }
            }
        }
    }


    for (int lobe_index = 0; lobe_index < num_lobes; lobe_index++)
    {
        lobes[lobe_index].weight = distribution.integrate(likelihoods[lobe_index]);

        /*
        assert(!std::isnan(lobes[lobe_index].concentration) &&
               !std::isnan(lobes[lobe_index].weight.R) &&
               !std::isnan(lobes[lobe_index].weight.G) &&
               !std::isnan(lobes[lobe_index].weight.B));
               */

//        std::cout << "weight: " << lobes[lobe_index].weight.energy() << std::endl;
    }

    for (std::size_t i = 0; i < lobes.size() && false; ++i)
    {
         // distribution.print();
         std::cout << "Mises_Fisher_spherical_function::generate_from_cube(): lobe " << i << ": " << lobes[i] <<
                      " lobe max: " << lobes[i].evaluate(lobes[i].mean_dir) <<
                      " cube tot: " << distribution.calc_total_energy() <<
                      " cube max: " << distribution.find_maximum() <<
                      std::endl;
    }

    return lobes;
}


//template <class Data>
//Mises_Fisher_spherical_function<Data> * Mises_Fisher_spherical_function<Data>::test()
//{
//    vector3d_t mean_dir(0.0f, 0.5f, 0.5f);
//    mean_dir.normalize();
//    float const concentration = 30;
//    float const weight = 1.0f;

//    Mises_fisher_lobe<float> lobe_0(mean_dir, concentration, weight);
//    Mises_fisher_lobe<float> lobe_1(vector3d_t(1.0f, 0.0f, 0.0f), 10.0f, weight);

//    std::cout << "test lobe 0: " << lobe_0 << std::endl;
//    std::cout << "test lobe 1: " << lobe_1 << std::endl;

//    Cube_spherical_function<float> * distribution   = Cube_spherical_function<float>::create_from_mises_fisher(lobe_0, 8, false, 2.0f);
//    Cube_spherical_function<float> * distribution_1 = Cube_spherical_function<float>::create_from_mises_fisher(lobe_1, 8, false, 1.0f);

//    // distribution->add(distribution_1);


//    Mises_Fisher_spherical_function * mf_sf = new Mises_Fisher_spherical_function;
//    mf_sf->generate_from_cube_spherical_function(distribution, 10);

//    std::cout << "lobes:" << std::endl;
//    for (unsigned int i = 0; i < mf_sf->get_lobes().size(); ++i)
//    {
//        std::cout << i << " " << mf_sf->get_lobes()[i] << std::endl;
//    }

//    return mf_sf;
//}





__END_YAFRAY


#endif // MISES_FISHER_FITTING_H
