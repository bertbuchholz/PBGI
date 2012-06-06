

#include <utilities/CubeRasterBuffer.h>

__BEGIN_YAFRAY

template <class Data>
void Mises_fisher_lobe<Data>::fit(Cube_raster_buffer<Data> const& likelihoods)
{
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
    mean_dir = newMean;
    concentration = newConcentration;
}

__END_YAFRAY
