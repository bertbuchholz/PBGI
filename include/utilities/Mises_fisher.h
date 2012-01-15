#ifndef MISES_FISHER_H
#define MISES_FISHER_H

#include <core_api/vector3d.h>
#include <utilities/CubeRasterBuffer.h>

__BEGIN_YAFRAY

struct Mises_fisher_lobe
{
    vector3d_t mean_dir;
    float concentration;
    float weight;

    Mises_fisher_lobe() {}

    Mises_fisher_lobe(
                vector3d_t mean_dir_,
                float concentration_,
                float weight_) :
        mean_dir(mean_dir_),
        concentration(concentration_),
        weight(weight_)
    { }

    float evaluate(vector3d_t const& dir) const
    {
        return weight * concentration / (2.0f * M_PI) * std::exp(concentration * (mean_dir * dir - 1.0f));
    }

    void randomize(random_t & my_random)
    {
        weight = 1;
        mean_dir = vector3d_t(my_random(), my_random(), my_random()) - vector3d_t(0.5f);
        mean_dir.normalize();
        concentration = 10;//0.1f + 2 * MathTools::frand();
    }

    void fit(Cube_raster_buffer const& likelihoods)
    {
        std::vector<Cube_cell> const& cells = likelihoods.get_cube_cells();

        const int size = cells.size();

        float likelihoodSum = 0.0f;
        for(int k = 0; k < size; k++)
        {
            Cube_cell const& c = cells[k];
            likelihoodSum += likelihoods.get_color(c)[0];
        }

        // Clamp to minimum
        const float minLikelihoodSum = 0.001f;
        likelihoodSum = std::max(0.1f, likelihoodSum);
        if (likelihoodSum == minLikelihoodSum)
        {
            std::cout << "Likelihood sum loss" << std::endl;
        }

        // Weight
        const float newWeight = likelihoodSum / float(size);

        // Mean
        vector3d_t weightedDirectionSum(0.0f);

        for(int k = 0; k < size; k++)
        {
            Cube_cell const& c = cells[k];
            vector3d_t const direction = likelihoods.get_cell_direction(c);

            const float likelihood = likelihoods.get_color(c)[0];
            weightedDirectionSum += likelihood * direction;
        }
        weightedDirectionSum /= likelihoodSum;

        vector3d_t newMean = weightedDirectionSum;
        newMean.normalize();

        // Concentration
        const float weightedDirectionSumLength = weightedDirectionSum.length();

        const float newConcentration =
                (3.0f * weightedDirectionSumLength - weightedDirectionSumLength * weightedDirectionSumLength * weightedDirectionSumLength) /
                (1.0f - weightedDirectionSumLength * weightedDirectionSumLength);

        // std::cout << concentration << std::endl;

        // Copy
        //weight = newWeight;
        mean_dir = newMean;
        concentration  = newConcentration;
    }

    friend std::ostream & operator <<(std::ostream & s, Mises_fisher_lobe const& l)
    {
        s << "mean: " << l.mean_dir << " k: " << l.concentration << " weight: " << l.weight;
        return s;
    }
};

__END_YAFRAY

#endif // MISES_FISHER_H
