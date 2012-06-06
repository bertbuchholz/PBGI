#ifndef MISES_FISHER_H
#define MISES_FISHER_H

#include <core_api/vector3d.h>
#include <core_api/color.h>
#include <utilities/mcqmc.h>


__BEGIN_YAFRAY

// template <class Data>
// class Cube_raster_buffer;

template <class Output_data> inline
Output_data convert_data(color_t const& color);

template <> inline
color_t convert_data<color_t>(color_t const& color)
{
    return color;
}

template <> inline
float convert_data<float>(color_t const& color)
{
    return color.energy();
}

template <class Output_data> inline
Output_data convert_data(float const& f);

template <> inline
color_t convert_data<color_t>(float const& f)
{
    return color_t(f);
}

template <> inline
float convert_data<float>(float const& f)
{
    return f;
}



inline
float convert_to_float(color_t const& color)
{
    return color.energy();
}

inline
float convert_to_float(float const f)
{
    return f;
}


inline
float normalize(float const& c)
{
    return c;
}

inline
color_t normalize(color_t const& c)
{
    color_t new_c = c / c.energy();
    return new_c;
}

template <class Data = float>
struct Mises_fisher_lobe
{
    vector3d_t mean_dir;
    float concentration;
    Data weight;

#ifdef PBGI_DEBUG
    vector3d_t pos;
#endif


    Mises_fisher_lobe() {}

    Mises_fisher_lobe(
                vector3d_t mean_dir_,
                float concentration_,
                Data weight_) :
        mean_dir(mean_dir_),
        concentration(concentration_),
        weight(weight_)
    { }

    Data evaluate(vector3d_t const& dir) const
    {
        return weight * (concentration / (2.0f * M_PI) * fExp(concentration * (mean_dir * dir - 1.0f)));
    }

    void randomize(random_t & my_random)
    {
        weight = Data(1.0f);
        mean_dir = vector3d_t(my_random(), my_random(), my_random()) - vector3d_t(0.5f);
        mean_dir.normalize();
        concentration = 1;//0.1f + 2 * MathTools::frand();
    }


    friend std::ostream & operator <<(std::ostream & s, Mises_fisher_lobe const& l)
    {
        s << "mean: " << l.mean_dir << " k: " << l.concentration << " weight: " << l.weight;
#ifdef PBGI_DEBUG
        s << " pos: " << l.pos;
#endif
        return s;
    }
};

/*
__END_YAFRAY
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
*/

__END_YAFRAY

#endif // MISES_FISHER_H
