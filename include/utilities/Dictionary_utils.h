#ifndef DICTIONARY_UTILS_H
#define DICTIONARY_UTILS_H

#include <vector>
#include <core_api/color.h>

__BEGIN_YAFRAY


template <class T>
std::vector<float> conv_to_vector(T const& t);

template <>
inline
std::vector<float> conv_to_vector<color_t>(color_t const& color)
{
    std::vector<float> result;
    result.push_back(color.R);
    result.push_back(color.G);
    result.push_back(color.B);
    return result;
}

template <>
inline
std::vector<float> conv_to_vector<float>(float const& f)
{
    std::vector<float> result;
    result.push_back(f);
    return result;
}


//template <class T>
//inline
//T conv_from_vector(std::vector<float> const& data, std::size_t & i);


template <class T>
inline
std::vector<T> conv_from_vector(std::vector<float> const& data);

template <>
inline
std::vector<color_t> conv_from_vector<color_t>(std::vector<float> const& data)
{
    std::vector<color_t> result;

    for (std::size_t i = 0; i < data.size(); i += 3)
    {
        color_t color;
        color.R = data[i];
        color.G = data[i + 1];
        color.B = data[i + 2];

        result.push_back(color);
    }

    return result;
}

template <>
inline
std::vector<float> conv_from_vector<float>(std::vector<float> const& data)
{
    std::vector<float> result;

    for (std::size_t i = 0; i < data.size(); ++i)
    {
        float f = data[i];
        result.push_back(f);
    }

    return result;
}



__END_YAFRAY

#endif // DICTIONARY_UTILS_H
