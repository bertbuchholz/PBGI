#ifndef SPHERICAL_HARMONICS_H
#define SPHERICAL_HARMONICS_H

#include <yafray_config.h>
#include <utilities/mcqmc.h>
#include <cmath>

__BEGIN_YAFRAY

template <class Point, class Color>
class GiSphericalHarmonics
{
    public:
    GiSphericalHarmonics(int b = 3)
    {
        bands = b;

        sh_color_coefficients = std::vector<Color>(bands * bands, Color(0.0f));
        sh_area_coefficients = std::vector<float>(bands * bands, 0.0f);
    }

    void calc_coefficients_random(Point const& normal, Color const& color, Color const& energy, float const area)
    {
        sh_color_coefficients = std::vector<Color>(bands * bands, Color(0.0f));
        sh_area_coefficients = std::vector<float>(bands * bands, 0.0f);

        random_t random;

        float res_u = 20.0f;
        float res_v = 20.0f;

        for (int u = 0; u < res_u; ++u)
        {
            for (int v = 0; v < res_v; ++v)
            {
                float x = (u + random()) / res_u;
                float y = (v + random()) / res_v;

                float theta = 2.0f * std::acos(std::sqrt(1.0f - x));
                float phi = 2.0f * M_PI * y;

                Point dir(std::sin(theta) * std::cos(phi),
                          std::sin(theta) * std::sin(phi),
                          std::cos(theta));

                for (int l = 0; l < bands; ++l)
                {
                    for(int m = -l; m <= l; ++m)
                    {
                        int sh_index = l * (l + 1) + m;

                        float Y_l_m = SH(l, m, theta, phi);

                        float area_coeff = area * std::max(dir * normal, 0.0f) * Y_l_m;

                        sh_color_coefficients[sh_index] += color * energy * area_coeff;
                        sh_area_coefficients[sh_index]  += area_coeff;
                    }
                }
            }
        }

        float inv = 4.0f * M_PI / (res_u * res_v);

        for (int sh_index = 0; sh_index < bands * bands; ++sh_index)
        {
            sh_color_coefficients[sh_index] *= inv;
            sh_area_coefficients[sh_index]  *= inv;
        }

    }

    void test(Point const& normal, float const area)
    {
        random_t random;

        float res_u = 20.0f;
        float res_v = 20.0f;

        float test_area_real = 0.0f;
        float test_area_sh = 0.0f;

        for (int u = 0; u < res_u; ++u)
        {
            for (int v = 0; v < res_v; ++v)
            {
                float x = (u + random()) / res_u;
                float y = (v + random()) / res_v;

                float theta = 2.0f * std::acos(std::sqrt(1.0f - x));
                float phi = 2.0f * M_PI * y;

                Point dir(std::sin(theta) * std::cos(phi),
                          std::sin(theta) * std::sin(phi),
                          std::cos(theta));

                test_area_real += area * std::max(dir * normal, 0.0f);

                test_area_sh += get_sh_area(theta, phi);
            }
        }

        float inv = 4.0f * M_PI / (res_u * res_v);

        test_area_real *= inv;
        test_area_sh *= inv;

        std::cout << "random, test area real: " << test_area_real << " sh: " << test_area_sh << std::endl;
    }


    float get_sh_area(Point const& dir) const
    {
        float theta = std::acos(dir.z);
        float phi = std::atan2(dir.y, dir.x);

        return get_sh_area(theta, phi);
    }

    float get_sh_area(float const theta, float const phi) const
    {
        float result = 0.0f;

        for (int l = 0; l < bands; ++l)
        {
            for(int m = -l; m <= l; ++m)
            {
                int sh_index = l * (l + 1) + m;

                float Y_l_m = SH(l, m, theta, phi);

                result += sh_area_coefficients[sh_index] * Y_l_m;
            }
        }

        return result;
    }

    Color get_sh_color(Point const& dir) const
    {
        float theta = std::acos(dir.z);
        float phi = std::atan2(dir.y, dir.x);

        return get_sh_color(theta, phi);
    }

    Color get_sh_color(float const theta, float const phi) const
    {
        Color result(0.0f);

        for (int l = 0; l < bands; ++l)
        {
            for(int m = -l; m <= l; ++m)
            {
                int sh_index = l * (l + 1) + m;

                float Y_l_m = SH(l, m, theta, phi);

                result += sh_color_coefficients[sh_index] * Y_l_m;
            }
        }

        return result;
    }

    Color const& get_color_coefficient(int index) const
    {
        return sh_color_coefficients[index];
    }

    float get_area_coefficient(int index) const
    {
        return sh_area_coefficients[index];
    }

    friend std::ostream & operator<<(std::ostream & s, GiSphericalHarmonics const& p)
    {
        for (int i = 0; i < 9; ++i)
        {
            s <<
                 p.sh_color_coefficients[i].R << " " <<
                 p.sh_color_coefficients[i].G << " " <<
                 p.sh_color_coefficients[i].B << " ";
        }

        s << std::endl;

        for (int i = 0; i < 9; ++i)
        {
            s << p.sh_area_coefficients[i] << " ";
        }

        return s;
    }

    int factorial(int num) const
    {
        int result = 1;
        for (int i = 1; i <= num; ++i)
        {
            result *= i;
        }

        return result;
    }

    float P(int const l, int const m, float x) const
    {
        float pmm = 1.0f;
        if (m > 0)
        {
            float somx2 = std::sqrt((1.0f - x) * (1.0f + x));
            float fact = 1.0f;
            for (int i = 1; i <= m; ++i)
            {
                pmm *= (-fact) * somx2;
                fact += 2.0f;
            }
        }

        if (l == m) return pmm;

        float pmmp1 = x * (2.0f * m + 1.0f) * pmm;

        if (l == m + 1) return pmmp1;

        float pll = 0.0f;

        for (int ll = m + 2; ll <= l; ++ll)
        {
            pll = ( (2.0f * ll - 1.0f) * x * pmmp1 - (ll + m - 1.0f) * pmm ) / (ll - m);
            pmm = pmmp1;
            pmmp1 = pll;
        }

        return pll;
    }

    float K(int const l, int const m) const
    {
        float temp = ((2.0f * l + 1.0f) * factorial(l - m)) / (4.0f * M_PI * factorial(l + m));
        return std::sqrt(temp);
    }

    float SH(int const l, int const m, float const theta, float const phi) const
    {
        const float sqrt2 = std::sqrt(2.0f);
        if (m == 0) return K(l, 0) * P(l, m, std::cos(theta));
        else if (m > 0) return sqrt2 * K(l, m) * std::cos(m * phi) * P(l, m, std::cos(theta));
        else return sqrt2 * K(l, -m) * std::sin(-m * phi) * P(l, -m, std::cos(theta));
    }

    friend std::istream & operator>>(std::istream & s, GiSphericalHarmonics & p)
    {
        for (int i = 0; i < 9; ++i)
        {
            s >>
                 p.sh_color_coefficients[i].R >>
                 p.sh_color_coefficients[i].G >>
                 p.sh_color_coefficients[i].B;
        }

        for (int i = 0; i < 9; ++i)
        {
            s >> p.sh_area_coefficients[i];
        }

        return s;
    }

    friend GiSphericalHarmonics operator+(GiSphericalHarmonics const& lhs, GiSphericalHarmonics const& rhs)
    {
        assert(lhs.bands == rhs.bands);

        GiSphericalHarmonics result(lhs.bands);

        for (int j = 0; j < lhs.bands; ++j)
        {
            result.sh_area_coefficients[j] = lhs.sh_area_coefficients[j] + rhs.sh_area_coefficients[j];
            result.sh_color_coefficients[j] = lhs.sh_color_coefficients[j] + rhs.sh_color_coefficients[j];
        }

        return result;
    }

    private:
    int bands;

    std::vector<Color> sh_color_coefficients;
    std::vector<float> sh_area_coefficients;
};



__END_YAFRAY

#endif // SPHERICAL_HARMONICS_H
