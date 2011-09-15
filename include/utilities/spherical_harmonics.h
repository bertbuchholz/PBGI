#ifndef SPHERICAL_HARMONICS_H
#define SPHERICAL_HARMONICS_H

#include <cmath>
//#include <stdint.h>

#include <boost/serialization/export.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/serialization.hpp>

#include <yafray_config.h>
#include <utilities/mcqmc.h>
#include <core_api/scene.h>
#include "CubeRasterBuffer.h"

__BEGIN_YAFRAY

static unsigned int factorial_table[] =
{
1,
1,
2,
6,
24,
120,
720,
5040,
40320,
362880,
3628800,
39916800,
479001600 };

/*,
6227020800,
87178291200,
1307674368000,
20922789888000,
355687428096000,
6402373705728000,
121645100408832000,
2432902008176640000,
51090942171709440000,
1124000727777607680000,
25852016738884976640000,
620448401733239439360000,
15511210043330985984000000,
403291461126605635584000000,
10888869450418352160768000000,
304888344611713860501504000000,
8841761993739701954543616000000,
265252859812191058636308480000000,
8222838654177922817725562880000000,
263130836933693530167218012160000000,
8683317618811886495518194401280000000 };
*/

template <class Color, class Vector>
Color estimate_light_sample_no_shadow_test(renderState_t &state, ray_t const& lightRay, Color const& light_color, const surfacePoint_t &surfel, const Vector &wo)
{
    color_t col(0.f);
    const material_t *material = surfel.material;

    // handle lights with delta distribution, e.g. point and directional lights
    Color surfCol = material->eval(state, surfel, wo, lightRay.dir, BSDF_ALL);
    // color_t transmitCol = scene->volIntegrator->transmittance(state, lightRay);
    col = surfCol * light_color * std::fabs(surfel.N*lightRay.dir); // * transmitCol;

    return col;
}

template <class Point, class Color>
class Spherical_function
{
public:
    virtual ~Spherical_function() {}

    virtual void calc_coefficients_random(renderState_t & state, surfacePoint_t const& surfel, float const area, ray_t const& lightRay, color_t const& light_color) = 0;

    virtual float get_area(Point const& dir) const  = 0;
    virtual Color get_color(Point const& dir) const = 0;

    virtual void  add(Spherical_function const* s) = 0;
    virtual void  normalize(float const factor) = 0;

    virtual Spherical_function * clone() = 0;

    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive & /*ar*/, const unsigned int /*version*/)
    {
    }
};

template <class Point, class Color>
class GiSphericalHarmonics : public Spherical_function<Point, Color>
{
    public:
    typedef float (GiSphericalHarmonics::*sh_coefficient_func)(Point const& dir) const; // sh_coefficient_func;

    GiSphericalHarmonics() :
        bands(3),
        exact(false)
    {
        init();
    }

    GiSphericalHarmonics(bool exact, int b)
    {
        this->exact = exact;

        if (exact)
        {
            bands = b;
        }
        else
        {
            bands = 3;
        }

        init();
    }

    void init()
    {
        sh_color_coefficients = std::vector<Color>(bands * bands, Color(0.0f));
        sh_area_coefficients = std::vector<float>(bands * bands, 0.0f);

        // index function: sh_index = l * (l + 1) + m
        // l   m   index
        // 0   0   0
        // 1  -1   1
        // 1   0   2
        // 1   1   3
        // 2  -2   4
        // 2  -1   5
        // 2   0   6
        // 2   1   7
        // 2   2   8

        sh_functions[0] = &GiSphericalHarmonics::SH_precomputed_0;
        sh_functions[1] = &GiSphericalHarmonics::SH_precomputed_1;
        sh_functions[2] = &GiSphericalHarmonics::SH_precomputed_2;
        sh_functions[3] = &GiSphericalHarmonics::SH_precomputed_3;
        sh_functions[4] = &GiSphericalHarmonics::SH_precomputed_4;
        sh_functions[5] = &GiSphericalHarmonics::SH_precomputed_5;
        sh_functions[6] = &GiSphericalHarmonics::SH_precomputed_6;
        sh_functions[7] = &GiSphericalHarmonics::SH_precomputed_7;
        sh_functions[8] = &GiSphericalHarmonics::SH_precomputed_8;
    }

    virtual void calc_coefficients_random(renderState_t & state, surfacePoint_t const& surfel, float const area, ray_t const& lightRay, color_t const& light_color)
    {
        sh_color_coefficients = std::vector<Color>(bands * bands, Color(0.0f));
        sh_area_coefficients = std::vector<float>(bands * bands, 0.0f);

        random_t random;

        vector3d_t const& normal = surfel.N;

        float res_u = 10.0f;
        float res_v = 20.0f;

        for (int u = 0; u < res_u; ++u)
        {
            for (int v = 0; v < res_v; ++v)
            {
                float x = (u + random()) / res_u;
                float y = (v + random()) / res_v;
                //                float x = u / res_u;
                //                float y = v / res_v;

                float theta = 2.0f * std::acos(std::sqrt(1.0f - x));
                float phi = 2.0f * M_PI * y;

                Point dir(std::sin(theta) * std::cos(phi),
                          std::sin(theta) * std::sin(phi),
                          std::cos(theta));

                // std::cout << "---" << std::endl;

                float const cos_dir_normal_abs = std::abs(dir * normal);
                vector3d_t const& wo = dir;
                color_t reflected_color = estimate_light_sample_no_shadow_test(state, lightRay, light_color, surfel, wo);

                if (exact)
                {
                    for (int l = 0; l < bands; ++l)
                    {
                        for(int m = -l; m <= l; ++m)
                        {
                            int sh_index = l * (l + 1) + m;

                            float Y_l_m = SH(l, m, theta, phi);

                            // std::cout << Y_l_m << " " << Y_index << std::endl;

                            // float base_coeff = std::max(dir * normal, 0.0f) * Y_l_m;
                            float base_coeff = Y_l_m;

                            // sh_color_coefficients[sh_index] += color * energy * base_coeff;
                            sh_color_coefficients[sh_index] += base_coeff * reflected_color;
                            sh_area_coefficients[sh_index]  += base_coeff * area * cos_dir_normal_abs;
                            // sh_area_coefficients[sh_index]  += base_coeff * area;
                        }
                    }
                }
                else
                {
                    for (int sh_index = 0; sh_index < 9; ++sh_index)
                    {
                        float const Y_l_m = SH_precomputed(sh_index, dir);

                        // float const base_coeff = std::max(dir * normal, 0.0f) * Y_l_m;
                        float const base_coeff = Y_l_m;

                        // sh_color_coefficients[sh_index] += color * energy * base_coeff;
                        sh_color_coefficients[sh_index] += base_coeff * reflected_color;
                        sh_area_coefficients[sh_index]  += base_coeff * area * cos_dir_normal_abs;
                        // sh_area_coefficients[sh_index]  += base_coeff * area;
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


                if (exact)
                {
                    test_area_sh += get_sh_area(theta, phi);
                }
                else
                {
                    test_area_sh += get_sh_area(dir);
                }

                test_area_real += area * std::max(dir * normal, 0.0f);
            }
        }

        float inv = 4.0f * M_PI / (res_u * res_v);

        test_area_real *= inv;
        test_area_sh *= inv;

        std::cout << "random, test area real: " << test_area_real << " sh: " << test_area_sh << std::endl;
    }

    virtual float get_area(Point const& dir) const
    {
        return get_sh_area(dir);
    }

    virtual Color get_color(Point const& dir) const
    {
        return get_sh_color(dir);
    }

    virtual void add(Spherical_function<Point, Color> const* s)
    {
        GiSphericalHarmonics const* sph_harmonics = dynamic_cast<GiSphericalHarmonics const*>(s);

        *this = *this + *sph_harmonics;
    }

    void normalize(float const factor)
    {
        normalize_color(factor);
    }

    void normalize_color(float const factor)
    {
        for (int j = 0; j < bands * bands; ++j)
        {
            sh_color_coefficients[j] *= factor;
        }
    }

    Spherical_function<Point, Color> * clone()
    {
        return new GiSphericalHarmonics(*this);
    }

    friend GiSphericalHarmonics operator+(GiSphericalHarmonics const& lhs, GiSphericalHarmonics const& rhs)
    {
        assert(lhs.bands == rhs.bands && lhs.exact == rhs.exact);

        GiSphericalHarmonics result(lhs.exact, lhs.bands);

        for (int j = 0; j < lhs.bands * lhs.bands; ++j)
        {
            result.sh_area_coefficients[j] = lhs.sh_area_coefficients[j] + rhs.sh_area_coefficients[j];
            result.sh_color_coefficients[j] = lhs.sh_color_coefficients[j] + rhs.sh_color_coefficients[j];
        }

        return result;
    }

    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & boost::serialization::base_object<Spherical_function<Point, Color> >(*this);

        ar & bands;

        ar & sh_color_coefficients;
        ar & sh_area_coefficients;

        ar & exact;
    }

private:
    float get_sh_area(Point const& dir) const
    {
        if (!exact)
        {
            float result = 0.0f;

            for (int sh_index = 0; sh_index < 9; ++sh_index)
            {
                float Y_l_m = SH_precomputed(sh_index, dir);

                result += sh_area_coefficients[sh_index] * Y_l_m;

            }

            return result;
        }

        float theta = std::acos(dir.z);
        float phi = std::atan2(dir.y, dir.x);
        return get_sh_area(theta, phi);
    }

    float get_sh_area(float const theta, float const phi) const
    {
        if (!exact)
            std::cout << "shouldn't be used!" << std::endl;

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
        if (!exact)
        {
            Color result(0.0f);

            for (int sh_index = 0; sh_index < 9; ++sh_index)
            {
                float Y_l_m = SH_precomputed(sh_index, dir);

                result += sh_color_coefficients[sh_index] * Y_l_m;

            }

            return result;
        }

        float theta = std::acos(dir.z);
        float phi = std::atan2(dir.y, dir.x);

        return get_sh_color(theta, phi);
    }

    Color get_sh_color(float const theta, float const phi) const
    {
        if (!exact)
            std::cout << "shouldn't be used!" << std::endl;

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

    int factorial(int num) const
    {
        if (num < 13) return factorial_table[num];

        std::cout << "factorial(): num > 12 ..." << std::endl;

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
        if (m == 0) return K(l, 0) * P(l, m, std::cos(theta));
        else if (m > 0) return M_SQRT2 * K(l, m) * std::cos(m * phi) * P(l, m, std::cos(theta));
        else return M_SQRT2 * K(l, -m) * std::sin(-m * phi) * P(l, -m, std::cos(theta));
    }

    float SH_precomputed(int const index, Point const& dir) const
    {
        sh_coefficient_func f = sh_functions[index];

        return (this->*(f))(dir);
    }

    inline float square(float x) const { return x * x; }

    // http://en.wikipedia.org/wiki/Table_of_spherical_harmonics

    float SH_precomputed_0(Point const& /* dir */) const
    {
        return 0.28209;
    }

    float SH_precomputed_1(Point const& dir) const
    {
        return 0.34549 * dir[0];
    }

    float SH_precomputed_2(Point const& dir) const
    {
        return 0.34549 * dir[1];
    }

    float SH_precomputed_3(Point const& dir) const
    {
        return 0.34549 * dir[2];
    }

    float SH_precomputed_7(Point const& dir) const
    {
        return 0.31539 * (2 * square(dir[2]) - square(dir[0]) - square(dir[1]));
    }

    float SH_precomputed_5(Point const& dir) const
    {
        return 1.0925 * dir[1] * dir[2];
    }

    float SH_precomputed_4(Point const& dir) const
    {
        return 1.0925 * dir[2] * dir[0];
    }

    float SH_precomputed_6(Point const& dir) const
    {
        return 1.0925 * dir[0] * dir[1];
    }

    float SH_precomputed_8(Point const& dir) const
    {
        return 0.54627 * (square(dir[0]) - square(dir[1]));
    }

    int bands;

    std::vector<Color> sh_color_coefficients;
    std::vector<float> sh_area_coefficients;

    sh_coefficient_func sh_functions[9];

    bool exact;
};


template <class Point, class Color>
class Cube_spherical_function : public Spherical_function<Point, Color>
{
public:
    Cube_spherical_function()
    {
        // init the buffers!
        color_buffer.setup_simple(8);
        area_buffer.setup_simple(8);
    }

    void calc_coefficients_random(renderState_t & state, surfacePoint_t const& surfel, float const area, ray_t const& lightRay, color_t const& light_color)
    {
        Point const& normal = surfel.N;

        std::vector<Cube_cell> const& cube_cells = color_buffer.get_cube_cells();

        for (unsigned int i = 0; i < cube_cells.size(); ++i)
        {
            Point const& dir = color_buffer.get_cell_center(cube_cells[i]);

            float const cos_dir_normal_abs = std::abs(dir * normal);
            vector3d_t const& wo = dir;
            color_t reflected_color = estimate_light_sample_no_shadow_test(state, lightRay, light_color, surfel, wo);

            color_buffer.set_color(cube_cells[i], reflected_color);
            area_buffer .set_color(cube_cells[i], color_t(area * cos_dir_normal_abs));
            // area_buffer .set_color(cube_cells[i], color_t(cos_dir_normal_abs));
        }
    }

    float get_area(Point const& dir) const
    {
        Cube_cell c;
        bool ok = area_buffer.get_cell(dir, c);
        assert(ok);

        Color area = area_buffer.get_color(c);

        return area[0];
    }

    Color get_color(Point const& dir) const
    {
        Cube_cell c;
        bool ok = color_buffer.get_cell(dir, c);
        assert(ok);

        return color_buffer.get_color(c);
    }

    void add(Spherical_function<Point, Color> const* s)
    {
        Cube_spherical_function const* csf = dynamic_cast<Cube_spherical_function<Point, Color> const*>(s);

        assert(color_buffer.get_resolution() == csf->color_buffer.get_resolution());

        std::vector<Cube_cell> const& cube_cells = color_buffer.get_cube_cells();

        for (unsigned int i = 0; i < cube_cells.size(); ++i)
        {
            Cube_cell const& c = cube_cells[i];

            Color new_color = color_buffer.get_color(c) + csf->color_buffer.get_color(c);
            float new_area  = area_buffer.get_color(c)[0] + csf->area_buffer.get_color(c)[0];

            color_buffer.set_color(c, new_color);
            area_buffer .set_color(c, new_area);
        }
    }

    void normalize(float const factor)
    {
        std::vector<Cube_cell> const& cube_cells = color_buffer.get_cube_cells();

        for (unsigned int i = 0; i < cube_cells.size(); ++i)
        {
            Cube_cell const& c = cube_cells[i];

            Color new_color = color_buffer.get_color(c) * factor;

            color_buffer.set_color(c, new_color);
        }
    }

    Spherical_function<Point, Color> * clone()
    {
        assert(false);
        return NULL;
    }

    Cube_raster_buffer const& get_color_buffer() const
    {
        return color_buffer;
    }

    Cube_raster_buffer const& get_area_buffer() const
    {
        return area_buffer;
    }

private:
    Cube_raster_buffer color_buffer;
    Cube_raster_buffer area_buffer;
};

__END_YAFRAY

#endif // SPHERICAL_HARMONICS_H
