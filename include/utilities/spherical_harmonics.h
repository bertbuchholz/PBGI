#ifndef SPHERICAL_HARMONICS_H
#define SPHERICAL_HARMONICS_H

#include <cmath>
//#include <stdint.h>

#include <boost/serialization/export.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/serialization.hpp>

#include <yafray_config.h>
#include <utilities/mcqmc.h>
#include <utilities/sample_utils.h>
#include <core_api/scene.h>
#include <core_api/surface.h>
#include <core_api/material.h>
#include <utilities/CubeRasterBuffer.h>
#include <utilities/Mises_fisher.h>
#include <utilities/Quaternion_Matrix.h>


#include <integrators/kmeans.hpp>

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

template <class Color, class Point>
Color estimate_light_sample_no_shadow_test(renderState_t &state, ray_t const& lightRay, Color const& light_color, const surfacePoint_t &surfel, const Point &wo)
{
    Color col(0.f);
    const material_t *material = surfel.material;

    // handle lights with delta distribution, e.g. point and directional lights
    Color surfCol = material->eval(state, surfel, wo, lightRay.dir, BSDF_ALL);
    // color_t transmitCol = scene->volIntegrator->transmittance(state, lightRay);
    col = surfCol * light_color * std::max(0.0f, surfel.N*lightRay.dir); // * transmitCol;

    return col;
}

typedef vector3d_t Point;
typedef color_t Color;

class GiSphericalHarmonics;
class Cube_spherical_function;
class Mises_Fisher_spherical_function;
class Indexed_Cube_Spherical_function;

class Spherical_function_visitor
{
public:
    virtual void visit(GiSphericalHarmonics * sf) = 0;
    virtual void visit(Cube_spherical_function * sf) = 0;
    virtual void visit(Mises_Fisher_spherical_function * sf) = 0;
    virtual void visit(Indexed_Cube_Spherical_function * sf) = 0;
    virtual void identify() { std::cout << "Spherical_function_visitor" << std::endl; }
};


class Spherical_function
{
public:
    virtual ~Spherical_function() {}

    virtual void calc_coefficients_random(renderState_t & state, surfacePoint_t const& surfel, float const area, ray_t const& lightRay, color_t const& light_color) = 0;

    virtual float get_area(Point const& dir) const  = 0;
    virtual Color get_color(Point const& dir) const = 0;

    virtual void add(Spherical_function const* s) = 0;
    virtual void normalize(float const factor) = 0;

    virtual void accept(Spherical_function_visitor * visitor) = 0;

    virtual std::vector< std::vector<float> > components_to_vectors() const { return std::vector< std::vector<float> >(); }
    virtual void from_component_vectors(std::vector< std::vector<float> > const& /* data */) {}

    virtual Spherical_function * clone() = 0;

    /*
    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive & , const unsigned int )
    {
    }
    */
};

class GiSphericalHarmonics : public Spherical_function
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

    virtual void add(Spherical_function const* s)
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

    Spherical_function * clone()
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

    /*
    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & boost::serialization::base_object<Spherical_function>(*this);

        ar & bands;

        ar & sh_color_coefficients;
        ar & sh_area_coefficients;

        ar & exact;
    }
    */


    virtual void accept(Spherical_function_visitor * visitor)
    {
        visitor->visit(this);
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


class Cube_spherical_function : public Spherical_function
{
public:
    Cube_spherical_function(int const resolution, bool use_jitter)
    {
        // init the buffers!
        color_buffer.setup_surfel_buffer(resolution);
        area_buffer. setup_surfel_buffer(resolution);

        _use_jitter = use_jitter;

        if (_use_jitter)
        {
            //_jitter_x = 4.0f * (rand() / float(RAND_MAX) - 0.5f);
            //_jitter_y = 4.0f * (rand() / float(RAND_MAX) - 0.5f);

            // vector3d_t axis(0.0f, 0.0f, 1.0f);
            vector3d_t axis(2.0f * rand() / float(RAND_MAX) - 1.0f, 2.0f * rand() / float(RAND_MAX) - 1.0f, 2.0f * rand() / float(RAND_MAX) - 1.0f);
            axis.normalize();
            // float angle = 2.0f * M_PI / 8.0f;
            float angle = rand() / float(RAND_MAX) * 2.0f * M_PI / 40.0f;

            Quaternion q_pos(axis, angle);
            _jitter_rotation_positive = q_pos.to_matrix();

            Quaternion q_neg(axis, -angle);
            _jitter_rotation_negative = q_neg.to_matrix();
        }
    }

    // only one sample on each cube cell (center) and the surfel (center)
    void sample_single(renderState_t & state, surfacePoint_t const& surfel, float const area, ray_t const& lightRay, color_t const& light_color)
    {
        Point const& normal = surfel.N;

        std::vector<Cube_cell> const& cube_cells = color_buffer.get_cube_cells();

        for (unsigned int i = 0; i < cube_cells.size(); ++i)
        {
            Cube_cell const& cell = cube_cells[i];

            Point dir = color_buffer.get_cell_center(cell);
            dir.normalize();

            float const cos_dir_normal_abs = std::abs(dir * normal);
            assert(cos_dir_normal_abs <= 1.0f);

            vector3d_t const& wo = dir;
            color_t reflected_color = estimate_light_sample_no_shadow_test(state, lightRay, light_color, surfel, wo);

            color_buffer.set_color(cell, reflected_color);
            area_buffer .set_color(cell, color_t(area * cos_dir_normal_abs));
        }
    }

    // one sample on the surfel (center) and many on each cube cell
    void sample_cube_cell_area(renderState_t & state, surfacePoint_t const& surfel, float const area, ray_t const& lightRay, color_t const& light_color)
    {
        Point const& normal = surfel.N;

        std::vector<Cube_cell> const& cube_cells = color_buffer.get_cube_cells();

        int const num_samples_per_cell = 8;
        int const combined_sample_count = num_samples_per_cell * cube_cells.size();

        float const cell_width_2 = 1.0f / color_buffer.get_resolution();

        int combined_samples_n = 0;

        random_t my_random;

        color_t emission = surfel.material->emission(state, surfel, vector3d_t());

        for (unsigned int i = 0; i < cube_cells.size(); ++i)
        {
            Cube_cell const& cell = cube_cells[i];

            color_t reflected_color(0.0f);
            float averaged_area = 0.0f;

            int const axis_0 = Cube_raster_buffer::get_corresponding_axis(0, cell.plane) / 2;
            int const axis_1 = Cube_raster_buffer::get_corresponding_axis(1, cell.plane) / 2;

            for (int j = 0; j  < num_samples_per_cell; ++j)
            {
                float s1, s2;

                // hammersley_2(combined_samples_n, combined_sample_count, s1, s2);
                s1 = my_random();
                s2 = my_random();

                vector3d_t offset(0.0f);
                offset[axis_0] = (2.0f * s1 - 1.0f) * cell_width_2;
                offset[axis_1] = (2.0f * s2 - 1.0f) * cell_width_2;

                Point dir;
                dir = color_buffer.get_cell_center(cell) + offset;

                if (_use_jitter)
                {
                    dir = _jitter_rotation_negative * dir;
                }

                dir.normalize();

                float const cos_dir_normal_abs = std::abs(dir * normal);
                assert(cos_dir_normal_abs <= 1.00001f);

                vector3d_t const& wo = dir;
                reflected_color += estimate_light_sample_no_shadow_test(state, lightRay, light_color, surfel, wo);

                averaged_area += area * cos_dir_normal_abs;

                ++combined_samples_n;
            }

            reflected_color *= 1.0f / float(num_samples_per_cell);
            averaged_area   /= float(num_samples_per_cell);

            reflected_color += emission;

            color_buffer.set_color(cell, reflected_color);
            area_buffer .set_color(cell, color_t(averaged_area));
        }

        assert(combined_samples_n == combined_sample_count);
    }



    // one sample on the surfel (center) and many on each cube cell
    void sample_cube_cell_area(Mises_fisher_lobe const& lobe)
    {
        Point const& normal = Point(0, 0, 1);
        float const area = 1.0f;

        std::vector<Cube_cell> const& cube_cells = color_buffer.get_cube_cells();

        int const num_samples_per_cell = 8;
        int const combined_sample_count = num_samples_per_cell * cube_cells.size();

        float const cell_width_2 = 1.0f / color_buffer.get_resolution();

        int combined_samples_n = 0;

        random_t my_random;

        for (unsigned int i = 0; i < cube_cells.size(); ++i)
        {
            Cube_cell const& cell = cube_cells[i];

            color_t reflected_color(0.0f);
            float averaged_area = 0.0f;

            int const axis_0 = Cube_raster_buffer::get_corresponding_axis(0, cell.plane) / 2;
            int const axis_1 = Cube_raster_buffer::get_corresponding_axis(1, cell.plane) / 2;

            for (int j = 0; j  < num_samples_per_cell; ++j)
            {
                float s1, s2;

                // hammersley_2(combined_samples_n, combined_sample_count, s1, s2);
                s1 = my_random();
                s2 = my_random();

                vector3d_t offset(0.0f);
                offset[axis_0] = (2.0f * s1 - 1.0f) * cell_width_2;
                offset[axis_1] = (2.0f * s2 - 1.0f) * cell_width_2;

                Point dir;
                dir = color_buffer.get_cell_center(cell) + offset;

                if (_use_jitter)
                {
                    dir = _jitter_rotation_negative * dir;
                }

                dir.normalize();

                float const cos_dir_normal_abs = std::abs(dir * normal);
                assert(cos_dir_normal_abs <= 1.00001f);

                reflected_color += lobe.evaluate(dir);

                averaged_area += area * cos_dir_normal_abs;

                ++combined_samples_n;
            }

            reflected_color *= 1.0f / float(num_samples_per_cell);
            averaged_area   /= float(num_samples_per_cell);

            color_buffer.set_color(cell, reflected_color);
            area_buffer .set_color(cell, color_t(averaged_area));
        }

        assert(combined_samples_n == combined_sample_count);
    }



    void calc_coefficients_random(renderState_t & state, surfacePoint_t const& surfel, float const area, ray_t const& lightRay, color_t const& light_color)
    {
        // sample_single(state, surfel, area, lightRay, light_color);
        sample_cube_cell_area(state, surfel, area, lightRay, light_color);
    }

    float get_area(Point const& dir) const
    {
        /*
        Cube_cell c;
        bool ok = area_buffer.get_cell(dir, c);
        assert(ok);
        */

        Point dir_tmp = dir;

        if (_use_jitter)
        {
            dir_tmp = _jitter_rotation_positive * dir_tmp;
            // c = get_positive_jittered_cell(c);
        }

        Cube_cell c = area_buffer.get_cell(dir_tmp);

        Color area  = area_buffer.get_color_interpolated(c);
        // Color area = area_buffer.get_color(c);

        return area[0];
    }

    Color get_color(Point const& dir) const
    {
        /*
        Cube_cell c;
        bool ok = color_buffer.get_cell(dir, c);
        assert(ok);
        */

        Point dir_tmp = dir;

        if (_use_jitter)
        {
            dir_tmp = _jitter_rotation_positive * dir_tmp;
            // c = get_positive_jittered_cell(c);
        }

        Cube_cell c = color_buffer.get_cell(dir_tmp);

        return color_buffer.get_color_interpolated(c);
        // return color_buffer.get_color(c);
    }

    void add(Spherical_function const* s)
    {
        Cube_spherical_function const* csf = dynamic_cast<Cube_spherical_function const*>(s);

        assert(color_buffer.get_resolution() == csf->color_buffer.get_resolution());

        std::vector<Cube_cell> const& cube_cells = color_buffer.get_cube_cells();

        //Matrix_3x3 rot_mat = csf->_jitter_rotation_negative;
        //rot_mat = rot_mat * _jitter_rotation_negative;

        for (unsigned int i = 0; i < cube_cells.size(); ++i)
        {
            Cube_cell const& c = cube_cells[i];

            Color new_color;
            float new_area;

            if (_use_jitter)
            {
                Point dir = color_buffer.get_cell_direction(c);
                dir = _jitter_rotation_negative * dir; // correct the direction back into the original, common coord system

                new_color = color_buffer.get_color(c)    + csf->get_color(dir);
                new_area  = area_buffer .get_color(c)[0] + csf->get_area (dir);
            }
            else
            {
                new_color = color_buffer.get_color(c)    + csf->color_buffer.get_color(c);
                new_area  = area_buffer .get_color(c)[0] + csf->area_buffer .get_color(c)[0];
            }

            color_buffer.set_color(c, new_color);
            area_buffer .set_color(c, new_area);

            if (!is_in_range(0.9f, 1.1f, area_buffer.get_color(c)[0] / new_area))
            {
                std::cout << "add failed: " << area_buffer.get_color(c)[0] << " " << new_area << std::endl;
                assert(false);
            }
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

    Spherical_function * clone()
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

    Cube_raster_buffer & get_color_buffer()
    {
        return color_buffer;
    }

    Cube_raster_buffer & get_area_buffer()
    {
        return area_buffer;
    }

    virtual void accept(Spherical_function_visitor * visitor)
    {
        visitor->visit(this);
    }

    virtual std::vector< std::vector<float> > components_to_vectors() const
    {
        std::vector< std::vector<float> > result(6);

        for (unsigned int i = 0; i < 6; ++i)
        {
            std::vector<float> color_component = color_buffer.component_to_vector(i, false);
            std::vector<float> area_component  = area_buffer .component_to_vector(i, true);

            assert(color_component.size() == area_component.size() * 3);

            std::copy(color_component.begin(), color_component.end(), std::back_inserter(result[i]));
            std::copy(area_component .begin(), area_component .end(), std::back_inserter(result[i]));
        }

        return result;
    }

    virtual void from_component_vectors(std::vector< std::vector<float> > const& data)
    {
        assert(data.size() == 6);

        int const entries = (color_buffer.get_resolution() * color_buffer.get_resolution() * 3 + area_buffer.get_resolution() * area_buffer.get_resolution());
        assert(entries == int(data[0].size()));

        for (unsigned int i = 0; i < 6; ++i)
        {
            if (data[i].size() == 0) continue;

            std::vector<float> area_vector;
            std::vector<float> color_vector;

            int const color_end_index = data[i].size() / 4 * 3;
            // int const color_end_index = 0;

            std::copy(data[i].begin(), data[i].begin() + color_end_index, std::back_inserter(color_vector));
            std::copy(data[i].begin() + color_end_index, data[i].end(),   std::back_inserter(area_vector));

            assert(color_vector.size() == area_vector.size() * 3);

            color_buffer.from_component_vector(i, color_vector, false);
            area_buffer. from_component_vector(i, area_vector,  true);
        }
    }

    // return a buffer filled with given color
    static Cube_spherical_function * create_filled_function(Color const& color, const float area)
    {
        Cube_spherical_function * sf = new Cube_spherical_function(4, false);

        std::vector<Cube_cell> const& cube_cells = sf->color_buffer.get_cube_cells();

        for (unsigned int i = 0; i < cube_cells.size(); ++i)
        {
            Cube_cell const& c = cube_cells[i];
            // sf->color_buffer.set_color(c, color);
            sf->area_buffer. set_color(c, color);
        }

        return sf;
    }

    // return a buffer where a single cell on component 0 is filled with given color
    static Cube_spherical_function * create_single_cell_function(Color const& color, const float area, float const r)
    {
        Cube_spherical_function * sf = new Cube_spherical_function(4, false);

        std::vector<Cube_cell> const& cube_cells = sf->color_buffer.get_cube_cells();

        int random_cell_index = r * cube_cells.size() / 6;

        Cube_cell const& c = cube_cells[random_cell_index];
        // sf->color_buffer.set_color(c, color);
        sf->area_buffer. set_color(c, area);

        return sf;
    }

    static Cube_spherical_function * create_area_function(const float area, Point const& normal)
    {
        Cube_spherical_function * sf = new Cube_spherical_function(8, false);

        std::vector<Cube_cell> const& cube_cells = sf->area_buffer.get_cube_cells();

        std::vector<float> original_values;

        for (unsigned int i = 0; i < cube_cells.size(); ++i)
        {
            Cube_cell const& cell = cube_cells[i];

            Point const dir = sf->area_buffer.get_cell_direction(cell);

            float const cos_dir_normal_abs = std::abs(dir * normal);
            assert(cos_dir_normal_abs <= 1.0f);

            sf->area_buffer.set_color(cell, color_t(area * cos_dir_normal_abs));

            original_values.push_back(area * cos_dir_normal_abs);
        }

        return sf;
    }

    static void test()
    {
        Point normal(0.0f, 0.0f, 1.0f);
        Cube_spherical_function * csf = create_area_function(1.0f, normal);

        float area;
        float orig_area;
        Point p;

        p = Point(0.0f, 0.0f, 1.0f);
        area = csf->get_area(p.normalize());
        orig_area = std::abs(p * normal);
        p = Point(1.0f, 0.0f, 0.0f);
        area = csf->get_area(p.normalize());
        orig_area = std::abs(p * normal);
        p = Point(1.0f, 1.0f, 1.0f);
        area = csf->get_area(p.normalize());
        orig_area = std::abs(p * normal);
        p = Point(1.0f, 0.9f, 1.0f);
        area = csf->get_area(p.normalize());
        orig_area = std::abs(p * normal);
        p = Point(1.0f, 0.8f, 1.0f);
        area = csf->get_area(p.normalize());
        orig_area = std::abs(p * normal);
        p = Point(1.0f, 0.7f, 1.0f);
        area = csf->get_area(p.normalize());
        orig_area = std::abs(p * normal);
        p = Point(1.0f, 0.0f, 1.0f);
        area = csf->get_area(p.normalize());
        orig_area = std::abs(p * normal);
    }

    static Cube_spherical_function * create_marked_area_cube()
    {
        Cube_spherical_function * sf = new Cube_spherical_function(8, false);

        {
            Point dir_right(1.0f, 0.5f, 0.0f);

            Cube_cell c_right = sf->area_buffer.get_cell(dir_right);
            sf->color_buffer.set_color(c_right, Color(1.0f, 0.0f, 0.0f));

            Point dir_up(1.0f, 0.0f, 0.5f);
            Cube_cell c_up = sf->area_buffer.get_cell(dir_up);
            sf->color_buffer.set_color(c_up, Color(0.0f, 1.0f, 0.0f));
        }

        {
            Point dir_right(0.5f, 1.0f, 0.0f);

            Cube_cell c_right = sf->area_buffer.get_cell(dir_right);
            sf->color_buffer.set_color(c_right, Color(1.0f, 0.0f, 0.0f));

            Point dir_up(0.0f, 1.0f, 0.5f);
            Cube_cell c_up = sf->area_buffer.get_cell(dir_up);
            sf->color_buffer.set_color(c_up, Color(0.0f, 1.0f, 0.0f));
        }

        {
            Point dir_right(0.5f, 0.0f, 1.0f);

            Cube_cell c_right = sf->area_buffer.get_cell(dir_right);
            sf->color_buffer.set_color(c_right, Color(1.0f, 0.0f, 0.0f));

            Point dir_up(0.0f, 0.5f, 1.0f);
            Cube_cell c_up = sf->area_buffer.get_cell(dir_up);
            sf->color_buffer.set_color(c_up, Color(0.0f, 1.0f, 0.0f));
        }

        {
            Point dir_right(-1.0f, 0.5f, 0.0f);

            Cube_cell c_right = sf->area_buffer.get_cell(dir_right);
            sf->color_buffer.set_color(c_right, Color(1.0f, 0.0f, 0.0f));

            Point dir_up(-1.0f, 0.0f, 0.5f);
            Cube_cell c_up = sf->area_buffer.get_cell(dir_up);
            sf->color_buffer.set_color(c_up, Color(0.0f, 1.0f, 0.0f));
        }

        {
            Point dir_right(0.5f, -1.0f, 0.0f);

            Cube_cell c_right = sf->area_buffer.get_cell(dir_right);
            sf->color_buffer.set_color(c_right, Color(1.0f, 0.0f, 0.0f));

            Point dir_up(0.0f, -1.0f, 0.5f);
            Cube_cell c_up = sf->area_buffer.get_cell(dir_up);
            sf->color_buffer.set_color(c_up, Color(0.0f, 1.0f, 0.0f));
        }

        {
            Point dir_right(0.5f, 0.0f, -1.0f);

            Cube_cell c_right = sf->area_buffer.get_cell(dir_right);
            sf->color_buffer.set_color(c_right, Color(1.0f, 0.0f, 0.0f));

            Point dir_up(0.0f, 0.5f, -1.0f);
            Cube_cell c_up = sf->area_buffer.get_cell(dir_up);
            sf->color_buffer.set_color(c_up, Color(0.0f, 1.0f, 0.0f));
        }

        return sf;
    }

    static Cube_spherical_function * create_from_mises_fisher(Mises_fisher_lobe const& lobe, int const resolution, bool const use_jitter, float const weight)
    {
        Cube_spherical_function * sf = new Cube_spherical_function(resolution, use_jitter);
        // std::vector<Cube_cell> const& cells = sf->get_area_buffer().get_cube_cells();

        sf->sample_cube_cell_area(lobe);

        /*
        for (unsigned int i = 0; i < cells.size(); ++i)
        {
            Cube_cell const& c = cells[i];

            vector3d_t const dir = sf->get_area_buffer().get_cell_direction(c);
            float const value = lobe.evaluate(dir);

            sf->get_area_buffer().set_color(c, color_t(value * weight));
        }
        */

        return sf;
    }

    int get_resolution() const
    {
        return color_buffer.get_resolution();
    }

private:
    Cube_raster_buffer color_buffer;
    Cube_raster_buffer area_buffer;

    bool _use_jitter;
    float _jitter_x, _jitter_y;
    Matrix_3x3 _jitter_rotation_positive;
    Matrix_3x3 _jitter_rotation_negative;
};


// how to replace the data in the final tree?

class Indexed_Cube_Spherical_function : public Spherical_function
{
public:
    Indexed_Cube_Spherical_function(std::vector< std::vector<float> > const& unpacked_dictionary, std::vector<Spherical_function*> const* dictionary, Cube_spherical_function const* csf) :
        _dictionary(dictionary)
    {
        // use the csf to find for each side the corresponding entry and fill _index_per_side
        std::vector< std::vector<float> > components = csf->components_to_vectors();

        Distance_function distance_fnc;

        for (unsigned int i = 0; i < components.size(); ++i)
        {
            float lowest_distance = 1e10f;
            int closest_dictionary_index;

            std::vector<float> const& component = components[i];

            for (unsigned int j = 0; j < unpacked_dictionary.size(); ++j)
            {
                assert(component.size() == unpacked_dictionary[j].size());

                float const distance = distance_fnc(unpacked_dictionary[j], component);

                if (distance < lowest_distance)
                {
                    lowest_distance = distance;
                    closest_dictionary_index = j;
                }
            }

            _index_per_side.push_back(closest_dictionary_index);
        }
    }


    virtual ~Indexed_Cube_Spherical_function() {}

    virtual void calc_coefficients_random(renderState_t & , surfacePoint_t const& , float const , ray_t const& , color_t const& )
    {
        std::cout << "Indexed_Cube_Spherical_function::calc_coefficients_random() shouldn't be used" << std::endl;
    }

    virtual float get_area(Point const& dir) const
    {
        Cube_spherical_function const* csf = static_cast<Cube_spherical_function*>((*_dictionary)[0]);
        Cube_cell c = csf->get_area_buffer().get_cell(dir);
        //return color_buffer.get_color_interpolated(c);

        int const dictionary_index_of_side = _index_per_side[c.plane];
        Cube_spherical_function const* dictionary_entry = static_cast<Cube_spherical_function*>((*_dictionary)[dictionary_index_of_side]);

        Abstract_frame_buffer * component = dictionary_entry->get_area_buffer().get_component(0);
        return component->get_color(c.pos[0], c.pos[1])[0];
    }

    virtual Color get_color(Point const& dir) const
    {
        Cube_spherical_function const* csf = static_cast<Cube_spherical_function*>((*_dictionary)[0]);
        Cube_cell c = csf->get_color_buffer().get_cell(dir);
        //return color_buffer.get_color_interpolated(c);

        int const dictionary_index_of_side = _index_per_side[c.plane];
        Cube_spherical_function const* dictionary_entry = static_cast<Cube_spherical_function*>((*_dictionary)[dictionary_index_of_side]);

        Abstract_frame_buffer * component = dictionary_entry->get_color_buffer().get_component(0);
        return component->get_color(c.pos[0], c.pos[1]);
    }

    Cube_raster_buffer const& get_color_buffer(int const component) const
    {
        assert(component >= 0 && component < 6);

        int const dictionary_index_of_side = _index_per_side[component];
        Cube_spherical_function const* dictionary_entry = static_cast<Cube_spherical_function*>((*_dictionary)[dictionary_index_of_side]);
        return dictionary_entry->get_color_buffer();
    }

    Cube_raster_buffer const& get_area_buffer(int const component) const
    {
        assert(component >= 0 && component < 6);

        int const dictionary_index_of_side = _index_per_side[component];
        Cube_spherical_function const* dictionary_entry = static_cast<Cube_spherical_function*>((*_dictionary)[dictionary_index_of_side]);
        return dictionary_entry->get_area_buffer();
    }

    virtual void add(Spherical_function const* /* s */) { std::cout << "Indexed_Cube_Spherical_function::add() shouldn't be used" << std::endl; }
    virtual void normalize(float const /* factor */) { std::cout << "Indexed_Cube_Spherical_function::normalize() shouldn't be used" << std::endl; }

    virtual void accept(Spherical_function_visitor * visitor)
    {
        visitor->visit(this);
    }

    virtual std::vector< std::vector<float> > components_to_vectors() const
    {
        std::vector< std::vector<float> > result(6);

        for (unsigned int i = 0; i < 6; ++i)
        {
            int const dict_index = _index_per_side[i];
            std::vector<float> component  = (*_dictionary)[dict_index]->components_to_vectors()[0];
            std::copy(component.begin(), component.end(), std::back_inserter(result[i]));
        }

        return result;
    }

    virtual void from_component_vectors(std::vector< std::vector<float> > const& /* data */) {}

    //virtual std::vector<float> into_vector() { return std::vector<float>(); }
    //virtual void from_vector(std::vector<float> const& data) { }

    int get_component_index(int const component) const
    {
        return _index_per_side[component];
    }

    virtual Spherical_function * clone() { return NULL; }

    int get_resolution() const
    {
        Cube_spherical_function const* dictionary_entry = static_cast<Cube_spherical_function*>((*_dictionary)[0]);
        return dictionary_entry->get_resolution();
    }

private:
    std::vector<Spherical_function*> const* _dictionary;
    std::vector<int> _index_per_side;
};



class Mises_Fisher_spherical_function : public Spherical_function
{
public:

    virtual void calc_coefficients_random(renderState_t & state, surfacePoint_t const& surfel, float const area, ray_t const& lightRay, color_t const& light_color)
    {
        std::cout << "Mises_Fisher_spherical_function::calc_coefficients_random(): not used" << std::endl;
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

                test_area_sh += get_area(dir);

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
        float result = 0.0f;

        for (unsigned int i = 0; i < _area_lobes.size(); ++i)
        {
            result += _area_lobes[i].evaluate(dir);
        }

        // result /= float(_area_lobes.size());

        return result;
    }

    virtual Color get_color(Point const& dir) const
    {
        float result = 0.0f;

        for (unsigned int i = 0; i < _color_lobes.size(); ++i)
        {
            result += _color_lobes[i].evaluate(dir);
        }

        // result /= float(_color_lobes.size());

        return color_t(result);
    }

    virtual void add(Spherical_function const* /* s */)
    {
        std::cout << "Mises_Fisher_spherical_function::add(): Should not be used" << std::endl;
        assert(false);

        // Mises_Fisher_spherical_function const* vmf = dynamic_cast<Mises_Fisher_spherical_function const*>(s);

    }

    void normalize(float const /* factor */)
    {
        std::cout << "Mises_Fisher_spherical_function::normalize(): Should not be used" << std::endl;
        assert(false);
    }

    void generate_from_cube_spherical_function(Cube_spherical_function const* csf, int const num_area_lobes, int const num_color_lobes)
    {
        _area_lobes  = generate_from_cube(csf->get_area_buffer(),  num_area_lobes);
        _color_lobes = generate_from_cube(csf->get_color_buffer(), num_color_lobes);
    }

    std::vector<Mises_fisher_lobe> generate_from_cube(Cube_raster_buffer const& distribution, int const lobeCount)
    {
        std::vector<Cube_cell> const& cells = distribution.get_cube_cells();

        const int size = cells.size();

        const int num_iterations = 30;

        std::vector<Mises_fisher_lobe> lobes(lobeCount);

        std::vector<Cube_raster_buffer> likelihoods(lobeCount);
        for (int i = 0; i < lobeCount; i++)
        {
            likelihoods[i].setup_surfel_buffer(distribution.get_resolution());
        }

        Cube_raster_buffer totalLikelihoods;
        totalLikelihoods.setup_surfel_buffer(distribution.get_resolution());

        // Random start guess
        random_t my_random;
        for (int i = 0; i < lobeCount; i++)
        {
            lobes[i].randomize(my_random);
            std::cout << "Mises_Fisher_spherical_function::generate_from_cube(): " << i << " " << lobes[i] << std::endl;
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
                for(int k = 0; k < lobeCount; k++)
                {
                    likelihoods[k].set_color(c, color_t(lobes[k].evaluate(direction)));
                }
            }

            for (int j = 0; j < size; j++)
            {
                // 2.) Sum over all lobes
                Cube_cell const& c = cells[j];

                totalLikelihoods.set_color(c, color_t(0.0f));
                for(int k = 0; k < lobeCount; k++)
                {
                    float current_likelihood = totalLikelihoods.get_color(c)[0];
                    current_likelihood += likelihoods[k].get_color(c)[0];
                    totalLikelihoods.set_color(c, color_t(current_likelihood));
                }
            }

            for (int j = 0; j < size; j++)
            {
                // 3.) Divide by total likelihood
                Cube_cell const& c = cells[j];
                const float pixel_value = distribution.get_color(c)[0];// + 0.01f;
                const float solid_angle = distribution.get_solid_angle(c); // / (4.0f * M_PI);
                for(int k = 0; k < lobeCount; k++)
                {
                    likelihoods[k].set_color(c, pixel_value * solid_angle * likelihoods[k].get_color(c)[0] / totalLikelihoods.get_color(c)[0]);
                }
            }

            // M-step
            for (int lobe_index = 0; lobe_index < lobeCount; lobe_index++)
            {
                lobes[lobe_index].fit(likelihoods[lobe_index]);

                std::cout << "step: " << step << " " << lobes[lobe_index] << std::endl;
            }
        }

        // calculate lobe weight
        for (int cell_index = 0; cell_index < size; cell_index++)
        {
            // 2.) Sum over all lobes
            Cube_cell const& c = cells[cell_index];

            totalLikelihoods.set_color(c, color_t(0.0f));
            for(int k = 0; k < lobeCount; k++)
            {
                float current_likelihood = totalLikelihoods.get_color(c)[0];
                current_likelihood += likelihoods[k].get_color(c)[0];
                totalLikelihoods.set_color(c, color_t(current_likelihood));
            }
        }

        for (int lobe_index = 0; lobe_index < lobeCount; lobe_index++)
        {
            vector3d_t const& mean = lobes[lobe_index].mean_dir;

            if (mean.length() < 0.001f)
            {
                lobes[lobe_index].weight = 0.0f;
                continue;
            }

            Cube_cell c = likelihoods[lobe_index].get_cell(mean);
            float const likelihood = likelihoods[lobe_index].get_color_interpolated(c)[0];
            float const total_likelihood = totalLikelihoods.get_color_interpolated(c)[0];
            float const amplitude  = distribution.get_color_interpolated(c)[0];

            lobes[lobe_index].weight = (likelihood / total_likelihood) * amplitude;

            // std::cout << "weight: " << lobes[lobe_index].weight << std::endl;
        }

        return lobes;
    }


    Spherical_function * clone()
    {
        return new Mises_Fisher_spherical_function(*this);
    }

    void add_area_lobe(Mises_fisher_lobe const& lobe)
    {
        _area_lobes.push_back(lobe);
    }

    std::vector<Mises_fisher_lobe> const& get_area_lobes() const
    {
        return _area_lobes;
    }

    std::vector<Mises_fisher_lobe> const& get_color_lobes() const
    {
        return _color_lobes;
    }

    virtual void accept(Spherical_function_visitor * visitor)
    {
        visitor->visit(this);
    }

    static Mises_Fisher_spherical_function * test()
    {
        vector3d_t mean_dir(0.0f, 0.5f, 0.5f);
        mean_dir.normalize();
        float const concentration = 30;
        float const weight = 1.0f;

        Mises_fisher_lobe lobe_0(mean_dir, concentration, weight);
        Mises_fisher_lobe lobe_1(vector3d_t(1.0f, 0.0f, 0.0f), 10.0f, weight);

        std::cout << "test lobe 0: " << lobe_0 << std::endl;
        std::cout << "test lobe 1: " << lobe_1 << std::endl;

        Cube_spherical_function * distribution   = Cube_spherical_function::create_from_mises_fisher(lobe_0, 8, false, 2.0f);
        Cube_spherical_function * distribution_1 = Cube_spherical_function::create_from_mises_fisher(lobe_1, 8, false, 1.0f);

        // distribution->add(distribution_1);


        Mises_Fisher_spherical_function * mf_sf = new Mises_Fisher_spherical_function;
        mf_sf->generate_from_cube_spherical_function(distribution, 10, 3);

        std::cout << "area lobes:" << std::endl;
        for (unsigned int i = 0; i < mf_sf->get_area_lobes().size(); ++i)
        {
            std::cout << i << " " << mf_sf->get_area_lobes()[i] << std::endl;
        }

        std::cout << "color lobes:" << std::endl;
        for (unsigned int i = 0; i < mf_sf->get_color_lobes().size(); ++i)
        {
            std::cout << i << " " << mf_sf->get_color_lobes()[i] << std::endl;
        }

        return mf_sf;
    }

private:
    std::vector<Mises_fisher_lobe> _area_lobes;
    std::vector<Mises_fisher_lobe> _color_lobes;
};



class Spherical_function_factory
{
public:
    virtual Spherical_function * create() const = 0;
    virtual std::string get_name() const = 0;
};

class Spherical_harmonics_factory : public Spherical_function_factory
{
public:
    Spherical_function * create() const
    {
        return new GiSphericalHarmonics(true, 3);
    }

    std::string get_name() const
    {
        return "SH";
    }
};


class Cube_spherical_function_factory : public Spherical_function_factory
{
public:
    Cube_spherical_function_factory(int const resolution, bool use_jitter) :
        _resolution(resolution),
        _use_jitter(use_jitter)
    { }

    Spherical_function * create() const
    {
        return new Cube_spherical_function(_resolution, _use_jitter);
    }

    std::string get_name() const
    {
        return "Cube";
    }

private:
    int _resolution;
    bool _use_jitter;
};


struct Spherical_node_representation
{
    Spherical_function * area;
    Spherical_function * color;
};


__END_YAFRAY

#endif // SPHERICAL_HARMONICS_H
