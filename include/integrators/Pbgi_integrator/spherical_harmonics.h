#ifndef SPHERICAL_HARMONICS_H
#define SPHERICAL_HARMONICS_H

#include <cmath>
//#include <stdint.h>


#include <tr1/functional>

#include <boost/optional/optional.hpp>
#include <boost/array.hpp>

#include <yafray_config.h>
#include <utilities/mcqmc.h>
#include <utilities/sample_utils.h>
#include <core_api/scene.h>
#include <core_api/surface.h>
#include <core_api/material.h>

#include <integrators/Pbgi_integrator/CubeRasterBuffer.h>

__BEGIN_YAFRAY

static float factorial_table[] =
{
    1.0f,
    1.0f,
    2.0f,
    6.0f,
    24.0f,
    120.0f,
    720.0f,
    5040.0f,
    40320.0f,
    362880.0f,
    3628800.0f,
    39916800.0f,
    479001600.0f,
    6227020800.0f,
    87178291200.0f,
    1307674368000.0f,
    20922789888000.0f,
    355687428096000.0f,
    6402373705728000.0f,
    121645100408832000.0f,
    2432902008176640000.0f,
    51090942171709440000.0f,
    1124000727777607680000.0f,
    25852016738884976640000.0f,
    620448401733239439360000.0f,
    15511210043330985984000000.0f,
    403291461126605635584000000.0f,
    10888869450418352160768000000.0f,
    304888344611713860501504000000.0f,
    8841761993739701954543616000000.0f,
    265252859812191058636308480000000.0f,
    8222838654177922817725562880000000.0f,
    263130836933693530167218012160000000.0f,
    8683317618811886495518194401280000000.0f };


inline
color_t estimate_light_sample_no_shadow_test(renderState_t &state, ray_t const& lightRay, color_t const& light_color, const surfacePoint_t &surfel, const vector3d_t &wo)
{
    const material_t *material = surfel.material;

    // handle lights with delta distribution, e.g. point and directional lights
    color_t surfCol = material->eval(state, surfel, wo, lightRay.dir, BSDF_ALL);
    // color_t transmitCol = scene->volIntegrator->transmittance(state, lightRay);
    color_t col = surfCol * light_color * std::max(0.0f, surfel.N*lightRay.dir); // * transmitCol;

    return col;
}

struct Cone
{
    vector3d_t dir;
    float cos_angle;
};

template <class Data>
class Abstract_spherical_function_estimator
{
public:
    virtual Data get_value(vector3d_t const& dir) const = 0;
};

template <class Data>
class Spherical_function_area_estimator : public Abstract_spherical_function_estimator<Data>
{
public:
    Spherical_function_area_estimator(vector3d_t const& normal, float const area) :
        _normal(normal),
        _area(area)
    {  }

    virtual Data get_value(vector3d_t const& dir) const
    {
        return std::max(0.0f, _normal * dir) * _area;
    }

private:
    vector3d_t _normal;
    float _area;
};

template <class Data>
class Spherical_function_light_color_estimator : public Abstract_spherical_function_estimator<Data>
{
public:
    Spherical_function_light_color_estimator(renderState_t & state, surfacePoint_t const& surfel, ray_t const& lightRay, color_t const& light_color, float const area) :
        _state(state),
        _surfel(surfel),
        _lightRay(lightRay),
        _light_color(light_color),
        _area(area)
    {  }

    virtual Data get_value(vector3d_t const& dir) const
    {
        return estimate_light_sample_no_shadow_test(_state, _lightRay, _light_color, _surfel, dir) * std::max(0.0f, _surfel.N * dir) * _area;
    }

private:
    renderState_t & _state;
    surfacePoint_t const& _surfel;
    ray_t const& _lightRay;
    color_t const& _light_color;
    float _area;
};

template <class Data>
class Spherical_function_emit_color_estimator : public Abstract_spherical_function_estimator<Data>
{
public:
    Spherical_function_emit_color_estimator(renderState_t & state, surfacePoint_t const& surfel, float const area) :
        _state(state),
        _surfel(surfel),
        _area(area)
    { }

    virtual Data get_value(vector3d_t const& dir) const
    {
        return _surfel.material->emit(_state, _surfel, dir) * _area;
    }

private:
    renderState_t & _state;
    surfacePoint_t const& _surfel;
    float _area;
};



inline float square(float x) { return x * x; }

template <class Data, int N_bands = 3>
class Spherical_harmonics
{
public:
    typedef std::tr1::function<float(vector3d_t const&)> Sh_coefficient_func;

    static const int Num_bands = N_bands;
    static const int Num_coeffs = Num_bands * Num_bands;

    Spherical_harmonics() :
        _exact(false)
    {
        init();
    }

    Spherical_harmonics(bool const exact)
    {
        if (N_bands != 3) assert(exact);

        _exact = exact;

        init();
    }

    void init()
    {
        _sh_coefficients.assign(Data(0.0f));
    }

    void calc_coefficients(Abstract_spherical_function_estimator<Data> const* estimator)
    {
        random_t random;

        float const res_u = 20.0f;
        float const res_v = 40.0f;

        for (int u = 0; u < res_u; ++u)
        {
            for (int v = 0; v < res_v; ++v)
            {
//                float x = (u + random()) / res_u;
//                float y = (v + random()) / res_v;
//                //                float x = u / res_u;
//                //                float y = v / res_v;

//                float theta = 2.0f * std::acos(std::sqrt(1.0f - x));
//                float phi = 2.0f * M_PI * y;


//                // dir corresponds to w_o, the outgoing direction
//                Point dir(std::sin(theta) * std::cos(phi),
//                          std::sin(theta) * std::sin(phi),
//                          std::cos(theta));

                float const v_tmp = 2.0f * M_PI * (v + random()) / res_v;
                float const u_tmp = 1.0f - 2.0f * (u + random()) / res_u;

                float const z = u_tmp;
                float const x = std::sqrt(1.0f - z * z) * std::cos(v_tmp);
                float const y = std::sqrt(1.0f - z * z) * std::sin(v_tmp);

                vector3d_t dir(x, y, z);

                float const theta = std::acos(z);
                float const phi = std::atan2(y, x);


                // std::cout << "---" << std::endl;

                Data const function_value = estimator->get_value(dir);

                if (_exact)
                {
                    for (int l = 0; l < Num_bands; ++l)
                    {
                        for(int m = -l; m <= l; ++m)
                        {
                            int sh_index = l * (l + 1) + m;

                            float Y_l_m = SH(l, m, theta, phi);
                            float base_coeff = Y_l_m;

                            _sh_coefficients[sh_index] += base_coeff * function_value;
                        }
                    }
                }
                else
                {
                    for (int sh_index = 0; sh_index < Num_coeffs; ++sh_index)
                    {
                        float const Y_l_m = SH_precomputed(sh_index, dir);
                        float const base_coeff = Y_l_m;

                        _sh_coefficients[sh_index]  += base_coeff * function_value;
                    }
                }
            }
        }

        float inv = 4.0f * M_PI / (res_u * res_v);

        for (int sh_index = 0; sh_index < Num_coeffs; ++sh_index)
        {
            _sh_coefficients[sh_index] *= inv;
        }
    }


    void get_precalculated_coefficients(Cube_raster_buffer< Spherical_harmonics<float> > const& normal_map,
                                        Data const& scale,
                                        vector3d_t const& surface_normal)
    {
        Cube_cell normal_cell = normal_map.get_cell(surface_normal);

        Spherical_harmonics<float> const& sph_harmonics = normal_map.get_data(normal_cell);

        boost::array<float, Num_coeffs> const& precalc_coefficients = sph_harmonics.get_coefficients();

        assert(precalc_coefficients.size() == _sh_coefficients.size());

        for (size_t i = 0; i < _sh_coefficients.size(); ++i)
        {
            _sh_coefficients[i] = scale * precalc_coefficients[i];
        }
    }


    Data get_value(vector3d_t const& dir) const
    {
        if (!_exact)
        {
            Data result = 0.0f;

            for (int sh_index = 0; sh_index < Num_coeffs; ++sh_index)
            {
                float Y_l_m = SH_precomputed(sh_index, dir);

                result += _sh_coefficients[sh_index] * Y_l_m;
            }

            return result;
        }

        float theta = fAcos(dir.z);
        float phi = std::atan2(dir.y, dir.x);
        return get_value(theta, phi);
    }

    void add(Spherical_harmonics<Data> const& s)
    {
        *this += s;
    }

    void operator+=(Spherical_harmonics const& sh)
    {
        assert(sh._exact == _exact);

        for (int i = 0; i < Num_coeffs; ++i)
        {
            _sh_coefficients[i] += sh._sh_coefficients[i];
        }
    }

    friend Spherical_harmonics operator+(Spherical_harmonics const& lhs, Spherical_harmonics const& rhs)
    {
        assert(lhs._exact == rhs._exact);

        Spherical_harmonics result(lhs._exact);

        for (int j = 0; j < Num_coeffs; ++j)
        {
            result._sh_coefficients[j] = lhs._sh_coefficients[j] + rhs._sh_coefficients[j];
        }

        return result;
    }

    friend Spherical_harmonics operator*(Spherical_harmonics<float> const& sh, Data const& scale)
    {
        Spherical_harmonics result(sh._exact);

        for (int j = 0; j < Num_coeffs; ++j)
        {
            result._sh_coefficients[j] = sh._sh_coefficients[j] * scale;
        }

        return result;
    }

    /*
    friend Spherical_harmonics<float> operator*(Spherical_harmonics const& sh, float const scale)
    {
        Spherical_harmonics result(sh._exact, sh._bands);

        for (int j = 0; j < sh._bands * sh._bands; ++j)
        {
            result._sh_coefficients[j] = sh._sh_coefficients[j] * scale;
        }

        return result;
    }

    friend Spherical_harmonics<color_t> operator*(Spherical_harmonics const& sh, color_t const scale)
    {
        Spherical_harmonics result(sh._exact, sh._bands);

        for (int j = 0; j < sh._bands * sh._bands; ++j)
        {
            result._sh_coefficients[j] = sh._sh_coefficients[j] * scale;
        }

        return result;
    }
    */

    /*
    std::vector< float > to_vector() const
    {
        std::vector<float> result;

        for (std::size_t i = 0; i < sh_coefficients.size(); ++i)
        {
            std::vector<float> value_vector = conv_to_vector(sh_coefficients[i]);

            for (std::size_t j = 0; j < value_vector.size(); ++j)
            {
                result.push_back(value_vector[j]);
            }
        }

        return result;
    }

    void from_vector(std::vector<float> const& data)
    {
        sh_coefficients = conv_from_vector<Data>(data);
        assert(int(sh_coefficients.size()) == bands * bands);
    }
    */



    static boost::array<Sh_coefficient_func, Num_coeffs> initialize_sh_functions()
    {
        boost::array<Sh_coefficient_func, Num_coeffs> sh_functions;

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

        sh_functions[0] = &Spherical_harmonics::SH_precomputed_0;
        sh_functions[1] = &Spherical_harmonics::SH_precomputed_1;
        sh_functions[2] = &Spherical_harmonics::SH_precomputed_2;
        sh_functions[3] = &Spherical_harmonics::SH_precomputed_3;
        sh_functions[4] = &Spherical_harmonics::SH_precomputed_4;
        sh_functions[5] = &Spherical_harmonics::SH_precomputed_5;
        sh_functions[6] = &Spherical_harmonics::SH_precomputed_6;
        sh_functions[7] = &Spherical_harmonics::SH_precomputed_7;
        sh_functions[8] = &Spherical_harmonics::SH_precomputed_8;

        return sh_functions;
    }


    Data const& get_coefficient(int index) const
    {
        return _sh_coefficients[index];
    }

    boost::array<Data, Num_coeffs> const& get_coefficients() const
    {
        return _sh_coefficients;
    }


private:
    Data get_value(float const theta, float const phi) const
    {
        if (!_exact)
            std::cout << "shouldn't be used!" << std::endl;

        Data result = 0.0f;

        for (int l = 0; l < Num_bands; ++l)
        {
            for(int m = -l; m <= l; ++m)
            {
                int sh_index = l * (l + 1) + m;

                float Y_l_m = SH(l, m, theta, phi);

                result += _sh_coefficients[sh_index] * Y_l_m;
            }
        }

        return result;
    }


    float factorial(int num) const
    {
        if (num < 34) return factorial_table[num];

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

    float SH_precomputed(int const index, vector3d_t const& dir) const
    {
        Sh_coefficient_func const& f = _sh_functions[index];

        return f(dir);
    }

    // http://en.wikipedia.org/wiki/Table_of_spherical_harmonics

    static float SH_precomputed_0(vector3d_t const& /* dir */)
    {
        return 0.28209;
    }

    static float SH_precomputed_1(vector3d_t const& dir)
    {
        return 0.34549 * dir[0];
    }

    static float SH_precomputed_2(vector3d_t const& dir)
    {
        return 0.34549 * dir[1];
    }

    static float SH_precomputed_3(vector3d_t const& dir)
    {
        return 0.34549 * dir[2];
    }

    static float SH_precomputed_7(vector3d_t const& dir)
    {
        return 0.31539 * (2 * square(dir[2]) - square(dir[0]) - square(dir[1]));
    }

    static float SH_precomputed_5(vector3d_t const& dir)
    {
        return 1.0925 * dir[1] * dir[2];
    }

    static float SH_precomputed_4(vector3d_t const& dir)
    {
        return 1.0925 * dir[2] * dir[0];
    }

    static float SH_precomputed_6(vector3d_t const& dir)
    {
        return 1.0925 * dir[0] * dir[1];
    }

    static float SH_precomputed_8(vector3d_t const& dir)
    {
        return 0.54627 * (square(dir[0]) - square(dir[1]));
    }

    bool _exact;

    boost::array<Data, Num_coeffs> _sh_coefficients;

    static boost::array<Sh_coefficient_func, Num_coeffs> _sh_functions;
};



template <class Data>
class Spherical_harmonics_factory
{
public:
    Spherical_harmonics_factory() :
        _exact(false)
    { }


    Spherical_harmonics_factory(bool const exact) :
        _exact(exact)
    { }

    Spherical_harmonics<Data> create() const
    {
        return Spherical_harmonics<Data>(_exact);
    }

private:
    bool _exact;
};


Spherical_harmonics<float> get_precalculated_sh(Cube_raster_buffer< Spherical_harmonics<float> > const& normal_map,
                                                vector3d_t const& surface_normal);




__END_YAFRAY

#endif // SPHERICAL_HARMONICS_H
