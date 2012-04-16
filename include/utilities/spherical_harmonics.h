#ifndef SPHERICAL_HARMONICS_H
#define SPHERICAL_HARMONICS_H

#include <cmath>
//#include <stdint.h>

#include <tr1/functional>

#include <boost/serialization/export.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/serialization.hpp>

#include <yafray_config.h>
#include <utilities/mcqmc.h>
#include <utilities/sample_utils.h>
#include <core_api/scene.h>
#include <core_api/surface.h>
#include <core_api/material.h>
#include <utilities/Quaternion_Matrix.h>
#include <utilities/Dictionary_utils.h>


#include <integrators/kmeans.hpp>

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


template <class Color, class Point>
Color estimate_light_sample_no_shadow_test(renderState_t &state, ray_t const& lightRay, Color const& light_color, const surfacePoint_t &surfel, const Point &wo)
{
    const material_t *material = surfel.material;

    // handle lights with delta distribution, e.g. point and directional lights
    Color surfCol = material->eval(state, surfel, wo, lightRay.dir, BSDF_ALL);
    // color_t transmitCol = scene->volIntegrator->transmittance(state, lightRay);
    Color col = surfCol * light_color * std::max(0.0f, surfel.N*lightRay.dir); // * transmitCol;

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
    std::vector<Cone> const& get_main_cones() const { return _main_cones; }

protected:
    std::vector<Cone> _main_cones;

};

template <class Data>
class Spherical_function_area_estimator : public Abstract_spherical_function_estimator<Data>
{
public:
    Spherical_function_area_estimator(vector3d_t const& normal, float const area) :
        _normal(normal),
        _area(area)
    {
        Cone c;
        // c.cos_angle = std::cos(M_PI / 2.0f);
        c.cos_angle = (M_PI / 2.0f);
        c.dir = normal;
        this->_main_cones.push_back(c);
    }

    virtual Data get_value(vector3d_t const& dir) const
    {
//        return std::abs(_normal * dir) * _area;
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
    Spherical_function_light_color_estimator(renderState_t & state, surfacePoint_t const& surfel, ray_t const& lightRay, color_t const& light_color) :
        _state(state),
        _surfel(surfel),
        _lightRay(lightRay),
        _light_color(light_color)
    {
        if (lightRay.dir * surfel.N <= 0.0f)
        {
            return;
        }

        Cone diffuse_cone;
        // diffuse_cone.cos_angle = std::cos(M_PI / 2.0f);
        diffuse_cone.cos_angle = (M_PI / 2.0f);
        diffuse_cone.dir = surfel.N;
        this->_main_cones.push_back(diffuse_cone);

        vector3d_t const reflected_dir = 2.0f * (lightRay.dir * surfel.N) * surfel.N - lightRay.dir;

        Cone specular_cone;
        float lobe_angle = M_PI / 2.0f;

        vector3d_t tangent, bitangent;
        createCS(reflected_dir, tangent, bitangent);

        std::string fuckyou;

        // FIXME: lobe_angle is okay, but get_value stays 0, maybe prob with the reflected_dir or the material?

        while (get_value(reflected_dir * std::cos(lobe_angle) + tangent * std::sin(lobe_angle)) < 0.1f && lobe_angle > 0.01f)
        {
            // std::cout << lobe_angle << " " << get_value(reflected_dir * std::cos(lobe_angle) + tangent * std::sin(lobe_angle)) << std::endl;
            lobe_angle /= 2.0f;
//            if (lobe_angle < 0.0001f)
//            {
//                std::cout << "ref: " << reflected_dir << " light: " << lightRay.dir << " N: " << surfel.N << std::endl;
//                assert(false);
//            }
        }

        lobe_angle += lobe_angle / 2.0f; // add some back to get to an area where the value is probably already a bit more than 0

        // specular_cone.cos_angle = std::cos(lobe_angle);
        specular_cone.cos_angle = (lobe_angle);
        specular_cone.dir = reflected_dir;
        this->_main_cones.push_back(specular_cone);
    }

    virtual Data get_value(vector3d_t const& dir) const
    {
        return estimate_light_sample_no_shadow_test(_state, _lightRay, _light_color, _surfel, dir);
    }

private:
    renderState_t & _state;
    surfacePoint_t const& _surfel;
    ray_t const& _lightRay;
    color_t const& _light_color;
};

template <class Data>
class Spherical_function_emit_color_estimator : public Abstract_spherical_function_estimator<Data>
{
public:
    Spherical_function_emit_color_estimator(renderState_t & state, surfacePoint_t const& surfel) :
        _state(state),
        _surfel(surfel)
    { }

    virtual Data get_value(vector3d_t const& dir) const
    {
        return _surfel.material->emission(_state, _surfel, dir);
    }

private:
    renderState_t & _state;
    surfacePoint_t const& _surfel;
};


typedef vector3d_t Point;
typedef color_t Color;

class Splat_cube_raster_buffer;
template <class Data>
class Cube_raster_buffer;

template <class Data>
class Spherical_harmonics;

template <class Data>
class Spherical_function;

template <class Data>
class Cube_spherical_function;

template <class Data>
class Mises_Fisher_spherical_function;

template <class Data>
class Indexed_Cube_Spherical_function;

template <class Data>
class Indexed_spherical_function;

template <class Data>
class Spherical_function_visitor
{
public:
    virtual void visit(Spherical_harmonics<Data> * sf) = 0;
    virtual void visit(Cube_spherical_function<Data> * sf) = 0;
    virtual void visit(Mises_Fisher_spherical_function<Data> * sf) = 0;
    virtual void visit(Indexed_spherical_function<Data> * sf) = 0;
    // virtual void visit(Indexed_Cube_Spherical_function<Data> * sf) = 0;
    virtual void identify() { std::cout << "Spherical_function_visitor" << std::endl; }
};


struct Spherical_node_representation
{
//    ~Spherical_node_representation();

    Spherical_function<float> * area;
    Spherical_function<color_t> * color;
};


template <class Data>
class Spherical_function
{
public:
    virtual ~Spherical_function() {}

    virtual void calc_coefficients_random(Abstract_spherical_function_estimator<Data> const* estimator) = 0;

    virtual Data get_value(Point const& dir) const  = 0;

    virtual void add(Spherical_function const* s) = 0;
    virtual void normalize(float const factor) = 0;

    virtual void add_light(Splat_cube_raster_buffer const& /* fb_in */) { std::cout << "Spherical_function::add_light(); " << std::endl; }

    virtual void accept(Spherical_function_visitor<Data> * visitor) = 0;

    virtual std::vector< float > to_vector() const { return std::vector<float>(); }
    virtual void from_vector(std::vector< float > const& /* data */) {}

    virtual void get_precalculated_coefficients(Cube_raster_buffer<Spherical_function<float> *> const& normal_map,
                                                Data const& scale,
                                                vector3d_t const& surface_normal) {}

    virtual std::vector< std::vector<float> > components_to_vectors() const { return std::vector< std::vector<float> >(); }
    virtual void from_component_vectors(std::vector< std::vector<float> > const& /* data */) {}

    virtual Spherical_function * clone() const = 0;

    /*
    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive & , const unsigned int )
    {
    }
    */
};

inline float square(float x) { return x * x; }

template <class Data>
class Spherical_harmonics : public Spherical_function<Data>
{
public:
    // typedef float (GiSphericalHarmonics::*Sh_coefficient_func)(Point const& dir) const; // sh_coefficient_func;
    typedef std::tr1::function<float(Point const&)> Sh_coefficient_func;

    Spherical_harmonics() :
        bands(3),
        exact(false)
    {
        init();
    }

    Spherical_harmonics(bool exact, int b)
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
        sh_coefficients = std::vector<Data>(bands * bands, Data(0.0f));
    }




    virtual void calc_coefficients_random(Abstract_spherical_function_estimator<Data> const* estimator)
    {
        random_t random;

        float const res_u = 10.0f;
        float const res_v = 20.0f;

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

                Point dir(x, y, z);

                float const theta = std::acos(z);
                float const phi = std::atan2(y, x);


                // std::cout << "---" << std::endl;

                Data const function_value = estimator->get_value(dir);

                if (exact)
                {
                    for (int l = 0; l < bands; ++l)
                    {
                        for(int m = -l; m <= l; ++m)
                        {
                            int sh_index = l * (l + 1) + m;

                            float Y_l_m = SH(l, m, theta, phi);
                            float base_coeff = Y_l_m;

                            sh_coefficients[sh_index] += base_coeff * function_value;
                        }
                    }
                }
                else
                {
                    for (int sh_index = 0; sh_index < 9; ++sh_index)
                    {
                        float const Y_l_m = SH_precomputed(sh_index, dir);
                        float const base_coeff = Y_l_m;

                        sh_coefficients[sh_index]  += base_coeff * function_value;
                    }
                }
            }
        }

        float inv = 4.0f * M_PI / (res_u * res_v);

        for (int sh_index = 0; sh_index < bands * bands; ++sh_index)
        {
            sh_coefficients[sh_index] *= inv;
        }
    }

    virtual void add_light(Splat_cube_raster_buffer const& fb_in);

    virtual void get_precalculated_coefficients(Cube_raster_buffer<Spherical_function<float> *> const& normal_map,
                                                Data const& scale,
                                                vector3d_t const& surface_normal);

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
                    test_area_sh += get_value(theta, phi);
                }
                else
                {
                    test_area_sh += get_value(dir);
                }

                test_area_real += area * std::max(dir * normal, 0.0f);
            }
        }

        float inv = 4.0f * M_PI / (res_u * res_v);

        test_area_real *= inv;
        test_area_sh *= inv;

        std::cout << "random, test area real: " << test_area_real << " sh: " << test_area_sh << std::endl;
    }

    virtual Data get_value(Point const& dir) const
    {
        if (!exact)
        {
            Data result = 0.0f;

            for (int sh_index = 0; sh_index < 9; ++sh_index)
            {
                float Y_l_m = SH_precomputed(sh_index, dir);

                result += sh_coefficients[sh_index] * Y_l_m;
            }

            return result;
        }

        float theta = std::acos(dir.z);
        float phi = std::atan2(dir.y, dir.x);
        return get_value(theta, phi);
    }

    virtual void add(Spherical_function<Data> const* s)
    {
        Spherical_harmonics const* sph_harmonics = dynamic_cast<Spherical_harmonics const*>(s);

        *this = *this + *sph_harmonics;
    }

    void normalize(float const factor)
    {
        for (int j = 0; j < bands * bands; ++j)
        {
            sh_coefficients[j] *= factor;
        }
    }


    Spherical_function<Data> * clone() const
    {
        return new Spherical_harmonics(*this);
    }

    friend Spherical_harmonics operator+(Spherical_harmonics const& lhs, Spherical_harmonics const& rhs)
    {
        assert(lhs.bands == rhs.bands && lhs.exact == rhs.exact);

        Spherical_harmonics result(lhs.exact, lhs.bands);

        for (int j = 0; j < lhs.bands * lhs.bands; ++j)
        {
            result.sh_coefficients[j] = lhs.sh_coefficients[j] + rhs.sh_coefficients[j];
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
        // assert(sh_coefficients.empty());

//        for (std::size_t i = 0; i < data.size(); )
//        {
//            Data value = conv_from_vector<Data>(data, i);
//            sh_coefficients.push_back(value);
//        }

        sh_coefficients = conv_from_vector<Data>(data);
        assert(int(sh_coefficients.size()) == bands * bands);
    }



    virtual void accept(Spherical_function_visitor<Data> * visitor)
    {
        visitor->visit(this);
    }


    static std::vector<Sh_coefficient_func> initialize_sh_functions()
    {
        std::vector<Sh_coefficient_func> sh_functions;

        sh_functions.resize(9);

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
        return sh_coefficients[index];
    }

private:
    Data get_value(float const theta, float const phi) const
    {
        if (!exact)
            std::cout << "shouldn't be used!" << std::endl;

        Data result = 0.0f;

        for (int l = 0; l < bands; ++l)
        {
            for(int m = -l; m <= l; ++m)
            {
                int sh_index = l * (l + 1) + m;

                float Y_l_m = SH(l, m, theta, phi);

                result += sh_coefficients[sh_index] * Y_l_m;
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

    float SH_precomputed(int const index, Point const& dir) const
    {
        Sh_coefficient_func const& f = _sh_functions[index];

        //return (this->*(f))(dir);
        return f(dir);
    }

    // http://en.wikipedia.org/wiki/Table_of_spherical_harmonics

    static float SH_precomputed_0(Point const& /* dir */)
    {
        return 0.28209;
    }

    static float SH_precomputed_1(Point const& dir)
    {
        return 0.34549 * dir[0];
    }

    static float SH_precomputed_2(Point const& dir)
    {
        return 0.34549 * dir[1];
    }

    static float SH_precomputed_3(Point const& dir)
    {
        return 0.34549 * dir[2];
    }

    static float SH_precomputed_7(Point const& dir)
    {
        return 0.31539 * (2 * square(dir[2]) - square(dir[0]) - square(dir[1]));
    }

    static float SH_precomputed_5(Point const& dir)
    {
        return 1.0925 * dir[1] * dir[2];
    }

    static float SH_precomputed_4(Point const& dir)
    {
        return 1.0925 * dir[2] * dir[0];
    }

    static float SH_precomputed_6(Point const& dir)
    {
        return 1.0925 * dir[0] * dir[1];
    }

    static float SH_precomputed_8(Point const& dir)
    {
        return 0.54627 * (square(dir[0]) - square(dir[1]));
    }

    int bands;
    bool exact;

    std::vector<Data> sh_coefficients;

    static std::vector<Sh_coefficient_func> _sh_functions;
};


template <class Data>
class Mises_fisher_lobe;

template <class Data>
class Cube_raster_buffer;

template <class Data>
class Mises_Fisher_spherical_function : public Spherical_function<Data>
{
public:

    virtual void calc_coefficients_random(Abstract_spherical_function_estimator<Data> const* estimator)
    {
        std::cout << "Mises_Fisher_spherical_function::calc_coefficients_random(): not used" << std::endl;
    }

    /*
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
    */

    virtual Data get_value(Point const& dir) const
    {
        Data result(0.0f);

        for (unsigned int i = 0; i < _lobes.size(); ++i)
        {
            result += _lobes[i].evaluate(dir);
        }

        // result /= float(_area_lobes.size());

        return result;
    }

    virtual void add(Spherical_function<Data> const* /* s */)
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

    void generate_from_cube_spherical_function(Cube_spherical_function<Data> const* csf, int const num_lobes);

    std::vector<Mises_fisher_lobe<Data> > generate_from_cube(Cube_raster_buffer<Data> const& distribution, int const num_lobes);

    Spherical_function<Data> * clone() const
    {
        return new Mises_Fisher_spherical_function(*this);
    }

    void add_lobe(Mises_fisher_lobe<Data> const& lobe)
    {
        _lobes.push_back(lobe);
    }

    std::vector<Mises_fisher_lobe<Data> > const& get_lobes() const
    {
        return _lobes;
    }

    virtual void accept(Spherical_function_visitor<Data> * visitor)
    {
        visitor->visit(this);
    }

    static Mises_Fisher_spherical_function * test();

private:
    std::vector<Mises_fisher_lobe<Data> > _lobes;
};



__END_YAFRAY
#include <utilities/Mises_fisher.h>
#include <utilities/CubeRasterBuffer.h>
__BEGIN_YAFRAY


template <class Data>
class Cube_spherical_function : public Spherical_function<Data>
{
public:
    Cube_spherical_function(int const resolution, bool use_jitter)
    {
        // init the buffers!
        buffer.setup_surfel_buffer(resolution);

        _use_jitter = use_jitter;

        if (_use_jitter)
        {   vector3d_t axis(2.0f * rand() / float(RAND_MAX) - 1.0f, 2.0f * rand() / float(RAND_MAX) - 1.0f, 2.0f * rand() / float(RAND_MAX) - 1.0f);
            axis.normalize();

            float angle = rand() / float(RAND_MAX) * 2.0f * M_PI / 40.0f;

            Quaternion q_pos(axis, angle);
            _jitter_rotation_positive = q_pos.to_matrix();

            Quaternion q_neg(axis, -angle);
            _jitter_rotation_negative = q_neg.to_matrix();
        }
    }

    // only one sample on each cube cell (center) and the surfel (center)
    void sample_single(Abstract_spherical_function_estimator<Data> const* estimator)
    {
        std::vector<Cube_cell> const& cube_cells = buffer.get_cube_cells();

        for (unsigned int i = 0; i < cube_cells.size(); ++i)
        {
            Cube_cell const& cell = cube_cells[i];

            Point dir = buffer.get_cell_center(cell);
            dir.normalize();

            Data value = estimator->get_value(dir);

            Data current_color = buffer.get_data(cell);
            value += current_color;

            buffer.set_data(cell, value);
        }
    }

    // one sample on the surfel (center) and many on each cube cell
    void sample_cube_cell_area(Abstract_spherical_function_estimator<Data> const* estimator)
    {
        std::vector<Cube_cell> const& cube_cells = buffer.get_cube_cells();

        int const num_samples_per_cell = 8;
        // int const combined_sample_count = num_samples_per_cell * cube_cells.size();

        float const cell_width_2 = 1.0f / buffer.get_resolution();

        int combined_samples_n = 0;

        random_t my_random;

        // FIXME: move into the estimator
        // color_t emission = surfel.material->emission(state, surfel, vector3d_t());

        bool use_cones = false;
        Cone specular_cone;

        if (estimator->get_main_cones().size() == 2)
        {
            use_cones = true;
            specular_cone = estimator->get_main_cones()[1];
        }

        for (unsigned int i = 0; i < cube_cells.size(); ++i)
        {
            _num_all_cells += 1;

            Cube_cell const& cell = cube_cells[i];

            vector3d_t const& cell_dir = buffer.get_cell_direction(cell);
            // if (use_cones && (cell_dir * specular_cone.dir) > specular_cone.cos_angle)
            if (use_cones && fAcos(cell_dir * specular_cone.dir) > specular_cone.cos_angle)
            {
                continue;
            }

            _num_used_cells += 1;

            Data value(0.0f);

            int const axis_0 = Cube_static_data::get_corresponding_axis(0, cell.plane) / 2;
            int const axis_1 = Cube_static_data::get_corresponding_axis(1, cell.plane) / 2;

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
                dir = buffer.get_cell_center(cell) + offset;

                if (_use_jitter)
                {
                    dir = _jitter_rotation_negative * dir;
                }

                dir.normalize();

                value += estimator->get_value(dir);

                ++combined_samples_n;
            }

            value *= 1.0f / float(num_samples_per_cell);

            Data current_value = buffer.get_data(cell);
            value += current_value;

            buffer.set_data(cell, value);
        }

        // assert(combined_samples_n == combined_sample_count);
    }



    // one sample on the surfel (center) and many on each cube cell
    void sample_cube_cell_area(Mises_fisher_lobe<Data> const& lobe)
    {
        std::vector<Cube_cell> const& cube_cells = buffer.get_cube_cells();

        int const num_samples_per_cell = 8;
        int const combined_sample_count = num_samples_per_cell * cube_cells.size();

        float const cell_width_2 = 1.0f / buffer.get_resolution();

        int combined_samples_n = 0;

        random_t my_random;

        for (unsigned int i = 0; i < cube_cells.size(); ++i)
        {
            Cube_cell const& cell = cube_cells[i];

            Data reflected_color(0.0f);

            int const axis_0 = Cube_static_data::get_corresponding_axis(0, cell.plane) / 2;
            int const axis_1 = Cube_static_data::get_corresponding_axis(1, cell.plane) / 2;

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
                dir = buffer.get_cell_center(cell) + offset;
                dir.normalize();

                if (_use_jitter)
                {
                    dir = _jitter_rotation_negative * dir;
                }


                reflected_color += lobe.evaluate(dir);

                ++combined_samples_n;
            }

            reflected_color *= 1.0f / float(num_samples_per_cell);

            Data current_color = buffer.get_data(cell);
            reflected_color += current_color;

            buffer.set_data(cell, reflected_color);

        }

        assert(combined_samples_n == combined_sample_count);
    }



    void calc_coefficients_random(Abstract_spherical_function_estimator<Data> const* estimator)
    // void calc_coefficients_random(renderState_t & state, surfacePoint_t const& surfel, float const area, ray_t const& lightRay, color_t const& light_color)
    {
        // sample_single(state, surfel, area, lightRay, light_color);
        sample_cube_cell_area(estimator);
    }

    Data get_value(Point const& dir) const
    {
        Point dir_tmp = dir;

        if (_use_jitter)
        {
            dir_tmp = _jitter_rotation_positive * dir_tmp;
            // c = get_positive_jittered_cell(c);
        }

        Cube_cell c = buffer.get_cell(dir_tmp);

        Data value  = buffer.get_data_interpolated(c);
        // Color area = area_buffer.get_color(c);

        return value;
    }



    void add(Spherical_function<Data> const* s)
    {
        Cube_spherical_function const* csf = dynamic_cast<Cube_spherical_function const*>(s);

        assert(buffer.get_resolution() == csf->buffer.get_resolution());

        std::vector<Cube_cell> const& cube_cells = buffer.get_cube_cells();

        for (unsigned int i = 0; i < cube_cells.size(); ++i)
        {
            Cube_cell const& c = cube_cells[i];

            Data new_value;

            if (_use_jitter)
            {
                Point dir = buffer.get_cell_direction(c);
                dir = _jitter_rotation_negative * dir; // correct the direction back into the original, common coord system

                new_value = buffer.get_data(c)    + csf->get_value(dir);
            }
            else
            {
                new_value = buffer.get_data(c)    + csf->buffer.get_data(c);
            }

            buffer.set_data(c, new_value);

            if (buffer.get_data(c) != new_value && !is_in_range(0.9f, 1.1f, buffer.get_data(c), new_value))
            {
                std::cout << "add failed: " << buffer.get_data(c) << " " << new_value << std::endl;
                assert(false);
            }
        }
    }

    void normalize(float const factor)
    {
        std::vector<Cube_cell> const& cube_cells = buffer.get_cube_cells();

        for (unsigned int i = 0; i < cube_cells.size(); ++i)
        {
            Cube_cell const& c = cube_cells[i];

            Data new_value = buffer.get_data(c) * factor;

            buffer.set_data(c, new_value);
        }
    }

    Spherical_function<Data> * clone() const
    {
        assert(false);
        return NULL;
    }

    Cube_raster_buffer<Data> const& get_buffer() const
    {
        return buffer;
    }

    Cube_raster_buffer<Data> & get_buffer()
    {
        return buffer;
    }

    virtual void accept(Spherical_function_visitor<Data> * visitor)
    {
        visitor->visit(this);
    }

    virtual std::vector< std::vector<float> > components_to_vectors() const
    {
        std::vector< std::vector<float> > result(6);

        for (unsigned int i = 0; i < 6; ++i)
        {
            std::vector<float> component = buffer.component_to_vector(i, false);
            std::copy(component.begin(), component.end(), std::back_inserter(result[i]));
        }

        return result;
    }

    virtual void from_component_vectors(std::vector< std::vector<float> > const& data)
    {
        assert(data.size() == 6);

        int const entries = (buffer.get_resolution() * buffer.get_resolution() * sizeof(Data) / sizeof(float));
        assert(entries == int(data[0].size()));

        for (unsigned int i = 0; i < 6; ++i)
        {
            if (data[i].size() == 0) continue;

            std::vector<float> value_vector;

            std::copy(data[i].begin(), data[i].end(), std::back_inserter(value_vector));

            buffer.from_component_vector(i, value_vector, false);
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


    static Cube_spherical_function * create_from_mises_fisher(Mises_fisher_lobe<Data> const& lobe, int const resolution, bool const use_jitter, float const weight)
    {
        Cube_spherical_function * sf = new Cube_spherical_function<Data>(resolution, use_jitter);
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
        return buffer.get_resolution();
    }

    static int _num_used_cells;
    static int _num_all_cells;

private:
    Cube_raster_buffer<Data> buffer;

    bool _use_jitter;
    float _jitter_x, _jitter_y;
    Matrix_3x3 _jitter_rotation_positive;
    Matrix_3x3 _jitter_rotation_negative;
};

template<class Data> int Cube_spherical_function<Data>::_num_used_cells = 0;
template<class Data> int Cube_spherical_function<Data>::_num_all_cells = 0;


// how to replace the data in the final tree?


template <class Data>
class Indexed_spherical_function : public Spherical_function<Data>
{
public:
    Indexed_spherical_function(std::vector< std::vector<float> > const& unpacked_dictionary,
                               std::vector<Spherical_function<Data> *> const* dictionary,
                               Spherical_function<Data> const* sf,
                               bool const keep_original = true) :
        _dictionary(dictionary),
        _original_sf(NULL)
    {
        Distance_function distance_fnc;

        float lowest_distance = 1e10f;
        int closest_dictionary_index = -1;

        std::vector<float> sf_vector = sf->to_vector();

        for (std::size_t i = 0; i < unpacked_dictionary.size(); ++i)
        {
            assert(sf_vector.size() == unpacked_dictionary[i].size());

            float const distance = distance_fnc(unpacked_dictionary[i], sf_vector);

            if (distance < lowest_distance)
            {
                lowest_distance = distance;
                closest_dictionary_index = i;
            }
        }

        assert(closest_dictionary_index >= 0);

        _index = closest_dictionary_index;

        if (keep_original)
        {
            _original_sf = sf->clone();
        }
    }

    Indexed_spherical_function(int const index,
                               std::vector<Spherical_function<Data> *> const* dictionary,
                               Spherical_function<Data> * sf = NULL) :
        _dictionary(dictionary),
        _original_sf(sf),
        _index(index)
    { }

    virtual ~Indexed_spherical_function() {}

    void calc_coefficients_random(Abstract_spherical_function_estimator<Data> const* estimator)
    {
        std::cout << "Indexed_spherical_function::calc_coefficients_random() shouldn't be used" << std::endl;
    }

    Data get_value(Point const& dir) const
    {
        return (*_dictionary)[_index]->get_value(dir);
    }

    virtual void add(Spherical_function<Data> const* /* s */) { std::cout << "Indexed_spherical_function::add() shouldn't be used" << std::endl; }
    virtual void normalize(float const /* factor */) { std::cout << "Indexed_spherical_function::normalize() shouldn't be used" << std::endl; }

    Spherical_function<Data> * get_indexed_sf()
    {
        return (*_dictionary)[_index];
    }

    // debugging
    Spherical_function<Data> * get_original_sf()
    {
        return _original_sf;
    }

    int get_index() const
    {
        return _index;
    }


    virtual void accept(Spherical_function_visitor<Data> * visitor)
    {
        visitor->visit(this);
    }

    virtual Spherical_function<Data> * clone() const { assert(false); return NULL; }

private:
    std::vector< Spherical_function<Data> *> const* _dictionary;
    Spherical_function<Data> * _original_sf; // debugging, comparison
    int _index;
};



#if 0
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

#endif




template <class Data>
class Spherical_function_factory
{
public:
    virtual Spherical_function<Data> * create() const = 0;
    virtual std::string get_name() const = 0;
    virtual std::string get_settings() const = 0;
};

template <class Data>
class Spherical_harmonics_factory : public Spherical_function_factory<Data>
{
public:
    Spherical_harmonics_factory(int const bands, bool const exact) :
        _bands(bands),
        _exact(exact)
    { }

    Spherical_function<Data> * create() const
    {
        return new Spherical_harmonics<Data>(_exact, _bands);
    }

    std::string get_name() const
    {
        return "SH";
    }

    std::string get_settings() const
    {
        std::stringstream ss;
        ss << "bands: " << _bands  << ", exact: " << _exact;
        return ss.str();
    }

private:
    int _bands;
    bool _exact;
};

template <class Data>
class Cube_spherical_function_factory : public Spherical_function_factory<Data>
{
public:
    Cube_spherical_function_factory(int const resolution, bool use_jitter) :
        _resolution(resolution),
        _use_jitter(use_jitter)
    { }

    Spherical_function<Data> * create() const
    {
        return new Cube_spherical_function<Data>(_resolution, _use_jitter);
    }

    std::string get_name() const
    {
        return "Cube";
    }

    std::string get_settings() const
    {
        std::stringstream ss;
        ss << "resolution: " << _resolution << ", jitter: " << _use_jitter;
        return ss.str();
    }

private:
    int _resolution;
    bool _use_jitter;
};


template <class Data>
void Mises_Fisher_spherical_function<Data>::generate_from_cube_spherical_function(Cube_spherical_function<Data> const* csf, int const num_lobes)
{
    // if (!csf->get_buffer().is_black())
    if (csf->get_buffer().calc_total_energy() > 0.001f)
    {
        _lobes = generate_from_cube(csf->get_buffer(), num_lobes);
    }
}



template <class Data>
std::vector<Mises_fisher_lobe<Data> > Mises_Fisher_spherical_function<Data>::generate_from_cube(Cube_raster_buffer<Data> const& distribution, int const num_lobes)
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
            // 3.) Divide by total likelihood
            Cube_cell const& c = cells[j];
            const Data pixel_value = distribution.get_data(c); // + 0.01f;
            const float solid_angle = distribution.get_solid_angle(c); // / (4.0f * M_PI);
            for (int k = 0; k < num_lobes; k++)
            {
                if (totalLikelihoods.get_data(c) > 0.0001f)
                {
                    likelihoods[k].set_data(c, pixel_value * solid_angle * likelihoods[k].get_data(c) / totalLikelihoods.get_data(c));
                }
            }
        }

        // M-step
        for (int lobe_index = 0; lobe_index < num_lobes; lobe_index++)
        {
            lobes[lobe_index].fit(likelihoods[lobe_index]);

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

    for (int lobe_index = 0; lobe_index < num_lobes; lobe_index++)
    {
        lobes[lobe_index].weight = distribution.integrate(likelihoods[lobe_index]);

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


template <class Data>
Mises_Fisher_spherical_function<Data> * Mises_Fisher_spherical_function<Data>::test()
{
    vector3d_t mean_dir(0.0f, 0.5f, 0.5f);
    mean_dir.normalize();
    float const concentration = 30;
    float const weight = 1.0f;

    Mises_fisher_lobe<float> lobe_0(mean_dir, concentration, weight);
    Mises_fisher_lobe<float> lobe_1(vector3d_t(1.0f, 0.0f, 0.0f), 10.0f, weight);

    std::cout << "test lobe 0: " << lobe_0 << std::endl;
    std::cout << "test lobe 1: " << lobe_1 << std::endl;

    Cube_spherical_function<float> * distribution   = Cube_spherical_function<float>::create_from_mises_fisher(lobe_0, 8, false, 2.0f);
    Cube_spherical_function<float> * distribution_1 = Cube_spherical_function<float>::create_from_mises_fisher(lobe_1, 8, false, 1.0f);

    // distribution->add(distribution_1);


    Mises_Fisher_spherical_function * mf_sf = new Mises_Fisher_spherical_function;
    mf_sf->generate_from_cube_spherical_function(distribution, 10);

    std::cout << "lobes:" << std::endl;
    for (unsigned int i = 0; i < mf_sf->get_lobes().size(); ++i)
    {
        std::cout << i << " " << mf_sf->get_lobes()[i] << std::endl;
    }

    return mf_sf;
}


template <class Data>
void Spherical_harmonics<Data>::add_light(Splat_cube_raster_buffer const& fb_in)
{
    std::vector<Cube_cell> const& cells = fb_in.get_cube_cells();

    std::vector<Data> new_sh_coefficients = std::vector<Data>(bands * bands, Data(0.0f));

    for (std::size_t i = 0; i < cells.size(); ++i)
    {
        vector3d_t const& cell_dir = fb_in.get_cell_direction(cells[i]);

        float const phi = std::atan2(cell_dir.y, cell_dir.x);
        float const theta = M_PI / 2.0f - std::atan(cell_dir.z / std::sqrt(cell_dir.x * cell_dir.x + cell_dir.y * cell_dir.y));

        Data const function_value = convert_data<Data>(fb_in.get_data(cells[i]));

        float const solid_angle = fb_in.get_solid_angle(cells[i]);

        if (exact)
        {
            for (int l = 0; l < bands; ++l)
            {
                for(int m = -l; m <= l; ++m)
                {
                    int sh_index = l * (l + 1) + m;

                    float Y_l_m = SH(l, m, theta, phi);
                    float base_coeff = Y_l_m;

                    new_sh_coefficients[sh_index] += base_coeff * function_value * solid_angle;
                }
            }
        }
        else
        {
            for (int sh_index = 0; sh_index < 9; ++sh_index)
            {
                float const Y_l_m = SH_precomputed(sh_index, cell_dir);
                float const base_coeff = Y_l_m;

                new_sh_coefficients[sh_index] += base_coeff * function_value * solid_angle;
            }
        }
    }


    // float inv = 4.0f * M_PI / (6 * fb_in.get_resolution() * fb_in.get_resolution());

    for (int sh_index = 0; sh_index < bands * bands; ++sh_index)
    {
        sh_coefficients[sh_index] += new_sh_coefficients[sh_index]; //  * inv;
    }
}

template <class Data>
void Spherical_harmonics<Data>::get_precalculated_coefficients(Cube_raster_buffer<Spherical_function<float> *> const& normal_map,
                                            Data const& scale,
                                            vector3d_t const& surface_normal)
{
    Cube_cell normal_cell = normal_map.get_cell(surface_normal);

    Spherical_function<float> * sf = normal_map.get_data(normal_cell);

    Spherical_harmonics<float> const* sph_harmonics = dynamic_cast<Spherical_harmonics<float> const*>(sf);

    for (size_t i = 0; i < sh_coefficients.size(); ++i)
    {
        sh_coefficients[i] = scale * sph_harmonics->get_coefficient(i);
    }
}



//Spherical_node_representation::~Spherical_node_representation()
//{
//    delete area;
//    delete color;
//}


__END_YAFRAY

#endif // SPHERICAL_HARMONICS_H
