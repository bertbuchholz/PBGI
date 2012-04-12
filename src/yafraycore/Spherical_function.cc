#include <utilities/spherical_harmonics.h>

__BEGIN_YAFRAY

//template <class Data>
//std::vector<typename GiSphericalHarmonics<Data>::Sh_coefficient_func> GiSphericalHarmonics<Data>::_sh_functions = GiSphericalHarmonics<Data>::initialize_sh_functions();

template <>
std::vector<GiSphericalHarmonics<color_t>::Sh_coefficient_func> GiSphericalHarmonics<color_t>::_sh_functions = GiSphericalHarmonics<color_t>::initialize_sh_functions();

template <>
std::vector<GiSphericalHarmonics<float>::Sh_coefficient_func> GiSphericalHarmonics<float>::_sh_functions = GiSphericalHarmonics<float>::initialize_sh_functions();


template <>
void GiSphericalHarmonics<color_t>::get_precalculated_coefficients(Cube_raster_buffer<Spherical_function<float> *> const& normal_map,
                                                                   vector3d_t const& surface_normal,
                                                                   color_t const& surface_color,
                                                                   vector3d_t const& light_dir,
                                                                   color_t const& light_color)
{
    Cube_cell normal_cell = normal_map.get_cell(surface_normal);

    Spherical_function<float> * sf = normal_map.get_data(normal_cell);

    GiSphericalHarmonics<float> const* sph_harmonics = dynamic_cast<GiSphericalHarmonics<float> const*>(sf);

    float const lambert_scale = light_dir * surface_normal;

    for (size_t i = 0; i < sh_coefficients.size(); ++i)
    {
        sh_coefficients[i] = surface_color * light_color * lambert_scale * sph_harmonics->get_coefficient(i);
    }
}


// TODO: for now, assume float means area, but this is nonesense, need to give something like the estimator stuff
template <>
void GiSphericalHarmonics<float>::get_precalculated_coefficients(Cube_raster_buffer<Spherical_function<float> *> const& normal_map,
                                                                 vector3d_t const& surface_normal,
                                                                 color_t const& surface_color,
                                                                 vector3d_t const& light_dir,
                                                                 color_t const& light_color)
{
    Cube_cell normal_cell = normal_map.get_cell(surface_normal);

    Spherical_function<float> * sf = normal_map.get_data(normal_cell);

    GiSphericalHarmonics<float> const* sph_harmonics = dynamic_cast<GiSphericalHarmonics<float> const*>(sf);

    for (size_t i = 0; i < sh_coefficients.size(); ++i)
    {
        sh_coefficients[i] = sph_harmonics->get_coefficient(i); // TODO: add area
    }
}


__END_YAFRAY
