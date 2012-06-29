#include <integrators/Pbgi_integrator/spherical_harmonics.h>

__BEGIN_YAFRAY

// NOTE: this specialization is for 3 bands (9 coeffs) only, not sure if functions can be partially specialzed

template <>
boost::array<Spherical_harmonics<color_t>::Sh_coefficient_func, 9> Spherical_harmonics<color_t>::_sh_functions = Spherical_harmonics<color_t>::initialize_sh_functions();

template <>
boost::array<Spherical_harmonics<float>::Sh_coefficient_func, 9> Spherical_harmonics<float>::_sh_functions = Spherical_harmonics<float>::initialize_sh_functions();

Spherical_harmonics<float> get_precalculated_sh(Cube_raster_buffer< Spherical_harmonics<float> > const& normal_map,
                                                          vector3d_t const& surface_normal)
{
    Cube_cell normal_cell = normal_map.get_cell(surface_normal);

    return normal_map.get_data(normal_cell);
}

__END_YAFRAY
