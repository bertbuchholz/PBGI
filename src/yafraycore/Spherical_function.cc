#include <utilities/spherical_harmonics.h>

__BEGIN_YAFRAY

template <>
std::vector<Spherical_harmonics<color_t>::Sh_coefficient_func> Spherical_harmonics<color_t>::_sh_functions = Spherical_harmonics<color_t>::initialize_sh_functions();

template <>
std::vector<Spherical_harmonics<float>::Sh_coefficient_func> Spherical_harmonics<float>::_sh_functions = Spherical_harmonics<float>::initialize_sh_functions();


__END_YAFRAY
