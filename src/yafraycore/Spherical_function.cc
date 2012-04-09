#include <utilities/spherical_harmonics.h>

__BEGIN_YAFRAY

//template <class Data>
//std::vector<typename GiSphericalHarmonics<Data>::Sh_coefficient_func> GiSphericalHarmonics<Data>::_sh_functions = GiSphericalHarmonics<Data>::initialize_sh_functions();

template <>
std::vector<GiSphericalHarmonics<color_t>::Sh_coefficient_func> GiSphericalHarmonics<color_t>::_sh_functions = GiSphericalHarmonics<color_t>::initialize_sh_functions();

template <>
std::vector<GiSphericalHarmonics<float>::Sh_coefficient_func> GiSphericalHarmonics<float>::_sh_functions = GiSphericalHarmonics<float>::initialize_sh_functions();


__END_YAFRAY
