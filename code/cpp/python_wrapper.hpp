#ifndef __python_wrapper_hpp__
#define __python_wrapper_hpp__

#include <iostream>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "physics.hpp"
#include "alp_spectrum.hpp"
#include "alp_decays.hpp"

using namespace pybind11::literals;

// The name of the module and other info
void module_info();
pybind11::dict constants_dict(
                              "alpha"_a=alpha, "hbar"_a=hbar, "c"_a=c0, "hbarc"_a=hbarc, "pc"_a=pc,
                              "d"_a=d, "distance_factor"_a=distance_factor, "renv"_a=r_env,
                              "ks2_payez15"_a=ks2_payez15, "t_eff_payez15"_a=t_eff_payez15, "snorm_payez15"_a=snorm_payez15,
                              "ks2_jaeckel17"_a=ks2_jaeckel17, "t_eff_jaeckel17"_a=t_eff_jaeckel17, "snorm_jaeckel17"_a=snorm_jaeckel17
                             );
#endif // defined __python_wrapper_hpp__
