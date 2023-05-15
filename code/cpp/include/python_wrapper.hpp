#ifndef __python_wrapper_hpp__
#define __python_wrapper_hpp__

#include <iostream>

#include <pybind11/pybind11.h>
#include <pybind11/embed.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "physics.hpp"
#include "alp_spectrum.hpp"
#include "alp_decays.hpp"

using namespace pybind11::literals;

// The name of the module and other info
void module_info();

const std::unordered_map<std::string, double> constants_umap {
   { "alpha", alpha }, { "hbar", hbar }, { "c", c0 }, { "hbarc", hbarc }, { "pc", pc },
   { "d", d }, { "distance_factor", distance_factor }, { "renv", r_env },
   { "ks2_payez15", ks2_payez15 }, { "t_eff_payez15", t_eff_payez15 }, { "snorm_payez15", snorm_payez15 },
   { "ks2_jaeckel17", ks2_jaeckel17 }, { "t_eff_jaeckel17", t_eff_jaeckel17 }, { "snorm_jaeckel17", snorm_jaeckel17 }
};

double get_umap_constant(std::string name);

#endif // defined __python_wrapper_hpp__
