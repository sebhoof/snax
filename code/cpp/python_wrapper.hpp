#ifndef __python_wrapper_hpp__
#define __python_wrapper_hpp__

#include <iostream>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "physics.hpp"
#include "alp_spectrum.hpp"
#include "alp_decays.hpp"

// The name of the module and other info
void module_info();

#endif // defined __python_wrapper_hpp__
