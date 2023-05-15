#include "python_wrapper.hpp"

using namespace pybind11::literals;

double py_gamma(double m, double ea) { return gamma(m, ea); };

PYBIND11_MODULE(pysnax, m)
{
  m.doc() = "Python wrapper for the SNax library, calculating axion and axion-like particle gamma-ray signatures from supernova SN1987A.";

  try { pybind11::module_::import("numpy"); } catch (...) { return; }

  m.def("module_info", &module_info, "Basic information about the module.");
  m.attr("constants") = constants_dict;
  m.def("get_constant", &get_umap_constant, "Return a constant used inside the module.", "name"_a);
  m.def("beta", &beta, "Beta function in special relativity.", "m"_a, "ea"_a);
  m.def("gamma", &py_gamma, "Gamma function in special relativity.", "m"_a, "ea"_a);
  m.def("tau", &tau, "ALP rest-frame lifetime (in eV^-1).", "m"_a, "g"_a);
  m.def("lam", &lam, "ALP lab-frame mean decay length (in m).", "m"_a, "g"_a, "ea"_a);
  m.def("sigma", &sigma, "Primakoff cross section for massive ALPs - including the exact limit of the massless case (in eV^-2).", "m"_a, "g"_a, "ea"_a, "ks2"_a=ks2_jaeckel17);
  m.def("spectrum0", &spectrum0, "Instantaneous emission spectrum for massless ALPs (in eV^-1).", "g"_a, "ea"_a, "ks2"_a=ks2_jaeckel17, "snorm"_a=snorm_jaeckel17, "t_eff"_a=t_eff_jaeckel17);
  m.def("spectrum", &spectrum, "Instantaneous emission spectrum for massive ALPs (in eV^-1).", "m"_a, "g"_a, "ea"_a, "ks2"_a=ks2_jaeckel17, "snorm"_a=snorm_jaeckel17, "t_eff"_a=t_eff_jaeckel17);
  m.def("temporal_spectrum", &temporal_spectrum, "Temporal emission spectrum for massive ALPs (in eV^-1 s^-1).", "m"_a, "g"_a, "ea"_a, "ta"_a, "norm"_a, "e0"_a, "p"_a, "ks2"_a=ks2_jaeckel17);
  m.def("calc_ea_min", &calc_ea_min, "Lower integration limit on E_a (in eV).", "m"_a, "ep"_a);
  m.def("calc_ea_max_weak", &calc_ea_max_weak, "Calculate a weak upper integration limit on E_a (in eV).", "m"_a, "ep"_a, "tp0"_a, "tem"_a = 0);
  m.def("inner_integrand_approx", &inner_integrand_approx, "Approximate integrand from the literature.", "ea"_a, "ep"_a, "tp1"_a, "m"_a, "g"_a, "ks2"_a=ks2_jaeckel17, "snorm"_a=snorm_jaeckel17, "t_eff"_a=t_eff_jaeckel17);
  m.def("ta_geometry", &ta_geometry, "Calculates the timescale associated with the geometrical condition (non-zero determinant) on the decay time.", "ea"_a, "ep"_a, "m"_a);
  m.def("ta_from_t", &ta_from_t, "t"_a, "ea"_a, "ep"_a, "m"_a);
  m.def("inner_integrand", &inner_integrand, "ea"_a, "ep"_a, "tp0"_a, "tp1"_a, "m"_a, "g"_a, "ks2"_a=ks2_jaeckel17, "snorm"_a=snorm_jaeckel17, "t_eff"_a=t_eff_jaeckel17);
  m.def("inner_integrand_with_tem", &inner_integrand_with_tem, "ea"_a, "ta"_a, "ep"_a, "tp0"_a, "tp1"_a, "m"_a, "g"_a, "norm"_a, "e0"_a, "p"_a, "ks2"_a=ks2_jaeckel17);
  m.def("acceptance_fraction", &acceptance_fraction, "eas"_a, "ebin"_a, "tbin"_a, "m"_a, "g"_a);
  m.def("acceptance_fraction_vegas", &acceptance_fraction_vegas, "ea"_a, "u_ldec"_a, "u_x0"_a, "ebin"_a, "tbin"_a, "m"_a, "g"_a);
}

void module_info()
{
  std::cout << "This is the SNAx library v2.0, calculating the gamma-ray signal predictions from ALPs after SN1987A since 2022!" << std::endl;
}

double get_umap_constant(std::string name)
{
  return constants_umap.at(name);
}