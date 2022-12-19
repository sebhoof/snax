#ifndef __alp_decays_hpp__
#define __alp_decays_hpp__

#include "physics.hpp"
#include "alp_spectrum.hpp"

double calc_ea_min(double m, double ep);
double calc_ea_max_weak(double m, double ep, double tp0, double tem=0);

double ta_from_t(double t, double ea, double ep, double m);
double inner_integrand_approx(double ea, double ep, double tp1, double m, double g, double ks2=ks2_jaeckel17, double snorm=snorm_jaeckel17, double t_eff=t_eff_jaeckel17);
double inner_integrand(double ea, double ep, double tp0, double tp1, double m, double g, double ks2=ks2_jaeckel17, double snorm=snorm_jaeckel17, double t_eff=t_eff_jaeckel17);
double inner_integrand_with_tem(double ea, double tem, double ep, double tp0, double tp1, double m, double g, double norm, double e0, double p, double ks2=ks2_jaeckel17);
std::vector<double> acceptance_fraction(std::vector<double> eas, std::pair<double,double> ebin, std::pair<double,double> tbin, double m, double g);
double acceptance_fraction_vegas(double ea, double u_ldec, double u_x0, std::pair<double,double> ebin, std::pair<double,double> tbin, double m, double g);

#endif // defined __alp_decays_hpp__