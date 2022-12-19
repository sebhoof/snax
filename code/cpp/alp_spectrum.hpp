#ifndef __alp_spectrum_hpp__
#define __alp_spectrum_hpp__

#include "physics.hpp"

// Supernova SN1987A constants

const double d = 51.4e3 * pc;         // Distance from Earth to SN1987A (in m)
const double d_in_s = d/c0;           // " (in s)
const double distance_factor = 1.0e-4/(4.0*pi*d*d); // (in cm^-2)
const double r_env = 3.0e10;          // SN1987a envelope below which photons are not released (in m)
const double r_env_in_s = r_env/c0;   // " (in s)

// Effective SN1987A properties, derived from arXiv:1410.3747
const double ks_payez15 = 1.73441270e7; // Effective SN1987A Debye screening scale (in eV)
const double ks2_payez15 = ks_payez15*ks_payez15;
const double t_eff_payez15 = 3.12857459e7; // Effective SN1987A temperature (in eV)
const double snorm_payez15 = 2.03595118e71; // Normalisation constant for the spectrum (in eV^-1)

// Alternative effective SN1987A properties given by arXiv:1702.02964
const double ks_jaeckel17 = 16.8e6;   // Effective SN1987A Debye screening scale (in eV)
const double ks2_jaeckel17 = ks_jaeckel17*ks_jaeckel17;
const double t_eff_jaeckel17 = 30.6e6; // Effective SN1987A temperature (in eV)
const double snorm_jaeckel17 = 2.54e71; // Normalisation constant for the spectrum (in eV^-1)

double sigma0(double g, double ea, double ks2=ks2_jaeckel17);
double sigma(double m, double g, double ea, double ks2=ks2_jaeckel17);

double spectrum0(double g, double ea, double ks2=ks2_jaeckel17, double snorm=snorm_jaeckel17, double t_eff=t_eff_jaeckel17);
double spectrum(double m, double g, double ea, double ks2=ks2_jaeckel17, double snorm=snorm_jaeckel17, double t_eff=t_eff_jaeckel17);
double temporal_spectrum(double m, double g, double ea, double tem, double norm, double e0, double p, double ks2=ks2_jaeckel17);

#endif // defined __alp_spectrum_hpp__