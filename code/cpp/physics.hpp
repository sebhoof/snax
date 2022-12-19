#ifndef __constants_hpp__
#define __constants_hpp__

#include <numbers>

const double pi = std::numbers::pi;
const double alpha = 1.0/137.035999084; // Fine structure constant alpha
const double hbar = 6.582119569e-16; // Planck's constant hbar (eVs)
const double c0 = 299792458.0; // Speed of light c (in m/s)
const double hbarc = hbar*c0; // bhar*c (in eVm)
const double pc = 3.08567758149137e16; // Parsec (in m)

double beta(double m, double ea);
double gamma(double m, double ea);
double tau(double m, double g);
double lam(double m, double g, double ea);

#endif // defined __constants_hpp__