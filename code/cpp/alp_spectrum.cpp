#include <cmath>
#include <iostream>

#include "alp_spectrum.hpp"

// Cross section for massless ALPs (in eV^-2)
double sigma0(double g, double ea, double ks2)
{
    const double prefactor = 0.125 * alpha;
    double x = 0.25*ks2/(ea*ea);
    double x1plog1px = (1.0 + x)*log1p(1/x) - 1.0;
    return prefactor * g*g * x1plog1px;
}

// Mass-dependent cross section (in eV^-2)
double sigma(double m, double g, double ea, double ks2)
{
    const double prefactor = 0.125 * alpha;
    double e2 = ea*ea;
    double m2 = m*m;
    double m4 = m2 * m2;
    // double eab = sqrt(e2 - m2);
    double b = beta(m, ea);
    double bp1 = b + 1.0;
    double yp = 2.0*e2*bp1 - m2;
    // Numerically stable way to compute 1 - sqrt(1 - x) = x / (1 + sqrt(1 - x))
    // With x = m2/e2, b = sqrt(1 - x), this leads to 2*e2*(1 - b) - m2 = 2*m2/(sqrt(-x + 1) + 1) - m2 = 2*m2/(b + 1) - m2
    double ym = 2.0*m2/bp1 - m2;
    if (ym < 0) { std::cerr << "Numerical instability in Primakoff cross section detected: ym = " << ym << std::endl; }

    double res = (1.0 + 0.25 * ks2 / e2 - 0.5 * m2 / e2) * log((yp + ks2) / (ym + ks2)) - b;

    double num = m4 + ks2 * yp;
    double denom = m4 + ks2 * ym;
    if (denom <= 0)
    {
        std::cerr << "Numerical instability in Primakoff cross section detected: denom = " << denom << std::endl;
    }
    else
    {
        res -= (0.25 * m4 / (ks2 * e2)) * log(num / denom);
    }
    return prefactor * g*g * res;
}

// Energy spectrum of massive or massless ALPs (in eV^-1) [Eq. (7) in arXiv:1702.02964]
double spectrum(double m, double g, double ea, double ks2, double snorm, double t_eff)
{   
    return snorm * ea*ea * sigma(m, g, ea, ks2) / (exp(ea/t_eff) - 1.0);
}

double spectrum0(double g, double ea, double ks2, double snorm, double t_eff)
{   
    return snorm * ea*ea * sigma0(g, ea, ks2) / (exp(ea/t_eff) - 1.0);
}

double temporal_spectrum(double m, double g, double ea, double tem, double norm, double e0, double p, double ks2)
{ 
   double er = 1e-6*ea/e0;
   return 1.0e46 * norm * pow(er, p) * exp(-(p+1.0)*er) * sigma(m, g, ea, ks2) / sigma0(1e-19, ea, ks2);
}