#include <iostream>
#include <cmath>
#include <complex>

#include "alp_spectrum.hpp"

// ALP-photon cross section (in eV^-2)
double sigma(double m, double g, double ea, double ks2)
{
    const double prefactor = 0.125 * alpha;
    if (m > ea) { return 0; }
    double res = prefactor * g*g;
    if (m > 0)
    {
        double e2 = ea*ea;
        double m2 = m*m;
        double m4 = m2 * m2;
        double b = beta(m, ea);
        double bp1 = b + 1.0;
        double yp = 2.0*e2*bp1 - m2;
        // Numerically stable way to compute 1 - sqrt(1 - x) = x / (1 + sqrt(1 - x))
        // With x = m2/e2, b = sqrt(1 - x), this leads to 2*e2*(1 - b) - m2 = 2*m2/(sqrt(-x + 1) + 1) - m2 = 2*m2/(b + 1) - m2
        double ym = 2.0*m2/bp1 - m2;
        if (ym < 0) { std::cerr << "Numerical instability in Primakoff cross section detected: ym = " << ym << std::endl; }
        double temp = (1.0 + 0.25 * ks2 / e2 - 0.5 * m2 / e2) * log((yp + ks2) / (ym + ks2)) - b;
        double num = m4 + ks2 * yp;
        double denom = m4 + ks2 * ym;
        if (denom <= 0)
        {
            std::cerr << "Numerical instability in Primakoff cross section detected: denom = " << denom << std::endl;
        }
        else
        {
        temp -= (0.25 * m4 / (ks2 * e2)) * log(num / denom);
        }    
        res *= temp;
    }
    else
    {
        double x = 0.25*ks2/(ea*ea);
        res *= (1.0 + x)*log1p(1/x) - 1.0;
    }
    return  res;
}

// Energy spectrum of massive or massless ALPs (in eV^-1) [Eq. (7) in arXiv:1702.02964]
double spectrum(double m, double g, double ea, double ks2, double snorm, double t_eff)
{   
    return snorm * ea*ea * sigma(m, g, ea, ks2) / (exp(ea/t_eff) - 1.0); // TODO: replace with expm1()
}

double spectrum0(double g, double ea, double ks2, double snorm, double t_eff)
{   
    return snorm * ea*ea * sigma(0, g, ea, ks2) / (exp(ea/t_eff) - 1.0); // TODO: replace with expm1()
}

double temporal_spectrum0(double g, double ea, double tem, double norm, double e0, double p, double ks2)
{ 
   double er = 1e-6*ea/e0;
   return 1.0e46 * norm * pow(er, p) * exp(-(p+1.0)*er);
}

double temporal_spectrum(double m, double g, double ea, double tem, double norm, double e0, double p, double ks2)
{ 
   double er = 1e-6*ea/e0;
   return 1.0e46 * norm * pow(er, p) * exp(-(p+1.0)*er) * sigma(m, g, ea, ks2) / sigma(0, g, ea, ks2);
}