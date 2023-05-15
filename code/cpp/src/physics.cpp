#include <cmath>

#include "physics.hpp"

// Beta and gamma factors in special relativity

double beta(double m, double ea) { return sqrt(1.0 - m*m/(ea*ea)); }

double gamma(double m, double ea) { return ea/m; }

// ALP frame lifetime (in eV^-1) and decay length (in lab system, in m)
double tau(double m, double g)
{
    double x = m*g;
    return 64.0*pi/(m*x*x);
}

double lam(double m, double g, double ea)
{
    double b = beta(m, ea);
    double gam = gamma(m, ea);
    double t = tau(m, g);
    return hbarc*b*gam*t;
};