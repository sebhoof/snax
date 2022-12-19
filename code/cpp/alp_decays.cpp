#include <cmath>
#include <iostream>
#include <random>

#include "alp_decays.hpp"

double calc_ea_min(double m, double ep) { return ep + 0.25*m*m/ep; }

double calc_ea_max_weak(double m, double ep, double tp0, double tem)
{
   double x = (2.0*d_in_s + tp0 - tem)*(tp0 - tem);
   x = d_in_s*d_in_s/x;
   return ep + 0.25*m*m*(1.0 + x)/ep;
}

double inner_integrand_approx(double ea, double ep, double tp1, double m, double g, double ks2, double snorm, double t_eff)
{
    double mtau_in_s = m*hbar*tau(m, g);
    double tcrit = 0.5*mtau_in_s/ep;
    double b = beta(m, ea);
    double t_env = 0.5*m*m*r_env_in_s/(ep*b*ea);
    double exp_sup = exp(-t_env/tcrit) - exp(-tp1/tcrit);
    // Return result in 1/cm^2 x 1 x 1/eV^2
    return 2.0 * distance_factor * exp_sup * spectrum(m, g, ea, ks2, snorm, t_eff) / ea;
}

double omsqrt1m(double x) { return x / (1.0 + sqrt(1.0 - x)); };

double ta_from_t(double t, double ea, double ep, double m)
{
    double delta = t + d_in_s;
    double denom = 2.0 * (1.0 - ep/ea);
    double omdiscr = 1.0 - d_in_s*d_in_s/(delta*delta);
    omdiscr *= 4.0*ep*(ea - ep)/(m*m);
    // Code relies on discriminat = 1 - omdiscr >= 0; the condition should be implemented in the integral limits
    // Make more precise computation of 1 - sqrt(1 - x) = x/(1 + sqrt(1 - x))
    return delta*omsqrt1m(omdiscr)/denom;
};

double t_from_ta(double ta, double ea, double ep, double m)
{
    double r = 0.5*m*m*ta/(ea*ep*d_in_s);
    double discr = 1.0 - (4.0*ep*(ea - ep)/m*m - 1.0)*r*r;
    return d_in_s*(r  - 1.0 + sqrt(discr));
}

double t_geometry(double ea, double ep, double m)
{
    double diff = ea - ep;
    double denom = diff - 0.25*m*m/ep;
    return d_in_s * (sqrt(diff/denom) - 1.0);
};

double ta_geometry(double ea, double ep, double m)
{
    double diff = ea - ep;
    double denom = diff*(diff - 0.25*m*m/ep);
    return 0.5*d_in_s * ea/sqrt(denom);
};

double inner_integrand(double ea, double ep, double tp0, double tp1, double m, double g, double ks2, double snorm, double t_eff)
{
    double b = beta(m, ea);
    double ta_env = r_env_in_s/b;
    double ta_tp0 = ta_from_t(tp0, ea, ep, m);
    double ta_min = fmax(ta_env, ta_tp0);
    double ta_geo = ta_geometry(ea, ep, m);
    double ta_tp1 = ta_from_t(tp1, ea, ep, m);
    double ta_max = fmin(ta_geo, ta_tp1);
    if (ta_max > ta_min)
    {
        double gam = gamma(m, ea);
        double gtau_in_s = gam*hbar*tau(m, g);
        double exp_int = exp(-ta_min/gtau_in_s) - exp(-ta_max/gtau_in_s);
        double erg_transform_x0_delta = 2.0/(b*ea); // in 1/eV
        // Return result in 1/cm^2 x 1 x 1/eV x 1/eV
        return distance_factor * spectrum(m, g, ea, ks2, snorm, t_eff) * erg_transform_x0_delta * exp_int;
    }
    return 0;
}

double inner_integrand_with_tem(double ea, double tem, double ep, double tp0, double tp1, double m, double g, double norm, double e0, double p, double ks2)
{
    double b = beta(m, ea);
    double ta_env = r_env_in_s/b;
    double ta_tp0ta = ta_from_t(tp0-tem, ea, ep, m);
    double ta_min = fmax(tem, fmax(ta_env, ta_tp0ta));
    double ta_geo = ta_geometry(ea, ep, m);
    double ta_tp1 = ta_from_t(tp1-tem, ea, ep, m);
    double ta_max = fmin(ta_geo, tp1);
    if (ta_max > ta_min)
    {
        double gam = gamma(m, ea);
        double gtau_in_s = gam*hbar*tau(m, g);
        double exp_int = exp(-ta_min/gtau_in_s) - exp(-ta_max/gtau_in_s);
        double erg_transform_x0_delta = 2.0/(b*ea); // in 1/eV
        // Return result in 1/cm^2 x 1 x 1/eV x 1/eV
        return distance_factor * temporal_spectrum(m, g, ea, tem, norm, e0, p, ks2) * erg_transform_x0_delta * exp_int;
    }
    return 0;
}

std::vector<double> acceptance_fraction(std::vector<double> eas, std::pair<double,double> ebin, std::pair<double,double> tbin, double m, double g)
{
    std::random_device rng;
    std::mt19937 gen(rng());
    std::uniform_real_distribution<double> unif_x0 (-1, 1);
    int n_trials = eas.size();
    int counts = 0;
    for (auto ea : eas)
    {
        std::exponential_distribution<double> exp_decay (1.0/lam(m, g, ea));
        double d_dec = exp_decay(gen);
        if (d_dec > r_env)
        {
            double b = beta(m, ea);
            // Draw random angles in the axion rest frame for x0 = cos(theta0) in U(-1,1)
            double x0 = unif_x0(gen);
            double bx0 = b * x0;
            // Transform photon angle and energy into the lab frame:
            double ep = 0.5 * ea * (1.0 + bx0);
            if ((ep >= ebin.first) and (ep <= ebin.second))
            {
                double x = (b + x0) / (1.0 + bx0);
                double determinant = d*d - d_dec*d_dec*(1.0 - x*x);
                if (determinant >= 0)
                {
                    double d_gamma = sqrt(determinant) - d_dec*x;
                    double t = (d_dec / b + d_gamma - d) / c0;
                    if ((t >= tbin.first) and (t <= tbin.second)) { counts += 1; }
                }
            }
        }
    }
    std::vector<double> res = { (2.0*counts)/n_trials, (2.0*sqrt(counts))/n_trials };
    return res;
}

double acceptance_fraction_vegas(double ea, double u_ldec, double u_x0, std::pair<double,double> ebin, std::pair<double,double> tbin, double m, double g)
{
    double d_dec = -lam(m, g, ea)*log(1.0 - u_ldec);
    if (d_dec > r_env)
    {
        double bea = sqrt(ea*ea - m*m);
        double x0 = -1.0 + 2.0*u_x0;
        double beax0 = bea * x0;
        // Transform photon angle and energy into the lab frame:
        double ep = 0.5 * (ea + beax0);
        if ((ep >= ebin.first) and (ep <= ebin.second))
        {
            double x = (bea + ea * x0) / (ea + beax0);
            double determinant = d*d - d_dec*d_dec*(1.0 - x*x);
            if (determinant >= 0)
            {
                double d_gamma = sqrt(determinant) - d_dec*x;
                double t = (ea * d_dec / bea + d_gamma - d) / c0;
                if ((t >= tbin.first) and (t <= tbin.second)) { return 2; }
            }
        }
    }

    return 0;
}