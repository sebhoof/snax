import os
import sys
script_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.abspath(script_dir+"/../lib/"))

import numpy as np
import pysnax as sn

from scipy.special import exprel
from scipy.interpolate import interp1d
from scipy.integrate import quad
from scipy.optimize import minimize_scalar

cnst = sn.constants

"""
Physical constants as used in the C++ code
"""
c = cnst['c']  # Speed of light (in m/s)

"""
Supernova 1987a constants
"""

d = cnst['d']  # Distance from Earth to SN1987A (in m)
distance_factor = cnst['distance_factor']
d_in_s = d/c
r_env = cnst['renv']  # SN1987a envelope below which photons are not released (in m)
r_env_in_s = r_env/c
ks = np.sqrt(cnst['ks2_jaeckel17'])  # Effective Debye screening scale (in eV)
spectral_norm = cnst['snorm_jaeckel17']  # Normalisation constant for the spectrum (in eV^-1)
t_eff = cnst['t_eff_jaeckel17']  # Effective temperature (in eV)

"""
Experimental parameters
"""

dt = 223.0  # Duration of the measurement (in s)
# Min. and max. values of the energy (the experimental window; in eV)
erg_window_lo = 25.0e6
erg_window_hi = 100.0e6

"""
Auxiliary functions
"""

# ALP decay length (in m) [Eq. (3) in arXiv:1702.02964]
def l_ALP(m, g, ea):
    x = 1.0e7 / m
    x2 = x * x
    y = 1.0e-19 / g
    return 4.0e13 * sn.beta(m, ea) * (ea / 1.0e8) * (x2 * x2) * (y * y)

# ALP energy spectrum (in eV^-1) [Eq. (7) in arXiv:1702.02964]
def spectrum(m, g, ea):
    return spectral_norm * ea * t_eff * sn.sigma(m, g, ea) / exprel(ea/t_eff)

def find_ea_peak(m):
   res = minimize_scalar(lambda e: -sn.spectrum(m, 1e-18, e), bracket=(m, max(2*m, 2*t_eff), max(10*m, 100*t_eff)))
   return res.x

def axion_fluence(m, g):
    ea_peak = find_ea_peak(m)
    pts0 = np.array([0.1, 1, 100])*ea_peak
    max_erg = max(100*ea_peak, 100*m)
    if (m < 1e4):
        res = quad(lambda e: spectral_norm*e*e*sn.sigma(0, g, e)/np.expm1(e/t_eff), m, max_erg, points=pts0)[0]
    else:
        res = quad(lambda e: sn.spectrum(m, g, e), m, max_erg, points=pts0)[0]
    return distance_factor * res

# Improved version of MC routines still available at https://github.com/marie-lecroq/ALP-fluence-calculation

# Number of particles of fixed m and g in the simulation
# N.B. Results sould be stable for number_alps > 1e7 (arXiv:1702.02964)
# Check the error of 'afrac' if necessary
number_alps = int(1e5)

# 'Naive' ALP fluence from SN1987A in (cm^-2)
def axion_fluence_sn1987a(m, g):
    max_erg = max(250*t_eff, 100*m)
    res = quad(lambda e: sn.spectrum(m, g, e), m, max_erg)[0]
    return distance_factor * res


# ALP emission spectrum data from Table 1, Payez et al. '15 [arXiv:1410.3747]
payez_data = np.array([
[5.000e-03, 6.240e-02, 3.520e+01, 2.250e+00],
[2.000e-01, 3.940e-01, 7.730e+01, 2.020e+00],
[5.000e-01, 8.050e-01, 9.850e+01, 2.065e+00],
[1.000e+00, 1.030e+00, 1.056e+02, 2.145e+00],
[1.500e+00, 1.100e+00, 1.076e+02, 2.190e+00],
[2.000e+00, 1.110e+00, 1.083e+02, 2.220e+00],
[3.000e+00, 1.070e+00, 1.078e+02, 2.280e+00],
[4.000e+00, 9.850e-01, 1.067e+02, 2.315e+00],
[5.000e+00, 8.820e-01, 1.053e+02, 2.340e+00],
[6.000e+00, 7.740e-01, 1.039e+02, 2.350e+00],
[7.000e+00, 6.660e-01, 1.024e+02, 2.355e+00],
[8.000e+00, 5.680e-01, 1.008e+02, 2.355e+00],
[9.000e+00, 4.850e-01, 9.940e+01, 2.350e+00],
[1.000e+01, 4.110e-01, 9.750e+01, 2.350e+00],
[1.100e+01, 3.490e-01, 9.580e+01, 2.350e+00],
[1.200e+01, 2.980e-01, 9.370e+01, 2.350e+00],
[1.300e+01, 2.530e-01, 9.160e+01, 2.350e+00],
[1.400e+01, 2.140e-01, 8.950e+01, 2.350e+00],
[1.500e+01, 1.780e-01, 8.720e+01, 2.355e+00],
[1.600e+01, 1.550e-01, 8.580e+01, 2.355e+00],
[1.700e+01, 1.300e-01, 8.290e+01, 2.370e+00],
[1.800e+01, 1.110e-01, 8.020e+01, 2.385e+00]
])

payez_c = interp1d(payez_data[:,0], 1e46*payez_data[:,1], kind='cubic', bounds_error=False, fill_value='extrapolate')
payez_e0 = interp1d(payez_data[:,0], 1e6*payez_data[:,2], kind='cubic', bounds_error=False, fill_value='extrapolate')
payez_beta = interp1d(payez_data[:,0], payez_data[:,3], kind='cubic', bounds_error=False, fill_value='extrapolate')

def spectrum_payez(m, g, ea, ta):
    x = g/1e-19
    y = ea/payez_e0(ta)
    b = payez_beta(ta)
    res = x*x * payez_c(ta) * pow(y, b) * np.exp(-(b+1.0)*y)
    if m > 1e4:
        ks2 = sn.constants['ks2_payez15']
        res *= sn.sigma(m, g, ea, ks2)/sn.sigma(0, g, ea, ks2)
    return res
