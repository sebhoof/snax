import os
import sys

import numpy as np

from scipy.integrate import quad
from scipy.interpolate import PchipInterpolator

import vegas

script_path = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, script_path+"/../code")

import pysnax as sn

from utils import run_mpi_job
from physics import find_ea_peak

cnst = sn.constants
teff = cnst["t_eff_jaeckel17"]
spectral_norm = cnst["snorm_jaeckel17"]
distance_factor = cnst["distance_factor"]
m0 = 1e4

def axion_fluence(m, g, pts0):
    max_erg = max(250*teff, 100*m)
    if (m < m0):
        res = quad(lambda e: spectral_norm*e*e*sn.sigma0(g, e)/np.expm1(e/teff), m, max_erg, points=[2.309*teff])[0]
    else:
        res = quad(lambda e: sn.spectrum(m, g, e), m, max_erg, points=pts0)[0]
    return distance_factor * res

def get_inv_cdf(ergs, cdf_vals):
   ergs_, cdf_vals_ = [ergs[0]], [cdf_vals[0]]
   i0 = 0
   for i in range(1, len(ergs)):
      if cdf_vals[i] > cdf_vals_[i0]:
         ergs_.append(ergs[i])
         cdf_vals_.append(cdf_vals[i])
         i0 += 1
   ergs_[-1] = ergs[-1]
   return PchipInterpolator(cdf_vals_, ergs_)

n_ergs0 = 2000
ea_peak = find_ea_peak(1)
ergs0 = np.linspace(1, 100*ea_peak, n_ergs0)
integrand = lambda ea: sn.spectrum(1, 1e-45, ea)
dcdf_vals = [quad(lambda ea: integrand(ea), ergs0[i], ergs0[i+1], epsabs=1e-10, epsrel=1e-10, limit=500)[0] for i in range(n_ergs0-1)]
cdf_vals0 = [0]
for i in range(n_ergs0-1):
    cdf_vals0.append( cdf_vals0[i] + dcdf_vals[i] )
cdf_vals0 = [c/cdf_vals0[-1] for c in cdf_vals0]
inv_cdf0 = get_inv_cdf(ergs0, cdf_vals0)

def job(m, g, ebin, tbin):
    neval = 5e8
    if (m > 1e3):
        neval = 1e7
    elif (m > 1e5):
        neval = 1e6
    ea_peak = find_ea_peak(m)
    pts = ea_peak*np.array([0.1, 1, 10])
    naive_one_photon_fluence = axion_fluence(m, g, pts0=pts)

    if (m > 1e4):
        n_ergs0 = 2000
        ergs0 = np.linspace(m, max(100*ea_peak, 100*m), n_ergs0)
        integrand = lambda ea: sn.spectrum(m, 1e-45, ea)
        dcdf_vals = [quad(lambda ea: integrand(ea), ergs0[i], ergs0[i+1], epsabs=1e-12, epsrel=1e-12, limit=500)[0] for i in range(n_ergs0-1)]
        cdf_vals0 = [0]
        for i in range(n_ergs0-1):
            cdf_vals0.append( cdf_vals0[i] + dcdf_vals[i] )
        cdf_vals0 = [c/cdf_vals0[-1] for c in cdf_vals0]
        inv_cdf = get_inv_cdf(ergs0, cdf_vals0)
    else:
        inv_cdf = inv_cdf0

    integrand = lambda u: sn.acceptance_fraction_vegas(inv_cdf(u[0]), u[1], u[2], ebin, tbin, m, g)
    integral = vegas.Integrator(3*[[0, 1]])
    result = integral(integrand, nitn=10, neval=neval, mpi=False)
    res = naive_one_photon_fluence * result.mean

    return [res]


ebins0 = [[25e6, 100e6]]
tbins0 = [0, 223] # Can only have length 2 for MC routines!
lggvals0 = np.linspace(-22.05, -12, 151)
lgmvals0 = np.linspace(1, 9, 51)

run_mpi_job(job, script_path+"/../results/sn1987a_alp_decays_mc_vegas", ebins0, tbins0, lgmvals0, lggvals0, save_temp_results=True, save_timing_info=False)
