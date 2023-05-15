import sys
import os
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.abspath(script_dir+"/../lib/"))

import numpy as np
import numpy.random as rd
import pysnax as sn

from scipy.integrate import quad
from scipy.interpolate import interp1d
from physics import axion_fluence_sn1987a, find_ea_peak

cnst = sn.constants
ks2_payez15 = cnst["ks2_payez15"]
snorm_payez15 = cnst["snorm_payez15"]
t_eff_payez15 = cnst["t_eff_payez15"]

# Compute the ALP decay photon fluence with the integral formalism, Eq. (2.8) in [arXiv:2212.09764]
def compute_decay_photon_fluence(m, g, ebin, tbins, ks2=ks2_payez15, snorm=snorm_payez15, t_eff=t_eff_payez15):
    ea_peak = find_ea_peak(m)
    pts0 = ea_peak*np.array([0.1, 1, 10])

    def f0(ea, ep, tp0, tp1):
        return sn.inner_integrand(ea, ep, tp0, tp1, m, g, ks2=ks2, snorm=snorm, t_eff=t_eff)

    def f1(ep, tp0, tp1):
        ea_min = sn.calc_ea_min(m, ep)
        ea_max = sn.calc_ea_max_weak(m, ep, tp0, 0)
        ea_max = min(ea_max, max(1000*ea_peak, 1000*m))
        return quad(f0, ea_min, ea_max, args=(ep, tp0, tp1), limit=2000, points=pts0)[0]

    results = [quad(f1, ebin[0], ebin[1], args=(tbins[i], tbins[i+1]), limit=2000)[0] for i in range(len(tbins)-1)]

    return np.array(results) # cm^-2

# Function to compute the inverse CDF for the m = 0 case.
# To be used together with MC rounties below.
def compute_massless_cdf(n_interp=2001):
    max_erg = 250*cnst['t_eff_jaeckel17']
    erg_vals = np.linspace(0.1, max_erg, num=n_interp)
    cdf_vals = []
    # N.B. It is okay to fix the coupling g to a convenient value as it just rescales the spectrum.
    norm = quad(lambda ea: sn.spectrum0(1e-30, ea), 0.1, max_erg)[0]
    cdf_vals = [quad(lambda ea: sn.spectrum0(1e-30, ea) / norm, 0.1, erg)[0] for erg in erg_vals]
    inv_cdf = interp1d(cdf_vals, erg_vals, kind="linear")
    return inv_cdf

# Compute the ALP decay photon fluence with direct, brute-force Monte Carlo simulations
# Improved version of MC routines, originally created for [arXiv:2205.13549], based on [arXiv:1702.02964]
# Previous code available at https://github.com/marie-lecroq/ALP-fluence-calculation
# Theoretically measured fluence given mass m/eV and ALP-photon coupling g/eV^-1 (in cm^-2)
def improved_expected_photon_fluence_from_mc(m, g, ebin, tbin, inv_cdf0, n_alps):
    naive_one_photon_fluence = axion_fluence_sn1987a(m, g)
    
    if (m < 1e4):
        inv_cdf = inv_cdf0
    else:
        new_max_erg = max(250*cnst['t_eff_jaeckel17'], 100*m)
        # Rescale the integral for numerical stability.
        # N.B. Can fix the value of g to a convenient value as it just rescales the spectrum.
        erg_pivot = 5.0 * m
        int_rescaling = sn.spectrum(m, 1.0e-20, erg_pivot)
        erg_vals = np.linspace(m, new_max_erg, num=1001)
        norm = quad(lambda ea: sn.spectrum(m, 1.0e-20, ea) / int_rescaling, m, new_max_erg)[0]
        norm = norm * int_rescaling
        cdf_vals = [quad(lambda ea: sn.spectrum(m, 1.0e-20, ea) / norm, m, erg)[0] for erg in erg_vals]
        inv_cdf = interp1d(cdf_vals, erg_vals, kind="linear")

    # Generate random numbers from [0,1] and transform to random energies via the inverse CDF.
    rng_ergs = rd.random_sample(n_alps)
    rng_ergs = inv_cdf(rng_ergs)
    # Consider only energies where an ALP of mass can be produced (for massless ALP approx. above).
    rng_ergs = rng_ergs[rng_ergs > m]
    afrac = sn.acceptance_fraction(rng_ergs, ebin, tbin, m, g)
    # N.B. afrac[1] contains the error estimate, which the user should also check
    res = naive_one_photon_fluence * afrac[0]

    return [res] # cm^-2