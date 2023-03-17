import os
import sys

import numpy as np

import numpy.random as rd
from scipy.integrate import quad
from scipy.interpolate import interp1d

script_path = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, script_path+"/../code")

import pysnax as sn

from utils import run_mpi_job

cnst = sn.constants
distance_factor = cnst["distance_factor"]
teff = cnst["t_eff_jaeckel17"]

# Improved version of MC routines still available at https://github.com/marie-lecroq/ALP-fluence-calculation

# Number of particles of fixed m and g in the simulation
# N.B.: according to arXiv:1702.02964, results are stable for number_alps > 1e7
number_alps = int(1e8)

# Axion fluence from SN1987A (= naive fluence for one photon per axion) in (cm^-2)
def axion_fluence(m, g):
    max_erg = max(250*teff, 100*m)
    res = quad(lambda e: sn.spectrum(m, g, e), m, max_erg)[0]
    return distance_factor * res

# Compute the CDF for the low-mass energy spectrum:
max_erg = 250*teff
erg_vals = np.linspace(0.1, max_erg, num=2001)
cdf_vals = []
# N.B. It is okay to fix the coupling g to a convenient value as it just rescales the spectrum.
norm = quad(lambda ea: sn.spectrum0(1e-30, ea), 0.1, max_erg)[0]
cdf_vals = [
    quad(lambda ea: sn.spectrum0(1e-30, ea) / norm, 0.1, erg)[0]
    for erg in erg_vals
]
inv_cdf_massless = interp1d(cdf_vals, erg_vals, kind="linear")


# Theoretically measured fluence given mass m/eV and ALP-photon coupling g/eV^-1 (in cm^-2)
def improved_expected_photon_fluence_from_mc(m, g, ebin, tbin):
    naive_one_photon_fluence = axion_fluence(m, g)
    
    if (m < 1e4):
        inv_cdf = inv_cdf_massless
    else:
        new_max_erg = max(max_erg, 100.0*m)
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
    rng_ergs = rd.random_sample(number_alps)
    rng_ergs = inv_cdf(rng_ergs)
    # Consider only energies where an ALP of mass can be produced (for massless ALP approx. above).
    rng_ergs = rng_ergs[rng_ergs > m]
    afrac = sn.acceptance_fraction(rng_ergs, ebin, tbin, m, g)
    res = naive_one_photon_fluence * afrac[0]

    return [res]

ebins0 = [[25e6, 100e6]]
tbins0 = [0, 223] # Can only have length 2 for MC routines!
lggvals0 = np.linspace(-22.05, -12, 51)
lgmvals0 = np.linspace(1, 9, 51)

run_mpi_job(improved_expected_photon_fluence_from_mc, script_path+"/../results/sn1987a_alp_decays_mc", ebins0, tbins0, lgmvals0, lggvals0, save_temp_results=True, save_timing_info=False)