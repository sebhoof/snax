import os
import time
import ctypes
import numpy as np
import numpy.random as rd
import matplotlib.pyplot as plt

from scipy import LowLevelCallable
from scipy.integrate import quad
from scipy.interpolate import interp1d

from physics import *
import pysnax as sn

# Number of particles of fixed m and g in the simulation
# N.B.: according to arXiv:1702.02964, results are stable for number_alps > 1e7
number_alps = int(1e8)

# Min. and max. values of the energy (for integrating the spectrum; in eV)
min_erg = 0.0
max_erg = 5.0e8

# Axion fluence from SN1987A (= naive fluence for one photon per axion) in (cm^-2)
def axion_fluence(m, g):
    res = quad(lambda e: spectrum(m, g, e), min_erg, max_erg)[0]
    return distance_factor * res

# Compute the CDF for the low-mass energy spectrum:
erg_vals = np.linspace(min_erg, max_erg, num=5001)
cdf_vals = []
# N.B. It is okay to fix the coupling g to a convenient value as it just rescales the spectrum.
norm = quad(lambda ea: spectrum(0.0, 1.0e-20, ea), min_erg, max_erg)[0]
cdf_vals = [
    quad(lambda ea: spectrum(0.0, 1.0e-20, ea) / norm, min_erg, erg)[0]
    for erg in erg_vals
]
inv_cdf_massless = interp1d(cdf_vals, erg_vals, kind="linear")


"""
Main functions
"""


# Theoretically measured fluence given mass m/eV and ALP-photon coupling g/eV^-1 (in cm^-2)
def expected_photon_fluence(m, g, verbose=0):
    start = time.time()
    naive_one_photon_fluence = axion_fluence(m, g)
    inv_cdf = inv_cdf_massless
    # Only recalculate CDF from spectrum for heavier axions;
    # otherwise the m = 0 version above is used.
    if m > 1.0e6:
        new_max_erg = max(max_erg, 100.0 * m)
        # Rescale the integral for numerical stability.
        # N.B. Can fix the value of g to a convenient value as it just rescales the spectrum.
        erg_pivot = 5.0 * m
        int_rescaling = spectrum(m, 1.0e-20, erg_pivot)
        erg_vals = np.linspace(m, new_max_erg, num=1001)
        norm = quad(
            lambda ea: spectrum(m, 1.0e-20, ea) / int_rescaling, m, new_max_erg
        )[0]
        norm = norm * int_rescaling
        cdf_vals = [
            quad(lambda ea: spectrum(m, 1.0e-20, ea) / norm, m, erg)[0]
            for erg in erg_vals
        ]
        inv_cdf = interp1d(cdf_vals, erg_vals, kind="linear")

    counts = 0

    # Generate random numbers from [0,1] and transform to random energies via the inverse CDF.
    rng_ergs = rd.random_sample(number_alps)
    rng_ergs = inv_cdf(rng_ergs)
    # Consider only energies where an ALP of mass can be produced (should be true for all).
    rng_ergs = rng_ergs[rng_ergs > m]

    for ea in rng_ergs:
        l_dec = l_ALP(m, g, ea)
        l1 = rd.exponential(scale=l_dec, size=None)
        if (l1 < d) and (l1 > r_env):
            bf = betafactor(m, ea)
            # Draw random angles in the axion rest frame;
            # cos(angle) from [-1,1] (i.e. angle in [0,pi]) and phi from [0,2pi]:
            phi = 2.0 * np.pi * rd.random_sample()
            cos_theta = 2.0 * rd.random_sample() - 1.0
            sin_theta = np.sin(np.arccos(cos_theta))
            cos_phi = np.cos(phi)
            # Photon 1. Transform angle and energy into the lab frame:
            ep1 = 0.5 * ea * (1.0 + bf * sin_theta * cos_phi)
            if (ep1 > erg_window_lo) and (ep1 < erg_window_hi):
                angle1 = np.arccos(
                    (bf + cos_phi * sin_theta)
                    / np.abs(1.0 + bf * cos_phi * sin_theta)
                )
                y1 = l1 * np.sin(angle1)
                L21 = -l1 * np.cos(angle1) + np.sqrt(d * d - y1 * y1)
                time1 = (l1 / bf + L21 - d) / c
                if (time1 < dt) and (time1 > 0):
                    counts += 1
            # Photon 2. Transform angle and energy into the lab frame:
            ep2 = 0.5 * ea * (1.0 - bf * sin_theta * cos_phi)
            if (ep2 > erg_window_lo) and (ep2 < erg_window_hi):
                angle2 = np.arccos(
                    (bf + cos_phi * sin_theta)
                    / np.abs(1.0 + bf * cos_phi * sin_theta)
                )
                y2 = l1 * np.sin(angle2)
                L22 = -l1 * np.cos(angle2) + np.sqrt(d * d - y2 * y2)
                time2 = (l1 / bf + L22 - d) / c
                if (time2 < dt) and (time2 > 0):
                    counts += 1
    lgm, lgg = np.log10(m), np.log10(g) + 9.0
    res = naive_one_photon_fluence * counts / float(number_alps)
    err = res / np.sqrt(counts) if (counts > 0) else np.nan
    if verbose > 1:
        print(
            "Calculation took {:.1f} mins. Result for {:.3f}, {:.3f} is {:.6e}.".format(
                (time.time() - start) / 60.0, lgm, lgg, res, err
            )
        )
    elif verbose > 0:
        print("{:.3f} {:.3f} {:.8e} {:.8e}".format(lgm, lgg, res, err))
    return np.array([lgm, lgg, res, err])


# Theoretically measured fluence given mass m/eV and ALP-photon coupling g/eV^-1 (in cm^-2)
def expected_photon_fluence_from_mc_fast(ebin, tbin, m, g, verbose=0):
    start = time.time()
    naive_one_photon_fluence = axion_fluence(m, g)
    
    new_max_erg = max(max_erg, 100.0 * m)
    # Rescale the integral for numerical stability.
    # N.B. Can fix the value of g to a convenient value as it just rescales the spectrum.
    erg_pivot = 5.0 * m
    int_rescaling = spectrum(m, 1.0e-20, erg_pivot)
    erg_vals = np.linspace(m, new_max_erg, num=1001)
    norm = quad(lambda ea: spectrum(m, 1.0e-20, ea) / int_rescaling, m, new_max_erg)[0]
    norm = norm * int_rescaling
    cdf_vals = [quad(lambda ea: spectrum(m, 1.0e-20, ea) / norm, m, erg)[0] for erg in erg_vals]
    inv_cdf = interp1d(cdf_vals, erg_vals, kind="linear")

    # Generate random numbers from [0,1] and transform to random energies via the inverse CDF.
    rng_ergs = rd.random_sample(number_alps)
    rng_ergs = inv_cdf(rng_ergs)
    #plt.hist(rng_ergs)
    #plt.show()
    # Consider only energies where an ALP of mass can be produced (should be true for all).
    rng_ergs = rng_ergs[rng_ergs > m]
    assert len(rng_ergs) == number_alps

    afrac = sn.acceptance_fraction(rng_ergs, ebin, tbin, m, g)

    lgm, lgg = np.log10(m), np.log10(g) + 9.0
    res = naive_one_photon_fluence * afrac[0]
    err = naive_one_photon_fluence * afrac[1]
    if verbose > 1:
        print(
            "Calculation took {:.1f} mins. Result for {:.3f}, {:.3f} is {:.6e}.".format(
                (time.time() - start) / 60.0, lgm, lgg, res, err
            )
        )
    elif verbose > 0:
        print("{:.3f} {:.3f} {:.8e} {:.8e}".format(lgm, lgg, res, err))
    return np.array([lgm, lgg, res, err])


# Alternative routine that saves the result for each parameter point to a local, unique file for
# each process. Not recommended for production runs; only for debugging and unstable systems.
def expected_photon_fluence_checkpoints(m, g, outpath, rank):
    lgm, lgg = np.log10(m), np.log10(g) + 9.0
    outfile = open(outpath + "/local_results_{}.txt".format(rank), "a")
    res = expected_photon_fluence(m, g)
    outfile.write("{:.3f} {:.3f} {:.10e}\n".format(lgm, lgg, res[2]))
    outfile.close()
    return res
