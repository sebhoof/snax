import os
import sys
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.abspath(script_dir+"/../"))
out_dir = script_dir+"/../../../results/"

import numpy as np

from utils import run_mpi_job
from alp_decay import compute_massless_cdf, improved_expected_photon_fluence_from_mc

# Number of particles of fixed m and g in the simulation
# N.B. Results sould be stable for number_alps > 1e7 (arXiv:1702.02964)
# Check the error of 'afrac' if necessary
n_alps = int(1e6)

# Compute the CDF for the low-mass energy spectrum:
inv_cdf0 = compute_massless_cdf()

ebins0 = [[25e6, 100e6]]
tbins0 = [0, 223] # Can only have length 2 for MC routines!
lggvals0 = np.linspace(-22.05, -12, 51)
lgmvals0 = np.linspace(1, 9, 51)

# Running the MPI job should take about ~ 0.6 CPUh (9 min on 4 cores)
kwargs = { 'inv_cdf0': inv_cdf0, 'n_alps': n_alps }
run_mpi_job(improved_expected_photon_fluence_from_mc, out_dir+"/sn1987a_alp_decays_mc", ebins0, tbins0, lgmvals0, lggvals0, save_temp_results=True, save_timing_info=False, **kwargs)