import os
import sys
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.abspath(script_dir+"/../"))
out_dir = script_dir+"/../../../results/"

import numpy as np

from utils import run_mpi_job
from alp_decay import compute_decay_photon_fluence

job = compute_decay_photon_fluence

ebins0 = [[10e6, 25e6], [25e6, 100e6]] 
tbins0 = [0, 223.232]
lggvals0 = np.linspace(-22, -12, 26)
lgmvals0 = np.linspace(1, 9, 26)

# Running the MPI job should take about ~ 0.4 CPUh (6 min on 4 cores)
run_mpi_job(job, out_dir+"sn1987a_alp_decays_quad", ebins0, tbins0, lgmvals0, lggvals0, save_temp_results=True, save_timing_info=False)