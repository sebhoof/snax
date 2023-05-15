import sys, os
script_path = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, script_path+"/../code")

import numpy as np

from utils import run_mpi_job
from alp_decay import compute_decay_photon_counts

job = compute_decay_photon_counts

ebins0 = [[10e6, 25e6], [25e6, 100e6]] 
tbins0 = [0, 223.232]
lggvals0 = np.linspace(-22, -12, 51)
lgmvals0 = np.linspace(1, 9, 51)

run_mpi_job(job, script_path+"/../results/sn1987a_alp_decays_quad", ebins0, tbins0, lgmvals0, lggvals0, save_temp_results=True, save_timing_info=False)

