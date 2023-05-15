import os
import sys
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.abspath(script_dir+"/../"))
out_dir = script_dir+"/../../../results/"

import numpy as np

from utils import run_mpi_job
from alp_conversion import compute_alp_conversion_photon_fluence

job = compute_alp_conversion_photon_fluence

ebins0 = [[10e6, 25e6], [25e6, 100e6]]
tbins0 = [max(0.005, min(2.048*i, 18.0)) for i in range(10)]
lgmvals0 = np.linspace(-11, -8, 11)
lggvals0 = np.linspace(-21, -19, 11)

# Running the MPI job should take about ~ 2.2 CPUh (30 min on 4 cores)
run_mpi_job(job, out_dir+"sn1987a_alp_conversions_quad", ebins0, tbins0, lgmvals0, lggvals0, save_temp_results=True, save_timing_info=False)