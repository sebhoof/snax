import os, sys
script_path = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, script_path+"/../code")

import numpy as np

from utils import run_mpi_job
from physics import compute_alp_conversion_photon_counts

job = compute_alp_conversion_photon_counts

ebins0 = [[10e6, 25e6], [25e6, 100e6]]
tbins0 = [max(0.005, min(2.048*i, 18.0)) for i in range(10)]
lgmvals0 = np.linspace(-11, -8, 151)
lggvals0 = np.linspace(-21, -19, 101)

run_mpi_job(job, script_path+"/../results/sn1987a_alp_conversions_quad", ebins0, tbins0, lgmvals0, lggvals0, save_temp_results=False, save_timing_info=False)
