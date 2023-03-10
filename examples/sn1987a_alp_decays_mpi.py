import sys, os
import numpy as np

from scipy.integrate import quad

script_path = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, script_path+"/../code")

import pysnax as sn

from utils import run_mpi_job
from physics import find_ea_peak


cnst = sn.constants
ks2_payez15 = cnst["ks2_payez15"]
t_eff_payez15 = cnst["t_eff_payez15"]
snorm_payez15 = cnst["snorm_payez15"]

def job(m, g, ebin, tbins):
    results = []
    ea_peak = find_ea_peak(m)
    pts0 = ea_peak*np.array([0.1, 1, 10])

    def f0(ea, ep, tp0, tp1):
        return sn.inner_integrand(ea, ep, tp0, tp1, m, g, ks2=ks2_payez15, snorm=snorm_payez15, t_eff=t_eff_payez15)

    def f1(ep, tp0, tp1):
        ea_min = sn.calc_ea_min(m, ep)
        ea_max = sn.calc_ea_max_weak(m, ep, tp0, 0)
        ea_max = min(ea_max, max(1000*ea_peak, 1000*m))
        return quad(f0, ea_min, ea_max, args=(ep, tp0, tp1), limit=2000, points=pts0)[0]

    for i in range(len(tbins)-1):
        integral = quad(f1, ebin[0], ebin[1], args=(tbins[i], tbins[i+1]), limit=2000)[0]
        results.append(integral)

    return results

nvals = 151
ebins0 = [[10e6, 25e6], [25e6, 100e6]] 
tbins0 = [0, 223.232]
lggvals0 = np.linspace(-22, -12, nvals)
lgmvals0 = np.linspace(1, 9, nvals)

run_mpi_job(job, script_path+"/../results/sn1987a_alp_decays_quad", ebins0, tbins0, lgmvals0, lggvals0, save_temp_results=False, save_timing_info=False)

