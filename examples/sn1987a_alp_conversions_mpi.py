import os
import sys
import numpy as np

from scipy.integrate import quad
from scipy.interpolate import interp1d

script_path = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, script_path+"/../code")

from gammaALPs.core import Source, ALP, ModuleList

import pysnax as sn

from utils import run_mpi_job
from physics import payez_e0, spectrum_payez

e_in_GeV = np.linspace(0.01, 0.1, 201)
src = Source(z=1.1555759101154975e-05, l=279.703427, b=-31.937066)
pa_in = np.diag([0, 0, 1])


def job(m, g, ebin, tbins):
    ml = ModuleList(ALP(m=1e9*m, g=1e20*g), src, pin=pa_in, EGeV=e_in_GeV, log_level='error')
    ml.add_propagation("GMF", 0, model='jansson12')
    px, py, pa = ml.run()
    p_int = interp1d(1e9*e_in_GeV, px[0]+py[0], kind='cubic', bounds_error=False, fill_value=0)
    lmspec = lambda e, t, g: spectrum_payez(e, t, g)*p_int(e)
    integrand = lambda t, g: quad(lmspec, ebin[0], ebin[1], args=(t, g), limit=1000, points=[payez_e0(t)])[0]
    pred =  [quad(integrand, tbins[i], tbins[i+1], args=(g), epsrel=1e-7, epsabs=1e-7, limit=1000)[0] for i in range(len(tbins)-1)]
    return pred


ebins0 = [[10e6, 25e6], [25e6, 100e6]]
tbins0 = [max(0.005, min(2.048*i, 18.0)) for i in range(10)]
lgmvals0 = np.linspace(-11, -8, 151)
lggvals0 = np.linspace(-21, -19, 101)

run_mpi_job(job, script_path+"/../results/sn1987a_alp_conversions_quad", ebins0, tbins0, lgmvals0, lggvals0, save_temp_results=False, save_timing_info=False)
