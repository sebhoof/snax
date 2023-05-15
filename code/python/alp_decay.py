import sys, os
import numpy as np

from scipy.integrate import quad

script_dir = os.path.dirname(os.path.realpath("__file__"))
sys.path.append(os.path.abspath(script_dir+"/../lib/"))

import pysnax as sn

from physics import find_ea_peak

cnst = sn.constants
ks2_payez15 = cnst["ks2_payez15"]
snorm_payez15 = cnst["snorm_payez15"]
t_eff_payez15 = cnst["t_eff_payez15"]


# m in eV, g in eV^-1
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

    return np.array(results)