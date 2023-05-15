import numpy as np

from scipy.integrate import quad
from scipy.interpolate import interp1d
from gammaALPs.core import Source, ALP, ModuleList
from physics import distance_factor, payez_e0, spectrum_payez

# Define the SN1987A source position, interpolation energies, initial state (= pure ALP)
src = Source(z=1.1555759101154975e-05, l=279.703427, b=-31.937066)
e_in_GeV = np.linspace(0.01, 0.1, 201)
pa_in = np.diag([0, 0, 1])

def compute_alp_conversion_photon_counts(m, g, ebin, tbins):
    # gammaALPs needs m [neV] and g [1e-11 GeV^-1]
    ml = ModuleList(ALP(m=1e9*m, g=1e20*g), src, pin=pa_in, EGeV=e_in_GeV, log_level='error')
    ml.add_propagation("GMF", 0, model='jansson12')
    px, py, _ = ml.run()
    p_int = interp1d(1e9*e_in_GeV, px[0]+py[0], kind='cubic', bounds_error=False, fill_value=0)
    lmspec = lambda e, t, g: distance_factor*spectrum_payez(0, g, e, t)*p_int(e)
    integrand = lambda t, g: quad(lmspec, ebin[0], ebin[1], args=(t, g), limit=1000, points=[payez_e0(t)])[0]
    pred =  [quad(integrand, tbins[i], tbins[i+1], args=(g), epsrel=1e-7, epsabs=1e-7, limit=1000)[0] for i in range(len(tbins)-1)]
    return pred