import numpy as np

from scipy.stats import poisson
from scipy.optimize import minimize

a_eff = [28, 115, 63] # cm^2

# Gaussian approximation of the likelihood for tp = [0, 223] s, ep = [25, 100] MeV
# as used in [arXiv:1702.02964]
sigma = np.sqrt(1393)/a_eff[2] # ~ 0.6/cm^2
def gaussian_loglike(x):
    y = x/sigma
    return -0.5*y*y

# Profile likelihoods

t_nu = 7*3600+35*60+41.374 # 7:35:41.374 UTC
dt = 2.048

def poissonian_loglike_bkg_single(b, t_off, data_off, bkg_factor=1):
   tmod = (t_off-t_nu)/dt
   ll = poisson.logpmf(k=data_off, mu=bkg_factor*(b[0]+b[1]*tmod))
   return np.sum(ll)

def poissonian_loglike_bkg(b, data_off, bkg_factor=1, n_ebins=2):
   tmod = (data_off[:,0]-t_nu)/dt
   ll = [poisson.logpmf(k=d, mu=bkg_factor*(b[2*i]+b[2*i+1]*tmod)) for i,d in enumerate(data_off[:,(-n_ebins):].T)]
   return np.sum(ll)

def poissonian_loglike(fluxes, data_off, data_on, bkg_factor=1, n_ebins=2):
    means = [np.mean(d) for d in data_off[:,(-n_ebins):].T]
    b0guess = sum([[mu/bkg_factor, 0.01] for mu in means], [])
    b0bounds = sum([[[max(1, (mu-2.0*np.sqrt(mu))/bkg_factor), (mu+2.0*np.sqrt(mu))/bkg_factor], [-0.1*bkg_factor, 0.1*bkg_factor]] for mu in means] ,[])
    tmod = (data_on[:,0]-t_nu)/dt
    ll_s = lambda b: np.sum([poisson.logpmf(k=d, mu=b[2*i]+b[2*i+1]*tmod+a_eff[-n_ebins+i]*fl) for i,(d,fl) in enumerate(zip(data_on[:,(-n_ebins):].T, fluxes))])
    ll_b = lambda b: poissonian_loglike_bkg(b, data_off, bkg_factor, n_ebins)
    ll = lambda b: -2.0*(ll_b(b) + ll_s(b))
    res = minimize(ll, x0=b0guess, bounds=b0bounds, method='Nelder-Mead')
    return -0.5*res.fun