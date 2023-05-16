# SNax &mdash; Supernovae and axion-like particle signatures

SNax is a code that calculates the signals of axion-like particles (ALPs) from supernovae (SNe).
It provides likelihood computations for ALP conversions in the Galactic magnetic field and ALP decays after SN 1987A.

The code was developed and improved by Sebastian Hoof and Lena Schulz.
Please cite ref. [[1]](#update) when using these.
The Monte Carlo (MC) routines for ALP decays were originally written by Marie Lecroq, Sebastian Hoof, and Csaba Bal&aacute;zs.
Please cite refs [[1,2]](#cosmoalp) when using these.


## Summary

Axion-like particles are produced in SNe and can decay into two photons or convert into photons while traveling through magnetic fields towards Earth.
By observing the absence of such photons from the direction of SN 1987A, we can establish limits on the ALP-photon coupling.

This code includes Monte Carlo (MC) routines for predicting the number of ALP decay photons, as described in ref. [[1]](#cosmoalp).

Additionally, in ref. [[2]](#update), we update the ALP-photon limits by digitizing data from refs [[3,4]](#data1) and developing new integration routines for ALP conversions and decays.


## Results

![ALP-photon limits](results/sn1987a_alp_limits_web.png)

The figure above shows the limits on the ALP-photon coupling (at 95% CL) obtained using our updated likelihood from ref. [[2]](#update).
The left panel displays the limit from ALP conversions in the Galactic magnetic field, while the right panel shows the limit from ALP decays.


## The code

The code implements optimized quadrature-based computation for ALP conversion (requires the `gammaALPs` package) and decays, as well as an older brute-force and improved Monte Carlo integration for ALP decays.
The code is developed for Python v3.9+ and C++-20-compliant compilers.

An simple example for how to use the code is included in the [snax_demo.ipynb](snax_demo.ipynb) Jupyter notebook.
Python scripts similar to those used in the paper are available in the [code/python/mpi_scripts/](code/python/mpi_scripts/) folder (require a correctly installed `mpi4py` package).


## Installation

1. Clone this repository: `git clone https://github.com/sebhoof/snax`.

2. Set up a build directory: `cd snax; mkdir build; cd build/`.

3. For macOS, install CMAKE using [Homebrew](https://brew.sh): `brew install cmake`.

4. Run `cmake ..` and then `make` to build everything. The C++ library and Python module `pysnax` should be generated in the newly created `code/lib` directory.

5. Ensure all required Python packages are installed: `python -m pip install numpy scipy mpi4py gammaALPs vegas`.

## References

If you use our code, please cite refs [[1,2]](#cosmoalp).

<a id="update">[1]</a> S. Hoof and L. Schulz, [&ldquo;*Updated constraints on axion-like particles from temporal information in supernova SN1987A gamma-ray data,*&rdquo;](https://doi.org/10.1088/1475-7516/2023/03/054) JCAP 03 (2023) 054 [[arXiv:2212.09764]](https://arxiv.org/abs/2212.09764).

<a id="cosmoalp">[2]</a> C. Bal&aacute;zs, S. Bloor, T. E. Gonzalo, *et al*. [&ldquo;*Cosmological constraints on decaying axion-like particles: a global analysis,*&rdquo;](https://doi.org/10.1088/1475-7516/2022/12/027) JCAP 12 (2022) 027 [[arXiv:2205.13549]](https://arxiv.org/abs/2205.13549).

<a id="data1">[3]</a> E. L. Chupp, W. T. Vestrand, and C. Reppin. [&ldquo;*Experimental Limits on the Radiative Decay of SN 1987A Neutrinos,*&rdquo;](https://doi.org/10.1103/PhysRevLett.62.505) Phys. Rev. Lett. **62**, 505 (1989).

<a id="data2">[4]</a> L. Oberauer, C. Hagner, G. Raffelt, and E. Rieger, [&ldquo;*Supernova bounds on neutrino radiative decays,*&rdquo;](https://doi.org/10.1016/0927-6505(93)90004-W) Astroparticle Physics **1**, 4 (1993), pp. 377&ndash;386.
