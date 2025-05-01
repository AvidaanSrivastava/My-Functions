# common_functions

Install using: pip install --upgrade --force-reinstall git+https://github.com/AvidaanSrivastava/My-Functions.git

Repository containing a python script which contains all useful python functions I've created. The following functions are in the script (more to be added as needed):

1. p_to_a: Get semi-major axis (a) from orbital period (P) and mass of the star (M_star)
2. calc_Teq: Get equilibrium temperature (Teq) and dayside temperature (Tday) for a planet given stellar temperature (T_star), semi-major axis (a), albedo (A), stellar radius (R_star) and heat distribution factor (q)
3. calc_vesc: Get escape velocity provided mass (M_planet) and radius (R_planet)
4. calc_insolation: Get insolation of an exoplanet relative to solar insolation given stellar temperature (T_star), stellar radius (R_star) and semi-major axis (a)
5. planck_func: Get specific intensity (freq space) from the Planck function given wavelength (wave) and temperature of source (T)
6. ESM: Get emission spectroscopic metric for an exoplanet to see if it is a good pick for JWST phase curves/secondary eclipse (Kempton et al. 2018: https://iopscience.iop.org/article/10.1088/1538-3873/aadf6f) given planet dayside temperature (T_day), stellar temperature (T_star), radius of planet (R_planet), radius of star (R_star) and magnitude of the star in K-band (k_mag)
7. calc_P_roche: Get Roche period and its ration to the orbital period as a measure of tidal disruption of a planet (Dai et al. 2024: https://iopscience.iop.org/article/10.3847/1538-3881/ad5a7d) given mass of the planet (M_planet), radius of the planet (R_planet) and orbital period (P)
8. etc_noiseless: Get number of measurements needed to obtain the RV precision gievn accuracey per point (svrad), expected semi-aplitude (kexp) and detection significance threshold (sig_thresh)
9. etc_noise: Get number of measurements needed to obtain the RV precision gievn accuracey per point (svrad), expected semi-aplitude (kexp), detection significance threshold (sig_thresh) and white noise (noise). The default noise level is 0.5-m/s and can be changed
10. rms: Get the root mean square of an array (arr)
11. sin_func: Get a sine wave given an array of timestamps (t), and constant amplitude (K) and phase (phi)
