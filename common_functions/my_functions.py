
## This script contains functions to calculate various astrophysical parameters


import numpy as np
from astropy.io import fits

def p_to_a(period, M_star, period_err=None, M_star_err=None):
    """
    Calculate the semi-major axis of a planet given its period and the mass of the star,
    and optionally propagate uncertainties.

    Parameters
    ----------
    period : float
        Period of the planet in days.
    M_star : float
        Mass of the star in solar masses.
    period_err : float, optional
        Uncertainty in the period (days).
    M_star_err : float, optional
        Uncertainty in the stellar mass (solar masses).

    Returns
    -------
    a : float
        Semi-major axis in meters.
    a_AU : float
        Semi-major axis in AU.
    a_err : float, optional
        Uncertainty in semi-major axis in meters (if errors provided).
    a_AU_err : float, optional
        Uncertainty in semi-major axis in AU (if errors provided).
    """
    period_sec = period * 86400   # From days to seconds
    G = 6.67430e-11
    M_star_kg = M_star * 1.989e30   # In kg
    a = (G * M_star_kg * period_sec**2 / (4 * np.pi**2))**(1/3)  # meters
    a_AU = a / 1.496e11  # AU

    if period_err is not None and M_star_err is not None:
        # Propagate errors: da/a = (1/3)*dM/M + (2/3)*dP/P
        rel_err = ((1/3)*(M_star_err/M_star))**2 + ((2/3)*(period_err/period))**2
        a_err = a * np.sqrt(rel_err)
        a_AU_err = a_AU * rel_err
        return a, a_AU, a_err, a_AU_err
    else:
        return a, a_AU


def calc_Teq(T_star, a, R_star, T_star_err=None, a_err=None, R_star_err=None, A=0, q=2/3):

    """ This function calculates the equilibrium temperature of a planet given its distance from the star,
    the star's temperature, the planet's albedo, and the star's radius.
    T_star: temperature of the star in Kelvin
    a: distance from the star in meters
    A: albedo of the planet
    R_star: radius of the star in solar radii
    q: heat redistribution (default is 2/3)
    T_eq: equilibrium temperature in Kelvin
    T_day: day side temperature in Kelvin
    """

    R_star = 6.96e8 * R_star    ## in meters

    x = R_star / (2 * a)
    Tday = T_star * (1 - A)**(1/4) * (x)**(1/2) * (q)**(1/4) * 2 ** .5
    Teq = T_star * (1 - A)**(1/4) * (x)**(1/2)

    if T_star_err is None and R_star_err is None and a_err is None:
       return Teq, Tday
    
    else:
        T_star = T_star
        R_star = R_star
        a = a

        Tday = T_star * (1 - A)**(1/4) * (x)**(1/2) * (q)**(1/4) * 2 ** .5
        Teq = T_star * (1 - A)**(1/4) * (x)**(1/2)

        # Calculate the error in the equilibrium temperature
        dTeq_T = 0.25 * Teq * (T_star_err / T_star)
        dTeq_R = 0.5 * Teq * (R_star_err / R_star)
        dTeq_a = -0.5 * Teq * (a_err / a)

        Teq_err = np.sqrt(dTeq_T**2 + dTeq_R**2 + dTeq_a**2)

        return Teq, Tday, Teq_err


def calc_vesc(M_planet, R_planet):

    """ This function calculates the escape velocity of a planet given its mass and radius.
    M_planet: mass of the planet in Earth masses
    R_planet: radius of the planet in Earth radii
    M: mass of the planet in Kg
    R: radius of the planet in meters
    G: gravitational constant in m^3 kg^-1 s^-2
    vesc: escape velocity in m/s
    """

    
    M = M_planet * 5.97e24     ## in Kg
    R = R_planet * 6.378e6       ## in meters
    G = 6.67430e-11

    vesc = np.sqrt(2 * G * M / R)

    return vesc


def calc_insolation(T_star, R_star, a, T_star_err=None, R_star_err=None, a_err=None):

    """ This function calculates the insolation of a planet given its distance from the star,
    the star's temperature, and the star's radius.
    T_star: temperature of the star in Kelvin
    R_star: radius of the star in solar radii
    a: distance from the star in AU
    S: insolation in units of S_earth
    """

    a = a 
    T_sun = 5778

    if T_star_err is None and R_star_err is None and a_err is None:
        # If no errors are provided, use the nominal values
        T_star = T_star
        R_star = R_star

        S = (R_star)**2 * (T_star/T_sun)**4 / a**2
        return S
    
    else:
        # If errors are provided, calculate the insolation with errors
        T_star = T_star
        R_star = R_star
        a = a

        S = (R_star)**2 * (T_star/T_sun)**4 / a**2

        # Calculate the error in the insolation
        dS_T = 4 * S * (T_star_err / T_star)
        dS_R = 2 * S * (R_star_err / R_star)
        dS_a = -2 * S * (a_err / a)

        S_err = np.sqrt(dS_T**2 + dS_R**2 + dS_a**2)

        return S, S_err



def planck_func(wave, T):

    """ This function calculates the specific intensity using Planck function for a given wavelength and temperature.
    wave: wavelength in meters
    T: temperature in Kelvin
    h: Planck's constant in J*s
    c: speed of light in m/s
    k: Boltzmann's constant in J/K
    nu: frequency in Hz
    B: Planck function in W/m^2/sr/m
    """
    
    h = 6.626e-34
    c = 3.0e8
    k = 1.38e-23
    nu = c / wave

    return 2 * h * nu**3 / c**2 / (np.exp(h * nu / k / T) - 1)


def ESM(T_day, T_star, R_planet, R_star, k_mag):

    """ This function calculates the Emission Spectroscopic Metric using the ESM formula.
    T_day: day side temperature in Kelvin
    T_star: temperature of the star in Kelvin
    R_planet: radius of the planet in Earth radii
    R_star: radius of the star in solar radii
    k_mag: magnitude of the star
    B_day: Planck function for the day side temperature
    B_star: Planck function for the star temperature
    esm: Emission Spectroscopic Metric in solar masses
    """

    B_day = planck_func(7.5*10**(-6), T_day)
    B_star = planck_func(7.5*10**(-6), T_star)

    esm = 4.29*10**6 * (R_planet * 6.278e6/(R_star*6.96e8))**2 * 10**(-1/5*k_mag) * B_day / B_star
    
    return esm


def calc_P_roche(M_planet, R_planet, P):

    """ This function calculates the Roche limit and the ratio of the Roche limit to the orbital period of a planet.
    M_planet: mass of the planet in Earth masses
    R_planet: radius of the planet in Earth radii
    P: period of the planet in days
    P_roche: Roche limit in hours
    P_ratio: ratio of the Roche limit to the period
    """

    P = P * 24       ## in hours

    M_planet = M_planet * 5.97e24   ## in Kg
    R_planet = R_planet * 6.378e6   ## in meters
    rho = M_planet / (4/3 * np.pi * R_planet**3) / 1000   ## in g/cm^3

    P_roche = 12.6 * rho ** (-.5)
    P_ratio = P_roche / P

    return P_roche, P_ratio


def etc_noiseless(svrad, kexp, sig_thresh):

    '''Using ESPRESSO ETC to calculate no. of measurements required to reach a specified RV precision.

    Parameters
    ----------  
    svrad : float
        RV precision from ETC (M2 model)

    kexp : float
        Expected semi amplitude of planet

    sig_thresh : float
        Desired RV precision (Eg: 10 for 10sigma)
    '''

    rv_precision = kexp / sig_thresh

    N = 2 / (rv_precision / svrad)**2

    return N, svrad


def etc_noise(svrad, kexp, sig_thresh, noise = .5):

    '''Using ESPRESSO ETC to calculate no. of measurements required to reach a specified RV precision.

    Parameters
    ----------  
    svrad : float
        RV precision from ETC (M2 model)

    kexp : float
        Expected semi amplitude of planet

    sig_thresh : float
        Desired RV precision (Eg: 10 for 10sigma)
    '''

    svrad_net = np.sqrt(svrad**2 + noise**2)
    print(svrad_net)

    rv_precision = kexp / sig_thresh

    N = 2 / (rv_precision / svrad_net)**2

    return N, svrad_net


def rms(arr):
    '''Calculates the RMS of an array.

    Parameters
    ----------
    arr : array
        Array of values

    Returns
    -------
    rms : float
        RMS of the array
    '''

    N = len(arr)
    if N == 0:
        return "ERROR: Empty array"
    elif N > 0:
        rms = np.sqrt(np.sum(arr**2) / N)

    return rms


def sin_func(t, K, phi):

    """ This function calculates the sine function for a given time, amplitude, and phase.
    t: time in seconds
    K: amplitude
    phi: phase in radians
    """

    return K * np.sin(t + phi)

def print_fits_header_names(file):
    """
    Function to get the header of a FITS file.
    :file: Full path of the file on the local machine
    """
    with fits.open(file) as hdul:
        for i in range(len(hdul)):
            print(i, hdul[i].name)


def derivative(y_col, x_col):

    """
    Function to calculate the derivative of a spectrum using the central difference method.
    :param y_col: The y values of the spectrum.
    :param x_col: The x values of the spectrum.
    :return: The derivative of the spectrum.
    """

    y_col = np.array(y_col).tolist()
    x_col = np.array(x_col).tolist()

    no_pixels = len(y_col)
    idx_list = list(map(int, np.linspace(1, no_pixels - 2, no_pixels - 2).tolist()))

    dy_dx_list = []
    for i in idx_list:

        pixel = i

        dy_dx = (y_col[pixel + 1] - y_col[pixel - 1]) / (x_col[pixel + 1] - x_col[pixel - 1])

        dy_dx_list.append(dy_dx)

    return dy_dx_list


def maroon_x_etc_sigRV(snr_peak, Teff):

    """
    Calculate the sigRV for a MAROON-X spectrum of M dwarfs with a given peak SNR and effective temperature.
    
    Parameters:
    snr_peak (float): The peak SNR of the spectrum.
    Teff (float): The effective temperature of the star.
    
    Returns:
    float: The calculated sigRV.
    """

    if Teff < 3500:
        sig_rv = 180 / (snr_peak ** 1.15) + 0.06
    
    elif Teff < 4000 and Teff >= 3500:
        sig_rv = 414 / (snr_peak ** 1.14) + 0.01

    
    return sig_rv


def plot_rv_phase_fold(time, rv, rv_error, period, tc, rv_fit, bin = False, nbins = 13):

    """
    Function to plot the phase folded RV data and the best fit model.

    Parameters
    ----------
    time : array
        Time in days.
    rv : array
        Radial velocity measurements.
    rv_error : array
        Radial velocity errors.
    period : float
        Orbital period of the planet in days.  
    tc : float
        Time of conjunction.
    rv_fit : float
        Best fit radial velocity model from juliet.
    bin : bool
        If True, bin the data.
    nbins : int
        Number of bins to use for the phase folded plot.
    """

    ## Get phases
    phases = juliet.utils.get_phases(time, period, tc)
    idx = np.argsort(phases)

    model_times = np.linspace(np.min(time) - 30, np.max(time) + 30, len(rv_fit))

    phases_kep = juliet.utils.get_phases(model_times, period, tc)
    idx_kep = np.argsort(phases_kep)

    ## Bin the RV data

    bins = np.linspace(-0.5, 0.5,nbins)
    increments = np.diff(bins)
    binned_rv = []
    binned_erv = []

    for i in range(len(bins)-1):

        bin_i = bins[i]
        this_bin_rv = []
        this_bin_erv = []

        for j in range(len(phases[idx])):

            phase_j = phases[idx][j]
            vel_j = rv[idx][j]
            evel_j = rv_error[idx][j]

            if (phase_j >= bin_i) & (phase_j < bins[i+1]):
                this_bin_rv.append(vel_j)
                this_bin_erv.append(evel_j ** 2)

        median_rv = np.median(this_bin_rv)
        std_rv = np.sqrt(np.sum(this_bin_erv)) / len(this_bin_erv)

        binned_rv.append(median_rv)
        binned_erv.append(std_rv)

    binned_rv = np.array(binned_rv)
    binned_erv = np.array(binned_erv)

    print(len(binned_rv))
    print(len(binned_erv))
    print(len(bins[:-1]))
    ## Get the best fit model

    offset = rv_fit[idx_kep][0]
    print('Offset:', offset)

    ## Plot the phase folded data
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.errorbar(phases[idx], rv[idx] - offset, yerr=rv_error[idx], fmt='o', label='RV data', color='blue', markersize=2)

    if bin == True:
        ax.errorbar(bins[:-1] + increments / 2, binned_rv - offset, yerr=binned_erv, fmt='o', label='Binned RV data', color='red', markersize=2)
    ax.plot(phases_kep[idx_kep], rv_fit[idx_kep] - offset, label='Best fit model', color='black')
    ax.set_xlabel('Phase')
    ax.set_ylabel('Radial Velocity (m/s)')
    ax.set_title('Phase Folded Radial Velocity')
    ax.legend()
    ax.grid()
    ax.set_xlim(-0.5, 0.5)
    plt.show()


def R_planet_error(delta, d_delta, r_star, d_r_star, r_planet):
    """
    Calculate the error in the planet radius based on the error in the transit depth and stellar radius.
    
    Parameters:
    delta (float): Transit depth
    d_delta (float): Error in transit depth
    r_star (float): Stellar radius
    d_r_star (float): Error in stellar radius
    r_planet (float): Calculated planet radius
    
    Returns:
    float: Error in planet radius
    """
    R_earth = 6378.1 * 10**3 # m
    r_planet = r_planet
    sigma_rp = np.sqrt((d_r_star / r_star)**2 + (0.5 * d_delta / delta)**2) * r_planet

    return sigma_rp # Return in R_earth
    

def mass_calc(K, Kerr, P, Perr, Mstar, Mstarerr):
    """
    Calculate the minimum mass (Msini) of a planet from RV semi-amplitude and its uncertainty.

    Parameters
    ----------
    K : float
        Radial velocity semi-amplitude (m/s).
    Kerr : float
        Uncertainty in the RV semi-amplitude (m/s).
    P : float
        Orbital period of the planet (days).
    Perr : float
        Uncertainty in the orbital period (days).
    Mstar : float
        Mass of the host star (in solar masses).
    Mstarerr : float
        Uncertainty in the stellar mass (in solar masses).

    Returns
    -------
    Mp : float
        Minimum mass of the planet (Earth masses).
    Mperr : float
        Uncertainty in the minimum mass (Earth masses).
    """

    # Calculate minimum mass (Msini) using radvel's Msini function (returns in Earth masses)
    Mp = radvel.utils.Msini(K=K, P=P, Mstar=Mstar, e=0)

    # Propagate uncertainties using partial derivatives (quadrature sum)
    Mperr = Mp * np.sqrt((Kerr / K)**2 + (2/3 * Mstarerr / Mstar)**2 + (1/3 * Perr / P)**2)

    return Mp, Mperr 


def sample_split_normal(mean, err_low, err_high, nsamples=100000):
    """
    Draw samples from an asymmetric (split) normal distribution.
    
    Parameters
    ----------
    mean : float
        Central value.
    err_low : float
        Lower 1-sigma error (positive number).
    err_high : float
        Upper 1-sigma error (positive number).
    nsamples : int
        Number of samples to draw.
    
    Returns
    -------
    samples : ndarray
        Array of samples.
    """
    # Split samples half/half
    n_low = nsamples // 2
    n_high = nsamples - n_low
    
    low_samples = np.random.normal(mean, err_low, n_low)
    high_samples = np.random.normal(mean, err_high, n_high)
    
    return np.concatenate([low_samples, high_samples])


def inclination_from_transit(b, b_err_low, b_err_high,
                                  aRs, aRs_err_low, aRs_err_high,
                                  nsamples=100000):
    """
    Compute orbital inclination (deg) with asymmetric uncertainties
    from impact parameter (b) and scaled semi-major axis (a/R*).
    
    Parameters
    ----------
    b : float
        Impact parameter.
    b_err_low, b_err_high : float
        Lower and upper 1-sigma uncertainties on b.
    aRs : float
        Scaled semi-major axis a/R*.
    aRs_err_low, aRs_err_high : float
        Lower and upper 1-sigma uncertainties on a/R*.
    nsamples : int
        Number of Monte Carlo samples.
    
    Returns
    -------
    inc_median : float
        Median inclination in degrees.
    inc_err_low : float
        Lower bound (16th percentile distance from median).
    inc_err_high : float
        Upper bound (84th percentile distance from median).
    """
    # Draw split-normal samples
    b_samples = sample_split_normal(b, b_err_low, b_err_high, nsamples)
    aRs_samples = sample_split_normal(aRs, aRs_err_low, aRs_err_high, nsamples)

    # Compute cos(i) = b / (a/R*)
    cosi = b_samples / aRs_samples
    cosi = np.clip(cosi, 0, 1)  # enforce physical bounds

    # Inclination samples
    inc_samples = np.degrees(np.arccos(cosi))

    # Credible interval
    inc_median = np.median(inc_samples)
    inc_lower = np.percentile(inc_samples, 16)
    inc_upper = np.percentile(inc_samples, 84)

    inc_err_low = inc_median - inc_lower
    inc_err_high = inc_upper - inc_median

    return inc_median, inc_err_low, inc_err_high


def planet_density(Mp, Rp, Mp_err=None, Rp_err=None):

    """ This function calculates the density of a planet given its mass and radius.
    Mp: mass of the planet in Earth masses
    Rp: radius of the planet in Earth radii
    rho: density of the planet in g/cm^3
    """

    if Mp_err is None and Rp_err is None:
        # If no errors are provided, use the nominal values
        Mp = Mp
        Rp = Rp

        rho = (Mp) / (Rp**3) * 5.51  # Density in g/cm^3
        return rho
    
    else:
        # If errors are provided, calculate the density with errors
        Mp = Mp
        Rp = Rp

        rho = (Mp) / (Rp**3) * 5.51  # Density in g/cm^3

        # Calculate the error in the density
        drho_M = rho * (Mp_err / Mp)
        drho_R = -3 * rho * (Rp_err / Rp)

        rho_err = np.sqrt(drho_M**2 + drho_R**2)

        return rho, rho_err



def calc_K_semi_amplitude(period, mass_planet, mass_star, e=0, i=90):
    """Calculate the expected RV semi-amplitude K in m/s.

    Args:
        period (float): Orbital period in days.
        mass_planet (float): Planet mass in Earth masses.
        mass_star (float, optional): Stellar mass in Solar masses.
        e (float, optional): Orbital eccentricity. Defaults to 0.
        i (float, optional): Orbital inclination in degrees. Defaults to 90.

    Returns:
        float: Expected RV semi-amplitude K in m/s.
    """

    M_earth = 5.972e24  # Mass of Earth in
    M_p = mass_planet * M_earth  # Convert planet mass to kg

    M_sun = 1.989e30  # Mass of Sun in kg
    M_s = mass_star * M_sun  # Convert stellar mass to kg

    P = period * 86400  # Convert period from days to seconds
    G = 6.67430e-11  # Gravitational constant in m^3 kg^-1 s^-2

    K = (2 * np.pi * G / P) ** (1 / 3) * M_p * np.sin(np.radians(i)) / (M_s + M_p) ** (2 / 3) / np.sqrt(1 - e ** 2)

    return K


def transit_duration(P_days, R_star_Rsun, R_planet_Rearth, M_star_Msun,
                     impact_parameter=0.0, inclination_deg=90.0):
    """
    Calculate total transit duration (T14) for a planet assuming circular orbit.

    Parameters
    ----------
    P_days : float
        Orbital period in days.
    R_star_Rsun : float
        Stellar radius in solar radii.
    R_planet_Rearth : float
        Planetary radius in Earth radii.
    M_star_Msun : float
        Stellar mass in solar masses.
    impact_parameter : float, optional
        Transit impact parameter (0 = central transit).
    inclination_deg : float, optional
        Orbital inclination in degrees (90 = edge-on).

    Returns
    -------
    T14_hours : float
        Total transit duration in hours.
    """
    # Convert to SI units
    P = P_days * day
    R_star = R_star_Rsun * R_sun
    R_planet = R_planet_Rearth * R_earth
    M_star = M_star_Msun * M_sun

    # Semi-major axis (Kepler's 3rd law)
    a = (G * M_star * P**2 / (4 * np.pi**2))**(1/3)

    # Convert inclination to radians
    i = np.radians(inclination_deg)

    # Compute argument for arcsin
    term = (R_star / a) * np.sqrt((1 + R_planet / R_star)**2 - impact_parameter**2) / np.sin(i)
    term = np.clip(term, -1, 1)  # ensure valid range

    # Transit duration (seconds)
    T14 = (P / np.pi) * np.arcsin(term)

    # Convert to hours
    return T14 / 3600.0


def TSM(Rp, Teq, Mp, Rs, jmag):

    '''
    
    Using Kempton et al. 2018 formula to calculate TSM
    
    Rp: Planet radius in Earth radii
    Teq: Planet equilibrium temperature in Kelvin
    Mp: Planet mass in Earth masses
    Rs: Stellar radius in Solar radii
    jmag: Stellar J magnitude
    
    '''

    tsm = 0.19 * (Rp**3 * Teq) / (Mp * Rs**2) * 10**(-1/5*jmag) 
    return tsm


## Getting HARPS 3.5 wavelength solution

def build_HARPS_wavelengths(hdr, array_size = (72, 4096)):
        """Compute the wavelength solution to this given spectra (EQ 4.1 of DRS manual)
        Convert from air wavelenbgths to vacuum
        """
        # degree of the polynomial
        d = hdr["HIERARCH ESO DRS CAL TH DEG LL"]
        # number of orders
        omax = hdr.get("HIERARCH ESO DRS CAL LOC NBO", array_size[0])
        xmax = array_size[1]

        # matrix X:
        #
        # axis 0: the entry corresponding to each coefficient
        # axis 1: each pixel number

        x = np.empty((d + 1, xmax), "int64")
        x[0].fill(1)  # x[0,*] = x^0 = 1,1,1,1,1,...
        x[1] = np.arange(xmax)

        for i in range(1, d):
            x[i + 1] = x[i] * x[1]

        # matrix A:
        #
        # axis 0: the different orders
        # axis 1: all coefficients for the given order

        A = np.reshape(
            [hdr["HIERARCH ESO DRS CAL TH COEFF LL" + str(i)] for i in range(omax * (d + 1))],
            (omax, d + 1),
        )  # slow 30 ms

        # the wavelengths for each order are a simple dot product between the coefficients and pixel-wise data (X)
        wavelengths = np.dot(A, x)

        vacuum_wavelengths = airtovac(wavelengths)
        return vacuum_wavelengths


import numpy as np


def airtovac(wave_air):
    """Taken from idl astrolib
    ;+
    ; NAME:
    ;       AIRTOVAC
    ; PURPOSE:
    ;       Convert air wavelengths to vacuum wavelengths
    ; EXPLANATION:
    ;       Wavelengths are corrected for the index of refraction of air under
    ;       standard conditions.  Wavelength values below 2000 A will not be
    ;       altered.  Uses relation of Ciddor (1996).
    ;
    ; CALLING SEQUENCE:
    ;       AIRTOVAC, WAVE_AIR, [ WAVE_VAC]
    ;
    ; INPUT/OUTPUT:
    ;       WAVE_AIR - Wavelength in Angstroms, scalar or vector
    ;               If this is the only parameter supplied, it will be updated on
    ;               output to contain double precision vacuum wavelength(s).
    ; OPTIONAL OUTPUT:
    ;        WAVE_VAC - Vacuum wavelength in Angstroms, same number of elements as
    ;                 WAVE_AIR, double precision
    ;
    ; EXAMPLE:
    ;       If the air wavelength is  W = 6056.125 (a Krypton line), then
    ;       AIRTOVAC, W yields an vacuum wavelength of W = 6057.8019
    ;
    ; METHOD:
    ;	Formula from Ciddor 1996, Applied Optics 62, 958
    ;
    ; NOTES:
    ;       Take care within 1 A of 2000 A.   Wavelengths below 2000 A *in air* are
    ;       not altered.
    ; REVISION HISTORY
    ;       Written W. Landsman                November 1991
    ;       Use Ciddor (1996) formula for better accuracy in the infrared
    ;           Added optional output vector, W Landsman Mar 2011
    ;       Iterate for better precision W.L./D. Schlegel  Mar 2011
    ;-
    """
    wave_vac = wave_air * 1.0
    g = wave_vac > 2000  # Only modify above 2000 A

    if np.sum(g):
        for iter in [0, 1]:
            if isinstance(g, np.ndarray):
                sigma2 = (1e4 / wave_vac[g]) ** 2.0  # Convert to wavenumber squared
                # Compute conversion factor
                fact = 1.0 + 5.792105e-2 / (238.0185 - sigma2) + 1.67917e-3 / (57.362 - sigma2)
                wave_vac[g] = wave_air[g] * fact  # Convert Wavelength
            else:  # scalar version
                sigma2 = (1e4 / wave_vac) ** 2.0  # Convert to wavenumber squared
                # Compute conversion factor
                fact = 1.0 + 5.792105e-2 / (238.0185 - sigma2) + 1.67917e-3 / (57.362 - sigma2)
                wave_vac = wave_air * fact  # Convert Wavelength

    return wave_vac

