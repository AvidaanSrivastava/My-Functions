
## This script contains functions to calculate various astrophysical parameters
## Owner: Avidaan Srivastava

import numpy as np

def p_to_a(period, M_star):

    """ This function calculates the semi-major axis of a planet given its period and the mass of the star.
    The formula used is Kepler's third law of planetary motion.
    period: period of the planet in days
    M_star: mass of the star in solar masses
    a: semi-major axis in meters
    a_AU: semi-major axis in AU
    """

    period = period * 86400   ## From days to seconds
    G = 6.67430e-11
    M_star = M_star * 1.989e30   ## In Kg
    a = (G * M_star * period**2 / (4*np.pi**2))**(1/3)     ## In meters
    a_AU = a / 1.496e11   ## In AU
    return a, a_AU


def calc_Teq(T_star, a, A, R_star, q=2/3):

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
    Tday = T_star * (1 - A)**(1/4) * (x)**(1/2) * (q)**(1/4)
    Teq = T_star * (1 - A)**(1/4) * (x)**(1/2)

    return Teq, Tday


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


def calc_insolation(T_star, R_star, a):

    """ This function calculates the insolation of a planet given its distance from the star,
    the star's temperature, and the star's radius.
    T_star: temperature of the star in Kelvin
    R_star: radius of the star in solar radii
    a: distance from the star in AU
    S: insolation in units of S_earth
    """

    a = a 
    T_sun = 5778

    S = (R_star)**2 * (T_star/T_sun)**4 / a**2

    return S


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

    B_day = planck_func(15*10**(-6), T_day)
    B_star = planck_func(15*10**(-6), T_star)

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


