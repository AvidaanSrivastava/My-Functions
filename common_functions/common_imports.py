# ========================
# Astropy and astronomy tools
# ========================
from astropy.io import fits, ascii
from astropy.table import Table
from astropy.time import Time
from astropy.timeseries import LombScargle
from astropy import constants as const
from astropy import units as u

# ========================
# Plotting
# ========================
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.gridspec import GridSpec
import matplotlib.gridspec as gridspec  # alias if needed
import seaborn as sns

# ========================
# Data & numerical packages
# ========================
import numpy as np
import pandas as pd
#from numpy import *  # Not recommended, but included if you rely on it
import pylab
import random
import time
import math
import cmath
import pickle
import shutil
import sys
import os
import glob
import string
import pprint
from copy import deepcopy

# ========================
# Scipy tools
# ========================
from scipy import interpolate, optimize
from scipy.linalg import solve, solve_banded
from scipy.optimize import curve_fit, least_squares, leastsq, minimize
from scipy.special import erf
from scipy.stats import median_abs_deviation as mad, norm
from scipy.signal import argrelextrema

# ========================
# Modeling & fitting
# ========================
from lmfit import Model
from lmfit.models import GaussianModel, ConstantModel

# ========================
# Astronomy-specific tools
# ========================
import radvel
import radvel.likelihood
from radvel.plot import orbit_plots, mcmc_plots
from dace_query.spectroscopy import Spectroscopy
import lightkurve as lk
import juliet
import emcee
import arviz as az
import corner
import george
from george import kernels
from george.gp import LinAlgError
from PyAstronomy.pyasl import foldAt

# ========================
# Progress bar
# ========================
from tqdm import tqdm

# ========================
# Collections
# ========================
from sortedcollections import SortedList

# ========================
# Custom tools
# ========================
from common_functions import *

# ========================
# Utility functions
# ========================
def mjd_to_date_formatter(x, pos):
    return Time(x, format='mjd').strftime('%Y-%m')
