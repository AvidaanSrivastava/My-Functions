## matplotlib.pyplot setting for paper-ready plots

import matplotlib.pyplot as plt
plt.rcParams.update({
    'figure.figsize': (5, 5),
    'figure.dpi': 300,
    'font.size': 10,
    'axes.titlesize': 10,
    'axes.labelsize': 10,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 7
})