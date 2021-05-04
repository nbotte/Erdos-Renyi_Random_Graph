import numpy as np
from scipy.stats import binom
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib.colors as colors
import seaborn as sns
from pathlib import Path
import math
import warnings
from typing import Iterable, Optional
from matplotlib.ticker import FormatStrFormatter

def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)


def subplots(*args, figsize=None, **kwargs):
    # figsize in cm instead of inches!
    if figsize is not None:
        return plt.subplots(*args, figsize=cm2inch(figsize), **kwargs)
    else:
        return plt.subplots(*args, **kwargs)


def get_rc_params(usetex=True, font_family='serif', latex_preamble=None, font_serif=None, font_sans_serif=None,
                  rc=None):
    if rc is None:
        rc = dict()
    rc["text.usetex"] = usetex
    rc["font.family"] = font_family
    rc["axes.edgecolor"] = 'darkgrey'
    if latex_preamble is not None:
        rc["text.latex.preamble"] = latex_preamble
    if font_serif is not None:
        rc["font.serif"] = font_serif
    if font_sans_serif is not None:
        rc["font.sans-serif"] = font_sans_serif
    return rc


def get_latex_preamble(use_libertine=True, use_fontenc=True, use_inputenc=True, additional_pkgs=tuple()):
    preamble = ''
    if use_libertine:
        preamble += r'\usepackage{libertine}' \
                    r'\usepackage[libertine]{newtxmath}'

    if use_fontenc:
        preamble += r'\usepackage[T1]{fontenc}'

    if use_inputenc:
        preamble += r'\usepackage[utf8]{inputenc}'

    for package in additional_pkgs:
        if isinstance(package, str):
            preamble += rf"\usepackage{{{package}}}"
        elif isinstance(package, dict):
            option = package["option"]
            name = package["name"]
            preamble += rf"\usepackage[{option}]{{{name}}}"
        else:
            raise ValueError(f"Unsupported type of package: {type(package)}!")
    return preamble


def func(x, a, b):
    return a*x**(-b)

def log_norm(x, mu, sigma):
    return 1/(np.sqrt(2*np.pi)*sigma*x)*np.exp(-((np.log(x)-
   mu)**2)/(2*sigma**2))

def exp(x, a, b):
    return np.exp(-a*x**(-b))

def power_exp(x, b, c):
    return x**(-b)*np.exp(c*x)

deg_distr = np.loadtxt('Degree_distribution_real_network_PGP.txt')

N = 10
p = 0.01

k = np.arange(0, 20, 1)

y = deg_distr[:,1]/np.trapz(deg_distr[:,0], deg_distr[:,1])

popt, pcov = curve_fit(log_norm, deg_distr[:,0][1::], np.abs(y[1::]))
popt1, pcov1 = curve_fit(func, deg_distr[:,0][1::20], np.abs(y[1::20]))
popt4, pcov4 = curve_fit(func, deg_distr[:,0][20::], np.abs(y[20::]))
popt2, pcov2 = curve_fit(exp, deg_distr[:,0][1::], np.abs(y[1::]))
popt3, pcov3 = curve_fit(power_exp, deg_distr[:,0][1::], np.abs(y[1::]))

print(popt1)
print(popt4)

x = np.linspace(1, 15, 500)
x1 = np.linspace(20, 150, 1000)
x2 = np.linspace(1, 150, 1000)

rc_params = get_rc_params(latex_preamble=get_latex_preamble(use_libertine=True))
sns.set_theme(font="DejaVu Serif", rc=rc_params, style="ticks", context="paper")

fig, ax = subplots(figsize=(8, 7))

ax.plot(deg_distr[:,0][1::], np.abs(y)[1::], 'o', markersize='3')
#ax.plot(x2, log_norm(x2, *popt), 'r--')
ax.plot(x2, func(x2, *popt1), 'y--')
ax.plot(x2, func(x2, *popt4), 'g--')
#ax.plot(x, exp(x, *popt2), 'g--')
#ax.plot(x, power_exp(x, *popt3), 'k--')

ax.set_xlabel(r"$k$", fontsize=10)
ax.set_ylabel(r"$P(k)$", fontsize=10)

ax.set_yscale('log')
ax.set_xscale('log')

ax.tick_params(labelsize=10, color='darkgrey')
#legend = ax.legend(loc='upper right')
#legend.get_frame().set_linewidth(0.0)
#ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.tight_layout()
plt.savefig('PGP_degree_distribution_8x7.png', dpi=500)
plt.show()
