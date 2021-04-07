import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import seaborn as sns
from pathlib import Path
import math
import warnings
from typing import Iterable, Optional
from matplotlib.ticker import FormatStrFormatter
from scipy.optimize import curve_fit

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

def func(x, a, b, c):
    return a*np.exp(b*x) + c

def func_trip(x, a, b, c, d):
    return a*x**3 + b*x**2 + c*x + d

echo_chambers1 = np.loadtxt("Echo_chamber_SBM_PR_01_0001_10x100_T=0.txt")
echo_chambers2 = np.loadtxt("Echo_chamber_SBM_PR_007-0004_10x100_all_res_stubb=0.txt")
echo_chambers3 = np.loadtxt("Echo_chamber_SBM_PR_005-0005_10x100_all_res_stubb=0.txt")
echo_chambers4 = np.loadtxt("Echo_chamber_SBM_PR_003-0008_10x100_all_res_stubb=0.txt")
echo_chambers5 = np.loadtxt("Echo_chamber_ER_PR_001_all_res_stubb=0.txt")

echo_chambers11 = np.loadtxt("Echo_chamber_SBM_REC_01-0001_10x100_probT=0.txt")
echo_chambers21 = np.loadtxt("Echo_chamber_SBM_REC_007-0004_10x100_fracRes=0_stubb=0.txt")
echo_chambers31 = np.loadtxt("Echo_chamber_SBM_REC_005-0005_10x100_fracRes=0_stubb=0.txt")
echo_chambers41 = np.loadtxt("Echo_chamber_SBM_REC_003-0008_10x100_all_res_stubb=0.txt")
echo_chambers51 = np.loadtxt("Echo_chamber_real_network_ER_REC_001_fracRes=0_stubb=0.txt")

echo = np.zeros(5)

echo[0] = np.mean(echo_chambers1[:,2])
echo[1] = np.mean(echo_chambers2[:,2])
echo[2] = np.mean(echo_chambers3[:,2])
echo[3] = np.mean(echo_chambers4[:,2])
echo[4] = np.mean(echo_chambers5[:,2])

echoREC = np.zeros(5)

echoREC[0] = np.mean(echo_chambers11[:,2])
echoREC[1] = np.mean(echo_chambers21[:,2])
echoREC[2] = np.mean(echo_chambers31[:,2])
echoREC[3] = np.mean(echo_chambers41[:,2])
echoREC[4] = np.mean(echo_chambers51[:,2])


err = np.zeros(5)
err[0] = 11.14
err[1] = 3.79
err[2] = 1.46
err[3] = 2.73
err[4] = 1.69

errREC = np.zeros(5)
errREC[0] = 1.60
errREC[1] = 1.93
errREC[2] = 1.04
errREC[3] = 1.58
errREC[4] = 1.44

'''echo1 = np.zeros(5)

echo1[0] = echo_chambers11[:,2][0]
echo1[1] = echo_chambers21[:,2][0]
echo1[2] = echo_chambers31[:,2][0]
echo1[3] = echo_chambers41[:,2][0]
echo1[4] = echo_chambers51[:,2][0]'''

mod = np.zeros(5)

mod[0] = 0.910383
mod[1] = 0.6
mod[2] = 0.47
mod[3] = 0.2
mod[4] = 0.

'''mod1 = np.zeros(5)

mod1[0] = 0.910383
mod1[1] = 0.830371
mod1[2] = 0.651065
mod1[3] = 0.469297
mod1[4] = 0.286703'''

popt, pcov = curve_fit(func_trip, mod, echo)
print(popt)
popt1, pcov1 = curve_fit(func_trip, mod, echoREC)
print(popt1)
xf = np.linspace(0,0.92,50)

rc_params = get_rc_params(latex_preamble=get_latex_preamble(use_libertine=True))
sns.set_theme(font="DejaVu Serif", rc=rc_params, style="whitegrid", context="paper")

fig, ax = subplots(figsize=(12, 10))

ax.errorbar(mod, echo, err/10, xerr=None, c='b', ls='', marker='o', label = r'Data points, PR')
ax.errorbar(mod, echoREC, errREC/10, xerr=None, c='r', ls='', marker='o', label = r'Data points, REC')

ax.plot(xf, func_trip(xf, *popt), 'b--')
ax.plot(xf, func_trip(xf, *popt1), 'r--')


#ax.plot(t[::10], y[::10], 'k--')
#plt.plot(t, fractions_25_av[:,0], label='Stubborness = 0.25')
#plt.plot(t, fractions_1_av[:,1], label='Stubborness = 1')
ax.set_xlabel(r"Modularity $Q$", fontsize=10)
ax.set_ylabel(r"Echo chamber size", fontsize=10)
#ax.set_yscale('log')
ax.tick_params(labelsize=10)
legend = ax.legend(loc='best')
legend.get_frame().set_linewidth(0.0)
#ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
#ax.set_ylim(-5, 200)
#ax.set_title('WS', fontsize=10)
plt.tight_layout()
plt.savefig("echo_chamber_vs_modularity_edge_prob_same_deg_PR_REC_12x10.png", dpi=500)
plt.show()
