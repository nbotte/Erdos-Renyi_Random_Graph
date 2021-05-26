import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
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

hist = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_real_network_lastfm_PR_8003_fracRes=0_stubb=0.txt")
hist1 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM-WS_PR_8003-4-0025-00001-powerlaw_fracRes=0_stubb=0_av10x5.txt")
hist2 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM-WS_PR_8003-4-0025-00001_53x151_fracRes=0_stubb=0_av10x5.txt")
#hist3 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_SBM_commOp0=01_other=50-50_PR_003-0008_10x100_T=0_random=55-45.txt")
#hist4 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_WS_PR_10-006_T=0.txt")
#hist5 = np.loadtxt("Hist_500_and_0_fraction_friends_opinion1_WS_PR_10-001_T=0_test.txt")
# probably need to change size of t and y back to 10!
t = np.zeros(10)
t1 = np.zeros(11)
y = np.zeros(10)

for i in range(10):
    t[i] = i/10
    y[i] = 1

for i in range(11):
    t1[i] = i/11


rc_params = get_rc_params(latex_preamble=get_latex_preamble(use_libertine=True))
sns.set_theme(font="DejaVu Serif", rc=rc_params, style="ticks", context="paper")

fig, ax = subplots(figsize=(8, 7))

ax.plot(t, hist[:,2], linestyle='solid', label=r"Last.fm")
ax.plot(t, hist1[:,2], linestyle='solid', label=r"SBM-WS1")
ax.plot(t, hist2[:,2], linestyle='solid', label=r"SBM-WS2")
#ax.plot(t, hist2[:,2], linestyle='solid', label=r"Comm.")
#ax.plot(t1, hist4[:,2], linestyle='solid', label=r"WS, $\left<cc\right> = 0.557 \pm 0.007$")
#ax.plot(t, hist5[:,2], linestyle='solid', label=r"WS, $\left<cc\right> = 0.648 \pm 0.003$")
ax.plot(t, y, 'k--')


#ax.plot(t[::10], y[::10], 'k--')
#plt.plot(t, fractions_25_av[:,0], label='Stubborness = 0.25')
#plt.plot(t, fractions_1_av[:,1], label='Stubborness = 1')
ax.set_xlabel(r"$\left<P_{\text{B}}^{\text{nn}}\right>$", fontsize=10)
ax.set_ylabel(r"$F_N(\left<P_{\text{B}}^{\text{nn}}\right>)$", fontsize=10)
#ax.set_yscale('log')
ax.tick_params(labelsize=10, color='darkgrey')
legend = ax.legend(loc="upper right")
legend.get_frame().set_linewidth(0.0)
#ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.set_ylim(0.2, 6.5)
#ax.set_title('WS', fontsize=10)
plt.tight_layout()
plt.savefig("echo_chambers_SBM_powerlaw_vs_regular_PR.png", dpi=500)
plt.show()
