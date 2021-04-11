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

fractions1 = np.loadtxt('Fraction_of_opinions_active_01_av_good_init_WS_PR_10000-6-001_fracRes=0_stubb=0_20-80.txt')
#fractions2 = np.loadtxt('Fraction_of_opinions_active_01_av_good_init_SBM_PR_003-0008_10x100_fracRes=0_stubb=0_20-80.txt')
#fractions3 = np.loadtxt('Fraction_of_opinions_active_01_av_good_init_SBM_REC_01-0001_10x100_fracRes=0_stubb=0_20-80.txt')
#fractions4 = np.loadtxt('Fraction_of_opinions_active_01_av_good_init_SBM_REC_003-0008_10x100_fracRes=0_stubb=0_20-80.txt')

#fractions_sameComm02 = np.loadtxt('Fraction_of_opinions_SBM_active_01_av_good_init_commOp0=01_other=44-56_PR_01-0001_10x100_T=0_sameComm.txt')
#fractions_PR = np.loadtxt('Fraction_of_opinions_Clustered_Cluster5_01-0001_50_50_no_stubb_paper8_active_01_av_good_init_PR.txt')
#fractions_50_av = np.loadtxt('Fraction_of_opinions_1_20_80_50_stubb_50_bern_050_av.txt')
#fractions_25_av = np.loadtxt('Fraction_of_opinions_1_20_80_50_stubb_25_bern_050_av.txt')
#fractions_1_av = np.loadtxt('Fraction_of_opinions_1_70_30_all_stubb_100_bern_5_av.txt')

t = np.zeros(500)
y = np.zeros(500)

for i in range(len(t)):
    t[i] = i
    y[i] = 0.2

rc_params = get_rc_params(latex_preamble=get_latex_preamble(use_libertine=True))
sns.set_theme(font="DejaVu Serif", rc=rc_params, style="whitegrid", context="paper")

fig, ax = subplots(figsize=(8, 7))

ax.errorbar(t[::50], fractions1[:,0][::50], np.sqrt(fractions1[:,2][::50]), label=r'WS, $N=10000,\ K = 6$')
#ax.errorbar(t[::50], fractions2[:,0][::50], np.sqrt(fractions2[:,2][::50]), label=r'Low mod, PR')
#ax.errorbar(t[::50], fractions3[:,0][::50], np.sqrt(fractions3[:,2][::50]), label=r'High mod, REC')
#ax.errorbar(t[::50], fractions4[:,0][::50], np.sqrt(fractions4[:,2][::50]), label=r'Low mod, REC')
#plt.plot(t[::10], fractions_sameComm02[:,0][::10], label='0.4 community opinion 0 (same comm.), other 50/50; REC')
ax.plot(t[::10], y[::10], 'k--')
#plt.plot(t, fractions_25_av[:,0], label='Stubborness = 0.25')
#plt.plot(t, fractions_1_av[:,1], label='Stubborness = 1')
ax.set_xlabel(r"Time $t$", fontsize=10)
ax.set_ylabel(r"Fraction of opinion 0 ($P_0(t)$)", fontsize=10)
ax.tick_params(labelsize=10)
legend = ax.legend(loc='upper right')
legend.get_frame().set_linewidth(0.0)
ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.set_ylim(-0.05, 0.5)
#ax.set_title('WS', fontsize=10)
plt.tight_layout()
plt.savefig("fraction_of_opinions_WS_10000-6_PR_20-80_8x7.png", dpi=500)
plt.show()

"""fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 4))
st = fig.suptitle("Stochastic block model (100 x 10), 10 x 10 averaged")
#axes[0].plot(t, y)
#axes[0].plot(t[::10], fractions_av0[:,0][::10], label='Opinion 0')
axes[0].plot(t[::10], fractions_REC[:,1][::10], label='Opinion 1, REC')
axes[0].legend(loc='best')
axes[0].set_xlabel("Timesteps t")
axes[0].set_ylabel("Opinion fraction inside a community, REC")
axes[0].set_ylim(0.39, 0.55)
axes[0].title.set_text("Opinion fraction inside a community vs time, REC, 50/50\n"r"$p_{cl}$ = 0.1, $p_{add}$ = 0.001")

axes[1].plot(t[::10], fractions_PR[:,1][::10], label='Opinion 1, PR')
#axes[1].plot(t[::10], fractions_av3[:,1][::10], label='Standard deviation')
#axes[1].plot(t, y)
axes[1].legend(loc='best')
axes[1].set_xlabel("Timesteps t")
axes[1].set_ylabel("Opinion fraction inside a community, PR")
axes[1].set_ylim(0.40, 0.45)
axes[1].title.set_text("Opinion fraction inside a community vs time, PR, 50/50\n"r"$p_{cl}$ = 0.1, $p_{add}$ = 0.001")
# shift subplots down:
st.set_y(0.03)
fig.subplots_adjust(bottom=0.1)

fig.tight_layout()
fig.savefig('Fraction_of_opinions_Clustered_Cluster5_01-0001_50_50_no_stubb_paper8_active_01_av_good_init_100x10_REC_PR.png')
plt.show()"""
