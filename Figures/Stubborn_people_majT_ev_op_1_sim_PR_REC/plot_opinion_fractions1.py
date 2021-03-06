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

fractions1 = np.loadtxt('Fraction_of_opinions_active_01_av_good_init_SBM_PR_01-0001_10x100_majT=0_one_sim.txt')
fractions2 = np.loadtxt('Fraction_of_opinions_active_01_av_good_init_SBM_PR_003-0008_10x100_majT=0_one_sim.txt')
fractions3 = np.loadtxt('Fraction_of_opinions_active_01_av_good_init_SBM_REC_01-0001_10x100_majT=0_one_sim.txt')
fractions4 = np.loadtxt('Fraction_of_opinions_active_01_av_good_init_SBM_REC_003-0008_10x100_majT=0_one_sim.txt')
#fractions5 = np.loadtxt('Fraction_of_opinions_active_01_av_good_init_SBM_PR_01-0001_10x100_fracRes=04_stubb=1_20-80.txt')
#fractions6 = np.loadtxt('Fraction_of_opinions_active_01_av_good_init_SBM_PR_01-0001_10x100_fracRes=08_stubb=1_20-80.txt')

'''fractions7 = np.loadtxt('Fraction_of_opinion0_comm_SBM_active_01_av_good_init_commOp0=01_other=50-50_PR_01-0001_10x100_T=0_random=55-45.txt')
fractions8 = np.loadtxt('Fraction_of_opinion0_comm_SBM_active_01_av_good_init_commOp0=02_other=50-50_PR_01-0001_10x100_T=0_random=60-40.txt')
fractions9 = np.loadtxt('Fraction_of_opinion0_comm_SBM_active_01_av_good_init_commOp0=03_other=50-50_PR_01-0001_10x100_T=0_random=65-35.txt')
fractions10 = np.loadtxt('Fraction_of_opinion0_comm_SBM_active_01_av_good_init_commOp0=04_other=50-50_PR_01-0001_10x100_T=0_random=70-30.txt')
fractions11 = np.loadtxt('Fraction_of_opinion0_comm_SBM_active_01_av_good_init_commOp0=05_other=50-50_PR_01-0001_10x100_T=0_random=75-25.txt')
fractions12 = np.loadtxt('Fraction_of_opinion0_comm_SBM_active_01_av_good_init_commOp0=06_other=50-50_PR_01-0001_10x100_T=0_random=80-20.txt')

fractions13 = np.loadtxt('Fraction_of_opinion0_comm_SBM_active_01_av_good_init_commOp0=01_other=50-50_PR_003-0008_10x100_T=0.txt')
fractions14 = np.loadtxt('Fraction_of_opinion0_comm_SBM_active_01_av_good_init_commOp0=02_other=50-50_PR_003-0008_10x100_T=0.txt')
fractions15 = np.loadtxt('Fraction_of_opinion0_comm_SBM_active_01_av_good_init_commOp0=03_other=50-50_PR_003-0008_10x100_T=0.txt')
fractions16 = np.loadtxt('Fraction_of_opinion0_comm_SBM_active_01_av_good_init_commOp0=04_other=50-50_PR_003-0008_10x100_T=0.txt')
fractions17 = np.loadtxt('Fraction_of_opinion0_comm_SBM_active_01_av_good_init_commOp0=05_other=50-50_PR_003-0008_10x100_T=0.txt')
fractions18 = np.loadtxt('Fraction_of_opinion0_comm_SBM_active_01_av_good_init_commOp0=06_other=50-50_PR_003-0008_10x100_T=0.txt')

fractions19 = np.loadtxt('Fraction_of_opinion0_comm_SBM_active_01_av_good_init_commOp0=01_other=50-50_PR_003-0008_10x100_T=0_random=55-45.txt')
fractions20 = np.loadtxt('Fraction_of_opinion0_comm_SBM_active_01_av_good_init_commOp0=02_other=50-50_PR_003-0008_10x100_T=0_random=60-40.txt')
fractions21 = np.loadtxt('Fraction_of_opinion0_comm_SBM_active_01_av_good_init_commOp0=03_other=50-50_PR_003-0008_10x100_T=0_random=65-35.txt')
fractions22 = np.loadtxt('Fraction_of_opinion0_comm_SBM_active_01_av_good_init_commOp0=04_other=50-50_PR_003-0008_10x100_T=0_random=70-30.txt')
fractions23 = np.loadtxt('Fraction_of_opinion0_comm_SBM_active_01_av_good_init_commOp0=05_other=50-50_PR_003-0008_10x100_T=0_random=75-25.txt')
fractions24 = np.loadtxt('Fraction_of_opinion0_comm_SBM_active_01_av_good_init_commOp0=06_other=50-50_PR_003-0008_10x100_T=0_random=80-20.txt')'''

t = np.zeros(500)
y = np.zeros(500)

for i in range(len(t)):
    t[i] = i
    y[i] = 0.5

'''x = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]

begin = np.zeros(6)
begin[0] = np.mean([np.mean(fractions1[:,0]), np.mean(fractions7[:,0]), np.mean(fractions13[:,0]), np.mean(fractions19[:,0])])
begin[1] = np.mean([np.mean(fractions2[:,0]), np.mean(fractions8[:,0]), np.mean(fractions14[:,0]), np.mean(fractions20[:,0])])
begin[2] = np.mean([np.mean(fractions3[:,0]), np.mean(fractions9[:,0]), np.mean(fractions15[:,0]), np.mean(fractions21[:,0])])
begin[3] = np.mean([np.mean(fractions4[:,0]), np.mean(fractions10[:,0]), np.mean(fractions16[:,0]), np.mean(fractions22[:,0])])
begin[4] = np.mean([np.mean(fractions5[:,0]), np.mean(fractions11[:,0]), np.mean(fractions17[:,0]), np.mean(fractions23[:,0])])
begin[5] = np.mean([np.mean(fractions6[:,0]), np.mean(fractions12[:,0]), np.mean(fractions18[:,0]), np.mean(fractions24[:,0])])

high_mod_comm = np.zeros(6)
high_mod_comm[0] = np.mean(fractions1[:,1])
high_mod_comm[1] = np.mean(fractions2[:,1])
high_mod_comm[2] = np.mean(fractions3[:,1])
high_mod_comm[3] = np.mean(fractions4[:,1])
high_mod_comm[4] = np.mean(fractions5[:,1])
high_mod_comm[5] = np.mean(fractions6[:,1])

high_mod_rand = np.zeros(6)
high_mod_rand[0] = np.mean(fractions7[:,1])
high_mod_rand[1] = np.mean(fractions8[:,1])
high_mod_rand[2] = np.mean(fractions9[:,1])
high_mod_rand[3] = np.mean(fractions10[:,1])
high_mod_rand[4] = np.mean(fractions11[:,1])
high_mod_rand[5] = np.mean(fractions12[:,1])

low_mod_comm = np.zeros(6)
low_mod_comm[0] = np.mean(fractions13[:,1])
low_mod_comm[1] = np.mean(fractions14[:,1])
low_mod_comm[2] = np.mean(fractions15[:,1])
low_mod_comm[3] = np.mean(fractions16[:,1])
low_mod_comm[4] = np.mean(fractions17[:,1])
low_mod_comm[5] = np.mean(fractions18[:,1])

low_mod_rand = np.zeros(6)
low_mod_rand[0] = np.mean(fractions19[:,1])
low_mod_rand[1] = np.mean(fractions20[:,1])
low_mod_rand[2] = np.mean(fractions21[:,1])
low_mod_rand[3] = np.mean(fractions22[:,1])
low_mod_rand[4] = np.mean(fractions23[:,1])
low_mod_rand[5] = np.mean(fractions24[:,1])'''


rc_params = get_rc_params(latex_preamble=get_latex_preamble(use_libertine=True))
sns.set_theme(font="DejaVu Serif", rc=rc_params, style="ticks", context="paper")

fig, ax = subplots(figsize=(8, 7))

ax.errorbar(t[::50], fractions1[:,0][::50], np.sqrt(fractions1[:,2][::50]), label=r'High $Q$')
ax.errorbar(t[::50], fractions3[:,1][::50], np.sqrt(fractions3[:,2][::50]), label=r'High $Q$')
ax.errorbar(t[::50], fractions2[:,1][::50], np.sqrt(fractions2[:,2][::50]), label=r'Low $Q$')
#ax.errorbar(t[::50], fractions4[:,0][::50], np.sqrt(fractions4[:,2][::50]), label=r'Low $Q$, REC')

'''ax.plot(x, begin, 'k', label=r'Average fraction opinion 0 at $t=0$')
ax.plot(x, high_mod_rand, 'cornflowerblue', linestyle='--', label=r'High mod, random distributed')
ax.plot(x, high_mod_comm, 'peru', linestyle='--', label=r'High mod, comm. distributed')
ax.plot(x, low_mod_rand, 'forestgreen', linestyle='--', label=r'Low mod, random distributed')
ax.plot(x, low_mod_comm, 'firebrick', linestyle='--', label=r'Low mod, comm. distributed')'''

ax.errorbar(t[::50], fractions4[:,1][::50], np.sqrt(fractions1[:,2][::50]), label=r'Low $Q$')
#ax.errorbar(t[::50], fractions5[:,0][::50], np.sqrt(fractions2[:,2][::50]), label=r'High $Q$')
#ax.errorbar(t[::50], fractions6[:,0][::50], np.sqrt(fractions3[:,2][::50]), label=r'High $Q$')

#ax.errorbar(t[::50], fractions7[:,0][::50], np.sqrt(fractions1[:,2][::50]), label=r'fracRes $= 0.8$')
#ax.errorbar(t[::50], fractions8[:,0][::50], np.sqrt(fractions2[:,2][::50]), label=r'REF, fracRes $= 0.8$')
#ax.errorbar(t[::50], fractions9[:,0][::50], np.sqrt(fractions3[:,2][::50]), label=r'fracRes $= 0.8$')
#plt.plot(t[::10], fractions_sameComm02[:,0][::10], label='0.4 community opinion 0 (same comm.), other 50/50; REC')
ax.plot(t[::10], y[::10], 'k--')
#plt.plot(t, fractions_25_av[:,0], label='Stubborness = 0.25')
#plt.plot(t, fractions_1_av[:,1], label='Stubborness = 1')
ax.set_xlabel(r"$t$", fontsize=10)
ax.set_ylabel(r"$P_\text{A}(t)$", fontsize=10)
ax.tick_params(labelsize=10, color='darkgrey')
h, l = ax.get_legend_handles_labels()
ph = [plt.plot([], marker="", ls="")[0]]*2
handles = ph + h
labels = [r"PR:", r"REC:"] + l
legend = ax.legend(handles, labels, ncol=3, handletextpad=0.6, borderpad=0.3, handlelength=1.5)
legend.get_frame().set_linewidth(0.0)
ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
#ax.set_ylim(-0.05, 0.4)
ax.set_ylim(0.47, 1.2)
#ax.set_title('WS', fontsize=10)
plt.tight_layout()
plt.savefig("fraction_of_opinions_SBM_majT=0_one_sim_8x7.png", dpi=500)
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
