
# See LV-Monovalue/Debug3 to add other numerical/analytical solvers

import time

import matplotlib.colorbar
import numpy as np
from tqdm import tqdm
from Browser import *


# TESTING

def plotDiagram(ax, res, l1, l2, lbl1, lbl2, plotType="Phase"):
    cb = np.empty((3, 4), dtype=matplotlib.colorbar.Colorbar)
    N_xticks = int(len(l1)/3)
    N_yticks = int(len(l2)/10)
    for i in range(3):
        for j in range(4):
            if plotType == "Phase":
                plt.subplots_adjust(wspace=1.)
                sh = ax[i, j].imshow(res[:, :, 4*i+j], aspect=len(l1)/len(l2))
                cb[i, j] = plt.colorbar(mappable=sh, ax=ax[i, j])
                ax[i, j].set_xticks(range(0, len(l1), N_xticks), ["{:10.0f}".format(s1) for s1 in l1[::N_xticks]], ha="right")
                ax[i, j].set_yticks(range(0, len(l2), N_yticks), ["{:10.1f}".format(s2) for s2 in l2[::N_yticks]])
            elif plotType == "Scatter":
                plt.subplots_adjust(hspace=0.5)
                ax[i, j].scatter(l1, res[:, 4*i+j], marker="+", alpha=0.1)
            ax[i, j].set_xlabel(lbl1)
            ax[i, j].set_ylabel(lbl2, labelpad=-20)
    #cb[0, 0].mappable.set_clim(0., 1.)
    ax[0, 0].set_title("p_surviving")
    ax[0, 1].set_title("(X/K)_mean")
    ax[0, 2].set_title("(X/K)_std")
    ax[0, 3].set_title("(X/K)_max")
    ax[1, 0].set_title("R_mean_X")
    ax[1, 1].set_title("K_mean_X")
    ax[1, 2].set_title("A_mean_X")
    ax[1, 3].set_title("C_mean_X")
    ax[2, 0].set_title("R_std_X")
    ax[2, 1].set_title("K_std_X")
    ax[2, 2].set_title("A_std_X")
    ax[2, 3].set_title("C_std_X")
    pass


def makePhaseDiagram(Niter, gs, ss, f):
    tic = time.perf_counter()
    res = np.zeros((len(ss), len(gs), 12))
    for _ in tqdm(range(Niter)):
        for iis, s in enumerate(ss):
            E = makeEasy_E(np.array([2., 2., 0., 2.]), np.array([0.5, 0.5, 1., 0.5]), [(s, f)])
            sols = solveGlobal(E, optPrint=False)
            for iig, g in enumerate(gs):
                sol, P = findSol_g(sols, g)
                res[iis, iig, :] += macroStats(E, P, sol, g)
    fig, ax = plt.subplots(3, 4)
    plotDiagram(ax, res/Niter, np.log10(gs), ss, "Log10(Gamma)", "Sigma")
    print(time.perf_counter()-tic)
    plt.show()


def makePhaseDiagram2(Niter, s1s, s2s, f1, f2, g):
    tic = time.perf_counter()
    res = np.zeros((len(s1s), len(s2s), 12))
    for _ in tqdm(range(Niter)):
        for is1, s1 in enumerate(s1s):
            for is2, s2 in enumerate(s2s):
                E = makeEasy_E(np.array([2., 2., 0., 2.]), np.array([0.5, 0.5, 1., 0.5]), [(s1, f1), (s2, f2)])
                sols = solveGlobal(E, optPrint=False)
                sol, P = findSol_g(sols, g)
                res[is1, is2, :] += macroStats(E, P, sol, g)
    resPrint = np.transpose(res / Niter, (1, 0, 2))
    fig, ax = plt.subplots(3, 4)
    plotDiagram(ax, resPrint, s1s, s2s, "Sigma_1", "Sigma_2")
    print(time.perf_counter()-tic)
    plt.show()


def makeIndepDiagram(Niter, gs):
    tic = time.perf_counter()
    fig, ax = plt.subplots(3, 4)
    for _ in tqdm(range(Niter)):
        E = makeEasy_E(np.array([2., 2., 0., 2.]), np.array([0.5, 0.5, 1., 0.5]), [])
        sols = solveGlobal(E, optPrint=False)
        resPrint = np.zeros((len(gs), 12))
        for ig, g in enumerate(gs):
            sol, P = findSol_g(sols, g)
            resPrint[ig] = macroStats(E, P, sol, g)
        plotDiagram(ax, resPrint, gs, [], "Gamma", "", plotType="Scatter")
    print(time.perf_counter()-tic)
    plt.show()


#makePhaseDiagram(1000, np.logspace(-2, 4., 101), np.linspace(-1., 1., 101), ("a", "c"))  # (100, 101, 101) -> 12mn
#makePhaseDiagram2(10, np.linspace(-1., 1., 21), np.linspace(-1., 1., 21), ("r", "K"), ("a", "c"), 5.)
#makeIndepDiagram(20, np.linspace(0., 1000., 101))


"""
tic = time.perf_counter()
#for i in range(1000):
i=14
print(i)
aa = np.random.uniform(-1., 1., 100)
np.random.seed(i)
E = makeEasy_E(np.array([2., 2., 0., 2.]), np.array([0.5, 0.5, 1., 0.5]), [])
sols = solveGlobal(E, True)
print(time.perf_counter()-tic)
"""

E = makeEasy_E(np.array([2., 2., 0., 2.]), np.array([0.5, 0.5, 1., 0.5]), [])
sols = solveGlobal(E, optPrint=False)
reconstructGlobal(E, sols, plotAnaSol)
plt.plot([-3.2, np.log10(2*sols[-1].g_i)], [0., 0.], c="black")
plt.xlabel("Log10(gamma)", fontsize=18)
plt.ylabel("X/K", fontsize=18)
plt.show()