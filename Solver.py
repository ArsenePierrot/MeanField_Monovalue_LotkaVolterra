import numpy as np

from Initialization import *


# ANALYTICAL TRINOME SOLVE
def trinomeG(E, P):
    n, r, K, a, c = unravel(E)
    sumE = lambda x: np.sum(K[P] * c[P] ** 2 / r[P] * x[P]) / n
    v1 = np.full(n, 1.)
    A0 = K
    A1 = K*c/r * (sumE(r/c) * a - sumE(a*r/c))
    A2 = K * (sumE(a**2) * sumE(v1) - sumE(a)**2) + \
         K*c*a/r * (sumE(a) * sumE(r/c) - sumE(v1) * sumE(a*r/c)) - \
         K*c/r * (sumE(a**2) * sumE(r/c) - sumE(a) * sumE(a*r/c))
    D2 = sumE(a ** 2) * sumE(v1) - sumE(a) ** 2
    return A2, A1, A0, D2


def solveTrinome(a2, a1, a0, g_0):
    x1, x2 = np.zeros_like(a2), np.zeros_like(a2)
    S = (np.abs(a2)/(a1**2) > 1E-10)
    x1[~S], x2[~S] = -a0[~S]/a1[~S], np.inf
    sr_d = np.sqrt(a1**2 - 4 * a2*a0)
    x1[S], x2[S] = (-a1[S] - sr_d[S]) / (2 * a2[S]), (-a1[S] + sr_d[S]) / (2 * a2[S])
    x1, x2 = np.minimum(x1, x2), np.maximum(x1, x2)

    res = np.full(len(a2), np.inf)
    S1, S2 = (x2 > g_0) & (x1 <= g_0), (x2 > g_0) & (x1 > g_0)  # both False if sr_d = nan
    res[S1], res[S2] = x2[S1], x1[S2]
    res[res > 1/ee] = np.inf

    # Final Check
    Sf = (res < np.full(len(a2), np.inf))
    v = a2[Sf] * res[Sf]**2 + a1[Sf] * res[Sf] + a0[Sf]
    if np.any(np.abs(v) > ee):
        print("WARNING! solveTrinome did a poor job finding roots.",
              "Err = " + str(np.max(np.abs(v))),
              "g_f = " + str(np.min(res)), sep="\n")
    return res


# GLOBAL SOLVING
class localSol:
    def __init__(self, g_i, g_f, idx):
        self.g_i = g_i
        self.g_f = g_f
        self.idx = idx


def next_gInterval(E, g_i, P_i):
    A2, A1, A0, D2 = trinomeG(E, P_i)
    g_f_candidates = solveTrinome(A2, A1, A0, g_i)
    g_f, idx = np.nanmin(g_f_candidates), np.nanargmin(g_f_candidates)
    return g_f, idx


def solveGlobal(E, optPrint):
    g_i, P_i, i_last = 0., np.full(len(E), True), None
    sols = []
    while g_i < np.inf:
        if optPrint:
            print("#"*30, "N = " + np.count_nonzero(P_i).__str__(), "g_i = " + g_i.__str__(), "i_last = " + i_last.__str__(), sep="\n")
        g_f, idx = next_gInterval(E, g_i, P_i)
        sols.append(localSol(g_i, g_f, idx))
        if not((i_last == idx) and (g_f < g_i + ee)):
            P_i[idx] = not(P_i[idx])
            if P_i[idx] and g_f != np.inf:
                print("WARNING! Resurrection of " + str(idx) + " at t=" + str(g_f))
        g_i, i_last = g_f, idx
    return sols