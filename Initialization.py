import numpy as np


# GLOBAL CONSTANTS
np.set_printoptions(linewidth=100000)
np.random.seed(0)
ee = 1E-6


# INITIALIZATION
N_g = 50
N = 100
pick = np.random.randint(N, size=(20,))


def makeEasy_E(m, std, sigma_l):
    vv = np.diag(std**2)
    dic = {"r": 0, "K": 1, "a": 2, "c": 3}
    for s, p in sigma_l:
        i, j = dic[p[0]], dic[p[1]]
        vv[i, j] = s * std[i] * std[j]
        vv[j, i] = vv[i, j]
    E = np.random.multivariate_normal(mean=m, cov=vv, size=N)
    E[:, 0], E[:, 1], E[:, 3] = np.maximum(E[:, 0], 0.01), np.maximum(E[:, 1], 0.01), np.maximum(E[:, 3], 0.01)
    return E


def unravel(E):
    return len(E), E[:, 0], E[:, 1], E[:, 2], E[:, 3]