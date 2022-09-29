from Solver import *
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider


# GLOBAL ANALYSIS
def reconstructGlobal(E, sols, fun):
    P_c, i_last = np.full(len(E), True), None
    for sol in sols:
        fun(E, P_c, sol)
        if not((i_last == sol.idx) and (sol.g_f < sol.g_i + ee)):
            P_c[sol.idx] = not(P_c[sol.idx])
        i_last = sol.idx
    pass


def findSol_g(sols, g):
    P_c, i, i_last = np.full(N, True), 0, None
    while sols[i].g_f < g:
        if not((i_last == sols[i].idx) and (sols[i].g_f < sols[i].g_i + ee)):
            P_c[sols[i].idx] = not(P_c[sols[i].idx])
        i_last = sols[i].idx
        i += 1
    return sols[i], P_c


def unravelSol(E, P, sol, clip_gf=True):
    A2, A1, A0, D2 = trinomeG(E, P)
    g_i, g_f = sol.g_i, sol.g_f
    if clip_gf and (g_f == np.inf):
        g_f = 2 * g_i
    x = lambda g: (A2 * g ** 2 + A1 * g + A0 * 1.) / (D2 * g ** 2 + 1.)
    return A2, A1, A0, D2, g_i, g_f, x


# BROWSING ACTIONS

def classicPrint(E, P, sol):
    n, r, K, a, c = unravel(E)
    A2, A1, A0, D2, gi, gf, x = unravelSol(E, P, sol, clip_gf=False)
    print("########################", sol,
          "g_i = " + str(gi), "g_f = " + str(gf),
          "P_c = " + str(P), "N_c = " + str(np.count_nonzero(P)),
          "x_i = " + str(x(gi)), "x_f = " + str(x(gf)),
          "A_2 = " + str(A2), sep="\n")


def plotAnaSol(E, P, sol):
    n, r, K, a, c = unravel(E)
    A2, A1, A0, D2, gi, gf, x = unravelSol(E, P, sol, clip_gf=True)
    gs = np.array([np.linspace(gi, gf, N_g)]).T
    plt.plot(np.log10(gs), (x(gs)/K), color=None, alpha=0.2)
    #plt.plot(gs, (x(gs)/K)[:, pick], color="red", alpha=0.2)
    plt.plot(np.log10(gs[-1]), 0., "kx")
    pass

def mean_std_w(values, weights):
    average = np.average(values, weights=weights)
    variance = np.average((values-average)**2, weights=weights)
    return average, np.sqrt(variance)


def macroStats(E, P, sol, g):
    n, r, K, a, c = unravel(E)
    A2, A1, A0, D2, gi, gf, x = unravelSol(E, P, sol)
    p, xw = np.mean((x(g) > 0.)), np.maximum(x(g), 0.)
    xSKm, xSKs = mean_std_w(xw/K, xw)
    rm, rs = mean_std_w(r, xw)
    Km, Ks = mean_std_w(K, xw)
    aam, aas = mean_std_w(a, xw)
    cm, cs = mean_std_w(c, xw)
    return np.array([
        p, xSKm, xSKs, np.max(xw/K), #p, np.mean(xw/K)/p, np.std(xw/K)/np.sqrt(p), np.max(xw/K),
        rm, Km, aam, cm,
        rs, Ks, aas, cs])


##############################################################################################################

# PLOT SLIDER
def plotSlider(E, sols):
    fig = plt.figure()
    ax = fig.add_subplot()
    t_slider = Slider(plt.axes([0.1, 0.02, 0.25, 0.02], facecolor="lightyellow"), "Gamma", -10., np.log(1.9 * sols[-1].g_i), 0.)

    def update(val):
        ax.clear()
        g_t = np.exp(t_slider.val)
        P_c, i = np.full(N, True), 0
        while sols[i].g_f < g_t:
            P_c[sols[i].idx] = not(P_c[sols[i].idx])
            i += 1
        sol_c = sols[i]
        n, r, K, a, c = unravel(E)
        A2, A1, A0, D2, gi, gf, x = unravelSol(E, P_c, sol_c)
        ax.scatter(r, c, color="red", alpha=np.clip(x(g_t)/K, 0., 1.), marker="o")
        fig.canvas.draw_idle()

    t_slider.on_changed(update)
    plt.show()