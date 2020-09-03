import sys
import numpy as np
import matplotlib.pyplot as plt


def analyze(filenames):

    fig, ax = plt.subplots(1, 1)

    maxN = 0

    for filename in filenames:
        Neval, I, err = np.loadtxt(filename, unpack=True, skiprows=1,
                                   usecols=[0,1,2])
        with open(filename, "r") as f:
            name = f.readline()

        ln = ax.plot(Neval, err, label=name)

        if "adapt" in name:    
            err_est = np.loadtxt(filename, unpack=True, skiprows=1,
                                 usecols=[3])
            ax.plot(Neval, err_est, color=ln[-1].get_color(), ls=':')

        maxN = max(maxN, Neval.max())

    n = np.geomspace(1.0, maxN, 100)

    ax.plot(n, np.power(n, -2), lw=2, ls='--', color='grey', alpha=0.8,
            label=r"$N^{-2}$")
    ax.plot(n, np.power(n, -4), lw=2, ls='--', color='grey', alpha=0.8,
            label=r"$N^{-4}$")
    ax.plot(n, np.power(n, -6), lw=2, ls='--', color='grey', alpha=0.8,
            label=r"$N^{-6}$")

    ax.set(xlabel=r"$N_{eval}$", xscale='log',
           ylabel=r"$L_1$ error", yscale='log', ylim=(3.0e-17, None))
    ax.legend()

    fig.tight_layout()
    figname = "convergence.pdf"
    print("Saving", figname)
    fig.savefig(figname)
    
    
if __name__ == "__main__":

    analyze(sys.argv[1:])

    plt.show()
