import sys
import numpy as np
import matplotlib.pyplot as plt


def analyze(filenames):

    fig, ax = plt.subplots(1, 1)

    for filename in filenames:
        Neval, I, err = np.loadtxt(filename, unpack=True, skiprows=1)
        with open(filename, "r") as f:
            name = f.readline()
        ax.plot(Neval, err, label=name)

    ax.set(xlabel=r"$N_{eval}$", xscale='log',
           ylabel=r"$L_1$ error", yscale='log')
    ax.legend()

    fig.tight_layout()
    figname = "convergence.pdf"
    print("Saving", figname)
    fig.savefig(figname)
    
    
if __name__ == "__main__":

    analyze(sys.argv[1:])

    plt.show()
