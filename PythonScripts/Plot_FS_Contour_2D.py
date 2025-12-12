import numpy as np
import matplotlib.pyplot as plt
import argparse


def plot_FS(file, nbands, nk=0):
    data = np.loadtxt(file, delimiter=",")
    if nk==0:
        nk = np.unique(data[:, 0]).shape[0]
    fig, ax = plt.subplots()
    # colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    for iband in range(nbands):
        cs = ax.contour(
            data[:, 0].reshape(nk, nk),
            data[:, 1].reshape(nk, nk),
            data[:, 3 + iband].reshape(nk, nk),
            levels=[0],
            # colors=colors[iband],
        )
        plt.clabel(cs, inline=1, fontsize=10, fmt=f"band {iband}")
    ax.set_aspect("equal")
    ax.set(xlabel=r"$k_x/\pi$", ylabel=r"$k_y/\pi$")
    ax.grid()
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot Fermi surface contour from ek.txt"
    )
    parser.add_argument("--file", dest="file", action="store", default="Gap.jsn")
    parser.add_argument("--nbands", type=int, dest="nbands", action="store", default=1)
    parser.add_argument("--nk", type=int, dest="nk", action="store", default=1)

    input_args = parser.parse_args()
    file = input_args.file
    nbands = input_args.nbands
    nk = input_args.nk
    if nk > 1:
        plot_FS(file, nbands, nk)
    else:
        plot_FS(file, nbands)
