from matplotlib.pyplot import *
from numpy import *
import numpy as np
import argparse
import pandas as pd


def build_bStar(b1, b2):
    #  Build 8 rec. lattice vectors that surround K0=(0,0)
    b1 = np.pi * np.array(b1, dtype=float)[0:2]
    b2 = np.pi * np.array(b2, dtype=float)[0:2]
    b = np.zeros((8, 2))
    b[0] = b1
    b[1] = b2
    b[2] = -b1
    b[3] = -b2
    b[4] = b1 + b2
    b[5] = b1 - b2
    b[6] = -b1 + b2
    b[7] = -b1 - b2
    bStar = b
    return bStar


def shift21BZ(Kin, bStar):
    b = bStar
    #  Now check if K-point is closest to K0=(0,0) and, if not, subtract from it the b-vector that it is closest to
    K = Kin.copy()
    for Kvec in K[:, 0:2]:
        dist0 = np.linalg.norm(Kvec)
        distClosest = dist0
        ibClosest = -1  # corresponds to K=(0,0)
        for ib in range(8):
            distb = np.linalg.norm(Kvec - b[ib])
            if distb < distClosest:
                distClosest = distb
                ibClosest = ib
        if ibClosest > -1:
            Kvec -= b[ibClosest]
    return K


# def eigen(file = "Gap.jsn",firstBZ=True,returnEvec=0, ncols=5, nrows=2, evList=list(range(0,10))):
def eigen(file, firstBZ, bStar, returnEvec, ncols, nrows, evList):
    data = eval(open(file).read())

    U = np.array(data["U"])
    Up = np.array(data["Up"])
    J = np.array(data["J"])
    Jp = np.array(data["Jp"])

    print("U,Up,J,Jp=", U, Up, J, Jp)

    e = np.array(data["Eigenvalues"])
    eSort = np.argsort(e)[0, ::-1]
    e = e[0, eSort][0:20]
    evec = np.array(data["Gap"])[eSort, :]
    kf = np.array(data["kfPoints"])
    # nkf = kf.shape[0]
    if firstBZ:
        kf = shift21BZ(kf, bStar)
        # for i in range(nkf):
        #     if kf[i, 0] > pi:
        #         kf[i, 0] -= 2.0 * pi
        #     if kf[i, 1] > pi:
        #         kf[i, 1] -= 2.0 * pi

    # kz0 = abs(kf[:,2]-0.0) < 1.0e-3
    nk = kf.shape[0]
    nkz = unique(kf[:, 2]).shape[0]
    print(nkz, "different kz values!")

    dd = pd.DataFrame([[kf[i, 0], kf[i, 1], evec[i, evList[0]], evList[0]] for i in range(nk)], columns=["kx", "ky", "evec", "index"])
    for j in range(1, len(evList)):
        dd2 = pd.DataFrame([[kf[i, 0], kf[i, 1], evec[i, evList[j]], evList[j]] for i in range(nk)], columns=["kx", "ky", "evec", "index"])
        dd = pd.concat([dd, dd2], ignore_index=True)

    dd.to_csv(file+"_Plot_Gap.jsn_eigenvectors.csv", index=False)

    # ncols = 5; nrows=2
    if (nrows > 1) | (ncols > 1):
        f, ax = subplots(ncols=ncols, nrows=nrows, figsize=(4 * ncols, 4 * nrows))
        for kz in range(nkz):
            index = 0
            if (ncols > 1) & (nrows > 1):
                for i in range(ncols):
                    for j in range(nrows):
                        # cc = (
                        #     evec[evList[index], :] / abs(evec[evList[index], :]).max()
                        # )  # normalize
                        evec[evList[index], :] /= abs(evec[evList[index], :]).max()
                        cc = evec[evList[index], :]
                        col = ["darkgreen"] * nk
                        for ik in range(0, nk):
                            if sign(cc[ik]) < 0:
                                col[ik] = "orange"
                        ax[j, i].scatter(
                            kf[:, 0] / pi,
                            kf[:, 1] / pi,
                            c=col,
                            # c=cc,
                            # cmap=get_cmap("RdBu_r"),
                            s=100 * abs(cc),
                            # s=50,
                            lw=0.2,
                            # vmin=-1,
                            # vmax=1,
                        )
                        ax[j, i].set_aspect("equal")
                        ax[j, i].grid(color="darkgrey")
                        ax[j, i].use_sticky_edges = False

                        # ax[j,i].margins(y=0.5, x=0.1)
                        ax[j, i].set_xlabel(r"$k_x/\pi$")
                        ax[j, i].set_ylabel(r"$k_y/\pi$")
                        ax[j, i].set_title(
                            r"$\lambda=$" + str(round(e[evList[index]], 4))
                        )
                        # ax[j, i].set(xticks=(-1, 0, 1), yticks=(-1, 0, 1))
                        index += 1
            else:
                for i in range(max([ncols, nrows])):
                    cc = (
                        evec[evList[index], :] / abs(evec[evList[index], :]).max()
                    )  # normalize
                    ax[i].scatter(
                        kf[:, 0] / pi,
                        kf[:, 1] / pi,
                        c=cc,
                        cmap=get_cmap("RdBu_r"),
                        s=50,
                        lw=0.2,
                        vmin=-1,
                        vmax=1,
                    )
                    ax[i].set_aspect("equal")
                    ax[i].grid(color="darkgrey")
                    ax[i].use_sticky_edges = False

                    # ax[i].margins(y=0.5, x=0.1)
                    ax[i].set_xlabel(r"$k_x/\pi$")
                    ax[i].set_ylabel(r"$k_y/\pi$")
                    ax[i].set_title(r"$\lambda=$" + str(round(e[evList[index]], 4)))
                    # ax[i].set(xticks=(-1, 0, 1), yticks=(-1, 0, 1))
                    index += 1
    else:  # just a single subplot
        f, ax = subplots(ncols=1, nrows=1, figsize=(6, 5))
        cc = evec[evList[0], :] / abs(evec[evList[0], :]).max()  # normalize
        cb = ax.scatter(
            kf[:, 0] / pi,
            kf[:, 1] / pi,
            c=cc,
            cmap=get_cmap("RdBu_r"),
            s=50,
            lw=0.2,
            vmin=-1,
            vmax=1,
        )
        f.colorbar(cb)
        ax.set_aspect("equal")
        ax.grid(color="darkgrey")
        ax.use_sticky_edges = False
        ax.set_xlabel(r"$k_x/\pi$")
        ax.set_ylabel(r"$k_y/\pi$")
        ax.set_title(r"$\lambda=$" + str(round(e[evList[0]], 4)))
        # ax.set(xticks=(-1, 0, 1), yticks=(-1, 0, 1))

    f.suptitle(r"U=" + str(U) + ", U'=" + str(Up) + ", J=" + str(J) + ", J'=" + str(Jp))
    f.tight_layout()
    show()

    return kf, evec[evList, :], f


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="plot leading eigenvectors")
    parser.add_argument("--file", dest="file", action="store", default="Gap.jsn")
    parser.add_argument(
        "--firstBZ", dest="firstBZ", action=argparse.BooleanOptionalAction
    )
    parser.add_argument("--ncols", dest="ncols", action="store", default=5)
    parser.add_argument("--nrows", dest="nrows", action="store", default=2)
    parser.add_argument(
        "--evList",
        dest="evList",
        nargs="+",
        type=int,
        action="store",
        default=list(range(0, 20)),
    )
    parser.add_argument(
        "--b1",
        type=float,
        nargs="+",
        default=[1, 0],
        help="reciprocal lattice vector b1 in units of pi",
    )
    parser.add_argument(
        "--b2",
        type=float,
        nargs="+",
        default=[0, 1],
        help="reciprocal lattice vector b2 in units of pi",
    )
    input_args = parser.parse_args()
    print(f"b1: {input_args.b1}")
    print(f"b2: {input_args.b2}")

    bStar = build_bStar(input_args.b1, input_args.b2)
    print("bStar=", bStar)

    kf, ev, fig = eigen(
        file=input_args.file,
        firstBZ=input_args.firstBZ,
        bStar=bStar,
        returnEvec=input_args.evList,
        ncols=int(input_args.ncols),
        nrows=int(input_args.nrows),
        evList=input_args.evList,
    )
    print(f"{kf.shape[0]} points on the Fermi surface")
    # for i in range(kf.shape[0]):
    #     print(kf[i, 0], " , ", kf[i, 1], " , ", ev[:, i])
