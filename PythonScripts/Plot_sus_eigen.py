# %%
import numpy as np
import matplotlib.pyplot as plt
import argparse
import seaborn as sns
import pandas as pd
from plotnine import *

script_name = "Plot_sus_eigen.py"

def plot_sus_eigen(file, nOrb):

    data = np.loadtxt(
        file,
        delimiter=",",
        skiprows=2,
    )

    nq = data.__len__()
    m = nOrb * nOrb
    qVec = np.zeros((nq, 3))
    chi0 = np.zeros((nq, m, m), dtype="complex")
    chi0Phys = np.zeros((nq), dtype="complex")
    step = 3
    for iq in range(nq):
        qVec[iq, :] = data[iq, 0:3]
        chiReal = data[iq, 4 : step + m * m + 1].reshape(m, m)
        chiImag = data[iq, step + m * m + 1 : step + 2 * m * m + 1].reshape(m, m)
        chi0[iq, :, :] = chiReal + 1j * chiImag
        chi0Phys[iq] = data[iq, step + 2 * m * m + 1]

    # Now look at diagonal elements only, i.e. chi_{l1,l1, l2,l2}, since only those contribute to the physical susceptibility
    chiDiag = chi0.reshape(nq, nOrb, nOrb, nOrb, nOrb)
    chiPhys = np.zeros((nq, nOrb, nOrb), dtype="complex")

    l = np.arange(nOrb)
    chiPhys[:, l[:, None], l] = chiDiag[:, l[:, None], l[:, None], l, l]

    lambdaSDiag = np.zeros((nq))
    lambdaSDiagSub = np.zeros((nq))
    evecSDiag = np.zeros((nq, nOrb), dtype="complex")
    evecSDiagSub = np.zeros((nq, nOrb), dtype="complex")

    # Vectorized Hermitian check (optional: can be commented out for performance)
    for iq in range(nq):
        if not np.allclose(chiPhys[iq], chiPhys[iq].conj().T):
            print(f"chiPhys not hermitian for iq={iq}")

    # Vectorized eigenvalue decomposition
    w, v = np.linalg.eigh(chiPhys)  # w: (nq, nOrb), v: (nq, nOrb, nOrb)
    lambdaSDiag = np.real(w[:, -1])
    lambdaSDiagSub = np.real(w[:, -2])
    evecSDiag = v[:, :, -1]
    evecSDiagSub = v[:, :, -2]

    # fig, ax = plt.subplots()
    # ax.plot(lambdaSDiag)
    # ax.set(
    #     xticks=[
    #         0,
    #         nq // 3,
    #         2 * nq // 3,
    #         3 * nq // 3,
    #     ]
    # )
    # # Gamma -> M -> X -> Gamma
    # ax.set(xticklabels=[r"$\Gamma$", "M", "X", r"$\Gamma$"])
    # ax.set(ylabel=r"$\lambda_S$", xlabel="$q$")
    # ax.set(xlabel="$q$")


    data = pd.DataFrame({"q": np.arange(nq), "lambda_S": lambdaSDiag})
    g = (
        ggplot(data, aes(x="q", y="lambda_S"))
        + geom_line()
        + theme_tufte(base_size=14, base_family="Arial", ticks=True)
        + theme(panel_grid_major=element_line(color="Lightgrey"),
            axis_ticks_length=10)
        + labs(x="q", y=r"$\lambda_S(q)$")
        + scale_x_continuous(
            breaks=[0, nq // 3, 2 * nq // 3, 3 * nq // 3],
            labels=[r"$\Gamma$", "M", "X", r"$\Gamma$"])

        )
    g.show()

    g.save(filename=file+"_"+script_name+"_lambdaS.png", dpi=300)
    data.to_csv(file+"_"+script_name+"_lambdaS.csv", index=False)

    imax = np.argmax(lambdaSDiag)
    print("Maximum eigenvalue at q =", qVec[imax], "with leading eigenvalues =", lambdaSDiag[imax], lambdaSDiagSub[imax])
    evecMax0 = evecSDiag[imax]
    evecMax1 = evecSDiagSub[imax]
    species = ("Layer 1", "Layer 2")
    weight_counts0 = {
        r"$d_{3z^2-r^2}$": evecMax0.real[::2],
        r"$d_{x^2-y^2}$": evecMax0.real[1::2],
    }
    weight_counts1 = {
        r"$d_{3z^2-r^2}$": evecMax1.real[::2],
        r"$d_{x^2-y^2}$": evecMax1.real[1::2],
    }

    data2 = pd.DataFrame(weight_counts0, index=species).reset_index().melt(id_vars='index')
    data3 = pd.DataFrame(weight_counts1, index=species).reset_index().melt(id_vars='index')

    g2 = (
        ggplot(data2, aes(x='index', y='value', fill='variable'))
        + geom_bar(stat='identity', position='stack', width=0.5)
        + theme_tufte(base_size=14, base_family="Arial", ticks=True)
        + theme(panel_grid_major=element_line(color="Lightgrey"),
            axis_ticks_length=10)
        + labs(x="", y=r"$\phi_\ell$ for in-plane $q_{max}$", fill="Orbital")
    )
    g2.show()

    g2.save(filename=file+"_"+script_name+"_evecMax.png", dpi=300)
    data2.to_csv(file+"_"+script_name+"_evecMax.csv", index=False)

    g3 = (
        ggplot(data3, aes(x='index', y='value', fill='variable'))
        + geom_bar(stat='identity', position='stack', width=0.5)
        + theme_tufte(base_size=14, base_family="Arial", ticks=True)
        + theme(panel_grid_major=element_line(color="Lightgrey"),
            axis_ticks_length=10)
        + labs(x="", y=r"$\phi_\ell$ for in-plane $q_{max}$", fill="Orbital")
    )
    g3.show()

    g3.save(filename=file+"_"+script_name+"_evecMaxSub.png", dpi=300)
    data3.to_csv(file+"_"+script_name+"_evecMaxSub.csv", index=False)

    # fig, ax = plt.subplots()
    # bottom = np.zeros(2)
    # for boolean, weight_count in weight_counts.items():
    #     p = ax.bar(species, weight_count, 0.5, label=boolean, bottom=bottom)
    #     bottom += weight_count
    # ax.legend()
    # ax.set_title(r"$\phi_\ell$ for in-plane $q_{max}$")
    # plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="plot leading eigenvectors")
    parser.add_argument("--file", dest="file", action="store", default="chi1Full.txt")
    parser.add_argument("--nOrb", dest="nOrb", type=int, action="store", default=1)

    input_args = parser.parse_args()

    plot_sus_eigen(input_args.file, input_args.nOrb)

