import numpy as np
import matplotlib.pyplot as plt
import argparse
import pandas as pd

def plot_Gammakkp(file, ind0 = 0):
    # Analyze momentum structure of pairing interaction Gamma(k,k')
    data = eval(open(file).read())
    # data.keys()
    gamma = np.array(data["GammaPP"])
    kf = np.array(data["kfPoints"])

    dd = pd.DataFrame([[kf[i, 0], kf[i, 1], gamma[ind0, i]] for i in range(kf.shape[0])], columns=["kx", "ky", "Gamma"])
    dd.to_csv(file + "_Plot_Gammakkp_k0.py_index" + str(ind0) + ".csv", index=False)

     # Plot Gamma(k0,k)

    fig, ax = plt.subplots()
    c = ax.scatter(x=kf[:, 0], y=kf[:, 1], c=gamma[ind0, :])
    plt.scatter(kf[ind0, 0], kf[ind0, 1], color="k", s=100)
    ax.set_aspect("equal")
    ax.set_xticks([0, np.pi, 2 * np.pi], ["0", r"$\pi$", r"$2\pi$"])
    ax.set_yticks([0, np.pi, 2 * np.pi], ["0", r"$\pi$", r"$2\pi$"])
    ax.set_title(r"$\Gamma^{pp}(k_0, k)$")
    fig.colorbar(c)
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="plot leading eigenvectors")
    parser.add_argument("--file", dest="file", action="store", default="chi1Full.txt")
    parser.add_argument("--index", dest="index", type=int, action="store", default=0)

    input_args = parser.parse_args()

    plot_Gammakkp(input_args.file, input_args.index)


