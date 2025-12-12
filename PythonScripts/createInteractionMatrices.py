import numpy as np
import argparse

def setupUMatrices(nOrb, nsites, orbToSite, U, Up, J, Jp):
    Us = np.zeros((nOrb**2, nOrb**2))
    Uc = np.zeros((nOrb**2, nOrb**2))
    for l1 in range(nOrb):
        for l2 in range(nOrb):
            if orbToSite[l1] != orbToSite[l2]:
                continue

            ind1 = l2 + l1 * nOrb
            ind2 = l1 + l1 * nOrb
            ind3 = l2 + l2 * nOrb
            ind4 = l1 + l2 * nOrb

            if l1 == l2:
                Us[ind1, ind1] = U[orbToSite[l1]]
                Uc[ind1, ind1] = -U[orbToSite[l1]]
            else:
                Us[ind1, ind1] = Up[orbToSite[l1]]
                Uc[ind1, ind1] = Up[orbToSite[l1]] - 2*J[orbToSite[l1]]
                Us[ind2, ind3] = J[orbToSite[l1]]
                Uc[ind2, ind3] = -2*Up[orbToSite[l1]] + J[orbToSite[l1]]
                Us[ind1, ind4] = Jp[orbToSite[l1]]
                Uc[ind1, ind4] = -Jp[orbToSite[l1]]
    return Us, Uc

def createUMatrices(fileID):

    # Edit the following material dependent parameters
    nOrb = 10
    nSites = 5

    orbToSite = np.zeros((nOrb), dtype='int')
    orbToSite[0] = 0
    orbToSite[1] = 0
    orbToSite[2] = 1
    orbToSite[3] = 1
    orbToSite[4] = 2
    orbToSite[5] = 2
    orbToSite[6] = 3
    orbToSite[7] = 3
    orbToSite[8] = 4
    orbToSite[9] = 4

    U = np.zeros((nSites))
    Up = np.zeros((nSites))
    J = np.zeros((nSites))
    Jp = np.zeros((nSites))

    U[0:2] = 0.64
    Up[0:2] = U[0:2] / 2
    J[0:2] = U[0:2] / 4
    Jp[0:2] = U[0:2] / 4

    U[2:5] = 0.8
    Up[2:5] = U[2:5] / 2
    J[2:5] = U[2:5] / 4
    Jp[2:5] = U[2:5] / 4

    #############################

    Us, Uc = setupUMatrices(nOrb, nSites, orbToSite, U, Up, J, Jp)

    np.savetxt("Us_"+fileID+".txt", Us, fmt='%.4e', delimiter=",")
    np.savetxt("Uc_"+fileID+".txt", Uc, fmt='%.4e', delimiter=",")

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="create and save RPA interaction matrices for multi-site problems")
    parser.add_argument('--fileID', dest="fileID", action="store", default="v0")
    input_args = parser.parse_args()

    createUMatrices(fileID = input_args.fileID)
