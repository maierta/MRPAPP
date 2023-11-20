# %%
from numpy import *

# %%
data = loadtxt("./hopping_BaFeAs_5Orb.csv", delimiter=",")
data.shape

# %%
orb1 = (data[:, 3] - 1).astype(int)
orb2 = (data[:, 4] - 1).astype(int)
dx = data[:, 0]
dy = data[:, 1]
dz = data[:, 2]
tr = data[:, 5]
multiplicity = array([1]*orb1.shape[0]).astype(int)
# %%

# %%
dataDCA = ndarray((orb1.shape[0], 8), dtype=object)
dataDCA[:, 0] = orb1
dataDCA[:, 1] = orb2
dataDCA[:, 2] = dx
dataDCA[:, 3] = dy
dataDCA[:, 4] = dz
dataDCA[:, 5] = tr
dataDCA[:, 6] = 0.0
dataDCA[:, 7] = multiplicity
dataDCA
dataDCAC = dataDCA[abs(dataDCA[:, 5]) >= 1.0e-5]
dataDCAC.shape
savetxt("t_ij_BaFeAs.txt", dataDCAC, delimiter=",", fmt=['%i', '%i', '%i', '%i', '%i', '%.5e', '%.5e', '%i'])

# %% Now generate U_ij matrix. This has a linear dimension of 2 * mOrbitals with 2 for the spin
# the first nOrbitals are for spin up, the second nOrbitals elements for spin down
# So the intra-orbital Hubbard U should appear in the spn-offidagonal, but orbital diagonal elements

# %%
U = 2.0
Uij = zeros((10, 10))
for i in range(5):
    Uij[i,i+5] = U
    Uij[5+i,i] = U

count_nonzero(Uij-Uij.T)
# %%
savetxt("U_ij_BaFeAs.txt", Uij, delimiter=",", fmt = "%5.2f")
