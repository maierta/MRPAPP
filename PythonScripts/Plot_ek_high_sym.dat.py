#! python
import sys
from numpy import *
import matplotlib as mpl
import matplotlib.pyplot as plt

plt.style.use("ggplot")
mpl.rcParams['lines.markeredgewidth']=0.0
mpl.rcParams['font.size']=16

data=loadtxt("./ek_high_sym.dat")
nk = data.shape[0]

fig,ax=plt.subplots()

ax.plot(range(0,nk),data[:,3:data.shape[1]],'-')
plt.tick_params(bottom='off',left='off',top='off',right='off')
plt.margins(0.05)
ax.set_xticks((0,nk/3,2*nk/3,nk))
ax.set_xticklabels(("$\Gamma$","X","M","$\Gamma$"))

ax.set_xlabel("$k$")
ax.set_ylabel("$E_n(k)$")



plt.show()

