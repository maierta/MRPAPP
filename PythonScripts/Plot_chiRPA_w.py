# %%

# Plot of ek_high_sym.dat (bands along high symm. directions in BZ)
from numpy import *
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
import matplotlib
import matplotlib as mpl

import sys
import imp
import matplotlib.gridspec as gridspec
from matplotlib.collections import LineCollection
import matplotlib.gridspec as gridspec


def plotChiRPAw(fileID = ""):

	print("fileID=",fileID)


	from matplotlib.colors import LinearSegmentedColormap
	colors=['#348ABD','#E24A33','#4431AC']

	plt.style.use(['ggplot']); mpl.rcParams['font.size'] = 12

	# gs = gridspec.GridSpec(1, 10)
	fig=plt.figure(figsize=(6,4))
	ax1=fig.add_subplot(111)

	data = np.loadtxt("chiRPA_"+fileID+".txt",delimiter=",")
	omega = data[:,3]
	chi = data[:,5]
	q = data[0,0:3]/pi
	qstr = str([round(i,2) for i in list(q)])


	ax1.plot(omega,chi,'-',lw=2,color=colors[0])

	ax1.set(ylabel=r"$\chi(q,\omega)$",xlabel="$\omega$")

	plt.savefig("chiRPAw_"+fileID+".pdf")

	# show()
# %%

if __name__ == "__main__":
    fileID = sys.argv[1]
    plotChiRPAw(fileID=fileID)



