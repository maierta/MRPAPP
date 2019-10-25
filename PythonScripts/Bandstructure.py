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


def plotBandstructure(fileID = ""):

	print("fileID=",fileID)


	from matplotlib.colors import LinearSegmentedColormap
	colors=['#348ABD','#E24A33','#4431AC']
	cmap = LinearSegmentedColormap.from_list('my_list',colors,N=256)


	# nbands = 2

	def plot_colourline(x,y,c,ax):
	    c = cmap(c)
	    for i in np.arange(len(x)-1):
	        ax.plot([x[i],x[i+1]], [y[i],y[i+1]], c=c[i], lw=3)
	    return



	plt.style.use(['ggplot']); mpl.rcParams['font.size'] = 12

	gs = gridspec.GridSpec(1, 10)
	fig=plt.figure(figsize=(12,4))

	ax1=fig.add_subplot(gs[0, 0:6])

	ax2=fig.add_subplot(gs[0, 6:10])
	 
	# fig=plt.figure(figsize=(8,6))

	# ax1=fig.add_subplot(111)


	# fig,ax=subplots(nrows=2,ncols=2)

	data = np.loadtxt("ek_high_sym_"+fileID+".txt")

	nbands = int(0.5*(data.shape[1]-3))
	nk = data.shape[0]

	ek = zeros((nk,nbands))
	cc = zeros((nk,nbands))

	ek = data[:,3:3+nbands+1]
	# cc = data[:-1,5:5+nbands+1]**2
	# ek2 = data[:,4] ; c2 = data[:-1,6]**2

	x=np.arange(0,ek.shape[0])

	for ib in range(nbands):
		ax1.plot(x,data[:,3+ib],'-',lw=2,color=colors[ib])

	# plot_colourline(x,ek1,c1,ax1)
	# plot_colourline(x,ek2,c2,ax1)


	ax1.axhline(y=0,color='grey',ls='--')
	# ax1.set(xlim=(x.min(),x.max()),ylim=(ek1.min()-0.5,ek2.max()+0.5))
	nkSeg = int(x.shape[0]/3)
	ax1.set(xticks=[0,nkSeg,2*nkSeg,3*nkSeg])
	ax1.set(xticklabels=[r"$\Gamma$","X","M",r"$\Gamma$"])
	ax1.set(ylabel=r"$E(k)$",xlabel="$k$")



	data2=np.loadtxt("ek_"+fileID+".txt")
	nklin = int(sqrt(data2.shape[0]))
	eps = 2/nklin

	for ib in range(nbands):
		ax2.contour(data2[:,3+ib].reshape(nklin,nklin),levels=[0],extent=(-1,1-eps,-1,1-eps),colors=colors[ib])

	ax2.set(xlim=(-1.1,1.1),ylim=(-1.1,1.1))
	ax2.set_xticks((-1,0,1))
	ax2.set_yticks((-1,0,1))
	ax2.set(xlabel=r"$k_x/\pi$",ylabel=r"$k_y/\pi$")
	ax2.set_aspect("equal")

	plt.subplots_adjust(left=0.1, bottom=0.2, right=1.0, top=0.9, wspace=1.5, hspace=0)

	plt.savefig("Bandstructure_"+fileID+".pdf")

	# show()
# %%

if __name__ == "__main__":
    fileID = sys.argv[1]
    plotBandstructure(fileID=fileID)



