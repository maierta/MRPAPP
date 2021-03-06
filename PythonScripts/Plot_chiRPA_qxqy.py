#! python
import sys
from numpy import *
from matplotlib.pyplot import *



def plotchi(file="chiRPA.txt"):

	print("Plotting file",file)
	imag = False # Plot imaginary part of chi; otherwise real part
	column = 4
	# column = 6
	if imag: column=5

	# file = "chiRPA_U1.0_Up0.5_J0.25_Jp0.25.txt"
	# file = "chiRPA.txt"

	data=loadtxt(file,delimiter=",")
	from mpl_toolkits.mplot3d import Axes3D

	nk = unique(data[:,0]).shape[0]
	nkz = unique(data[:,2]).shape[0]
	#x,y = mgrid[0:1:nk*1j,0:1:nk*1j]
	x = data[:,0].reshape(nk,nk)/pi
	y = data[:,1].reshape(nk,nk)/pi

	omega = data[0,3]

	ax=figure().gca(projection='3d')
	z=data[:,column]
	ax.plot_surface(x,y,z.reshape(nk,nk),rstride=1,cstride=1,linewidth=1.0,alpha=0.75,cmap='jet')
	zmax=ax.get_zlim()[1]
	ax.set_zlim(0,zmax)
	ax.contour(x,y,z.reshape(nk,nk),zdir='z',offset=0,cmap='jet')

	# ax.set_aspect(0.75)

	ax.set_xlabel(r"$q_x/\pi$",fontsize=13,labelpad=10)
	ax.set_ylabel(r"$q_y/\pi$",fontsize=13,labelpad=10)
	if imag:
	    ax.set_title(r"$\chi\prime\prime(q,\omega)$ for $\omega=$"+str(round(omega,5)),fontsize=13)
	else:
	    ax.set_title(r"$\chi\prime(q,\omega)$ for $\omega=$"+str(round(omega,5)),fontsize=13)

	# ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
	# ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
	# ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))

	# ax.set_title(file)

	show()

	fig2 = figure(figsize=(8,6))
	ax2=fig2.add_subplot(111)

	cset = ax2.imshow(z.reshape(nk,nk),extent=(x.min(),x.max(),y.min(),y.max()),interpolation="bilinear",origin='lower',cmap='jet')
	fig2.colorbar(cset)
	ax2.set_xlabel(r"$q_x/\pi$",fontsize=13,labelpad=10)
	ax2.set_ylabel(r"$q_y/\pi$",fontsize=13,labelpad=10)

	if imag:
	    ax2.set_title(r"$\chi\prime\prime(q,\omega)$ for $\omega=$"+str(round(omega,5)))
	else:
	    ax2.set_title(r"$\chi\prime(q,\omega)$ for $\omega=$"+str(round(omega,5)))

	show()

if __name__ == "__main__":
    if len(sys.argv)>1:
    	file = sys.argv[1]
    	plotchi(file)
    else:
    	plotchi()
    
