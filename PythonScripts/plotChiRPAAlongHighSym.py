import os
import sys
import time
import imp


from numpy import *
from matplotlib.pyplot import *
import matplotlib
matplotlib.rcParams.update({'font.size': 14})
matplotlib.rcParams.update({'lines.linewidth': 2})


def plotchi(U,Up,J,Jp):

	sys.path.append("/Users/t7m/Dropbox/C++Code/")

	import calcAndPlotRPAChiQ as c
	imp.reload(c)

	# U=1.0; Up=0.5; J=0.25; Jp=0.25
	# U=1.0; Up=0.8; J=0.1; Jp=0.1
	# U=1.3; Up=1.05; J=0.125; Jp=0.125
	# U=1.2; Up=0.6; J=0.3; Jp=0.3
	# U=1.0; Up=0.8; J=0.1; Jp=0.1
	# U=1.0; Up=0.; J=0.; Jp=0.
	# U=0.8; Up=0.8; J=0.4; Jp=0.4

	chi = c.calcChiQRPA("chi0Full.txt",U=U,Up=Up,J=J,Jp=Jp,sublattice=0)
	nqx = int(sqrt(chi.shape[0]))

	chiPlotDiag = zeros((3*nqx))
	chiPlotOffD = zeros((3*nqx))

	chiPlotDiag[0:nqx] = 2.*chi[chi[:,1]==0][:,4]
	chiPlotOffD[0:nqx] = chi[chi[:,1]==0][:,6]

	chiPlotDiag[nqx:2*nqx] = 2.*chi[chi[:,0]==3.1415926536][:,4]
	chiPlotOffD[nqx:2*nqx] = chi[chi[:,0]==3.1415926536][:,6]

	chiPlotDiag[2*nqx:3*nqx] = 2.*chi[chi[:,0]==chi[:,1]][::-1][:,4]
	chiPlotOffD[2*nqx:3*nqx] = chi[chi[:,0]==chi[:,1]][::-1][:,6]

	if (True):
		fig, ax = subplots(figsize=(10,6))
		ax.plot(chiPlotDiag,label=r"$\chi_{intra}$")
		ax.plot(chiPlotOffD,label=r"$\chi_{inter}$")

		nq = chiPlotDiag.shape[0]

		ax.set(xticks=(0,int(nq/3),int(2*nq/3),int(nq)),xticklabels=(r"$\Gamma$","X","M",r"$\Gamma$"),ylabel=r"$\chi(q)$",xlabel="q")
		ax.axvline(x=0,linestyle="--",color="lightgray",lw=0.5)
		ax.axvline(x=int(nq/3),linestyle="--",color="lightgray",lw=0.5)
		ax.axvline(x=int(2*nq/3),linestyle="--",color="lightgray",lw=0.5)
		ax.axvline(x=int(nq),linestyle="--",color="lightgray",lw=0.5)
		ax.set_title("U="+str(U)+", U'="+str(Up)+", J="+str(J)+", J'="+str(Jp))
		ax.legend()
		show()

if __name__ == "__main__":
	U  = float(sys.argv[1])
	if len(sys.argv)>2:
		Up = float(sys.argv[2])
		J  = float(sys.argv[3])
		Jp = float(sys.argv[4])
	else:
		Up = U/2; J=U/4; Jp=U/4

	plotchi(U=U,Up=Up,J=J,Jp=Jp)
