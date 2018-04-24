from mayavi import mlab
from numpy import *

bands=loadtxt("ek.dat")

# nkx = unique(bands[:,0]).shape[0]
# nky = unique(bands[:,1]).shape[0]
nkz = unique(bands[:,2]).shape[0]
nkx = int(sqrt(bands.shape[0]/nkz))
nky=nkx
nbands = bands.shape[1]-3

x,y,z=mgrid[-1:1:1j*nkx,-1:1:1j*nky,-1.25:1.25:1j*nkz]

for i in range(3,nbands+3):
	FS = bands[:,i].min()*bands[:,i].max()
	if FS < 0:
		ek = bands[:,i].reshape(nkx,nky,nkz,order='F')
		mlab.contour3d(x,y,z,ek,contours=[0])

mlab.outline()





