# import commands
import os
import sys
import time

from numpy import *

import shutil
import os

from pylab import *
from matplotlib.pyplot import *
from mpl_toolkits.mplot3d import Axes3D

# data=eval(open("susOfQFull.dat").read())
# data=data['chi0OfQ']



def calcChiQRPA(file,U=0,Up=0,J=0,Jp=0,skip=2,sublattice=1,variable="q",model="",ax=[]):

    # print "Us.shape: ",Us.shape[0]
    # d = loadtxt(file,delimiter=",",skiprows=1)
    fin  = open(file,'r')
    fin.readline()
    string = fin.readline()
    d = [int(x) for x in string.strip().split(',')]
    nqx = d[0]
    nqy = d[1]
    nqz = d[2]
    nw  = d[3]
    print ("nqx,nqy,nqz,nw=",nqx,nqy,nqz,nw)
    data=loadtxt(file,delimiter=",",skiprows=skip)
    norb = int(sqrt(sqrt(0.5*(data.shape[1]-6))))
    print (norb, "orbitals.")
    m=norb*norb

    if (model=='SrRuO'): 
        Us = setupIntMatricesSrRuO(norb, U, Up, J, Jp)
    else:
        Us = setupIntMatrices(norb,U,Up,J,Jp,sublattice)

    nq=data.__len__()
    chi0q = zeros((m,m),dtype='complex')
    chiPhysQ=zeros((nq),dtype='complex')
    chiOffQ=zeros((nq),dtype='complex')
    chiReal = zeros((m,m))
    chiImag = zeros((m,m))
    step = 3
    if (nw>0): 
        step = 4
        omega = data[:,3]
    for iq in range(0,nq):
        chiReal = data[iq,step:m*m+step].reshape(m,m)
        chiImag = data[iq,m*m+step:2*m*m+step].reshape(m,m)
        # chiPhys = data[iq,2*m*m+3]
        chi0q[:,:] = chiReal[:,:] + 1j*chiImag[:,:]
        chi0q[:,:] = calcChiSRPA(Us,chi0q) 
        # Calculate physical chi
        if (model=='SrRuO'):
            chiPhysQ[iq] = calcChiPhysSrRuO(chi0q)
        else:
            for l1 in range(0,norb):
                ind1=l1+l1*norb
                for l2 in range(0,norb):
                    ind2=l2+l2*norb
                    # if (nw==0 or sum(omega)==0):
                    chiPhysQ[iq] += real(chi0q[ind1,ind2]) +1j*imag(chi0q[ind1,ind2])
                    # else:
                    # chiPhysQ[iq] += imag(chi0q[ind1,ind2])
            # Calculate off-diagonal chi
            for l1 in range(0,norb):
                for l1 in range(0,norb):
                    if l2!=l1:
                        ind1=l1+l2*norb; ind2=ind1
                        if (nw==0 or sum(omega)==0):
                            chiOffQ[iq] += real(chi0q[ind1,ind2]) +1j*imag(chi0q[ind1,ind2])

    # for iq in range(0,nq):
        # if (chiPhysQ[iq]<0): print "chi(q) negative!! q=",data[iq,0],",",data[iq,1],",",data[iq,2]
    chiPhysQ /= 2
    if sublattice==1: chiPhysQ /= 2
    if (variable=="q"):
        qx = data[:,0].reshape(nqz,nqx,nqy)[0,:,:]
        qy = data[:,1].reshape(nqz,nqx,nqy)[0,:,:]
        qz = data[:,2].reshape(nqz,nqx,nqy)[:,0,0]
        nqlin = int(sqrt(nq))

        omega = data[0,3]
        if omega!=0.0:
            chi  = imag(chiPhysQ)
            chi2 = imag(chiOffQ)
        else:
            chi  = real(chiPhysQ)
            chi2 = real(chiOffQ)

        # chi = chiPhysQ.reshape(nqlin,nqlin)
        # chi = vstack((chi,chi[0]))
        # qx = vstack((qx,qx[0]))
        # qx = hstack((qx,-qx[:,0].reshape(nqlin+1,1)))
        # qy = vstack((qy,-qy[0]))
        # qy = hstack((qy,qy[:,0].reshape(nqlin+1,1)))
        # chi = hstack((chi,chi[:,0].reshape(nqlin+1,1)))
        # chi=chi.reshape((nqlin+1)*(nqlin+1))

   # Plotting
        
        for i in range(0,qz.shape[0]):
            # if sublattice==1: plotChiRPAQ(ax,0.5*(qx-qy),0.5*(qx+qy),qz,chi,i,U,Up,J,Jp)
            # if sublattice==0: plotChiRPAQ(ax,qx,qy,qz,chi,i,U,Up,J,Jp)

            plotChiRPAQ(ax,qx,qy,qz,omega,chi,i,U,Up,J,Jp)
            plotChiRPAQ(ax,qx,qy,qz,omega,chi2,i,U,Up,J,Jp)
            show()

    return column_stack((data[:,0],data[:,1],data[:,2],data[:,3],real(chiPhysQ),imag(chiPhysQ),real(chiOffQ),imag(chiOffQ)))
        # return qx,qy,qz,chi
    # elif (variable=="omega"):
        # if (ax==[]): ax=figure().add_subplot(111)
        # ax.plot(omega,chiPhysQ)
        # ax.set_xlabel("$\omega$")
        # ax.set_ylabel("$\chi''(q,\omega)$")
        # return omega,chiPhysQ



def plotChiRPAQ(ax,qx,qy,qz,omega,chiPhysQ,qzIndex,U,Up,J,Jp):
    nqz=qz.shape[0]
    nqx=int(sqrt(chiPhysQ.shape[0]//nqz))

    nqy=nqx
    if ax == []:
        fig = figure(figsize=figaspect(0.75))
        ax=fig.gca(projection='3d')
    ax1=ax
    # ax1.plot_surface(qx/pi,qy/pi,chiPhysQ.reshape(nqz,nqx,nqy)[qzIndex,:,:],cmap=cm.jet,rstride=1,cstride=1)
    ax1.plot_surface(qx/pi,qy/pi,chiPhysQ.reshape(nqz,nqx,nqy)[qzIndex,:,:],rstride=1,cstride=1,alpha=0.75,cmap='jet',lw=0.25)
    zmax=ax1.get_zlim()[1]
    xmin=ax1.get_xlim()[0]
    xmax=ax1.get_xlim()[1]
    ymin=ax1.get_ylim()[0]
    ymax=ax1.get_ylim()[1]
    ax1.set_zlim(0,zmax)
    ax1.contour(qx/pi,qy/pi,chiPhysQ.reshape(nqz,nqx,nqy)[qzIndex,:,:],zdir='z',offset=0,cmap='jet')

    # ax1.set_aspect(0.75)


    # cset=ax1.contour(qx/pi,qy/pi,chiPhysQ.reshape(nqz,nqx,nqy)[qzIndex,:,:],zdir='x',offset=ymin,cmap='jet')
    # cset=ax1.contour(qx/pi,qy/pi,chiPhysQ.reshape(nqz,nqx,nqy)[qzIndex,:,:],zdir='y',offset=xmin,cmap='jet')
    ax1.set_title("U="+str(U)+", U'="+str(Up)+", J="+str(J)+", J'="+str(Jp)+"; $q_z=$"+str(qz[qzIndex]))

    # ax1.set_title("U="+str(U)+", J="+str(J))
    
    ax1.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax1.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax1.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))

    ax1.set_xlabel(r"$q_x/\pi$",labelpad=15)
    ax1.set_ylabel(r"$q_y/\pi$",labelpad=15)
    if omega==0:
        ax1.set_zlabel(r"$\chi'(q,\omega=0)$",labelpad=15)
    else:
        ax1.set_zlabel(r"$\chi''(q,\omega=$"+str(round(omega,4))+"$)$",labelpad=15)


    # show()


def setupIntMatrices(norb,U,Up,J,Jp,sublattice=0):
    print("norb=",norb)
    limit = int(norb)
    if sublattice==1: limit=int(norb/2)

    m=int(limit*limit)
    print("m=",m)
    gamma = zeros((m,m))
    n=limit

    # diagonal terms
    for l1 in range(0,n):
        for l2 in range(0,n):
            ind1=l2+l1*n
            if (l1==l2):
                gamma[ind1,ind1] = U
            else:
                gamma[ind1,ind1] = Up

    # off-diagonal terms
    for l1 in range(0,n):
        ind1 = l1+l1*n
        for l2 in range(0,n):
            ind2 = l2+l2*n
            if (l2!=l1):
                gamma[ind1,ind2] = J
    # pair-hopping terms
    for l1 in range(0,n):
        for l2 in range(0,n):
            if (l2!=l1):
                ind1=l2+l1*n
                ind2=l1+l2*n
                gamma[ind1,ind2] = Jp    

    if sublattice==1:
        gammaFull = zeros((norb*norb,norb*norb))
        for l1 in range(0,limit):
            for l2 in range(0,limit):
                for l3 in range(0,limit):
                    for l4 in range(0,limit):
                        ind1=l2+l1*limit;
                        ind2=l4+l3*limit;
                        ind3=l2+l1*norb;
                        ind4=l4+l3*norb;
                        gammaFull[ind3,ind4] = gamma[ind1,ind2]
                        ind3=l2+limit+(l1+limit)*norb;
                        ind4=l4+limit+(l3+limit)*norb;
                        gammaFull[ind3,ind4] = gamma[ind1,ind2]
    if sublattice==0:
        gammaFull = gamma

    return gammaFull

def setupIntMatricesSrRuO(norb,U,Up,J,Jp, component="zz"):
    if (norb!=6):
        sys.exit("Wrong number of orbitals for SrRuO setup!")
    msize = norb*norb
    gammaFull = zeros((msize,msize))

    indexToOrb = zeros((msize,2),dtype="int")
    for l1 in range(norb):
        for l2 in range(norb):
                    ind=l2+l1*norb
                    indexToOrb[ind,0] = l1
                    indexToOrb[ind,1] = l2

    spinOfEll = zeros((6),dtype="int"); orbOfEll = zeros((6),dtype='int')
    spinOfEll[0] = +1; orbOfEll[0] = 0;
    spinOfEll[1] = +1; orbOfEll[1] = 1; 
    spinOfEll[2] = -1; orbOfEll[2] = 2; 
    spinOfEll[3] = -1; orbOfEll[3] = 0; 
    spinOfEll[4] = -1; orbOfEll[4] = 1; 
    spinOfEll[5] = +1; orbOfEll[5] = 2; 

    for ind1 in range(msize):
        for ind2 in range(msize):
                l1 = indexToOrb[ind1,0]; l2 = indexToOrb[ind1,1]
                l3 = indexToOrb[ind2,0]; l4 = indexToOrb[ind2,1]
                s1 = spinOfEll[l1]; s2 = spinOfEll[l2]
                s3 = spinOfEll[l3]; s4 = spinOfEll[l4]
                o1 = orbOfEll[l1] ; o2 = orbOfEll[l2]
                o3 = orbOfEll[l3] ; o4 = orbOfEll[l4]

                # U-terms
                if (o1 == o2 & o1 == o3 & o1 == o4):
                    if (s1 == -s2 & s1 == s3 & s1 == -s4):  gammaFull[ind1,ind2] += U 
                    if (s1 == s2 & s1 == -s3 & s1 == -s4):  gammaFull[ind1,ind2] -= U
                # U'-terms
                if (o1 != o2 & o1 == o3 & o2 == o4):
                    if (s1 == s3 & s2 == s4): gammaFull[ind1,ind2] += Up
                if (o1 == o2 & o1 != o3 & o3 == o4):
                    if (s1 == s2 & s3 == s4): gammaFull[ind1,ind2] -= Up
                # J-terms
                if (o1 == o2 & o1 != o3 & o3 == o4):
                    if (s1 == s3 & s2 == s4): gammaFull[ind1,ind2] += J
                if (o1 != o2 & o1 == o3 & o2 == o4):
                    if (s1 == s2 & s3 == s4): gammaFull[ind1,ind2] -= J
                # J'-terms
                if (o1 != o2 & o1 == o4 & o2 == o3):
                    if (s1 == s3 & s2 == s4 & s1 != s2): gammaFull[ind1,ind2] += Jp
                    if (s1 == s2 & s3 == s4 & s1 != s3): gammaFull[ind1,ind2] -= Jp
 

    return gammaFull

def calcChiSRPA(Us,chi0):
    m=chi0.shape[0]
    chiRPA = 0.0*ones_like(chi0)
    Uchi0Plus1 = -dot(Us,chi0[:,:])+eye(m)
    denom = linalg.inv(Uchi0Plus1)
    chiRPA[:,:] = dot(chi0[:,:],denom)
    return chiRPA

def calcChiPhysSrRuO(chi0, component="zz"):

    norb = 6
    msize = 6*6
    spinOfEll = zeros((6),dtype="int"); orbOfEll = zeros((6),dtype='int')
    spinOfEll[0] = +1; orbOfEll[0] = 0;
    spinOfEll[1] = +1; orbOfEll[1] = 1; 
    spinOfEll[2] = -1; orbOfEll[2] = 2; 
    spinOfEll[3] = -1; orbOfEll[3] = 0; 
    spinOfEll[4] = -1; orbOfEll[4] = 1; 
    spinOfEll[5] = +1; orbOfEll[5] = 2; 

    indexToOrb = zeros((msize,2),dtype="int")
    for l1 in range(norb):
        for l2 in range(norb):
                    ind=l2+l1*norb
                    indexToOrb[ind,0] = l1
                    indexToOrb[ind,1] = l2

    chiPhys = 0.0
    for ind1 in range(msize):
        for ind2 in range(msize):
            l1 = indexToOrb[ind1,0]; l2 = indexToOrb[ind1,1]
            l3 = indexToOrb[ind2,0]; l4 = indexToOrb[ind2,1]
            s1 = spinOfEll[l1]; s2 = spinOfEll[l2]
            s3 = spinOfEll[l3]; s4 = spinOfEll[l4]
            o1 = orbOfEll[l1] ; o2 = orbOfEll[l2]
            o3 = orbOfEll[l3] ; o4 = orbOfEll[l4]
            if (o1 == o2 & o3 == o4):
                # if (component == "zz" & s1 == s2 & s3 == s4):
                if (s1 == s2 & s3 == s4):
                    chiPhys += chi0[ind1,ind2] * s1 * s3
                # elif (component == "+-" & s1 == s3 & s2 == s4 & s1 != s4):
                #     chiPhys += chi0[ind1,ind2]
    return 0.25*chiPhys;

# file = "/Projects/LaOFeAs/Papers/15_LiFeAS/DATA/Siggi5orbitalModel/3D/susOfQFull.txt"
# plotChiQRPA(file,33,33,7)
