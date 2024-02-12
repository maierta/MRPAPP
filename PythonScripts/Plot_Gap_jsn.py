from matplotlib.pyplot import *
from numpy import *
import argparse

# def eigen(file = "Gap.jsn",firstBZ=True,returnEvec=0, ncols=5, nrows=2, evList=list(range(0,10))):
def eigen(file,firstBZ,returnEvec, ncols, nrows, evList):
    data=eval(open(file).read())

    U  = array(data['U'])
    Up = array(data['Up'])
    J  = array(data['J'])
    Jp = array(data['Jp'])

    print("U,Up,J,Jp=",U,Up,J,Jp)

    e = array(data['Eigenvalues'])
    eSort = argsort(e)[0,::-1]
    e = e[0,eSort][0:20]
    # ratio1 = e[0]/e[0]
    # ratio2 = e[1]/e[0]
    # ratio3 = e[2]/e[0]
    # ratio4 = e[3]/e[0]
    # ratio5 = e[4]/e[0]

    # print ("Leading 5 eigenvectors: ",e)
    # print ("Ratios to first: ",ratio1,ratio2,ratio3,ratio4,ratio5)

    # evec = array(data['Eigenvectors'])[eSort,:]
    evec = array(data['Gap'])[eSort,:]
    kf = array(data['kfPoints'])
    nkf = kf.shape[0]
    if firstBZ:
        for i in range(nkf):
            if kf[i,0] > pi: kf[i,0] -= 2.*pi
            if kf[i,1] > pi: kf[i,1] -= 2.*pi


    # kz0 = abs(kf[:,2]-0.0) < 1.0e-3
    nkz = unique(kf[:,2]).shape[0]
    print(nkz,"different kz values!")

    # ncols = 5; nrows=2
    if (nrows > 1) | (ncols > 1):
        f, ax =subplots(ncols=ncols,nrows=nrows,figsize=(4*ncols,4*nrows))
        for kz in range(nkz):
            index = 0
            if (ncols>1) & (nrows>1):
                for i in range(ncols):
                    for j in range(nrows):
                        ax[j,i].scatter(kf[:,0]/pi,kf[:,1]/pi,c=evec[evList[index],:],cmap=get_cmap('RdBu_r'),s=50,lw=0.2)
                        ax[j,i].set_aspect('equal')
                        ax[j,i].grid(color='darkgrey')
                        ax[j,i].use_sticky_edges = False

                        # ax[j,i].margins(y=0.5, x=0.1)
                        ax[j,i].set_xlabel(r"$k_x/\pi$")
                        ax[j,i].set_ylabel(r"$k_y/\pi$")
                        ax[j,i].set_title(r'$\lambda=$'+str(round(e[evList[index]],4)))
                        index += 1
            else:
                for i in range(max(ncols, nrows)):
                    ax[i].scatter(kf[:,0]/pi,kf[:,1]/pi,c=evec[evList[index],:],cmap=get_cmap('RdBu_r'),s=50,lw=0.2)
                    ax[i].set_aspect('equal')
                    ax[i].grid(color='darkgrey')
                    ax[i].use_sticky_edges = False

                    # ax[i].margins(y=0.5, x=0.1)
                    ax[i].set_xlabel(r"$k_x/\pi$")
                    ax[i].set_ylabel(r"$k_y/\pi$")
                    ax[i].set_title(r'$\lambda=$'+str(round(e[evList[index]],4)))
                    index += 1
    else:
        f, ax =subplots(ncols=1,nrows=1,figsize=(4,4))
        ax.scatter(kf[:,0]/pi,kf[:,1]/pi,c=evec[evList[0],:],cmap=get_cmap('RdBu_r'),s=50,lw=0.2)
        ax.set_aspect('equal')
        ax.grid(color='darkgrey')
        ax.use_sticky_edges = False
        ax.set_xlabel(r"$k_x/\pi$")
        ax.set_ylabel(r"$k_y/\pi$")
        ax.set_title(r'$\lambda=$'+str(round(e[evList[0]],4)))

    f.suptitle(r"U="+str(U)+", U'="+str(Up)+", J="+str(J)+", J'="+str(Jp))
    # f.tight_layout()
    show()

    return kf,evec[returnEvec,:]

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="plot leading eigenvectors")
    parser.add_argument('--file', dest="file", action="store", default="Gap.jsn")
    parser.add_argument('--ncols', dest="ncols", action="store", default=5)
    parser.add_argument('--nrows', dest="nrows", action="store", default=2)
    parser.add_argument('--evList', dest="evList", nargs='+', type=int, action="store", default=list(range(0, 10)))
    input_args = parser.parse_args()
    kf, ev = eigen(file = input_args.file, firstBZ=True, returnEvec=0, ncols=int(input_args.ncols), nrows=int(input_args.nrows), evList=input_args.evList)
    for i in range(kf.shape[0]):
        print(kf[i, 0], " , ", kf[i, 1], " , ", ev[i])
