from matplotlib.pyplot import *
from numpy import *

def eigen(file = "Gap.jsn",firstBZ=True,returnEvec=0):
    data=eval(open(file).read())

    U  = array(data['U'])
    Up = array(data['Up'])
    J  = array(data['J'])
    Jp = array(data['Jp'])

    print("U,Up,J,Jp=",U,Up,J,Jp)

    e = array(data['Eigenvalues'])
    eSort = argsort(e)[0,::-1]
    e = e[0,eSort][0:10]
    ratio1 = e[0]/e[0]
    ratio2 = e[1]/e[0]
    ratio3 = e[2]/e[0]
    ratio4 = e[3]/e[0]
    ratio5 = e[4]/e[0]

    print ("Leading 5 eigenvectors: ",e)
    print ("Ratios to first: ",ratio1,ratio2,ratio3,ratio4,ratio5)

    evec = array(data['Eigenvectors'])[eSort,:]
    kf = array(data['kfPoints'])
    nkf = kf.shape[0]
    if firstBZ:
        for i in range(nkf):
            if kf[i,0] > pi: kf[i,0] -= 2.*pi
            if kf[i,1] > pi: kf[i,1] -= 2.*pi


    kz0 = abs(kf[:,2]-0.0) < 1.0e-3
    
    f, ((ax1,ax2,ax3,ax4,ax5),(ax6,ax7,ax8,ax9,ax10))=subplots(ncols=5,nrows=2,figsize=(20,8))
    ax1.scatter(kf[kz0,0],kf[kz0,1],c=evec[0,kz0],cmap=get_cmap('BrBG'),s=50,lw=0.2)
    ax2.scatter(kf[kz0,0],kf[kz0,1],c=evec[1,kz0],cmap=get_cmap('BrBG'),s=50,lw=0.2)
    ax3.scatter(kf[kz0,0],kf[kz0,1],c=evec[2,kz0],cmap=get_cmap('BrBG'),s=50,lw=0.2)
    ax4.scatter(kf[kz0,0],kf[kz0,1],c=evec[3,kz0],cmap=get_cmap('BrBG'),s=50,lw=0.2)
    ax5.scatter(kf[kz0,0],kf[kz0,1],c=evec[4,kz0],cmap=get_cmap('BrBG'),s=50,lw=0.2)
    ax6.scatter(kf[kz0,0],kf[kz0,1],c=evec[5,kz0],cmap=get_cmap('BrBG'),s=50,lw=0.2)
    ax7.scatter(kf[kz0,0],kf[kz0,1],c=evec[6,kz0],cmap=get_cmap('BrBG'),s=50,lw=0.2)
    ax8.scatter(kf[kz0,0],kf[kz0,1],c=evec[7,kz0],cmap=get_cmap('BrBG'),s=50,lw=0.2)
    ax9.scatter(kf[kz0,0],kf[kz0,1],c=evec[8,kz0],cmap=get_cmap('BrBG'),s=50,lw=0.2)
    ax10.scatter(kf[kz0,0],kf[kz0,1],c=evec[9,kz0],cmap=get_cmap('BrBG'),s=50,lw=0.2)
    
    i=0
    for ax in (ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9,ax10):
        ax.set_aspect('equal')
        ax.set_title(r'$\lambda=$'+str(round(e[i],4)))
        i+=1

    f.suptitle(r"U="+str(U)+", U'="+str(Up)+", J="+str(J)+", J'="+str(Jp))
    show()

    return kf,evec[returnEvec,kz0]

if __name__ == "__main__":
    eigen()