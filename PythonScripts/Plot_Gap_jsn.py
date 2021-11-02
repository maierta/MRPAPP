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
    
    ncols = 5; nrows=2
    f, ax =subplots(ncols=ncols,nrows=nrows,figsize=(4*ncols,4*nrows))
    for kz in range(nkz):
        index = 0
        for i in range(ncols):
            for j in range(nrows):
                ax[j,i].scatter(kf[:,0],kf[:,1],c=evec[index,:],cmap=get_cmap('BrBG'),s=50,lw=0.2)
                ax[j,i].set_aspect('equal')
                ax[j,i].set_title(r'$\lambda=$'+str(round(e[index],4)))
                index += 1
    

    f.suptitle(r"U="+str(U)+", U'="+str(Up)+", J="+str(J)+", J'="+str(Jp))
    show()

    return kf,evec[returnEvec,:]

if __name__ == "__main__":
    file = sys.argv[1]
    eigen(file)