import argparse
from numpy import *
from matplotlib.pyplot import *


style.use("bmh")

def plotchi(file, interpolate, column=8):

    print("Plotting file",file)
    imag = False # Plot imaginary part of chi; otherwise real part
    # column = 8
    # column = 6
    print("Showing column ", column)
    if imag: column=5

    # file = "chiRPA_U1.0_Up0.5_J0.25_Jp0.25.txt"
    # file = "chiRPA.txt"

    data=loadtxt(file,delimiter=",")

    nk = unique(data[:,0]).shape[0]
    # nkz = unique(data[:,2]).shape[0]
    #x,y = mgrid[0:1:nk*1j,0:1:nk*1j]
    x = data[:,0].reshape(nk,nk)/pi
    y = data[:,1].reshape(nk,nk)/pi

    omega = data[0,3]

    ax=figure().add_subplot(projection='3d')
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
        # ax.set_title(r"$\chi\prime(q,\omega)$ for $\omega=$"+str(round(omega,5)),fontsize=13)
        ax.set_zlabel(r"$\chi\prime(q,\omega=0)$",fontsize=13,labelpad=10)
    ax.set_facecolor("white")
    # ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    # ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    # ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))

    # ax.set_title(file)

    # show()

    fig2 = figure(figsize=(8,6))
    ax2=fig2.add_subplot(111)

    if interpolate:
        # cset = ax2.imshow(z.reshape(nk,nk),extent=(x.min(),x.max(),y.min(),y.max()),interpolation="spline16",origin='lower',cmap='jet')
        # cset = ax2.imshow(z.reshape(nk,nk),extent=(x.min(),x.max(),y.min(),y.max()),interpolation="sinc",origin='lower',cmap='jet')
        cset = ax2.imshow(z.reshape(nk,nk),extent=(x.min(),x.max(),y.min(),y.max()),interpolation="bessel",origin='lower',cmap='jet')
    else:
        cset = ax2.imshow(z.reshape(nk,nk),extent=(x.min(),x.max(),y.min(),y.max()),interpolation="none",origin='lower',cmap='jet')
    fig2.colorbar(cset)
    ax2.set_xlabel(r"$q_x/\pi$",fontsize=13,labelpad=10)
    ax2.set_ylabel(r"$q_y/\pi$",fontsize=13,labelpad=10)


    if imag:
        ax2.set_title(r"$\chi\prime\prime(q,\omega)$ for $\omega=$"+str(round(omega,5)))
    else:
        ax2.set_title(r"$\chi\prime(q,\omega)$ for $\omega=$"+str(round(omega,5)))

    show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="plot chi(q,w")
    parser.add_argument('--file', dest="file", action="store", default="chiRPA.txt")
    parser.add_argument('--interpolation', action=argparse.BooleanOptionalAction)
    parser.add_argument('--column', dest="column", action="store", default=4)
    input_args = parser.parse_args()

    plotchi(input_args.file, interpolate=input_args.interpolation, column=int(input_args.column))
    
