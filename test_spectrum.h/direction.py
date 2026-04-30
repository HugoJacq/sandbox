""" 

This script plots the surface elevation, surface currents and the
associated wave spectrum for differents directions 

How to use this script:
make direction.tst
python3 direction.py

Note: Make sure to have Ndir equal in both direction.c and direction.py
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os.path
import sys

# add libpy
dirname = os.path.dirname(__file__)
filename = os.path.join(dirname, '../libpy/')
sys.path.append( filename )
from fftlib import get_spec_1D, get_spec_2D, azimuthal_integral

# some parameters
L=200.
kp=10*np.pi/L
omega=np.sqrt(9.81*kp)
cphase = omega/kp
N_grid = 256 # 1024
Ndir = 4
x = np.linspace(-L/2,L/2,N_grid)
y = x

# files
file = "./direction/eta_C"
fileu = "./direction/u_C"
filev = "./direction/v_C"

# data from C code
raw_eta_C = np.fromfile(file)
eta_C = np.reshape(raw_eta_C, (N_grid,N_grid,Ndir))
u_C = np.reshape(np.fromfile(fileu), (N_grid,N_grid,Ndir))
v_C = np.reshape(np.fromfile(filev), (N_grid,N_grid,Ndir))

Fk = []
krs = []
for k in range(Ndir):
    print(f"direction is: {k}pi/4")

    """
    Computing the azimuthal integrated spectrum from eta
    """
    kr_C, F_eta_C = get_spec_1D(eta_C[:,:,k], eta_C[:,:,k], L/N_grid)
    krs.append(kr_C)
    Fk.append(F_eta_C)

    """
    Variance check (should all be equal !) 
    """

    kr = krs[k]
    F = Fk[k]
    dkr = kr[1] - kr[0]
    print("var(eta)=%f" %(np.nansum(F*dkr)))

    """
    Computing mean currents at z=eta
    """
    meanU = np.mean(u_C[:,:,k])
    meanV = np.mean(v_C[:,:,k])
    print(f"<U>={meanU}, <V>={meanV}")
    print(f"<dir>={np.arctan(meanV/meanU)*180/np.pi} (deg)")

    """
    Plotting the surface fields
    """

    fig, ax = plt.subplots(1,1,figsize = (7,5),constrained_layout=True,dpi=100)
    s=ax.pcolormesh(x/L, y/L, eta_C[:,:,k].T*kp, cmap='bwr')
    # ax.annotate("", xytext=(0, 0), xy=(meanU*10, meanV*10),
    #         arrowprops=dict(arrowstyle="->", lw=0.5))
    ax.set_xlabel(r'$x$/$L$')
    ax.set_ylabel(r'$y$/$L$')
    plt.colorbar(s,ax=ax,label=r'$\eta kp$')
    fig.savefig(f'eta_dir{k}.png')
    # U
    fig, ax = plt.subplots(1,1,figsize = (7,5),constrained_layout=True,dpi=100)
    s=ax.pcolormesh(x/L, y/L, u_C[:,:,k].T/cphase, cmap='bwr',  vmin=-0.3,vmax=0.3)
    # ax.annotate("", xytext=(0, 0), xy=(meanU*10, meanV*10),
    #         arrowprops=dict(arrowstyle="->", lw=0.5))
    ax.set_xlabel(r'$x$/$L$')
    ax.set_ylabel(r'$y$/$L$')
    plt.colorbar(s,ax=ax,label=r'u/$c_p$')
    fig.savefig(f'u_dir{k}.png')
    # V
    fig, ax = plt.subplots(1,1,figsize = (7,5),constrained_layout=True,dpi=100)
    s=ax.pcolormesh(x/L, y/L, v_C[:,:,k].T/cphase, cmap='bwr', vmin=-0.3,vmax=0.3)
    # ax.annotate("", xytext=(0, 0), xy=(meanU*10, meanV*10),
    #         arrowprops=dict(arrowstyle="->", lw=0.5))
    ax.set_xlabel(r'$x$/$L$')
    ax.set_ylabel(r'$y$/$L$')
    plt.colorbar(s,ax=ax,label=r'v/$c_p$')
    fig.savefig(f'v_dir{k}.png')

    """
    Plotting the 2D spectrum
    """
    kx,ky,spec2D = get_spec_2D(eta_C[:,:,k].T, eta_C[:,:,k].T, Delta=L/N_grid)
    fig, ax = plt.subplots(1,1,figsize = (7,5),constrained_layout=True,dpi=100)
    s=ax.pcolormesh(kx*2*np.pi/kp,ky*2*np.pi/kp,(spec2D/(2*np.pi))*kp**3,cmap='viridis', shading='nearest', 
                    norm=matplotlib.colors.LogNorm(vmin=1e-7, vmax=1e-2))
    ax.set_xlabel('kx/kp')
    ax.set_ylabel('ky/kp')
    ax.set_xlim([-5,5])
    ax.set_ylim([-5,5])
    plt.colorbar(s,ax=ax,label=r'$F(k_x,k_y) kp^{3}$')
    ax.set_title(rf'direction={k}$\pi$/4')
    fig.savefig(f'F_kxky_dir{k}.png')






