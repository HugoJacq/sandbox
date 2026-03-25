import numpy as np
import xarray as xr
import matplotlib.pyplot as plt


g = 9.81
L0 = 200.0 # m
kp = 10  * np.pi / L0
wp = np.sqrt(g*kp)
fp = wp/(2*np.pi)
Tp = 1/fp
N = 1024
L = 200.
dpi=200

# My data
path15 = './'
path30 = "./figs_N10_P0.02_L30/"
filedata = "data.nc"
file30 = path30 + filedata
file15 = path15 + filedata
Jpath = "data_Jiarong/"

file_grads = ['dudz.nc','dudy.nc','dudx.nc','dvdz.nc','dvdy.nc','dvdx.nc','dwdz.nc','dwdy.nc', 'dwdx.nc']
file_omeg = ['omegaxp.nc','omegayp.nc','omegazp.nc','enstrophy.nc']
file_diss = ['epsilon.nc']# We build a dict, easier to process
d_data = {"L15":xr.open_mfdataset([file15]+file_grads+file_omeg+file_diss),
          "L30":xr.open_mfdataset([file30] + [path30+file for file in (file_grads+file_omeg+file_diss)])}

if True:
    print('\nProfiles ')
    time_=120 # s
    clrs = ['b','g']
    # Jiarong's data
    Jux_lagr = np.loadtxt(Jpath+'2025_fig7a_layer.txt',skiprows=1,delimiter=",")
    Jens_lagr = np.loadtxt(Jpath+'2025_fig7b_layer.txt',skiprows=1,delimiter=",")
    Jdiss_lagr = np.loadtxt(Jpath+'2025_fig7c_layer.txt',skiprows=1,delimiter=",")

    fig, ax = plt.subplots(1,3,figsize = (9,3),constrained_layout=True,dpi=dpi)
    ax[0].plot(Jux_lagr[:,0], Jux_lagr[:,1], c='gray', label='paper', marker='x', markerfacecolor='None')
    ax[1].semilogx(Jens_lagr[:,0], Jens_lagr[:,1], c='gray', label='paper', marker='x', markerfacecolor='None')
    ax[2].semilogx(Jdiss_lagr[:,0], Jdiss_lagr[:,1], c='gray', label='paper', marker='x', markerfacecolor='None' )

    for k,name in enumerate(d_data.keys()):
        ds1=d_data[name].sel(time=time_)
        color = clrs[k]

        z_lagr = ds1.z.mean(['x','y'])
        ux_lagr = ds1['u.x'].mean(['x','y'])
        ens_lagr = ds1['enstrophy'].mean(['x','y'])
        diss_lagr = ds1['epsilon'].mean(['x','y'])

        # U
        ax[0].plot(ux_lagr, z_lagr, c=color, ls='-', label=name,  marker='s', markerfacecolor='None')
        ax[0].set_xlim([0,0.2])

        # Vorticity
        ax[1].semilogx(ens_lagr, z_lagr, c=color, ls='-', label=name,  marker='s', markerfacecolor='None')
        ax[1].set_xlim([1e-4,1])

        # Dissipation
        ax[2].semilogx(diss_lagr, z_lagr, c=color, ls='-', label=name, marker='s', markerfacecolor='None')
        ax[2].set_xlim([1e-3,10])

    ax[0].set_xlabel('<u> (m/s)')
    ax[1].set_xlabel('<enstrophy>')
    ax[2].set_xlabel('<diss>')
    ax[1].set_title(f'time = {time_} s ($w_p t=$'+str(int(wp*time_))+')')
    for axe in ax:
        axe.legend(frameon=False)
        axe.set_ylim([-10,0])
        axe.set_ylabel('z (m)')
    fig.savefig('comp_avg_profiles_L15_L30.svg')
