import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pathlib import Path
import colorcet as cc

# add libpy
import os.path
import sys

dirname = os.path.dirname(__file__)
filename = os.path.join(dirname, "../libpy/")
sys.path.append(filename)
from data_reader import read_bas_data


# General parameters
L0 = 200
H0 = 50
N = 1024
kp = 10 * np.pi / L0
lambdaP = 2 * np.pi / kp
Ndeux = 2e-6
g = 9.81
betaT = 2e-4
Ts = 20.0
hchunks = 256

# inpath = "N512_nl30_0.000002_Tinizl/"
# infile = "out.nc"
# inpath = "~/work/BasiLagrangian/models/N1024/"
# infile = "out954.nc"
inpath = "N1024_nl60_cuda/"
infile = "out551.nc"
outfile = "L2"
lsfx = -1
it = -1

# opening file
ds, grid = read_bas_data(
    inpath + infile, chunks={"level": -1, "x": hchunks, "y": hchunks}
)

ds = ds.isel(time=it)  # for this specific file I select a snapshot only

dx = (ds.x[1] - ds.x[0]).values
dy = dx
X = np.linspace(-L0 / 2, L0 / 2 + dx, N) / L0
Y = np.linspace(-L0 / 2, L0 / 2 + dy, N) / L0
Zm = ds.z.mean(dim=["x", "y"])

wet_z = ds.z.isel(level=-1).min().values
print("Depth always under the surface: %f (m)" % wet_z)


def triple_avg(ds, var, maindir="x"):
    u_m = ds[var].mean(dim=["x", "y"])
    u_w = ds[var].mean(dim=["y"])  # assuming that main direction is 'x'
    u_f = ds[var] - u_m - u_w
    return u_m, u_w, u_f


# Triple decomposition following Tsai 2023
um = triple_avg(ds, "u.x")
vm = triple_avg(ds, "u.y")
wm = triple_avg(ds, "u.z")
Tm = triple_avg(ds, "T")

varnames = ["T"]  # ["u.x", "u.y", "u.z", "T"]
nicenames = ["T", "u", "v", "w"]

for k, var in enumerate(varnames):  # ["u.x", "u.y", "u.z", "T"]
    um = triple_avg(ds, var)
    # profile
    fig, ax = plt.subplots(1, 1, figsize=(5, 7), constrained_layout=True, dpi=100)
    ax.plot(um[0], Zm, c="k")
    ax.set_xlabel(r"$\overline{%s}$" % nicenames[k])
    ax.set_ylabel("z")
    # original surface
    fig, ax = plt.subplots(1, 1, figsize=(9, 7), constrained_layout=True, dpi=100)
    s = ax.pcolormesh(
        X, Y, ds[var].isel(level=-1), cmap="Greys_r"
    )  # , vmin=-2, vmax=2)
    plt.colorbar(s, ax=ax, label=r"$\overline{%s}$" % nicenames[k])
    ax.set_xlabel("x/L0")
    ax.set_ylabel("y/L0")
    # fluctuations surface
    fig, ax = plt.subplots(1, 1, figsize=(9, 7), constrained_layout=True, dpi=100)
    s = ax.pcolormesh(X, Y, um[2].isel(level=-1), cmap="Greys_r")  # , vmin=-2, vmax=2)
    plt.colorbar(s, ax=ax, label=r"%s'" % nicenames[k])
    ax.set_xlabel("x/L0")
    ax.set_ylabel("y/L0")
    # wave
    fig, ax = plt.subplots(1, 1, figsize=(9, 4), constrained_layout=True, dpi=100)
    s = ax.pcolormesh(X, Zm, um[1], cmap="Greys_r")  # , vmin=-2, vmax=2)
    plt.colorbar(s, ax=ax, label=r"%s'" % nicenames[k])
    ax.set_xlabel("x/L0")
    ax.set_ylabel("y/L0")


"""
Notes: méthode tirée de Tsai 2023
* méthode prometteuse pour les vitesses
* ne fonctionne pas tel quelle pour la température (gradients H << gradient vertical)
"""

plt.show()
