import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

# add libpy
import os.path
import sys

dirname = os.path.dirname(__file__)
filename = os.path.join(dirname, "../libpy/")
sys.path.append(filename)
from tools import build_z_ds


def main():
    # General parameters
    L0 = 200
    H0 = 20
    N = 512
    g = 9.81
    kp = 10 * np.pi / L0
    lambdaP = 2 * np.pi / kp
    omegap = np.sqrt(g * kp)
    Tp = 2 * np.pi / omegap
    Ndeux = 2e-6
    betaT = 2e-4
    Ts = 20.0
    analysis = "Lagrange"  # Euler or Lagrange
    hchunks = 64

    inpath = "ml_breaking_strat.cuda/out.nc"
    inpath = "N512_40_RE40000_TiniZl_remap1.0_cpu/out.nc"
    ds = xr.open_dataset(
        inpath, chunks={"time": -1, "level": -1, "x": hchunks, "y": hchunks}
    )
    ds = build_z_ds(ds)
    nl = len(ds.level)
    nt = len(ds.time)
    dx = (ds.x[1] - ds.x[0]).values
    dy = dx
    X = np.linspace(-L0 / 2, L0 / 2 + dx, N) / L0
    Y = np.linspace(-L0 / 2, L0 / 2 + dy, N) / L0
    LEVELS = np.arange(0, nl)

    Tm = ds.T.mean(dim=["x", "y"]).compute()
    Zm = ds.z.mean(dim=["x", "y"]).compute()

    fig, ax = plt.subplots(1, 1, figsize=(5, 5), constrained_layout=True, dpi=100)
    for it in range(0, nt):
        if it == 0:
            alpha = 1
            c = "gray"
        else:
            alpha = (it + 1) / (2 * nt)
            c = "b"

        ax.plot((Tm.isel(time=it) - Ts) * 1000, Zm.isel(time=it), alpha=alpha, c=c)
    ax.set_ylim((-H0, 0))
    ax.set_xlim((-25, 1))
    ax.set_xlabel(r"$\overline{T}-T_s$ (mK)")
    ax.set_ylabel(r"$\overline{z}$")
    fig.savefig("leaks_T.pdf")
    plt.show()


if __name__ == "__main__":
    main()
