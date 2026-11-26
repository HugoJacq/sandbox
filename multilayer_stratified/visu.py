import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import matplotlib as mpl

# loading data
ds = xr.open_dataset("ml_convection.cuda/out.nc")
dsm = ds.mean(dim=["x", "y"])
nl = len(ds.level)
nt = len(ds.time)
N = len(ds.x)
h_ = 50.0
T0 = 20.0
dz = h_ / nl

# > T_profile
colors = mpl.colormaps["plasma"](np.linspace(0.1, 1, nt))
fig, ax = plt.subplots(1, 1, figsize=(7, 3), constrained_layout=True, dpi=100)
for it in range(nt):
    ax.plot(dsm.T.isel(time=it) / T0, ds.level * dz, c=colors[it])
ax.set_xlabel("T/T0")
ax.set_ylabel("")
fig.savefig("T_profile.png")


plt.show()
