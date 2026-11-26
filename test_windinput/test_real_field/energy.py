"""
Just testing the output of the case of Rui with small CFL
CFL_H = 0.05
CFL = 0.05

(he indicated CFL_H=0.01, but it is very long so I tried with 0.05 first)
It seems that the simulation is still unstable, energy grows ...
"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt


SCATTER = False
PLOT_KE = True
PLOT_GPE = True
CASE_PLOT_KEGPE = "hold"
skip = 1
Hskp = 0.247363
L0 = 200
kp = 2 * 5 * np.pi / L0
g = 9.81
P = 0.02
rho = 1000
NT0 = 30
T0 = 2.0 * np.pi / np.sqrt(g * kp)

E0th = rho * g * (Hskp / kp / 4) ** 2 * L0**2

data = {
    "decay": np.loadtxt("PM_decay.cuda/out", skiprows=1),
    # "hold": np.loadtxt("PM_seastate_hold.cuda_dev/out", skiprows=1),
    "hold": np.loadtxt("PM_seastate_hold.cuda_release/out", skiprows=1),
    # "constant": np.loadtxt("PM_seastate.cuda/out", skiprows=1),
}
colors = {"hold": "r", "decay": "k", "constant": "orange"}
colorke = "b"
colorgpe = "g"
lss = {"hold": "--", "decay": "-", "constant": "dashdot"}


time = {}
ke = {}
gpe = {}
E = {}
E0 = {}

for case in data.keys():
    time[case] = data[case][:, 0]
    ke[case] = data[case][:, 1]
    gpe[case] = data[case][:, 2]
    E[case] = ke[case] + gpe[case]
    E0[case] = E[case][0]

print("========================")
print("PM spectrum P = %f" % P)
print("kpHs initial = %f" % Hskp)
print("L0 = %f" % L0)
print("E0th/rho (equipartition Ek Ep) = %f" % (E0th / rho))
print("E0/rho (actual) = %f" % (E["hold"][0] / rho))
print("========================")


fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True, dpi=200)

ax.hlines(1.0, 0, NT0, colors="gray", alpha=0.7)


if SCATTER:
    for case in data.keys():
        ax.scatter(
            time[case][::skip],
            E[case][::skip] / E0,
            marker="x",
            c=colors[case],
            ls=lss[case],
            label=f"E ({case})",
        )
else:
    for case in data.keys():
        if PLOT_KE and case == CASE_PLOT_KEGPE:
            ax.plot(
                time[case][::skip],
                2 * ke[case][::skip] / E0[case],
                color=colorke,
                ls=lss[case],
            )
        if PLOT_GPE and case == CASE_PLOT_KEGPE:
            ax.plot(
                time[case][::skip],
                2 * gpe[case][::skip] / E0[case],
                color=colorgpe,
                ls=lss[case],
            )

        ax.plot(
            time[case][::skip],
            E[case][::skip] / E0[case],
            color=colors[case],
            ls=lss[case],
            label=f"E ({case})",
        )
ax.legend(loc="lower left")
ax.set_ylabel("E/E0")
ax.set_ylim((0.9, 1.1))
ax.set_xlabel("t/T0")
fig.savefig("energy.png")


# ==============================================
print("\n-> Computing current drift")

ds = xr.open_dataset(
    "PM_seastate_hold.cuda_release/out.nc",
    # "PM_seastate_hold.cuda_dev/out.nc",
    chunks={"time": 1, "x": -1, "y": -1},
)

Qu = (ds.h * ds["u.x"]).sum(dim=("level")) / ds.h.sum(dim=("level"))
Qu = Qu.mean(dim=["x", "y"])

fig, ax = plt.subplots(figsize=(3, 3), constrained_layout=True, dpi=200)
ax.plot(ds.time / T0, Qu, c="k")
ax.set_ylabel("mean drift (m/s)")
ax.set_xlabel("t/T0")
fig.savefig("mean_drift_u.png")
plt.show()
