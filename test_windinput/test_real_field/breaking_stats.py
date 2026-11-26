import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from scipy import ndimage as ndi


import os.path
import sys

dirname = os.path.dirname(__file__)
filename = os.path.join(dirname, "../../libpy/")
sys.path.append(filename)

from breaking import (
    get_bins,
    label_periodic,
    btrack,
    simple_mapping_parallel,
    track_table,
    length1,
    bbarycenter,
    speed_direction,
)

path = "/home/jacqhugo/basilisk/wiki/sandbox/hugoj/test_windinput/test_real_field/PM_seastate_hold.cuda/"
myfile = "out.nc"
seltime = slice(-4, -1, 1)
L0 = 200
kp = 5 * 2 * np.pi / L0
lvl_surf = -1
threshold = -2 * kp
hchunks = 256

# =========================================
print("Step 0: get the data")
# chunking along 'time' but no 'x' or 'y' as we compute XY diags and we want to // on time
ds = xr.open_dataset(path + myfile, chunks={"time": 1, "x": -1, "y": -1})

delta = (ds.x[1] - ds.x[0]).values
nt = len(ds.time)
bins = get_bins(kp)
bins_center = bins[1:] - (bins[2] - bins[1]) / 2

dmask = simple_mapping_parallel(
    ds.eta,
    ds["u.x"],
    ds["u.y"],
    lvl_surf=lvl_surf,
    delta=delta,
    bins=bins,
    threshold=threshold,
    output_file="crest_mask.nc",
)
