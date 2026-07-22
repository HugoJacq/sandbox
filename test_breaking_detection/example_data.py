"""
This script extract a snapshot of the first layer of a simulation of a breaking wave field

The fields extracted are: eta, u.x, u.y
"""

import xarray as xr
import matplotlib.pyplot as plt

path = '/home/jacqhugo/basilisk/wiki/sandbox/hugoj/reproducing_jiarongs_plots/N10_P0.02_L30/out.nc'
data = xr.open_dataset(path)

data = data.isel(level=-1, time=-1).drop_vars(['zb','h','w'])
data.encoding['unlimited_dims'] = {}
print(data.encoding)
print(data)
data.to_netcdf('example_file.nc')
