"""
# Breaking front detection for a multilayer simulation

This script tests the breaking detection method from Wu et al. 2023
"""


import xarray as xr


# add libpy
dirname = os.path.dirname(__file__)
filename = os.path.join(dirname, '../libpy/')
sys.path.append( filename )
from breaking import *




