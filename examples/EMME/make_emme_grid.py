import numpy as np
from mpl_toolkits.basemap import Basemap, shiftgrid
from scipy.interpolate import RegularGridInterpolator
import netCDF4

import pyroms
import xarray as xr

from bathy_smoother import bathy_smoothing
from bathy_smoother.bathy_tools import RoughnessMatrix

'''




USE wrf_to_roms.py INSTEAD





derived from make_YELLOW_grd_v1.py
'''

# Grid was generated by MATLAB, WRF2ROMS.m
roms_grid = '/project/k1090/osipovs/Data/COAWST/EMME/roms_model_grid.nc'
ds = xr.open_dataset(roms_grid)

# ETOPO bathymetry
# ds = xr.open_dataset('/home/osipovs/Data/NOAA/ETOPO/ETOPO1_Bed_g_gdal.grd')
nc = netCDF4.Dataset('/home/osipovs/Data/NOAA/ETOPO/ETOPO_2022_v1_30s_N90W180_bed.nc', 'r')
topo = nc.variables['z'][:]
lons = nc.variables['lon'][:]
lats = nc.variables['lat'][:]

# depth positive
topo = -topo

# fix minimum depth
hmin = 5
topo = np.where(topo < hmin, hmin, topo)

# interpolate new bathymetry
lon, lat = np.meshgrid(lons, lats)
my_interpolating_function = RegularGridInterpolator((lats, lons), topo, method='linear')
h = my_interpolating_function((ds.lat_rho, ds.lon_rho))

# ensure that depth is always deeper than hmin
h = np.where(h < hmin, hmin, h)

# set depth to hmin where masked
idx = np.where(ds.mask_rho == 0)
h[idx] = hmin

# save raw bathymetry
hraw = h.copy()

# check bathymetry roughness
RoughMat = RoughnessMatrix(h, ds.mask_rho)
print('Max Roughness value is: ', RoughMat.max())

# smooth the raw bathy using the direct iterative method from Martinho and Batteen (2006)
rx0_max = 0.35
h = bathy_smoothing.smoothing_Positive_rx0(ds.mask_rho, h, rx0_max)

# check bathymetry roughness again
RoughMat = RoughnessMatrix(h, ds.mask_rho)
print('Max Roughness value is: ', RoughMat.max())

ds['h'][:] = h  # update the grid


ds.to_netcdf(roms_grid + 'depth_updated')

#%%
hgrd = ds

# vertical coordinate
theta_b = 2
theta_s = 7.0
Tcline = 50
N = 30
vgrd = pyroms.vgrid.s_coordinate_4(h, theta_b, theta_s, Tcline, N, hraw=hraw)



# ROMS grid
grd_name = 'EMME'
grd = pyroms.grid.ROMS_Grid(grd_name, hgrd, vgrd)

# write grid to netcdf file
pyroms.grid.write_ROMS_grid(grd, filename=roms_grid + '_pyroms')