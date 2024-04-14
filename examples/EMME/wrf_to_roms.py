import numpy as np
import netCDF4
from mpl_toolkits.basemap import Basemap
from scipy.interpolate import interp2d, RegularGridInterpolator
import xarray as xr

import os
os.environ['PYROMS_GRIDID_FILE'] = '/home/osipovs/PycharmProjects/pyroms/examples/EMME/grids.txt'

import pyroms
import pyroms_toolbox
from bathy_smoother import *

'''
Generate ROMS grid from WRF grid. Prescribe Land Mask from Ocean Reanalysis GLORY
gogomamba
mamba activate py38_pyroms
'''

#%%
wrf_geo_fp = "/project/k1090/osipovs/Data/AirQuality/EMME/geo_em.d01.nc"
bathyfile = '/home/osipovs/Data/NOAA/ETOPO/ETOPO_2022_v1_30s_N90W180_bed.nc'
roms_grid_fp = '/project/k1090/osipovs/Data/COAWST/EMME/roms_model_grid.nc'

ncid = netCDF4.Dataset(wrf_geo_fp, "r")

rlat = ncid.variables['XLAT_M'][0,:,:]
rlon = ncid.variables['XLONG_M'][0,:,:]
latv = ncid.variables['XLAT_V'][0,1:-1,:]
lonv = ncid.variables['XLONG_V'][0,1:-1,:]
latu = ncid.variables['XLAT_U'][0,:,1:-1]
lonu = ncid.variables['XLONG_U'][0,:,1:-1]
f = ncid.variables['F'][:]
cosang = ncid.variables['COSALPHA'][:]
sinang = ncid.variables['SINALPHA'][:]
mask = ncid.variables['LANDMASK'][:]
f = ncid.variables['F'][:]
corner_lats = ncid.corner_lats
corner_lons = ncid.corner_lons

map_proj = ncid.MAP_PROJ
if (map_proj == 1):
   my_proj = 'lcc'
   lat_1 = ncid.TRUELAT1
   lat_2 = ncid.TRUELAT2
   lon_0 = ncid.STAND_LON
   llcrnrlon = corner_lons[12]
   llcrnrlat = corner_lats[12]
   urcrnrlon = corner_lons[14]
   urcrnrlat = corner_lats[14]

ncid.close()

map = Basemap(projection=my_proj, lat_1=lat_1, lat_2=lat_2, lon_0=lon_0, llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat, urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat, resolution='h')

# define the 4 corners of the grid
# first point is the top left corner then counter clock wise rotation
lon0=corner_lons[13] ; lat0=corner_lats[13]
lon1=corner_lons[12] ; lat1=corner_lats[12]
lon2=corner_lons[15] ; lat2=corner_lats[15]
lon3=corner_lons[14] ; lat3=corner_lats[14]

#generate the new grid
lonp=np.array([lon0, lon1, lon2, lon3])
latp=np.array([lat0, lat1, lat2, lat3])

# shift data so lons go from 0 to 360 instead of -180 to 180.
lonp = np.where(lonp < 0, lonp+360, lonp)

beta = np.array([1, 1, 1, 1])

Mp, Lp  = rlon.shape
# if you have problems with loading gridgen, then one solution is to install gridgen into the mamba env and to link it
# ln -sf /home/osipovs/apps/mambaforge/envs/py38_pyroms/lib/libgridgen.* .
# https://github.com/ESMG/pyroms/issues/26
hgrd = pyroms.grid.Gridgen(lonp, latp, beta, (Mp+1,Lp+1), proj=map)

lonv, latv = map(hgrd.x_vert, hgrd.y_vert, inverse=True)
hgrd = pyroms.grid.CGrid_geo(lonv, latv, map)
hgrd.lon_rho = np.where(hgrd.lon_rho < 0, hgrd.lon_rho+360, hgrd.lon_rho)

hgrd.mask_rho = (1.0 - mask)
#%% Override the land sea mask. Set the land mask from GLORYS
glorys_fp = '/project/k1090/osipovs/Data/ECMWF/GLORYS12V1/GLOBAL_MULTIYEAR_PHY_001_030/cmems_mod_glo_phy_my_0.083_P1D-m_202112/2017/01/mercatorglorys12v1_gl12_mean_20170101_R20170104.nc'
glorys_ds = xr.open_dataset(glorys_fp).isel(time=0)
glorys_ds['landmask'] = xr.where(glorys_ds.zos.isnull(), 1, 0)

# subset to fit into memory
# irange = (2300, 3300)
# jrange = (900, 1750)
glorys_ds = glorys_ds.isel(longitude=slice(2300, 3300), latitude=slice(900, 1750))
my_interpolating_function = RegularGridInterpolator((glorys_ds.latitude, glorys_ds.longitude), glorys_ds.landmask.to_numpy(), method='nearest')
landmask = my_interpolating_function((hgrd.lat_rho, hgrd.lon_rho))

hgrd.mask_rho = (1.0 - landmask)  # Convert to ROMS watermask.
#%% Post Processing: exclude Azov Sea as it is too shallow
hgrd.mask_rho[419:447, 165:203] = 0
print('Removed Azov Sea')

# debug / diag
# (glorys_ds.zos.isel(time=0)<0).any()
# glorys_ds.landmask.sum()

# import matplotlib.pyplot as plt
# plt.ion()
# plt.show()
# plt.clf()
# glorys_ds.landmask.plot()

#%% prescribe Bathemytry from ETOPO
ncid = netCDF4.Dataset(bathyfile, "r")

lons = ncid.variables['lon'][:]
lats = ncid.variables['lat'][:]
topo = ncid.variables['z'][:]
ncid.close()

# depth positive
topo = -topo

# fix minimum depth
hmin = 5
topo = np.where(topo < hmin, hmin, topo)

# interpolate new bathymetry
lon, lat = np.meshgrid(lons, lats)
my_interpolating_function = RegularGridInterpolator((lats, lons), topo, method='linear')
h = my_interpolating_function((hgrd.lat_rho, hgrd.lon_rho))

# ensure that depth is always deeper than hmin
h = np.where(h < hmin, hmin, h)

#%% Smooth the Bathymetry
# check bathymetry roughness
hgrd.mask_rho = np.reshape(hgrd.mask_rho, (Mp,Lp))
RoughMat = bathy_tools.RoughnessMatrix(h, hgrd.mask_rho)
print('Max Roughness value is: {}'.format(RoughMat.max()))
#h = creep.cslf(h, nan, -200., 200.)
h = np.where(np.isnan(h), 5500.0, h)
hraw = h.copy()

# smooth the raw bathy using the direct iterative method from Martinho and Batteen
# (2006)
rx0_max = 0.35
h = bathy_smoothing.smoothing_Positive_rx0(hgrd.mask_rho, h, rx0_max)

# check bathymetry roughness again
RoughMat = bathy_tools.RoughnessMatrix(h, hgrd.mask_rho)
print('Max Roughness value is: {}'.format(RoughMat.max()))
h = pyroms_toolbox.shapiro_filter.shapiro2(h, 2)

#%% Export the grid to NC
hgrd.h = h
# define vertical grd
Vtrans = 2
theta_s = 7.0
theta_b = 0.1
Tcline = 250
N = 50
vgrd = pyroms.vgrid.s_coordinate_4(h, theta_b, theta_s, Tcline, N, hraw=hraw)
vgrd.h = h   # what the hell??

#ROMS grid
grd_name = 'EMME'
grd = pyroms.grid.ROMS_Grid(grd_name, hgrd, vgrd)

#write grid to netcdf file
pyroms.grid.write_ROMS_grid(grd, filename=roms_grid_fp)

