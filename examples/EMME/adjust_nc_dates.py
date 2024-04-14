import xarray as xr
import datetime as dt
import pandas as pd
'''
A bunch of scripts to manually adjust the time steps in NC
'''

data_dir = '/project/k1090/osipovs/Data/ECMWF/GLORYS12V1/GLOBAL_MULTIYEAR_PHY_001_030/cmems_mod_glo_phy_my_0.083_P1D-m_202112/2017/01/'
dst_dir='/project/k1090/osipovs/Data/COAWST/EMME/ROMS/IC_BC/'

#%% prep cold start
fp_in = dst_dir + 'mercatorglorys12v1_gl12_mean_20170101_R20170104_clim_EMME.nc'
fp_out = dst_dir + 'mercatorglorys12v1_gl12_mean_20170101_R20170104_clim_EMME_cold_start.nc'
ds = xr.open_dataset(fp_in)

keys = 'v,u,vbar,ubar,zeta'.split(',')
for key in keys:
    print('Processing {}'.format(key))
    ds[key][:] = 0

# date is centered midday. correct to 00:00
ds['ocean_time'] = ds['ocean_time'].copy(data=[dt.datetime(2017, 1, 1, 0, 0), ])

ds.to_netcdf(fp_out)

ds.to_netcdf('/home/osipovs/Temp/mercatorglorys12v1_gl12_mean_20170101_R20170104_clim_EMME_cold_start.nc')

#%% update the time in a cold spinup restart: 2017-01-02 -> 2017-01-01
ds = xr.open_dataset('/project/k1090/osipovs/Data/COAWST/EMME/ROMS/ocean_rst.nc.cold_spinup_1day')
ds = xr.open_dataset('/scratch/osipovs/Models/COAWST/run_emme_2017/cold_start/ocean_rst.nc')
ds['ocean_time'] = ds['ocean_time'].copy(data=[dt.datetime(2017, 1, 1, 0, 0), ])
ds.to_netcdf('/scratch/osipovs/Models/COAWST/run_emme_2017/cold_start/ocean_rst.nc.cold_spinup_1day')

#%% move first date 12 hours back to provide 2017-01-01 exactly
ds = xr.open_dataset('/project/k1090/osipovs/Data/COAWST/EMME/ROMS/IC_BC/mercatorglorys12v1_gl12_mean_2017_bdry_EMME.nc')

dates=[]
for d in ds['ocean_time'].values:
    print(d)
    dates += [pd.to_datetime(d),]
dates[0] = pd.to_datetime('2017-01-01')

ds['ocean_time'] = ds['ocean_time'].copy(data=dates)

# ds.to_netcdf('/project/k1090/osipovs/Data/COAWST/EMME/ROMS/IC_BC/mercatorglorys12v1_gl12_mean_201701_bdry_EMME.nc_time_fixed')
ds.to_netcdf('/home/osipovs/Temp/mercatorglorys12v1_gl12_mean_2017_bdry_EMME.nc_time_fixed')

