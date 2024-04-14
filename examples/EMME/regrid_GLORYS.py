import xarray as xr
import xesmf

def regrid_GLORYS(fld, method='nearest_s2d', irange=None, jrange=None):
    roms_grid_fp = '/project/k1090/osipovs/Data/COAWST/EMME/roms_model_grid.nc'
    coords = xr.open_dataset(roms_grid_fp)
    coords = coords.rename({'lon_rho': 'lon', 'lat_rho': 'lat'})
    glorys_grid_fp = '/project/k1090/osipovs/Data/ECMWF/GLORYS12V1/GLOBAL_MULTIYEAR_PHY_001_030/cmems_mod_glo_phy_my_0.083_P1D-m_202112/2017/01/mercatorglorys12v1_gl12_mean_20170101_R20170104.nc'
    gsource = xr.open_dataset(glorys_grid_fp)
    gsource = gsource.rename({'longitude': 'lon', 'latitude': 'lat'})

    if irange is not None:
        gsource = gsource.isel(lon=slice(irange[0],irange[1]))

    if jrange is not None:
        gsource = gsource.isel(lat=slice(jrange[0],jrange[1]))

    regrid = xesmf.Regridder(
        gsource,
        coords,
        method=method,
        periodic=False,
        #filename='regrid_t.nc',
        #reuse_weights=True
    )
    tdest = regrid(fld)
    return tdest
