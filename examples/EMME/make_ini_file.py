import subprocess
import os
import sys
import subprocess
import numpy as np
from datetime import datetime
import matplotlib
matplotlib.use('Agg')

import os
os.environ['PYROMS_GRIDID_FILE'] = '/home/osipovs/PycharmProjects/pyroms/examples/EMME/grids.txt'

if 'ESMFMKFILE' not in os.environ:  # os.environ.get('READTHEDOCS') and
    # RTD doesn't activate the env, and esmpy depends on a env var set there
    # We assume the `os` package is in {ENV}/lib/pythonX.X/os.py
    # See conda-forge/esmf-feedstock#91 and readthedocs/readthedocs.org#4067
    print('fixing ESMFMKFILE env variable')
    from pathlib import Path
    os.environ['ESMFMKFILE'] = str(Path(os.__file__).parent.parent / 'esmf.mk')

import pyroms
import pyroms_toolbox
from remap_clm import remap_clm
from remap_clm_uv import remap_clm_uv

lst_year = sys.argv[1:]
lst_year = [2017, ]

'''
GLORYS12V1 reanalysis page 
https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030/description

Download using python:
cd /project/k1090/osipovs/Data/ECMWF/GLORYS12V1
mamba activate cmc-beta
copernicus-marine get -i cmems_mod_glo_phy_my_0.083_P1D-m --filter *2017/01*
'''

data_dir = '/project/k1090/osipovs/Data/ECMWF/GLORYS12V1/GLOBAL_MULTIYEAR_PHY_001_030/cmems_mod_glo_phy_my_0.083_P1D-m_202112/2017/01/'
dst_dir='/project/k1090/osipovs/Data/COAWST/EMME/ROMS/IC_BC/'

lst_file = []

for year in lst_year:
    year = str(year)
    # lst = subprocess.getoutput('ls ' + data_dir + '*GLBy*' + year + '*')
    lst = subprocess.getoutput('ls ' + data_dir + '*0101*.nc')
#   lst = subprocess.getoutput('ls ' + data_dir + '*GLBy*' + year + '*12.nc')
    lst = lst.split()
    lst_file = lst_file + lst

print('Build CLM file from the following file list: {}'.format('\n'.join(lst_file)))
print(' ')

irange = (2300, 3300)
jrange = (900, 1750)

src_grd_file=lst_file[0]
src_grd = pyroms_toolbox.Grid_GLORYS.get_nc_Grid_GLORYS(src_grd_file, irange=irange, jrange=jrange)
roms_grid_fp = '/project/k1090/osipovs/Data/COAWST/EMME/roms_model_grid.nc'
# dst_grd = pyroms.grid.get_ROMS_grid('EMME', grid_file=roms_grid_fp, hist_file=roms_grid_fp)
dst_grd = pyroms.grid.get_ROMS_grid('EMME')
#%%
for file in lst_file:
# remap
    zeta = remap_clm(file, 'zos', src_grd, dst_grd, dst_dir=dst_dir, irange=irange, jrange=jrange)
    dst_grd = pyroms.grid.get_ROMS_grid('EMME', zeta=zeta)  # , grid_file=roms_grid_fp, hist_file=roms_grid_fp
    remap_clm(file, 'thetao', src_grd, dst_grd, dst_dir=dst_dir, irange=irange, jrange=jrange)
    remap_clm(file, 'so', src_grd, dst_grd, dst_dir=dst_dir, irange=irange, jrange=jrange)
    remap_clm_uv(file, src_grd, dst_grd, dst_dir=dst_dir, irange=irange, jrange=jrange)

# merge file
    clim_file = dst_dir + file.rsplit('/')[-1][:-3] + '_clim_' + dst_grd.name + '.nc'

    out_file = dst_dir + file.rsplit('/')[-1][:-3] + '_zos_clim_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-O', out_file, clim_file)
    print(command)
    subprocess.check_call(command)
    os.remove(out_file)
    out_file = dst_dir + file.rsplit('/')[-1][:-3] + '_thetao_clim_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-A', out_file, clim_file)
    print(command)
    subprocess.check_call(command)
    os.remove(out_file)
    out_file = dst_dir + file.rsplit('/')[-1][:-3] + '_so_clim_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-A', out_file, clim_file)
    print(command)
    subprocess.check_call(command)
    os.remove(out_file)
    out_file = dst_dir + file.rsplit('/')[-1][:-3] + '_u_clim_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-A', out_file, clim_file)
    print(command)
    subprocess.check_call(command)
    os.remove(out_file)
    out_file = dst_dir + file.rsplit('/')[-1][:-3] + '_v_clim_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-A', out_file, clim_file)
    print(command)
    subprocess.check_call(command)
    os.remove(out_file)

#%%

print("DONE")

