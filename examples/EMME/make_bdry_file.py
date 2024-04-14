import matplotlib
matplotlib.use('Agg')

import subprocess
import os
import sys
import numpy as np
from multiprocessing import Pool

import pyroms
import pyroms_toolbox

import os
os.environ['PYROMS_GRIDID_FILE'] = '/home/osipovs/PycharmProjects/pyroms/examples/EMME/grids.txt'

if 'ESMFMKFILE' not in os.environ:  # os.environ.get('READTHEDOCS') and
    # RTD doesn't activate the env, and esmpy depends on a env var set there
    # We assume the `os` package is in {ENV}/lib/pythonX.X/os.py
    # See conda-forge/esmf-feedstock#91 and readthedocs/readthedocs.org#4067
    print('fixing ESMFMKFILE env variable')
    from pathlib import Path
    os.environ['ESMFMKFILE'] = str(Path(os.__file__).parent.parent / 'esmf.mk')

from remap_bdry import remap_bdry
from remap_bdry_uv import remap_bdry_uv

lst_year = sys.argv[1:]
lst_year = [2017, ]

data_dir = '/project/k1090/osipovs/Data/ECMWF/GLORYS12V1/GLOBAL_MULTIYEAR_PHY_001_030/cmems_mod_glo_phy_my_0.083_P1D-m_202112/2017/*/'
dst_dir = '/project/k1090/osipovs/Data/COAWST/EMME/ROMS/IC_BC/'

#%%

lst_file = []

for year in lst_year:
    # lst = subprocess.getoutput('ls ' + data_dir + 'GLORYS_REANALYSIS_' + year + '-01-*')
    lst = subprocess.getoutput('ls ' + data_dir + '*')
    lst = lst.split()
    lst_file = lst_file + lst

print('Build OBC file from the following file list:\n{}'.format('\n'.join(lst_file)))
print(' ')

irange = (2300, 3300)
jrange = (900, 1750)

# I have processed Jan already
lst_file = lst_file[31:]

src_grd_file=lst_file[0]
src_grd = pyroms_toolbox.Grid_GLORYS.get_nc_Grid_GLORYS(src_grd_file, irange=irange, jrange=jrange)
roms_grid_fp = '/project/k1090/osipovs/Data/COAWST/EMME/roms_model_grid.nc'
dst_grd = pyroms.grid.get_ROMS_grid('EMME') #, grid_file=roms_grid_fp, hist_file=roms_grid_fp)

#%%
def do_file(file):
    zeta = remap_bdry(file, 'zos', src_grd, dst_grd, dst_dir=dst_dir, irange=irange, jrange=jrange)
    dst_grd2 = pyroms.grid.get_ROMS_grid('EMME', zeta=zeta)  # , grid_file=roms_grid_fp, hist_file=roms_grid_fp,
    remap_bdry(file, 'thetao', src_grd, dst_grd2, dst_dir=dst_dir, irange=irange, jrange=jrange)
    remap_bdry(file, 'so', src_grd, dst_grd2, dst_dir=dst_dir, irange=irange, jrange=jrange)
    remap_bdry_uv(file, src_grd, dst_grd2, dst_dir=dst_dir, irange=irange, jrange=jrange)

    # merge file
    bdry_file = dst_dir + file.rsplit('/')[-1][:-3] + '_bdry_' + dst_grd.name + '.nc'

    out_file = dst_dir + file.rsplit('/')[-1][:-3] + '_zos_bdry_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-O', out_file, bdry_file)
    subprocess.check_call(command)
    os.remove(out_file)
    out_file = dst_dir + file.rsplit('/')[-1][:-3] + '_thetao_bdry_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-A', out_file, bdry_file)
    subprocess.check_call(command)
    os.remove(out_file)
    out_file = dst_dir + file.rsplit('/')[-1][:-3] + '_so_bdry_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-A', out_file, bdry_file)
    subprocess.check_call(command)
    os.remove(out_file)
    out_file = dst_dir + file.rsplit('/')[-1][:-3] + '_u_bdry_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-A', out_file, bdry_file)
    subprocess.check_call(command)
    os.remove(out_file)
    out_file = dst_dir + file.rsplit('/')[-1][:-3] + '_v_bdry_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-A', out_file, bdry_file)
    subprocess.check_call(command)
    os.remove(out_file)

#%%
# for file in lst_file:  # serial
#     do_file(file)

processes = 4
p = Pool(processes)
results = p.map(do_file, lst_file)
