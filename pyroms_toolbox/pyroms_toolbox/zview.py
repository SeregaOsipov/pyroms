import numpy as np
import matplotlib.pyplot as plt
from  matplotlib import cm, colors
from mpl_toolkits.basemap import Basemap
import pyroms
import pyroms_toolbox


def zview(var, tindex, depth, grid, filename=None, \
          cmin=None, cmax=None, clev=None, clb_format='%.2f', \
          fill=False, contour=False, d=4, range=None, fts=None, \
          title=None, clb=True, pal=None, proj='merc', \
          fill_land=False, outfile=None):
    """
    map = zview(var, tindex, depth, grid, {optional switch})

    optional switch:
      - filename         if defined, load the variable from file
      - cmin             set color minimum limit
      - cmax             set color maximum limit
      - clev             set the number of color step
      - fill             use contourf instead of pcolor
      - contour          overlay contour (request fill=True)
      - range            set axis limit
      - fts              set font size (default: 12)
      - title            add title to the plot
      - clb              add colorbar (defaul: True)
      - pal              set color map (default: cm.jet)
      - proj             set projection type (default: merc)
      - fill_land        fill land masked area with gray (defaul: True)
      - outfile          if defined, write figure to file

    plot a constant-z slice of variable var. If filename is provided,
    var must be a string and the variable will be load from the file.
    grid can be a grid object or a gridid. In the later case, the grid
    object correponding to the provided gridid will be loaded.
    If proj is not None, return a Basemap object to be used with quiver
    for example.
    """

    # get grid
    if type(grid).__name__ == 'ROMS_Grid':
        grd = grid
    else:
        grd = pyroms.grid.get_ROMS_grid(grid)


    # get variable
    if filename == None:
        var = var
    else:
        data = pyroms.io.Dataset(filename)

        var = data.variables[var]

    Np, Mp, Lp = grd.vgrid.z_r[0,:].shape

    if tindex == -1:
        assert len(var.shape) == 3, 'var must be 3D (no time dependency).'
        N, M, L = var.shape
    else:
        assert len(var.shape) == 4, 'var must be 4D (time plus space).'
        K, N, M, L = var.shape

    # determine where on the C-grid these variable lies
    if N == Np and M == Mp and L == Lp:
        Cpos='rho'
        mask = grd.hgrid.mask_rho

    if N == Np and M == Mp and L == Lp-1:
        Cpos='u'
        mask = grd.hgrid.mask_u

    if N == Np and M == Mp-1 and L == Lp:
        Cpos='v'
        mask = grd.hgrid.mask_v

    # get constante-z slice
    if tindex == -1:
        var = var[:,:,:]
    else:
        var = var[tindex,:,:,:]

    depth = -abs(depth)
    if fill == True:
        zslice, lon, lat = pyroms.tools.zslice(var, depth, grd, \
                                  Cpos=Cpos)
    else:
        zslice, lon, lat = pyroms.tools.zslice(var, depth, grd, \
                                  Cpos=Cpos, vert=True)

    # plot
    if cmin is None:
        cmin = zslice.min()
    else:
        cmin = float(cmin)

    if cmax is None:
        cmax = zslice.max()
    else:
        cmax = float(cmax)

    if clev is None:
        clev = 100.
    else:
        clev = float(clev)

    dc = (cmax - cmin)/clev ; vc = np.arange(cmin,cmax+dc,dc)

    if pal is None:
        pal = cm.jet
    else:
        pal = pal

    if fts is None:
        fts = 12
    else:
        fts = fts

    #pal.set_over('w', 1.0)
    #pal.set_under('w', 1.0)
    #pal.set_bad('w', 1.0)

    pal_norm = colors.BoundaryNorm(vc,ncolors=256, clip = False)

    if range is None:
        lon_min = lon.min()
        lon_max = lon.max()
        lon_0 = (lon_min + lon_max) / 2.
        lat_min = lat.min()
        lat_max = lat.max()
        lat_0 = (lat_min + lat_max) / 2.
    else:
        lon_min = range[0]
        lon_max = range[1]
        lon_0 = (lon_min + lon_max) / 2.
        lat_min = range[2]
        lat_max = range[3]
        lat_0 = (lat_min + lat_max) / 2.

    # clear figure
    #plt.clf()

    # Arctic 4
    if proj is 'stere':
        map = Basemap(projection='stere', lat_0=90, lon_0=180, llcrnrlon=-210,
	        llcrnrlat=40, urcrnrlon=-50, urcrnrlat=50, resolution='h')
        x, y = list(map(lon,lat))
    elif proj is not None:
        map = Basemap(projection=proj, llcrnrlon=lon_min, llcrnrlat=lat_min, \
                 urcrnrlon=lon_max, urcrnrlat=lat_max, lat_0=lat_0, lon_0=lon_0, \
                 resolution='h', area_thresh=5.)
        #map = pyroms.utility.get_grid_proj(grd, type=proj)
        x, y = list(map(lon,lat))

    if fill_land and proj is not None:
        # fill land and draw coastlines
        map.drawcoastlines()
        map.fillcontinents(color='grey')
    else:
        if proj is not None:
            Basemap.pcolor(map, x, y, mask, vmin=-2, cmap=cm.gray)
            pyroms_toolbox.plot_coast_line(grd, map)
        else:
            plt.pcolor(lon, lat, mask, vmin=-2, cmap=cm.gray)
            pyroms_toolbox.plot_coast_line(grd)

    if fill:
        if proj is not None:
            cf = Basemap.contourf(map, x, y, zslice, vc, cmap = pal, \
                                  norm = pal_norm)
        else:
            cf = plt.contourf(lon, lat, zslice, vc, cmap = pal, \
                              norm = pal_norm)
    else:
        if proj is not None:
            cf = Basemap.pcolor(map, x, y, zslice, cmap = pal, norm = pal_norm)
        else:
            cf = plt.pcolor(lon, lat, zslice, cmap = pal, norm = pal_norm)

    if clb:
        clb = plt.colorbar(cf, fraction=0.075,format=clb_format)
        for t in clb.ax.get_yticklabels():
            t.set_fontsize(fts)

    if contour:
        if not fill:
            raise Warning('Please run again with fill=True to overlay contour.')
        else:
            if proj is not None:
                Basemap.contour(map, x, y, zslice, vc[::d], colors='k', linewidths=0.5, linestyles='solid')
            else:
                plt.contour(lon, lat, zslice, vc[::d], colors='k', linewidths=0.5, linestyles='solid')

    if proj is None and range is not None:
        plt.axis(range)


    if title is not None:
            plt.title(title, fontsize=fts+4)

    if proj is not None:
        map.drawmeridians(np.arange(lon_min,lon_max, (lon_max-lon_min)/5.001), \
                          labels=[0,0,0,1], fmt='%.1f')
        map.drawparallels(np.arange(lat_min,lat_max, (lat_max-lat_min)/5.001), \
                          labels=[1,0,0,0], fmt='%.1f')

    if outfile is not None:
        if outfile.find('.png') != -1 or outfile.find('.svg') != -1 or \
           outfile.find('.eps') != -1:
            print('Write figure to file', outfile)
            plt.savefig(outfile, dpi=200, facecolor='w', edgecolor='w', \
                        orientation='portrait')
        else:
            print('Unrecognized file extension. Please use .png, .svg or .eps file extension.')


    if proj is None:
        return
    else:
        return map
