"""Some tools to generate AZMP ResDocs

Contains:
- build_bottomT_climato(infile, year_lims)

----------

Atlantic Zone Monitoring Program @NAFC:
https://azmp-nl.github.io/

"""

__author__ = 'Frederic.Cyr@dfo-mpo.gc.ca'
__version__ = '0.1'

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import warnings
#from sys import version_info
# read/write tools
import netCDF4
import h5py
import xarray as xr
import pandas as pd
# maps
os.environ['PROJ_LIB'] = '/home/cyrf0006/anaconda3/share/proj'
#from mpl_toolkits.basemap import Basemap
import cartopy.crs as ccrs
import cartopy.feature as cfeature
# interpolation tools
from scipy.interpolate import griddata
from scipy.interpolate import interp1d
from scipy.interpolate.interpnd import _ndim_coords_from_arrays
from scipy.spatial import cKDTree
# Shaping tools
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from shapely.ops import unary_union
## AZMP custom imports
sys.path.append('/home/jcoyne/Documents/AZMP-NL_python-toolbox/python-toolbox/azmp_modules')
import azmp_utils as azu
## for scorecards
import unicodedata
from matplotlib.colors import from_levels_and_colors
# For shapefiles
import shapefile 
import cmocean as cmo


def is_number(s):
    '''
    Used for differentiate numbers from letters in scorecards.
    https://www.pythoncentral.io/how-to-check-if-a-string-is-a-number-in-python-including-unicode/
    
    '''
    try:
        float(s)
        return True
    except ValueError:
        pass 
    try:
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass 
    return False

def add_path(PATH):
    """ Since this module uses (e.g.) bathymetry data, the path to the file must be specified if not already permanent.

    Usage ex: (for bathymetry data)
    import azmp_report_tools as az_r
    az_r.add_path('/home/cyrf0006/data/GEBCO/')
    az_r.add_path('/home/cyrf0006/data/dev_database/')

    ** Turns out to be a useless function...
    """

    sys.path.append(PATH)  # or .insert(0, YOUR_PATH) may give higher priority


def bottom_temperature(season, year, zmin=0, zmax=1000, dz=5, clim_fill=True, netcdf_path='/home/cyrf0006/data/dev_database/netCDF/', climato_file=''):
    """
    Bottom temperature maps for AZMP ResDoc

    This function is adapted from azmp_bottomT.py
    To generate bottom climato:
    import numpy as np
    import azmp_utils as azu
    dc = .1
    lonLims = [-60, -43] # fish_hab region
    latLims = [39, 56]
    lonLims = [-60, -45] # FC AZMP report region
    latLims = [42, 56]
    lonLims = [-63, -45] # include 2H in above
    latLims = [42, 58]
    lon_reg = np.arange(lonLims[0]+dc/2, lonLims[1]-dc/2, dc)
    lat_reg = np.arange(latLims[0]+dc/2, latLims[1]-dc/2, dc)
    azu.get_bottomT_climato('/home/cyrf0006/data/dev_database/netCDF/*.nc', lon_reg, lat_reg, season='fall', h5_outputfile='Tbot_climato_fall_0.10.h5') 

    * see: /home/cyrf0006/AZMP/state_reports/bottomT

    usage ex:
    >> import azmp_report_tools as azrt
    >> azrt.bottom_temperature(season='spring', year='2019'):

    ** For NAFO SA4:
    >> azrt.bottom_temperature(season='summer', year='2019', climato_file='Tbot_climato_SA4_summer_0.10.h5'):

    Frederic.Cyr@dfo-mpo.gc.ca - January 2020
    """

    if len(climato_file) == 0:
        if season=='spring':
            climato_file='operation_files/Tbot_climato_spring_0.10.h5'
        elif season=='fall':
            climato_file='operation_files/Tbot_climato_fall_0.10.h5'
        elif season=='summer':
            climato_file='operation_files/Tbot_climato_summer_0.10.h6'

    year_file = netcdf_path + str(year) + '.nc'


    ## ---- Load Climato data ---- ##
    print('Load ' + climato_file)
    h5f = h5py.File(climato_file, 'r')
    Tbot_climato = h5f['Tbot'][:]
    lon_reg = h5f['lon_reg'][:]
    lat_reg = h5f['lat_reg'][:]
    Zitp = h5f['Zitp'][:]
    h5f.close()

    ## ---- Derive some parameters ---- ##
    lon_0 = np.round(np.mean(lon_reg))
    lat_0 = np.round(np.mean(lat_reg))
    lonLims = [lon_reg[0], lon_reg[-1]]
    latLims = [lat_reg[0], lat_reg[-1]]
    lon_grid, lat_grid = np.meshgrid(lon_reg,lat_reg)
    dc = np.diff(lon_reg[0:2])

    ## ---- NAFO divisions ---- ##
    nafo_div = azu.get_nafo_divisions()

    ## ---- Get CTD data --- ##
    print('Get ' + year_file)
    ds = xr.open_dataset(year_file)

    # Selection of a subset region
    ds = ds.where((ds.longitude>lonLims[0]) & (ds.longitude<lonLims[1]), drop=True)
    ds = ds.where((ds.latitude>latLims[0]) & (ds.latitude<latLims[1]), drop=True)
    # Select time (save several options here)
    if season == 'summer':
        #ds = ds.sel(time=ds['time.season']=='JJA')
        ds = ds.sel(time=((ds['time.month']>=6)) & ((ds['time.month']<=9)))
    elif season == 'spring':
        #ds = ds.sel(time=ds['time.season']=='MAM')
        ds = ds.sel(time=((ds['time.month']>=4)) & ((ds['time.month']<=6)))
    elif season == 'fall':
        #ds = ds.sel(time=ds['time.season']=='SON')
        ds = ds.sel(time=((ds['time.month']>=9)) & ((ds['time.month']<=12)))
    else:
        print('!! no season specified, used them all! !!')

    # Isolate for temperatures of interest
    ds = ds.sel(level=ds['level']<zmax)
    da_temp = ds['temperature']
    lons = np.array(ds.longitude)
    lats = np.array(ds.latitude)
    #To Pandas Dataframe
    df_temp = da_temp.to_pandas()
    # Remove empty columns
    idx_empty_rows = df_temp.isnull().all(1).values.nonzero()[0]
    df_temp = df_temp.dropna(axis=0,how='all')
    lons = np.delete(lons,idx_empty_rows)
    lats = np.delete(lats,idx_empty_rows)
    print(' -> Done!')

    ## --- fill 3D cube --- ##
    print('Fill regular cube')
    z = df_temp.columns.values
    V = np.full((lat_reg.size, lon_reg.size, z.size), np.nan)

    # Aggregate on regular grid
    for i, xx in enumerate(lon_reg):
        for j, yy in enumerate(lat_reg):
            idx = np.where((lons>=xx-dc/2) & (lons<xx+dc/2) & (lats>=yy-dc/2) & (lats<yy+dc/2))
            tmp = np.array(df_temp.iloc[idx].mean(axis=0))
            idx_good = np.argwhere((~np.isnan(tmp)) & (tmp<30))
            if np.size(idx_good)==1:
                V[j,i,:] = np.array(df_temp.iloc[idx].mean(axis=0))
            elif np.size(idx_good)>1: # vertical interpolation between pts
                interp = interp1d(np.squeeze(z[idx_good]), np.squeeze(tmp[idx_good]))
                idx_interp = np.arange(int(idx_good[0]),int(idx_good[-1]+1))
                V[j,i,idx_interp] = interp(z[idx_interp]) # interpolate only where possible (1st to last good idx)

    # horizontal interpolation at each depth
    lon_grid, lat_grid = np.meshgrid(lon_reg,lat_reg)
    lon_vec = np.reshape(lon_grid, lon_grid.size)
    lat_vec = np.reshape(lat_grid, lat_grid.size)
    for k, zz in enumerate(z):
        # Meshgrid 1D data (after removing NaNs)
        tmp_grid = V[:,:,k]
        tmp_vec = np.reshape(tmp_grid, tmp_grid.size)
        # griddata (after removing nans)
        idx_good = np.argwhere(~np.isnan(tmp_vec))
        if idx_good.size > 5: # will ignore depth where no data exist
            LN = np.squeeze(lon_vec[idx_good])
            LT = np.squeeze(lat_vec[idx_good])
            TT = np.squeeze(tmp_vec[idx_good])
            zi = griddata((LN, LT), TT, (lon_grid, lat_grid), method='linear')

            #Mask the array where the max interpolation distance is exceeded
            THRESHOLD = 2
            tree = cKDTree(np.array([LN,LT]).T)
            xi = _ndim_coords_from_arrays(tuple([lon_grid,lat_grid]))
            dists, indexes = tree.query(xi)
            zi[dists > THRESHOLD] = np.nan

            V[:,:,k] = zi
            
        else:
            continue
    print(' -> Done!')

    # mask using bathymetry (I don't think it is necessary, but make nice figures)
    for i, xx in enumerate(lon_reg):
        for j,yy in enumerate(lat_reg):
            if Zitp[j,i] > -10: # remove shallower than 10m
                V[j,i,:] = np.nan


    # getting bottom temperature
    print('Getting bottom Temp.')    
    Tbot = np.full([lat_reg.size,lon_reg.size], np.nan)
    for i, xx in enumerate(lon_reg):
        for j,yy in enumerate(lat_reg):
            bottom_depth = -Zitp[j,i] # minus to turn positive
            temp_vec = V[j,i,:]
            idx_good = np.squeeze(np.where(~np.isnan(temp_vec)))
            if idx_good.size:
                idx_closest = np.argmin(np.abs(bottom_depth-z[idx_good]))
            else:
                continue

            if np.abs([idx_closest] - bottom_depth) <= 20:
                Tbot[j,i] = temp_vec[idx_good[idx_closest]]
            elif np.abs(z[idx_closest] - bottom_depth) <= 50:
                #print('used data located [30,50]m from bottom')
                Tbot[j,i] = temp_vec[idx_good[idx_closest]]

    print(' -> Done!')


    # Use climatology to fill missing pixels
    if clim_fill:
        print('Fill NaNs with climatology')
        Tbot_orig = Tbot.copy()
        #Tbot_fill = Tbot_climato[Tbot_orig.isna()]
        #Tbot[Tbot_orig.isna()] = Tbot_climato[Tbot_orig.isna()]
        Tbot_fill = Tbot_climato[np.isnan(Tbot_orig)]
        Tbot[np.isnan(Tbot_orig)] = Tbot_climato[np.isnan(Tbot_orig)]



    # Mask data outside Nafo div.
    print('Mask according to NAFO division for ' + season)
    # Polygons
    polygon4R = Polygon(zip(nafo_div['4R']['lon'], nafo_div['4R']['lat']))
    polygon3K = Polygon(zip(nafo_div['3Kx']['lon'], nafo_div['3Kx']['lat']))
    polygon3L = Polygon(zip(nafo_div['3L']['lon'], nafo_div['3L']['lat']))
    polygon3N = Polygon(zip(nafo_div['3N']['lon'], nafo_div['3N']['lat']))
    polygon3O = Polygon(zip(nafo_div['3O']['lon'], nafo_div['3O']['lat']))
    polygon3Ps = Polygon(zip(nafo_div['3Ps']['lon'], nafo_div['3Ps']['lat']))
    polygon2J = Polygon(zip(nafo_div['2J']['lon'], nafo_div['2J']['lat']))
    polygon2H = Polygon(zip(nafo_div['2H']['lon'], nafo_div['2H']['lat']))

    # Contour of data to mask
    contour_mask = np.load('operation_files/100m_contour_labrador.npy')
    polygon_mask = Polygon(contour_mask)

    if season == 'spring':
        for i, xx in enumerate(lon_reg):
            for j,yy in enumerate(lat_reg):
                point = Point(lon_reg[i], lat_reg[j])
                if polygon3L.contains(point) | polygon3N.contains(point) | polygon3O.contains(point) | polygon3Ps.contains(point) | polygon4R.contains(point):
                    pass #nothing to do but cannot implement negative statement "if not" above
                else:
                    Tbot[j,i] = np.nan
                    Tbot_orig[j,i] = np.nan

    elif season == 'fall':
        for i, xx in enumerate(lon_reg):
            for j,yy in enumerate(lat_reg):
                point = Point(lon_reg[i], lat_reg[j])
                if polygon2H.contains(point) | polygon2J.contains(point) | polygon3K.contains(point) | polygon3L.contains(point) | polygon3N.contains(point) | polygon3O.contains(point):
                    pass #nothing to do but cannot implement negative statement "if not" above
                else:
                    Tbot[j,i] = np.nan ### <--------------------- Do mask the fall / OR / 
                    Tbot_orig[j,i] = np.nan
                    #Tbot[j,i] = np.nan ### <--------------------- Do not mask the fall!!!!!

                if polygon_mask.contains(point): # mask data near Labrador in fall
                    Tbot[j,i] = np.nan 
                    Tbot_orig[j,i] = np.nan

    elif season == 'summer':
        for i, xx in enumerate(lon_reg):
            for j,yy in enumerate(lat_reg):
                point = Point(lon_reg[i], lat_reg[j])
                # Just mask labrador
                if polygon_mask.contains(point): # mask data near Labrador in fall
                    Tbot[j,i] = np.nan 
                    Tbot_orig[j,i] = np.nan

    else:
        print('no division mask, all data taken')

    print(' -> Done!')

    # Temperature anomaly:
    anom = Tbot-Tbot_climato
    div_toplot = ['2H', '2J', '3K', '3L', '3N', '3O', '3Ps', '4R']

    #Set up the map
    land_10m = cfeature.NaturalEarthFeature('physical','land','10m',
    facecolor='tan')

    ## ---- Plot Anomaly ---- ##
    fig = plt.figure(figsize=(6,9))
    ax = plt.subplot(projection=ccrs.Mercator())

    #Plot the coastline
    ax.set_facecolor('white')
    ax.add_feature(land_10m,zorder=2)
    ax.set_extent([lonLims[0],lonLims[1],latLims[0],latLims[1]])

    #Plot the temperature anomalies
    levels = np.linspace(-3.5, 3.5, 8)
    xi, yi = np.meshgrid(lon_reg, lat_reg)
    c = ax.contourf(xi, yi, anom, levels, cmap=cmo.cm.balance, extend='both', transform=ccrs.PlateCarree())

    #Plot the bathymetry
    cc = ax.contour(xi, yi, -Zitp, [100, 500, 1000, 4000], colors='grey', linewidths=0.75, transform=ccrs.PlateCarree())
    ax.clabel(cc, inline=1, fontsize=7, fmt='%d')

    #Plot the titles
    if season=='fall':
        plt.title('Fall Bottom Temperature ' + year + ' Anomaly')
    elif season=='spring':
        plt.title('Spring Bottom Temperature ' + year + ' Anomaly')
    else:
        plt.title('Bottom Temperature ' + year + '  Anomaly')

    #Add gridlines
    gl = ax.gridlines(draw_labels=['bottom'], xlocs=[-50,-55,-60], ylocs=[40, 45, 50, 55, 60],
        dms=True, x_inline=False, y_inline=False, linestyle='--')
    gl.xlabel_style = {'size': 8}

    #Add a colourbar
    cax = fig.add_axes([0.15, 0.07, 0.73, 0.025])
    cb = plt.colorbar(c, cax=cax, orientation='horizontal')
    cb.set_label(r'$\rm T(^{\circ}C)$', fontsize=8, fontweight='normal')
    cb.ax.tick_params(labelsize=8)

    #Plot the divisions
    for div in div_toplot:
        div_lon, div_lat = nafo_div[div]['lon'], nafo_div[div]['lat']
        ax.plot(div_lon, div_lat, 'k', linewidth=1, transform=ccrs.PlateCarree())
        ax.text(np.mean(div_lon), np.mean(div_lat), div, fontsize=8, color='black', fontweight='bold', transform=ccrs.PlateCarree()) 


    # Save Figure
    outfile = 'bottom_temp_anomaly_' + season + '_' + year + '.png'
    fig.savefig(outfile, dpi=300)
    os.system('convert -trim ' + outfile + ' ' + outfile)
    # Save French Figure
    plt.sca(ax)
    if season=='fall':
        plt.title(u'Anomalie de température au fond - Automne ' + year)
    elif season=='spring':
        plt.title(u'Anomalie de température au fond - Printemp ' + year)
    else:
        plt.title(u'Anomalie de température au fond ' + year)
    outfile = 'bottom_temp_anomaly_' + season + '_' + year + '_FR.png'
    fig.savefig(outfile, dpi=300)
    os.system('convert -trim ' + outfile + ' ' + outfile)



    ## ---- Plot Temperature ---- ##


    #Set up the map
    land_10m = cfeature.NaturalEarthFeature('physical','land','10m',
    facecolor='tan')

    ## ---- Plot Anomaly ---- ##
    fig = plt.figure(figsize=(6,9))
    ax = plt.subplot(projection=ccrs.Mercator())

    #Plot the coastline
    ax.set_facecolor('white')
    ax.add_feature(land_10m,zorder=2)
    ax.set_extent([lonLims[0],lonLims[1],latLims[0],latLims[1]])

    #Plot the temperature 
    levels = np.linspace(-2,6,17)
    xi, yi = np.meshgrid(lon_reg, lat_reg)
    ax.contourf(xi, yi, Tbot, levels, cmap=cmo.cm.thermal, extend='both', transform=ccrs.PlateCarree(), alpha=0.3)
    c = ax.contourf(xi, yi, Tbot_orig, levels, cmap=cmo.cm.thermal, extend='both', transform=ccrs.PlateCarree())

    #Plot the bathymetry
    cc = ax.contour(xi, yi, -Zitp, [100, 500, 1000, 4000], colors='grey', linewidths=0.75, transform=ccrs.PlateCarree())
    ax.clabel(cc, inline=1, fontsize=7, fmt='%d')

    #Plot the titles
    if season=='fall':
        plt.title('Fall Bottom Temperature ' + year)
    elif season=='spring':
        plt.title('Spring Bottom Temperature ' + year)
    else:
        plt.title('Bottom Temperature ' + year)

    #Plot the cast locations
    ax.scatter(lons, lats, s=2, marker='.', c='k', transform=ccrs.PlateCarree())

    #Add gridlines
    gl = ax.gridlines(draw_labels=['bottom'], xlocs=[-50,-55,-60], ylocs=[40, 45, 50, 55, 60],
        dms=True, x_inline=False, y_inline=False, linestyle='--')
    gl.xlabel_style = {'size': 8}

    #Add a colourbar
    cax = fig.add_axes([0.15, 0.07, 0.73, 0.025])
    cb = plt.colorbar(c, cax=cax, orientation='horizontal')
    cb.set_label(r'$\rm T(^{\circ}C)$', fontsize=8, fontweight='normal')
    cb.ax.tick_params(labelsize=8)

    #Plot the divisions
    for div in div_toplot:
        div_lon, div_lat = nafo_div[div]['lon'], nafo_div[div]['lat']
        ax.plot(div_lon, div_lat, 'k', linewidth=1, transform=ccrs.PlateCarree())
        ax.text(np.mean(div_lon), np.mean(div_lat), div, fontsize=8, color='black', fontweight='bold', transform=ccrs.PlateCarree())


    # Save Figure
    outfile = 'bottom_temp_' + season + '_' + year + '.png'
    fig.savefig(outfile, dpi=300)
    os.system('convert -trim ' + outfile + ' ' + outfile)
    # Save French Figure
    plt.sca(ax)
    if season=='fall':
        plt.title(u'Température au fond - Automne ' + year )
    elif season=='spring':
        plt.title(u'Température au fond - Printemp ' + year )
    else:
        plt.title(u'Température au fond ' + year )
    outfile = 'bottom_temp_' + season + '_' + year + '_FR.png'
    fig.savefig(outfile, dpi=300)
    os.system('convert -trim ' + outfile + ' ' + outfile)



    ## ---- Plot Climato ---- ##




    #Set up the map
    land_10m = cfeature.NaturalEarthFeature('physical','land','10m',
    facecolor='tan')

    ## ---- Plot Anomaly ---- ##
    fig = plt.figure(figsize=(6,9))
    ax = plt.subplot(projection=ccrs.Mercator())

    #Plot the coastline
    ax.set_facecolor('white')
    ax.add_feature(land_10m,zorder=2)
    ax.set_extent([lonLims[0],lonLims[1],latLims[0],latLims[1]])

    #Plot the temperature 
    levels = np.linspace(-2,6,17)
    xi, yi = np.meshgrid(lon_reg, lat_reg)
    c = ax.contourf(xi, yi, Tbot_climato, levels, cmap=cmo.cm.thermal, extend='both', transform=ccrs.PlateCarree())

    #Plot the bathymetry
    cc = ax.contour(xi, yi, -Zitp, [100, 500, 1000, 4000], colors='grey', linewidths=0.75, transform=ccrs.PlateCarree())
    ax.clabel(cc, inline=1, fontsize=7, fmt='%d')

    #Plot the titles
    if season=='fall':
        plt.title('Fall Bottom Temperature Climatology')
    elif season=='spring':
        plt.title('Spring Bottom Temperature Climatology')
    else:
        plt.title('Bottom Temperature Climatology')

    #Add gridlines
    gl = ax.gridlines(draw_labels=['left','bottom'], xlocs=[-50,-55,-60], ylocs=[40, 45, 50, 55, 60],
        dms=True, x_inline=False, y_inline=False, linestyle='--')    
    gl.xlabel_style = {'size': 8}

    #Add a colourbar
    cax = fig.add_axes([0.15, 0.07, 0.73, 0.025])
    cb = plt.colorbar(c, cax=cax, orientation='horizontal')
    cb.set_label(r'$\rm T(^{\circ}C)$', fontsize=8, fontweight='normal')
    cb.ax.tick_params(labelsize=8)

    #Plot the divisions
    for div in div_toplot:
        div_lon, div_lat = nafo_div[div]['lon'], nafo_div[div]['lat']
        ax.plot(div_lon, div_lat, 'k', linewidth=1, transform=ccrs.PlateCarree())
        ax.text(np.mean(div_lon), np.mean(div_lat), div, fontsize=8, color='black', fontweight='bold', transform=ccrs.PlateCarree())


    # Save Figure
    outfile = 'bottom_temp_climato_' + season + '_' + year + '.png'
    fig.savefig(outfile, dpi=300)
    os.system('convert -trim ' + outfile + ' ' + outfile)
    # Save French Figure
    plt.sca(ax)
    if season=='fall':
        plt.title(u'Climatologie de température au fond - Automne ' + year )
    elif season=='spring':
        plt.title(u'Climatologie de température au fond - Printemps ' + year )
    else:
        plt.title(u'Climatologie de température au fond ' + year )
    outfile = 'bottom_temp_climato_' + season + '_' + year + '_FR.png'
    fig.savefig(outfile, dpi=300)
    os.system('convert -trim ' + outfile + ' ' + outfile)


    # Convert to a subplot
    os.system('montage bottom_temp_climato_' + season + '_' + year + '.png bottom_temp_' + season + '_' + year + '.png bottom_temp_anomaly_' + season + '_' + year + '.png  -tile 3x1 -geometry +10+10  -background white  bottomT_' + season + year + '.png')
    # in French
    os.system('montage bottom_temp_climato_' + season + '_' + year + '_FR.png bottom_temp_' + season + '_' + year + '_FR.png bottom_temp_anomaly_' + season + '_' + year + '_FR.png  -tile 3x1 -geometry +10+10  -background white  bottomT_' + season + year + '_FR.png')
    # Move to year folder
    os.system('cp bottomT_' + season + year + '.png bottomT_' + season + year + '_FR.png ' + year)


#### bottom_salinity
def bottom_salinity(season, year, zmin=0, zmax=1000, dz=5, proj='merc', clim_fill=True, netcdf_path='/home/cyrf0006/data/dev_database/netCDF/', climato_file=''):
    '''
    Bottom salinity maps for AZMP ResDoc

    This function is adapted from azmp_bottomS.py

    * process in : /home/cyrf0006/AZMP/state_reports/bottomT

    usage ex:
    >> import azmp_report_tools as azrt
    >> azrt.bottom_salinity(season='spring', year='2019'):

    Frederic.Cyr@dfo-mpo.gc.ca - January 2020

    '''

    '''
    #TEMPORARY, FOR TESTING PURPOSES ONLY
    season='fall' 
    year='2018'
    zmin=0
    zmax=1000
    clim_fill=True
    dz=5
    proj='merc'
    climato_file='Sbot_climato_'+season+'_0.10.h5'
    netcdf_path='/home/jcoyne/Documents/CASH/Combined_Data/CASTS_new-vertical_v2/'
    '''

    if len(climato_file) == 0:
        if season=='spring':
            climato_file='operation_files/Sbot_climato_spring_0.10.h5'
        elif season=='fall':
            climato_file='operation_files/Sbot_climato_fall_0.10.h5'
        elif season=='summer':
            climato_file='operation_files/Sbot_climato_summer_0.10.h6'

    year_file = netcdf_path + str(year) + '.nc'


    ## ---- Load Climato data ---- ##    
    print('Load ' + climato_file)
    h5f = h5py.File(climato_file, 'r')
    Sbot_climato = h5f['Sbot'][:]
    lon_reg = h5f['lon_reg'][:]
    lat_reg = h5f['lat_reg'][:]
    Zitp = h5f['Zitp'][:]
    h5f.close()

    ## ---- Derive some parameters ---- ##    
    lon_0 = np.round(np.mean(lon_reg))
    lat_0 = np.round(np.mean(lat_reg))
    lonLims = [lon_reg[0], lon_reg[-1]]
    latLims = [lat_reg[0], lat_reg[-1]]
    lon_grid, lat_grid = np.meshgrid(lon_reg,lat_reg)
    dc = np.diff(lon_reg[0:2])

    ## ---- NAFO divisions ---- ##
    nafo_div = azu.get_nafo_divisions()

    ## ---- Get CTD data --- ##
    print('Get ' + year_file)
    ds = xr.open_mfdataset(year_file)
    # Selection of a subset region
    ds = ds.where((ds.longitude>lonLims[0]) & (ds.longitude<lonLims[1]), drop=True)
    ds = ds.where((ds.latitude>latLims[0]) & (ds.latitude<latLims[1]), drop=True)
    # Select time (save several options here)
    if season == 'summer':
        ds = ds.sel(time=((ds['time.month']>=6)) & ((ds['time.month']<=9)))
    elif season == 'spring':
        ds = ds.sel(time=((ds['time.month']>=4)) & ((ds['time.month']<=6)))
    elif season == 'fall':
        ds = ds.sel(time=((ds['time.month']>=9)) & ((ds['time.month']<=12)))
    else:
        print('!! no season specified, used them all! !!')

    # Restrict max depth to zmax defined earlier
    ds = ds.sel(level=ds['level']<zmax)
    lons = np.array(ds.longitude)
    lats = np.array(ds.latitude)
    da_sal = ds['salinity']
    #To Pandas Dataframe
    df_sal = da_sal.to_pandas()
    # Remove empty columns & drop coordinates (for cast identification on map)
    idx_empty_rows = df_sal.isnull().all(1).values.nonzero()[0]
    df_sal = df_sal.dropna(axis=0,how='all')
    lons = np.delete(lons,idx_empty_rows)
    lats = np.delete(lats,idx_empty_rows)
    print(' -> Done!')

    ## --- fill 3D cube --- ##  
    print('Fill regular cube')
    z = df_sal.columns.values
    V = np.full((lat_reg.size, lon_reg.size, z.size), np.nan)

    # Aggregate on regular grid
    for i, xx in enumerate(lon_reg):
        for j, yy in enumerate(lat_reg):    
            idx = np.where((lons>=xx-dc/2) & (lons<xx+dc/2) & (lats>=yy-dc/2) & (lats<yy+dc/2))
            tmp = np.array(df_sal.iloc[idx].mean(axis=0))
            idx_good = np.argwhere(~np.isnan(tmp))
            if np.size(idx_good)==1:
                V[j,i,:] = np.array(df_sal.iloc[idx].mean(axis=0))
            elif np.size(idx_good)>1: # vertical interpolation between pts
                interp = interp1d(np.squeeze(z[idx_good]), np.squeeze(tmp[idx_good]))
                idx_interp = np.arange(int(idx_good[0]),int(idx_good[-1]+1))
                V[j,i,idx_interp] = interp(z[idx_interp]) # interpolate only where possible (1st to last good idx)

    # horizontal interpolation at each depth
    lon_grid, lat_grid = np.meshgrid(lon_reg,lat_reg)
    lon_vec = np.reshape(lon_grid, lon_grid.size)
    lat_vec = np.reshape(lat_grid, lat_grid.size)
    for k, zz in enumerate(z):
        # Meshgrid 1D data (after removing NaNs)
        tmp_grid = V[:,:,k]
        tmp_vec = np.reshape(tmp_grid, tmp_grid.size)
        # griddata (after removing nans)
        idx_good = np.argwhere(~np.isnan(tmp_vec))
        if idx_good.size>5: # will ignore depth where no data exist
            LN = np.squeeze(lon_vec[idx_good])
            LT = np.squeeze(lat_vec[idx_good])
            TT = np.squeeze(tmp_vec[idx_good])
            zi = griddata((LN, LT), TT, (lon_grid, lat_grid), method='linear')

            #Mask the array where the max interpolation distance is exceeded
            THRESHOLD = 2
            tree = cKDTree(np.array([LN,LT]).T)
            xi = _ndim_coords_from_arrays(tuple([lon_grid,lat_grid]))
            dists, indexes = tree.query(xi)
            zi[dists > THRESHOLD] = np.nan

            V[:,:,k] = zi
        else:
            continue
    print(' -> Done!')


    # mask using bathymetry (I don't think it is necessary, but make nice figures)
    for i, xx in enumerate(lon_reg):
        for j,yy in enumerate(lat_reg):
            if Zitp[j,i] > -10: # remove shallower than 10m
                V[j,i,:] = np.nan

    # getting bottom temperature
    print('Getting bottom Temp.')    
    Sbot = np.full([lat_reg.size,lon_reg.size], np.nan) 
    for i, xx in enumerate(lon_reg):
        for j,yy in enumerate(lat_reg):
            bottom_depth = -Zitp[j,i] # minus to turn positive
            temp_vec = V[j,i,:]
            idx_good = np.squeeze(np.where(~np.isnan(temp_vec)))
            if idx_good.size:
                idx_closest = np.argmin(np.abs(bottom_depth-z[idx_good]))
            else:
                continue

            if np.abs([idx_closest] - bottom_depth) <= 20:
                Sbot[j,i] = temp_vec[idx_good[idx_closest]]
            elif np.abs(z[idx_closest] - bottom_depth) <= 50:
                #print('used data located [30,50]m from bottom')
                Sbot[j,i] = temp_vec[idx_good[idx_closest]]

    print(' -> Done!')    

    # Use climatology to fill missing pixels
    if clim_fill:
        print('Fill NaNs with climatology')
        Sbot_orig = Sbot.copy()
        #Tbot_fill = Tbot_climato[Tbot_orig.isna()]
        #Tbot[Tbot_orig.isna()] = Tbot_climato[Tbot_orig.isna()]
        Sbot_fill = Sbot_climato[np.isnan(Sbot_orig)]
        Sbot[np.isnan(Sbot_orig)] = Sbot_climato[np.isnan(Sbot_orig)]            


    # Mask data outside Nafo div.
    print('Mask according to NAFO division for ' + season)
    # Polygons
    polygon4R = Polygon(zip(nafo_div['4R']['lon'], nafo_div['4R']['lat']))
    polygon3K = Polygon(zip(nafo_div['3Kx']['lon'], nafo_div['3Kx']['lat']))
    polygon3L = Polygon(zip(nafo_div['3L']['lon'], nafo_div['3L']['lat']))
    polygon3N = Polygon(zip(nafo_div['3N']['lon'], nafo_div['3N']['lat']))
    polygon3O = Polygon(zip(nafo_div['3O']['lon'], nafo_div['3O']['lat']))
    polygon3Ps = Polygon(zip(nafo_div['3Ps']['lon'], nafo_div['3Ps']['lat']))
    polygon2J = Polygon(zip(nafo_div['2J']['lon'], nafo_div['2J']['lat']))
    polygon2H = Polygon(zip(nafo_div['2H']['lon'], nafo_div['2H']['lat']))

    # Contour of data to mask
    contour_mask = np.load('operation_files/100m_contour_labrador.npy')
    polygon_mask = Polygon(contour_mask)

    if season == 'spring':
        for i, xx in enumerate(lon_reg):
            for j,yy in enumerate(lat_reg):
                point = Point(lon_reg[i], lat_reg[j])
                if polygon3L.contains(point) | polygon3N.contains(point) | polygon3O.contains(point) | polygon3Ps.contains(point) | polygon4R.contains(point):
                    pass #nothing to do but cannot implement negative statement "if not" above
                else:
                    Sbot[j,i] = np.nan
                    Sbot_orig[j,i] = np.nan

    elif season == 'fall':
        for i, xx in enumerate(lon_reg):
            for j,yy in enumerate(lat_reg):
                point = Point(lon_reg[i], lat_reg[j])
                if polygon2H.contains(point) | polygon2J.contains(point) | polygon3K.contains(point) | polygon3L.contains(point) | polygon3N.contains(point) | polygon3O.contains(point):
                    pass #nothing to do but cannot implement negative statement "if not" above
                else:
                    Sbot[j,i] = np.nan ### <--------------------- Do not mask the fall!!!!!
                    Sbot_orig[j,i] = np.nan

                if polygon_mask.contains(point): # mask data near Labrador in fall
                    Sbot[j,i] = np.nan 
                    Sbot_orig[j,i] = np.nan
    else:
        print('no division mask, all data taken')

    print(' -> Done!')    

    # Salinity anomaly:
    anom = Sbot-Sbot_climato
    div_toplot = ['2H', '2J', '3K', '3L', '3N', '3O', '3Ps', '4R']

    #Set up the map
    land_10m = cfeature.NaturalEarthFeature('physical','land','10m',
    facecolor='tan')

    ## ---- Plot Anomaly ---- ##
    fig = plt.figure(figsize=(6,9))
    ax = plt.subplot(projection=ccrs.Mercator())

    #Plot the coastline
    ax.set_facecolor('white')
    ax.add_feature(land_10m,zorder=2)
    ax.set_extent([lonLims[0],lonLims[1],latLims[0],latLims[1]])

    #Plot the salinity anomalies
    levels = np.array([-1, -.8, -.6, -.4, -.2, .2, .4, .6, .8, 1])
    xi, yi = np.meshgrid(lon_reg, lat_reg)
    c = ax.contourf(xi, yi, anom, levels, cmap=cmo.cm.balance, extend='both', transform=ccrs.PlateCarree())
    
    #Plot the bathymetry
    cc = ax.contour(xi, yi, -Zitp, [100, 500, 1000, 4000], colors='grey', linewidths=0.75, transform=ccrs.PlateCarree())
    ax.clabel(cc, inline=1, fontsize=7, fmt='%d')

    #Plot the titles
    if season=='fall':
        plt.title('Fall Bottom Salinity ' + year + ' Anomaly')
    elif season=='spring':
        plt.title('Spring Bottom Salinity ' + year + ' Anomaly')
    else:
        plt.title('Bottom Salinity ' + year + '  Anomaly')
    
    #Add gridlines
    gl = ax.gridlines(draw_labels=['bottom'], xlocs=[-50,-55,-60], ylocs=[40, 45, 50, 55, 60],
        dms=True, x_inline=False, y_inline=False, linestyle='--')    
    gl.xlabel_style = {'size': 8}

    #Add a colourbar
    cax = fig.add_axes([0.15, 0.07, 0.73, 0.025])
    cb = plt.colorbar(c, cax=cax, orientation='horizontal')
    cb.set_label(r'$\rm S$', fontsize=8, fontweight='normal')
    cb.ax.tick_params(labelsize=8)
    
    #Plot the divisions
    for div in div_toplot:
        div_lon, div_lat = nafo_div[div]['lon'], nafo_div[div]['lat']
        ax.plot(div_lon, div_lat, 'k', linewidth=1, transform=ccrs.PlateCarree())
        ax.text(np.mean(div_lon), np.mean(div_lat), div, fontsize=8, color='black', fontweight='bold', transform=ccrs.PlateCarree()) 


    # Save Figure
    outfile = 'bottom_sal_anomaly_' + season + '_' + year + '.png'
    fig.savefig(outfile, dpi=300)
    os.system('convert -trim ' + outfile + ' ' + outfile)
    # Save French Figure
    plt.sca(ax)
    if season=='fall':
        plt.title(u'Anomalie de salinité au fond - Automne ' + year )
    elif season=='spring':
        plt.title(u'Anomalie de salinité au fond - Printemp ' + year )
    else:
        plt.title(u'Anomalie de salinité au fond ' + year )
    outfile = 'bottom_sal_anomaly_' + season + '_' + year + '_FR.png'
    fig.savefig(outfile, dpi=300)
    os.system('convert -trim ' + outfile + ' ' + outfile)




    ## ---- Plot Salinity ---- ##


    #Set up the map
    land_10m = cfeature.NaturalEarthFeature('physical','land','10m',
    facecolor='tan')

    fig = plt.figure(figsize=(6,9))
    ax = plt.subplot(projection=ccrs.Mercator())

    #Plot the coastline
    ax.set_facecolor('white')
    ax.add_feature(land_10m,zorder=2)
    ax.set_extent([lonLims[0],lonLims[1],latLims[0],latLims[1]])

    #Plot the salinity anomalies
    levels = np.linspace(30, 36, 13)
    xi, yi = np.meshgrid(lon_reg, lat_reg)
    ax.contourf(xi, yi, Sbot, levels, cmap=cmo.cm.haline, extend='both', transform=ccrs.PlateCarree(), alpha=0.3)
    c = ax.contourf(xi, yi, Sbot_orig, levels, cmap=cmo.cm.haline, extend='both', transform=ccrs.PlateCarree())
    
    #Plot the bathymetry
    cc = ax.contour(xi, yi, -Zitp, [100, 500, 1000, 4000], colors='grey', linewidths=0.75, transform=ccrs.PlateCarree())
    ax.clabel(cc, inline=1, fontsize=7, fmt='%d')

    #Plot the titles
    if season=='fall':
        plt.title('Fall Bottom Salinity ' + year)
    elif season=='spring':
        plt.title('Spring Bottom Salinity ' + year)
    else:
        plt.title('Bottom Salinity ' + year)

    #Plot the cast locations
    ax.scatter(lons, lats, s=2, marker='.', c='k', transform=ccrs.PlateCarree())

    #Add gridlines
    gl = ax.gridlines(draw_labels=['bottom'], xlocs=[-50,-55,-60], ylocs=[40, 45, 50, 55, 60],
        dms=True, x_inline=False, y_inline=False, linestyle='--')    
    gl.xlabel_style = {'size': 8}

    #Add a colourbar
    cax = fig.add_axes([0.15, 0.07, 0.73, 0.025])
    cb = plt.colorbar(c, cax=cax, orientation='horizontal')
    cb.set_label(r'$\rm S$', fontsize=8, fontweight='normal')
    cb.ax.tick_params(labelsize=8)
    
    #Plot the divisions
    for div in div_toplot:
        div_lon, div_lat = nafo_div[div]['lon'], nafo_div[div]['lat']
        ax.plot(div_lon, div_lat, 'k', linewidth=1, transform=ccrs.PlateCarree())
        ax.text(np.mean(div_lon), np.mean(div_lat), div, fontsize=8, color='black', fontweight='bold', transform=ccrs.PlateCarree()) 


    # Save Figure
    outfile = 'bottom_sal_' + season + '_' + year + '.png'
    fig.savefig(outfile, dpi=300)
    os.system('convert -trim ' + outfile + ' ' + outfile)
    # Save French Figure
    plt.sca(ax)
    if season=='fall':
        plt.title(u'Salinité au fond - Automne ' + year )
    elif season=='spring':
        plt.title(u'Salinité au fond - Printemp ' + year )
    else:
        plt.title(u'Salinité au fond ' + year )
    outfile = 'bottom_sal_' + season + '_' + year + '_FR.png'
    fig.savefig(outfile, dpi=300)
    os.system('convert -trim ' + outfile + ' ' + outfile)


    ## ---- Plot Climato ---- ##


    #Set up the map
    land_10m = cfeature.NaturalEarthFeature('physical','land','10m',
    facecolor='tan')

    fig = plt.figure(figsize=(6,9))
    ax = plt.subplot(projection=ccrs.Mercator())

    #Plot the coastline
    ax.set_facecolor('white')
    ax.add_feature(land_10m,zorder=2)
    ax.set_extent([lonLims[0],lonLims[1],latLims[0],latLims[1]])

    #Plot the salinity anomalies
    levels = np.linspace(30, 36, 13)
    xi, yi = np.meshgrid(lon_reg, lat_reg)
    c = ax.contourf(xi, yi, Sbot_climato, levels, cmap=cmo.cm.haline, extend='both', transform=ccrs.PlateCarree())
    
    #Plot the bathymetry
    cc = ax.contour(xi, yi, -Zitp, [100, 500, 1000, 4000], colors='grey', linewidths=0.75, transform=ccrs.PlateCarree())
    ax.clabel(cc, inline=1, fontsize=7, fmt='%d')

    #Plot the titles
    if season=='fall':
        plt.title('Fall Bottom Salinity Climatology')
    elif season=='spring':
        plt.title('Spring Bottom Salinity Climatology')
    else:
        plt.title('Bottom Salinity Climatology')

    #Add gridlines
    gl = ax.gridlines(draw_labels=['left','bottom'], xlocs=[-50,-55,-60], ylocs=[40, 45, 50, 55, 60],
        dms=True, x_inline=False, y_inline=False, linestyle='--')    
    gl.xlabel_style = {'size': 8}

    #Add a colourbar
    cax = fig.add_axes([0.15, 0.07, 0.73, 0.025])
    cb = plt.colorbar(c, cax=cax, orientation='horizontal')
    cb.set_label(r'$\rm S$', fontsize=8, fontweight='normal')
    cb.ax.tick_params(labelsize=8)
    
    #Plot the divisions
    for div in div_toplot:
        div_lon, div_lat = nafo_div[div]['lon'], nafo_div[div]['lat']
        ax.plot(div_lon, div_lat, 'k', linewidth=1, transform=ccrs.PlateCarree())
        ax.text(np.mean(div_lon), np.mean(div_lat), div, fontsize=8, color='black', fontweight='bold', transform=ccrs.PlateCarree()) 


    # Save Figure
    outfile = 'bottom_sal_climato_' + season + '_' + year + '.png'
    fig.savefig(outfile, dpi=300)
    os.system('convert -trim ' + outfile + ' ' + outfile)
    # Save French Figure
    plt.sca(ax)
    if season=='fall':
        plt.title(u'Climatoligie de salinité au fond - Automne ' + year )
    elif season=='spring':
        plt.title(u'Climatologie de salinité au fond - Printemp ' + year )
    else:
        plt.title(u'Climatologie de salinité au fond ' + year )
    outfile = 'bottom_sal_climato_' + season + '_' + year + '_FR.png'
    fig.savefig(outfile, dpi=300)
    os.system('convert -trim ' + outfile + ' ' + outfile)


    # Convert to a subplot
    os.system('montage bottom_sal_climato_' + season + '_' + year + '.png bottom_sal_' + season + '_' + year + '.png bottom_sal_anomaly_' + season + '_' + year + '.png  -tile 3x1 -geometry +10+10  -background white  bottomS_' + season + year + '.png') 
    # French
    os.system('montage bottom_sal_climato_' + season + '_' + year + '_FR.png bottom_sal_' + season + '_' + year + '_FR.png bottom_sal_anomaly_' + season + '_' + year + '_FR.png  -tile 3x1 -geometry +10+10  -background white  bottomS_' + season + year + '_FR.png') 
    # Move to year folder
    os.system('cp bottomS_' + season + year + '.png bottomS_' + season + year + '_FR.png ' + year)


#### bottom_stats
def bottom_stats(years, season, proj='merc', plot=False, netcdf_path='/home/cyrf0006/data/dev_database/netCDF/', climato_file=''):

    '''
        Function bottom_stats() based on script azmp_bottom_stats.py

        See the latter on how to generate the bottom climatology.
        See it also for specific usage such as Plaice - COSEWIC analysis.

        usage example:
        >> import azmp_report_tools as azrt
        >> import numpy as np
        >> azrt.bottom_stats(years=np.arange(1980, 2020), season='fall')

        *** This needs to be improve because at the moment I need to comment the generation of .pkl file to not over-write when I change my map region.        
                
        For NAFO SA4:
        >> azrt.bottom_stats(years=np.arange(1980, 2020), season='summer', climato_file='Tbot_climato_SA4_summer_0.10.h5')
        
        For COSEWIC:
        >> azrt.bottom_stats(years=np.arange(1980, 2020), season='summer', climato_file='Tbot_climato_2GH_summer_0.10.h5')
        >> azrt.bottom_stats(years=np.arange(1980, 2020), season='summer', climato_file='Tbot_climato_SA45_summer_0.10.h5')

        
        Frederic.Cyr@dfo-mpo.gc.ca - January 2020
    '''

    '''
    #FOR TESTING PURPOSES ONLY
    years=np.arange(1980,2023)
    season='fall'
    proj='merc'
    plot=False
    climato_file=''
    netcdf_path='/home/jcoyne/Documents/CASH/Combined_Data/CASTS_new-vertical_v2/'
    '''

    #Finding the mean of an empty slice or double scalar errorys
    np.seterr(divide='ignore', invalid='ignore')
    warnings.simplefilter("ignore", category=RuntimeWarning)

    # load climato
    if len(climato_file) == 0:
        if season=='spring':
            climato_file='operation_files/Tbot_climato_spring_0.10.h5'
        elif season=='fall':
            climato_file='operation_files/Tbot_climato_fall_0.10.h5'
        elif season=='summer':
            climato_file='operation_files/Tbot_climato_summer_0.10.h6'
    else:
        print('Climato file provided')
               
    h5f = h5py.File(climato_file, 'r')
    Tbot_climato = h5f['Tbot'][:]
    lon_reg = h5f['lon_reg'][:]
    lat_reg = h5f['lat_reg'][:]
    Zitp = h5f['Zitp'][:]
    h5f.close()

    # Derive some map parameters
    lon_0 = np.round(np.mean(lon_reg))
    lat_0 = np.round(np.mean(lat_reg))
    lonLims = [lon_reg[0], lon_reg[-1]]
    latLims = [lat_reg[0], lat_reg[-1]]

    # NAFO divisions
    nafo_div = azu.get_nafo_divisions()
    polygon3L = Polygon(zip(nafo_div['3L']['lon'], nafo_div['3L']['lat']))
    polygon3N = Polygon(zip(nafo_div['3N']['lon'], nafo_div['3N']['lat']))
    polygon3O = Polygon(zip(nafo_div['3O']['lon'], nafo_div['3O']['lat']))
    shape = [polygon3L, polygon3N, polygon3O]
    shape_3LNO = unary_union(shape)
    shape_3M = Polygon(zip(nafo_div['3M']['lon'], nafo_div['3M']['lat']))
    shape_3Ps = Polygon(zip(nafo_div['3Ps']['lon'], nafo_div['3Ps']['lat']))
    shape_2G = Polygon(zip(nafo_div['2G']['lon'], nafo_div['2G']['lat']))
    shape_2H = Polygon(zip(nafo_div['2H']['lon'], nafo_div['2H']['lat']))
    shape_2J = Polygon(zip(nafo_div['2J']['lon'], nafo_div['2J']['lat']))
    shape_3K = Polygon(zip(nafo_div['3K']['lon'], nafo_div['3K']['lat']))
    shape_3L = Polygon(zip(nafo_div['3L']['lon'], nafo_div['3L']['lat']))
    shape_3O = Polygon(zip(nafo_div['3O']['lon'], nafo_div['3O']['lat']))
    shape = [shape_2J, shape_2H]
    shape_2HJ = unary_union(shape)
    shape = [shape_2G, shape_2H]
    shape_2GH = unary_union(shape)
    shape_4R = Polygon(zip(nafo_div['4R']['lon'], nafo_div['4R']['lat']))
    shape_4S = Polygon(zip(nafo_div['4S']['lon'], nafo_div['4S']['lat']))
    shape_4T = Polygon(zip(nafo_div['4T']['lon'], nafo_div['4T']['lat']))
    shape = [shape_4R, shape_4S]
    shape_4RS = unary_union(shape)
    shape = [shape_4R, shape_4S, shape_4T]
    shape_4RST = unary_union(shape)
    polygon4Vs = Polygon(zip(nafo_div['4Vs']['lon'], nafo_div['4Vs']['lat']))
    polygon4Vn = Polygon(zip(nafo_div['4Vn']['lon'], nafo_div['4Vn']['lat']))
    polygon4W = Polygon(zip(nafo_div['4W']['lon'], nafo_div['4W']['lat']))
    polygon4X = Polygon(zip(nafo_div['4X']['lon'], nafo_div['4X']['lat']))
    shape = [polygon4Vs.buffer(0), polygon4Vn.buffer(0), polygon4W.buffer(0), polygon4X.buffer(0)]
    shape_4VWX = unary_union(shape)
    shape_5Y = Polygon(zip(nafo_div['5Y']['lon'], nafo_div['5Y']['lat']))

    dict_stats_3LNO = {}
    dict_stats_3M = {}
    dict_stats_3Ps = {}
    dict_stats_3K = {}
    dict_stats_3L = {}
    dict_stats_3O = {}
    dict_stats_2G = {}
    dict_stats_2H = {}
    dict_stats_2J = {}
    dict_stats_2HJ = {}
    dict_stats_2GH = {}
    dict_stats_4R = {}
    dict_stats_4S = {}
    dict_stats_4RS = {}
    dict_stats_4RT = {}
    dict_stats_4RST = {}
    dict_stats_4T = {}
    dict_stats_4VWX = {}
    dict_stats_5Y = {}

    # Loop on years
    df_list = []
    for year in years:
        print(' ---- ' + str(year) + ' ---- ')
        year_file = netcdf_path + str(year) + '.nc'
        Tdict = azu.get_bottomT(year_file, season, climato_file)    
        Tbot = Tdict['Tbot']
        lons = Tdict['lons']
        lats = Tdict['lats']
        anom = Tbot-Tbot_climato


        # NAFO division stats    
        dict_stats_2GH[str(year)] = azu.polygon_temperature_stats(Tdict, shape_2GH)
        dict_stats_2G[str(year)] = azu.polygon_temperature_stats(Tdict, shape_2G)
        dict_stats_2H[str(year)] = azu.polygon_temperature_stats(Tdict, shape_2H)
        dict_stats_2J[str(year)] = azu.polygon_temperature_stats(Tdict, shape_2J)
        dict_stats_2HJ[str(year)] = azu.polygon_temperature_stats(Tdict, shape_2HJ)
        dict_stats_3LNO[str(year)] = azu.polygon_temperature_stats(Tdict, shape_3LNO)
        dict_stats_3M[str(year)] = azu.polygon_temperature_stats(Tdict, shape_3M)
        dict_stats_3Ps[str(year)] = azu.polygon_temperature_stats(Tdict, shape_3Ps)
        dict_stats_3K[str(year)] = azu.polygon_temperature_stats(Tdict, shape_3K)
        dict_stats_3L[str(year)] = azu.polygon_temperature_stats(Tdict, shape_3L)
        dict_stats_3O[str(year)] = azu.polygon_temperature_stats(Tdict, shape_3O)
        ## dict_stats_4R[str(year)] = azu.polygon_temperature_stats(Tdict, shape_4R)
        ## dict_stats_4S[str(year)] = azu.polygon_temperature_stats(Tdict, shape_4S)
        ## dict_stats_4RS[str(year)] = azu.polygon_temperature_stats(Tdict, shape_4RS)
        ## dict_stats_4RST[str(year)] = azu.polygon_temperature_stats(Tdict, shape_4RST)
        ## dict_stats_4T[str(year)] = azu.polygon_temperature_stats(Tdict, shape_4T)
        ## dict_stats_4VWX[str(year)] = azu.polygon_temperature_stats(Tdict, shape_4VWX)
        #dict_stats_5Y[str(year)] = azu.polygon_temperature_stats(Tdict, shape_5Y)

        # Append bottom temperature for multi-index export
        df = pd.DataFrame(index=lat_reg, columns=lon_reg)
        df.index.name='latitude'
        df.columns.name='longitude'
        df[:] = Tbot
        df_list.append(df)


        if plot:
            div_toplot = ['2H', '2J', '3K', '3L', '3N', '3O', '3Ps', '4R', '4Vn', '4Vs', '4W', '4X']
    
            # 1.1 - Plot Anomaly
            fig, ax = plt.subplots(nrows=1, ncols=1)
            m = Basemap(ax=ax, projection='merc',lon_0=lon_0,lat_0=lat_0, llcrnrlon=lonLims[0],llcrnrlat=latLims[0],urcrnrlon=lonLims[1],urcrnrlat=latLims[1], resolution= 'i')
            levels = np.linspace(-3.5, 3.5, 8)
            xi, yi = m(*np.meshgrid(lon_reg, lat_reg))
            c = m.contourf(xi, yi, anom, levels, cmap=plt.cm.RdBu_r, extend='both')
            cc = m.contour(xi, yi, -Zitp, [100, 500, 1000, 4000], colors='grey');
            plt.clabel(cc, inline=1, fontsize=10, fmt='%d')
            if season=='fall':
                plt.title('Fall Bottom Temperature Anomaly')
            elif season=='spring':
                plt.title('Spring Bottom Temperature Anomaly')
            else:
                plt.title('Bottom Temperature Anomaly')
            m.fillcontinents(color='tan');
            m.drawparallels([40, 45, 50, 55, 60], labels=[1,0,0,0], fontsize=12, fontweight='normal');
            m.drawmeridians([-60, -55, -50, -45], labels=[0,0,0,1], fontsize=12, fontweight='normal');
            cax = fig.add_axes([0.16, 0.05, 0.7, 0.025])
            cb = plt.colorbar(c, cax=cax, orientation='horizontal')
            cb.set_label(r'$\rm T(^{\circ}C)$', fontsize=12, fontweight='normal')
            #cax = plt.axes([0.85,0.15,0.04,0.7], facecolor='grey')
            #cb = plt.colorbar(c, cax=cax)
            #cb.set_label(r'$\rm T(^{\circ}C)$', fontsize=12, fontweight='normal')
            for div in div_toplot:
                div_lon, div_lat = m(nafo_div[div]['lon'], nafo_div[div]['lat'])
                m.plot(div_lon, div_lat, 'k', linewidth=2)
                ax.text(np.mean(div_lon), np.mean(div_lat), div, fontsize=12, color='black', fontweight='bold')    
            # Save Figure
            fig.set_size_inches(w=7, h=8)
            fig.set_dpi(200)
            outfile = 'bottom_temp_anomaly_' + season + '_' + str(year) + '.png'
            fig.savefig(outfile)

            # 1.2 - Plot Temperature
            fig, ax = plt.subplots(nrows=1, ncols=1)
            m = Basemap(ax=ax, projection='merc',lon_0=lon_0,lat_0=lat_0, llcrnrlon=lonLims[0],llcrnrlat=latLims[0],urcrnrlon=lonLims[1],urcrnrlat=latLims[1], resolution= 'i')
            levels = np.linspace(-2, 6, 9)
            xi, yi = m(*np.meshgrid(lon_reg, lat_reg))
            c = m.contourf(xi, yi, Tbot, levels, cmap=plt.cm.RdBu_r, extend='both')
            cc = m.contour(xi, yi, -Zitp, [100, 500, 1000, 4000], colors='grey');
            plt.clabel(cc, inline=1, fontsize=10, fmt='%d')
            if season=='fall':
                plt.title('Fall Bottom Temperature')
            elif season=='spring':
                plt.title('Spring Bottom Temperature')
            else:
                plt.title('Bottom Temperature')
            m.fillcontinents(color='tan');
            m.drawparallels([40, 45, 50, 55, 60], labels=[1,0,0,0], fontsize=12, fontweight='normal');
            m.drawmeridians([-60, -55, -50, -45], labels=[0,0,0,1], fontsize=12, fontweight='normal');
            x, y = m(lons, lats)
            m.scatter(x,y, s=50, marker='.',color='k')
            cax = fig.add_axes([0.16, 0.05, 0.7, 0.025])
            cb = plt.colorbar(c, cax=cax, orientation='horizontal')
            cb.set_label(r'$\rm T(^{\circ}C)$', fontsize=12, fontweight='normal')
            #cax = plt.axes([0.85,0.15,0.04,0.7], facecolor='grey')
            #cb = plt.colorbar(c, cax=cax)
            #cb.set_label(r'$\rm T(^{\circ}C)$', fontsize=12, fontweight='normal')
            for div in div_toplot:
                div_lon, div_lat = m(nafo_div[div]['lon'], nafo_div[div]['lat'])
                m.plot(div_lon, div_lat, 'k', linewidth=2)
                ax.text(np.mean(div_lon), np.mean(div_lat), div, fontsize=12, color='black', fontweight='bold')
            # Save Figure
            fig.set_size_inches(w=7, h=8)
            fig.set_dpi(200)
            outfile = 'bottom_temp_' + season + '_' + str(year) + '.png'
            fig.savefig(outfile)
            plt.close('all')

    df_2G = pd.DataFrame.from_dict(dict_stats_2G, orient='index')
    df_2H = pd.DataFrame.from_dict(dict_stats_2H, orient='index')
    df_2J = pd.DataFrame.from_dict(dict_stats_2J, orient='index')
    df_2HJ = pd.DataFrame.from_dict(dict_stats_2HJ, orient='index')
    df_2GH = pd.DataFrame.from_dict(dict_stats_2GH, orient='index')
    df_3Ps = pd.DataFrame.from_dict(dict_stats_3Ps, orient='index')
    df_3LNO = pd.DataFrame.from_dict(dict_stats_3LNO, orient='index')
    df_3M = pd.DataFrame.from_dict(dict_stats_3M, orient='index')
    df_3K = pd.DataFrame.from_dict(dict_stats_3K, orient='index')
    df_3L = pd.DataFrame.from_dict(dict_stats_3L, orient='index')
    df_3O = pd.DataFrame.from_dict(dict_stats_3O, orient='index')
    df_4R = pd.DataFrame.from_dict(dict_stats_4R, orient='index')
    df_4S = pd.DataFrame.from_dict(dict_stats_4S, orient='index')
    df_4RS = pd.DataFrame.from_dict(dict_stats_4RS, orient='index')
    df_4RST = pd.DataFrame.from_dict(dict_stats_4RST, orient='index')
    df_4T = pd.DataFrame.from_dict(dict_stats_4T, orient='index')
    df_4VWX = pd.DataFrame.from_dict(dict_stats_4VWX, orient='index')
    #df_5Y = pd.DataFrame.from_dict(dict_stats_5Y, orient='index')

    outname = 'stats_3Ps_' + season + '.pkl'
    df_3Ps.to_pickle(outname)
    outname = 'stats_3LNO_' + season + '.pkl'
    df_3LNO.to_pickle(outname)
    outname = 'stats_3M_' + season + '.pkl'
    df_3M.to_pickle(outname)
    outname = 'stats_3K_' + season + '.pkl'
    df_3K.to_pickle(outname)
    outname = 'stats_3L_' + season + '.pkl'
    df_3L.to_pickle(outname)
    outname = 'stats_3O_' + season + '.pkl'
    df_3O.to_pickle(outname)
    outname = 'stats_2G_' + season + '.pkl'
    df_2G.to_pickle(outname)
    outname = 'stats_2H_' + season + '.pkl'
    df_2H.to_pickle(outname)
    outname = 'stats_2J_' + season + '.pkl'
    df_2J.to_pickle(outname)
    outname = 'stats_2HJ_' + season + '.pkl'
    df_2HJ.to_pickle(outname)
    outname = 'stats_2GH_' + season + '.pkl'
    df_2GH.to_pickle(outname)
    ## outname = 'stats_4R_' + season + '.pkl'
    ## df_4R.to_pickle(outname)
    ## outname = 'stats_4S_' + season + '.pkl'
    ## df_4S.to_pickle(outname)
    ## outname = 'stats_4RS_' + season + '.pkl'
    ## df_4RS.to_pickle(outname)
    ## outname = 'stats_4RST_' + season + '.pkl'
    ## df_4RST.to_pickle(outname)
    ## outname = 'stats_4T_' + season + '.pkl'
    ## df_4T.to_pickle(outname)
    ## outname = 'stats_4VWX_' + season + '.pkl'
    ## df_4VWX.to_pickle(outname)
    #outname = 'stats_5Y_' + season + '.pkl'
    #df_5Y.to_pickle(outname)

    # Save in multi-index  dataFrame
    year_index = pd.Series(years)
    year_index.name='year'
    df_mindex = pd.concat(df_list,keys=year_index)
    df_mindex.to_pickle(season + '_bottom_temperature.pkl')

    #### bottom_stats
def sfa_bottom_stats(years, season, proj='merc', plot=False, netcdf_path='/home/cyrf0006/data/dev_database/netCDF/', sfas=[2,3,4], climato_file='', clim_fill=False, plot_biomass=False, plot_biomass_anomaly=False, var='temperature'):

    '''
    Function sfa_bottom_stats() is based bottom_stats(), but sfa instead of NAFO divs.

    Run in:
    /home/cyrf0006/AZMP/state_reports/bottomT/NSRF

    Climato obtained with:
    Tbot_dict = azu.get_bottomT_climato()

    usage example:
    import azmp_report_tools as azrt
    import numpy as np
    azrt.sfa_bottom_stats(years=np.arange(2006, 2022), season='summer', climato_file='Tbot_climato_NSRFx_summer_2006-2021.h5')

    *** This needs to be improve because at the moment I need to comment the generation of .pkl file to not over-write when I change my map region.        

    April 2022: modified to use climatology to fill missing pixels during a specific year
    azrt.sfa_bottom_stats(years=np.arange(2006, 2022), season='summer', climato_file='Tbot_climato_NSRFx_summer_2006-2021.h5', clim_fill=True)

    Nov. 2022: Modified to add Pandalus biomass on plot:
    azrt.sfa_bottom_stats(years=np.arange(2006, 2022), season='summer', climato_file='Tbot_climato_NSRFx_summer_2006-2021.h5', plot=True, clim_fill=True, plot_biomass=True)

    Feb. 2023: Modified to also compute botto salinity stats (now variable option)
        Note that variable Tbot now refers to temperature or salinity in the function
    azrt.sfa_bottom_stats(years=np.arange(2006, 2022), season='summer', climato_file='Sbot_climato_NSRFx_summer_2006-2021.h5', plot=True, clim_fill=True, plot_biomass=True ,var='salinity')

    Frederic.Cyr@dfo-mpo.gc.ca - January 2021
    '''

    '''
    #TEMPORARY, FOR TESTING PURPOSES
    years=np.arange(2006,2023)
    season='fall'
    climato_file='Tbot_climato_NSRF_'+season+'_0.10.h5'
    plot=True
    var='temperature'
    clim_fill=True
    plot_biomass=True
    '''

    # load climato
    if len(climato_file) == 0: # climato not provided (default)
        if var=='salinity':
            print('ERROR! Please provide climatology file')
            return
        if season == 'fall':
            climato_file = 'Tbot_climato_NSRF_fall.h5'
        elif season == 'spring':
            climato_file = 'Tbot_climato_NSRF_spring.h5'
        elif season == 'summer':
            climato_file = 'Tbot_climato_NSRF_summer.h5'
    else:
        print('Climato file provided')
               
    h5f = h5py.File(climato_file, 'r')
    if var == 'temperature':
        Tbot_climato = h5f['Tbot'][:]
    elif var == 'salinity':
        Tbot_climato = h5f['Sbot'][:]
    lon_reg = h5f['lon_reg'][:]
    lat_reg = h5f['lat_reg'][:]
    Zitp = h5f['Zitp'][:]
    h5f.close()

    # Derive some map parameters
    lon_0 = np.round(np.mean(lon_reg))
    lat_0 = np.round(np.mean(lat_reg))
    lonLims = [lon_reg[0], lon_reg[-1]]
    latLims = [lat_reg[0], lat_reg[-1]]

    ## get Pandalus biomass data
    df_biomass = pd.read_excel('/home/cyrf0006/research/Pandalus_project/NSRF_PBPM_biomass_standardized.xlsx')
    df_biomass.set_index('Year', inplace=True)
    df_biomass = df_biomass[['Latitude', 'Longitude', 'PB TB(kg/km2)', 'PM TB(kg/km2)']]
    df_biomass = df_biomass.replace(0.0, 1) # replace zeros by ones for log transform

    ## Get SFAs data
    myshp = open('/home/cyrf0006/github/AZMP-NL/utils/SFAs/SFAs_PANOMICS_Fall2020_shp/SFAs_PANOMICS_Fall2020.shp', 'rb')
    mydbf = open('/home/cyrf0006/github/AZMP-NL/utils/SFAs/SFAs_PANOMICS_Fall2020_shp/SFAs_PANOMICS_Fall2020.dbf', 'rb')
    r = shapefile.Reader(shp=myshp, dbf=mydbf, encoding = "ISO8859-1")
    records = r.records()
    shapes = r.shapes()
    
    # Fill dictionary with shapes
    shrimp_area = {}
    for idx, rec in enumerate(records):
        if rec[1] == 'Eastern Assessment Zone':
            shrimp_area['2'] = np.array(shapes[idx].points)
        elif rec[1] == 'Western Assessment Zone':
            shrimp_area['3'] = np.array(shapes[idx].points)
        else:
            shrimp_area[rec[0]] = np.array(shapes[idx].points)            

    ## sfa0 = Polygon(shrimp_area['0'])
    ## sfa1 = Polygon(shrimp_area['1'])
    sfa2 = Polygon(shrimp_area['2'])
    sfa3 = Polygon(shrimp_area['3'])
    sfa4 = Polygon(shrimp_area['4'])
    sfa5 = Polygon(shrimp_area['5'])
    sfa6 = Polygon(shrimp_area['6'])
    sfa7 = Polygon(shrimp_area['7'])

    ## dict_stats_sfa0 = {}
    ## dict_stats_sfa1 = {}
    dict_stats_sfa2 = {}
    dict_stats_sfa3 = {}
    dict_stats_sfa4 = {}
    dict_stats_sfa5 = {}
    dict_stats_sfa6 = {}
    dict_stats_sfa7 = {}
   
    # Loop on years
    df_list = []
    for year in years:
        print(' ---- ' + np.str(year) + ' ---- ')
        year_file = netcdf_path + np.str(year) + '.nc'
        if var == 'temperature':
            Tdict = azu.get_bottomT(year_file, season, climato_file)    
            Tbot = Tdict['Tbot']
        else:
            Tdict = azu.get_bottomS(year_file, season, climato_file)    
            Tbot = Tdict['Sbot']
        lons = Tdict['lons']
        lats = Tdict['lats']

        # Use climatology to fill missing pixels
        if clim_fill:
            print('Fill NaNs with climatology')
            Tbot_orig = Tbot.copy()
            #Tbot_fill = Tbot_climato[Tbot_orig.isna()]
            #Tbot[Tbot_orig.isna()] = Tbot_climato[Tbot_orig.isna()]
            Tbot_fill = Tbot_climato[np.isnan(Tbot_orig)]
            Tbot[np.isnan(Tbot_orig)] = Tbot_climato[np.isnan(Tbot_orig)]            
        # Calculation of anomaly
        anom = Tbot-Tbot_climato
        
        # NAFO division stats
        if season == 'summer':
            # To use different definition for each SFAs
            ## dict_stats_sfa2[np.str(year)] = azu.polygon_temperature_stats(Tdict, sfa2, nsrf=True)
            ## dict_stats_sfa3[np.str(year)] = azu.polygon_temperature_stats(Tdict, sfa3, nsrf=True)
            ## dict_stats_sfa4[np.str(year)] = azu.polygon_temperature_stats(Tdict, sfa4, nsrf=True)
            dict_stats_sfa2[np.str(year)] = azu.polygon_temperature_stats(Tdict, sfa2, var=var)
            dict_stats_sfa3[np.str(year)] = azu.polygon_temperature_stats(Tdict, sfa3, var=var)
            dict_stats_sfa4[np.str(year)] = azu.polygon_temperature_stats(Tdict, sfa4, var=var)
        elif season=='fall':
            dict_stats_sfa4[np.str(year)] = azu.polygon_temperature_stats(Tdict, sfa4, var=var)
            dict_stats_sfa5[np.str(year)] = azu.polygon_temperature_stats(Tdict, sfa5, var=var)
            dict_stats_sfa6[np.str(year)] = azu.polygon_temperature_stats(Tdict, sfa6, var=var)
            dict_stats_sfa7[np.str(year)] = azu.polygon_temperature_stats(Tdict, sfa7, var=var)
        elif season == 'spring':
            print('Not implemented, check this!')
            keyboard

        # Append bottom temperature for multi-index export
        df = pd.DataFrame(index=lat_reg, columns=lon_reg)
        df.index.name='latitude'
        df.columns.name='longitude'
        df[:] = Tbot
        df_list.append(df)

        if plot:
            div_toplot = [ '2', '3', '4']
            if var == 'temperature':
                levels = np.linspace(-2, 6, 9)
                levels_anom = np.linspace(-3.5, 3.5, 8)
                CMAP = cmo.cm.thermal
            elif var == 'salinity':
                levels = np.linspace(30, 35, 11)
                levels_anom = np.array([-1, -.8, -.6, -.4, -.2, .2, .4, .6, .8, 1])
                CMAP = cmo.cm.haline

            if year == years[0]: 
            # 1.0 - Plot Climatology
                fig, ax = plt.subplots(nrows=1, ncols=1)
                m = Basemap(ax=ax, projection='merc',lon_0=lon_0,lat_0=lat_0, llcrnrlon=lonLims[0],llcrnrlat=latLims[0],urcrnrlon=lonLims[1],urcrnrlat=latLims[1], resolution= 'i')
                xi, yi = m(*np.meshgrid(lon_reg, lat_reg))
                cclim = m.contourf(xi, yi, Tbot_climato, levels, cmap=CMAP, extend='both')
                cc = m.contour(xi, yi, -Zitp, [100, 500, 1000, 4000], colors='grey');
                plt.clabel(cc, inline=1, fontsize=10, fmt='%d')
                m.fillcontinents(color='tan');
                m.drawparallels([56, 58, 60, 62, 64, 66], labels=[1,0,0,0], fontsize=12, fontweight='normal');
                m.drawmeridians([-68, -66, -64, -62, -60, -58, -56], labels=[0,0,0,1], fontsize=12, fontweight='normal');
                if var == 'temperature':
                    if season=='fall':
                        plt.title('Fall Bottom Temperature')
                    elif season=='spring':
                        plt.title('Spring Bottom Temperature')
                    else:
                        plt.title('Bottom Temperature')
                elif var == 'salinity':
                    if season=='fall':
                        plt.title('Fall Bottom Salinity')
                    elif season=='spring':
                        plt.title('Spring Bottom Salinity')
                    else:
                        plt.title('Bottom Salinity')
                # add colorbar
                cax = fig.add_axes([0.25, 0.055, 0.5, 0.02])
                cb = plt.colorbar(cclim, cax=cax, orientation='horizontal')
                if var == 'temperature':
                    cb.set_label(r'$\rm T(^{\circ}C)$', fontsize=12, fontweight='normal')
                elif var == 'salinity':
                    cb.set_label('S', fontsize=12, fontweight='normal')
                # plot SFA divisons
                for div in div_toplot:
                    div_lon, div_lat = m(shrimp_area[div][:,0], shrimp_area[div][:,1])
                    m.plot(div_lon, div_lat, 'k', linewidth=2)
                    ax.text(np.mean(div_lon), np.mean(div_lat), 'SFA'+div, fontsize=12, color='black', fontweight='bold')
 
                # Save Figure
                fig.set_size_inches(w=7, h=8)
                fig.set_dpi(200)
                outfile = 'bottom_clim_' + season  + '.png'
                fig.savefig(outfile)
                plt.close('all')
                os.system('convert -trim ' + outfile + ' ' + outfile)

            # 1.1 - Plot Anomaly
            fig, ax = plt.subplots(nrows=1, ncols=1)
            m = Basemap(ax=ax, projection='merc',lon_0=lon_0,lat_0=lat_0, llcrnrlon=lonLims[0],llcrnrlat=latLims[0],urcrnrlon=lonLims[1],urcrnrlat=latLims[1], resolution= 'i')
            xi, yi = m(*np.meshgrid(lon_reg, lat_reg))
            c = m.contourf(xi, yi, anom, levels_anom, cmap=plt.cm.RdBu_r, extend='both')
            cc = m.contour(xi, yi, -Zitp, [100, 500, 1000, 4000], colors='grey');
            plt.clabel(cc, inline=1, fontsize=10, fmt='%d')
            m.fillcontinents(color='tan');
            m.drawparallels([56, 58, 60, 62, 64, 66], labels=[1,0,0,0], fontsize=12, fontweight='normal');
            m.drawmeridians([-68, -66, -64, -62, -60, -58, -56], labels=[0,0,0,1], fontsize=12, fontweight='normal');
            if var == 'temperature':
                if season=='fall':
                    plt.title('Fall Bottom Temperature Anomaly')
                elif season=='spring':
                    plt.title('Spring Bottom Temperature Anomaly')
                else:
                    plt.title('Bottom Temperature Anomaly')
            elif var == 'salinity':
                if season=='fall':
                    plt.title('Fall Bottom Salinity Anomaly')
                elif season=='spring':
                    plt.title('Spring Bottom Salinity Anomaly')
                else:
                    plt.title('Bottom Salinity Anomaly')
            # add colorbar
            cax = fig.add_axes([0.25, 0.055, 0.5, 0.02])
            cb = plt.colorbar(c, cax=cax, orientation='horizontal')
            if var == 'temperature':
                cb.set_label(r'$\rm T(^{\circ}C)$', fontsize=12, fontweight='normal')
            elif var == 'salinity':
                cb.set_label('S', fontsize=12, fontweight='normal')
            # plot SFA divisons
            for div in div_toplot:
                div_lon, div_lat = m(shrimp_area[div][:,0], shrimp_area[div][:,1])
                m.plot(div_lon, div_lat, 'k', linewidth=2)
                ax.text(np.mean(div_lon), np.mean(div_lat), 'SFA'+div, fontsize=12, color='black', fontweight='bold')
            # Add Biomass
            if plot_biomass_anomaly:
                biomass = df_biomass[df_biomass.index==int(year)]
                xb, yb = m(biomass.Longitude.values, biomass.Latitude.values)
                for idx, tmp in enumerate(xb):
                    #m.scatter(xb[idx], yb[idx], s=np.log10(biomass['PB TB(kg/km2)'].iloc[idx])*30, c='orange', alpha=.5, zorder=200)
                    #m.scatter(xb[idx], yb[idx], s=np.log10(biomass['PM TB(kg/km2)'].iloc[idx])*30, c='slategray', alpha=.5, zorder=200)
                    m.scatter(xb[idx], yb[idx], s=biomass['PB TB(kg/km2)'].iloc[idx]/50, c='orange', alpha=.5)
                    m.scatter(xb[idx], yb[idx], s=biomass['PM TB(kg/km2)'].iloc[idx]/50, c='darkgray', alpha=.5)
                # add legend
                xbl, ybl = m(-69.5, 66.8)
                xml, yml = m(-69.5, 66.5)
                m.scatter(xbl, ybl, s=60, c='orange', alpha=.7, zorder=200)
                m.scatter(xml, yml, s=60, c='slategray', alpha=.7, zorder=200)
                #ax.text(xbl, ybl, r'  P. Borealis (100 $\rm kg\,km^{-2}$)', zorder=200, verticalalignment='center')
                #ax.text(xml, yml, r'  P. Montagui (100 $\rm kg\,km^{-2}$)', zorder=200, verticalalignment='center')
                ax.text(xbl, ybl, r'  P. Borealis (3000 $\rm kg\,km^{-2}$)', zorder=200, verticalalignment='center')
                ax.text(xml, yml, r'  P. Montagui (3000 $\rm kg\,km^{-2}$)', zorder=200, verticalalignment='center')
            # Save Figure
            fig.set_size_inches(w=7, h=8)
            fig.set_dpi(200)
            outfile1 = 'bottom_' + var + '_anomaly_' + season + '_' + np.str(year) + '.png'
            fig.savefig(outfile1)
            os.system('convert -trim ' + outfile1 + ' ' + outfile1)

            # 1.2 - Plot Temperature
            fig, ax = plt.subplots(nrows=1, ncols=1)
            m = Basemap(ax=ax, projection='merc',lon_0=lon_0,lat_0=lat_0, llcrnrlon=lonLims[0],llcrnrlat=latLims[0],urcrnrlon=lonLims[1],urcrnrlat=latLims[1], resolution= 'i')
            xi, yi = m(*np.meshgrid(lon_reg, lat_reg))
            if clim_fill:
                #c = m.contourf(xi, yi, Tbot, levels, cmap=plt.cm.RdBu_r, extend='both', alpha=.3)
                #cfill = m.contourf(xi, yi, Tbot_orig, levels, cmap=plt.cm.RdBu_r, extend='both', alpha=1)
                c = m.contourf(xi, yi, Tbot, levels, cmap=CMAP, extend='both', alpha=.3)
                cfill = m.contourf(xi, yi, Tbot_orig, levels, cmap=CMAP, extend='both', alpha=.8)
            else:
                cfill = m.contourf(xi, yi, Tbot, levels, cmap=CMAP, extend='both')
            cc = m.contour(xi, yi, -Zitp, [100, 500, 1000, 4000], colors='grey');
            plt.clabel(cc, inline=1, fontsize=10, fmt='%d')
            m.fillcontinents(color='tan');
            m.drawparallels([56, 58, 60, 62, 64, 66], labels=[1,0,0,0], fontsize=12, fontweight='normal');
            m.drawmeridians([-68, -66, -64, -62, -60, -58, -56], labels=[0,0,0,1], fontsize=12, fontweight='normal');
            x, y = m(lons, lats)
            m.scatter(x,y, s=50, marker='.',color='k')
            if var == 'temperature':
                if season=='fall':
                    plt.title('Fall Bottom Temperature ' + str(year))
                elif season=='spring':
                    plt.title('Spring Bottom Temperature ' + str(year))
                else:
                    plt.title('Bottom Temperature ' + str(year))
            elif var == 'salinity':
                if season=='fall':
                    plt.title('Fall Bottom Salinity ' + str(year))
                elif season=='spring':
                    plt.title('Spring Bottom Salinity ' + str(year))
                else:
                    plt.title('Bottom Salinity ' + str(year))
            # add colorbar
            cax = fig.add_axes([0.25, 0.055, 0.5, 0.02])
            cb = plt.colorbar(cclim, cax=cax, orientation='horizontal')
            if var == 'temperature':
                cb.set_label(r'$\rm T(^{\circ}C)$', fontsize=12, fontweight='normal')
            elif var == 'salinity':
                cb.set_label('S', fontsize=12, fontweight='normal')                    
            # plot SFA divisons
            for div in div_toplot:
                div_lon, div_lat = m(shrimp_area[div][:,0], shrimp_area[div][:,1])
                m.plot(div_lon, div_lat, 'k', linewidth=2)
                ax.text(np.mean(div_lon), np.mean(div_lat), 'SFA'+div, fontsize=12, color='black', fontweight='bold')
            # Add Biomass
            if plot_biomass:
                biomass = df_biomass[df_biomass.index==int(year)]
                xb, yb = m(biomass.Longitude.values, biomass.Latitude.values)
                for idx, tmp in enumerate(xb):
                    #m.scatter(xb[idx], yb[idx], s=np.log10(biomass['PB TB(kg/km2)'].iloc[idx])*30, c='red', alpha=.5, zorder=200)
                    #m.scatter(xb[idx], yb[idx], s=np.log10(biomass['PM TB(kg/km2)'].iloc[idx])*30, c='slategray', alpha=.5, zorder=200)
                    m.scatter(xb[idx], yb[idx], s=biomass['PB TB(kg/km2)'].iloc[idx]/50, c='red', alpha=.5)
                    m.scatter(xb[idx], yb[idx], s=biomass['PM TB(kg/km2)'].iloc[idx]/50, c='dimgray', alpha=.5)
                # add legend
                xbl, ybl = m(-69.5, 66.8)
                xml, yml = m(-69.5, 66.5)
                m.scatter(xbl, ybl, s=60, c='red', alpha=.7, zorder=200)
                m.scatter(xml, yml, s=60, c='slategray', alpha=.7, zorder=200)
                #ax.text(xbl, ybl, r'  P. Borealis (100 $\rm kg\,km^{-2}$)', zorder=200, verticalalignment='center')
                #ax.text(xml, yml, r'  P. Montagui (100 $\rm kg\,km^{-2}$)', zorder=200, verticalalignment='center')
                ax.text(xbl, ybl, r'  P. Borealis (3000 $\rm kg\,km^{-2}$)', zorder=200, verticalalignment='center')
                ax.text(xml, yml, r'  P. Montagui (3000 $\rm kg\,km^{-2}$)', zorder=200, verticalalignment='center')
            # Save Figure
            fig.set_size_inches(w=7, h=8)
            fig.set_dpi(200)
            if clim_fill:
                outfile2 = 'bottom_' + var + '_filled_' + season + '_' + np.str(year) + '.png'
            else:
                outfile2 = 'bottom_' + var + '_' + season + '_' + np.str(year) + '.png'
            fig.savefig(outfile2)
            plt.close('all')
            os.system('convert -trim ' + outfile2 + ' ' + outfile2)

            # Montage figure
            if var == 'temperature':
                montage_out = 'sfa_bottomT_' + str(year) + '.png'
            elif var == 'salinity':
                montage_out = 'sfa_bottomS_' + str(year) + '.png'
            os.system('montage ' + outfile + ' ' + outfile2 + ' ' + outfile1 + ' ' + '-tile 3x1 -geometry +10+10  -background white ' + montage_out)
            os.system('rm ' + outfile2 + ' ' + outfile1) 
          
            # clear memory
            del Tdict, Tbot, anom

    ## df_sfa0 = pd.DataFrame.from_dict(dict_stats_sfa0, orient='index')
    ## df_sfa1 = pd.DataFrame.from_dict(dict_stats_sfa1, orient='index')
    df_sfa2 = pd.DataFrame.from_dict(dict_stats_sfa2, orient='index')
    df_sfa3 = pd.DataFrame.from_dict(dict_stats_sfa3, orient='index')
    df_sfa4 = pd.DataFrame.from_dict(dict_stats_sfa4, orient='index')
    df_sfa5 = pd.DataFrame.from_dict(dict_stats_sfa5, orient='index')
    df_sfa6 = pd.DataFrame.from_dict(dict_stats_sfa6, orient='index')
    df_sfa7 = pd.DataFrame.from_dict(dict_stats_sfa7, orient='index')

    ## outname = 'stats_sfa0_' + season + '.pkl' # Some are missing here
    ## df_sfa0.to_pickle(outname)
    ## outname = 'stats_sfa1_' + season + '.pkl'
    ## df_sfa1.to_pickle(outname)
    if var == 'salinity':
        season = season + '_' + var
    outname = 'stats_sfa2_' + season + '.pkl'
    df_sfa2.to_pickle(outname)
    outname = 'stats_sfa3_' + season + '.pkl'
    df_sfa3.to_pickle(outname)   
    outname = 'stats_sfa4_' + season + '.pkl'
    df_sfa4.to_pickle(outname)
    outname = 'stats_sfa5_' + season + '.pkl'
    df_sfa5.to_pickle(outname)   
    outname = 'stats_sfa6_' + season + '.pkl'
    df_sfa6.to_pickle(outname)
    outname = 'stats_sfa7_' + season + '.pkl'
    df_sfa7.to_pickle(outname)   

    # Save in multi-index  dataFrame (still useful?)
    #year_index = pd.Series(years)
    #year_index.name='year'
    #df_mindex = pd.concat(df_list,keys=year_index)
    #df_mindex.to_pickle(season + '_sfa_bottom_temperature.pkl')


def bottom_scorecards(years, clim_year=[1991, 2020]):


    '''
    To generate AZMP score cards for bottom temperature
    Uses pickled object generated by bottom_stats()

    usage example:


    bottom_scorecards(years=[1980, 2019], clim_year=[1981, 2010]):

    '''

    #### ------------- For fall ---------------- ####
    # 0.
    infile = 'stats_2H_fall.pkl'
    df = pd.read_pickle(infile)
    df.index = pd.to_datetime(df.index) # update index to datetime
    df = df[(df.index.year>=years[0]) & (df.index.year<=years[-1])]
    # Flag bad years (no or weak sampling):
    bad_years = np.array([1980, 1982, 1984, 1985, 1986, 1987, 1988, 1989, 1990, 1992, 1993, 1994, 1995, 1996, 2000, 2002, 2003, 2005, 2007, 2009, 2022])
    for i in bad_years:
        df[df.index.year==i]=np.nan
    year_list = df.index.year.astype('str')
    year_list = [i[2:4] for i in year_list] # 2-digit year
    df['area_colder0'] = df['area_colder0']/1000 # In 1000km
    df['area_colder1'] = df['area_colder1']/1000 # In 1000km
    df['area_warmer2'] = df['area_warmer2']/1000
    df_clim = df[(df.index.year>=clim_year[0]) & (df.index.year<=clim_year[1])]
    std_anom = (df-df_clim.mean(axis=0))/df_clim.std(axis=0)
    # Save in .csv for future use
    std_anom.to_csv('bottomT_stn_anom_2H_fall.csv', sep=',', float_format='%0.3f')
    # Transpose for table
    std_anom = std_anom.T
    std_anom['MEAN'] = df_clim.mean(axis=0)
    std_anom['SD'] = df_clim.std(axis=0)
    std_anom = std_anom.reindex(['Tmean', 'Tmean_sha200', 'area_warmer2', 'area_colder1'])
    std_anom = std_anom.rename({'Tmean': r'$\rm T_{bot}$', 'Tmean_sha200': r'$\rm T_{bot_{<200m}}$', 'area_warmer2': r'$\rm Area_{>2^{\circ}C}$', 'area_colder1': r'$\rm Area_{<1^{\circ}C}$'})
    std_anom.rename(columns={'MEAN': r'$\rm \overline{x}$', 'SD': r'sd'}, inplace=True)
    
    # Get text values +  cell color
    year_list.append(r'$\rm \overline{x}$') # add 2 extra columns
    year_list.append(r'sd')   
    vals = np.around(std_anom.values,1)
    vals[vals==-0.] = 0.
    vals_color = vals.copy()
    vals_color[-1,] = vals_color[-1,]*-1 # Reverse last row colorscale
    vals_color[:,-1] = 0 # No color to last two columns (mean and STD)
    vals_color[:,-2] = 0
    #vals_color[(vals_color<0.5) & (vals_color>-.5)] = 0.

    # Build the colormap
    vmin = -3.49
    vmax = 3.49
    midpoint = 0
    levels = np.linspace(vmin, vmax, 15)
    midp = np.mean(np.c_[levels[:-1], levels[1:]], axis=1)
    colvals = np.interp(midp, [vmin, midpoint, vmax], [-1, 0., 1])
    normal = plt.Normalize(-3.49, 3.49)
    reds = plt.cm.Reds(np.linspace(0,1, num=7))
    blues = plt.cm.Blues_r(np.linspace(0,1, num=7))
    whites = [(1,1,1,1)]*2
    colors = np.vstack((blues[0:-1,:], whites, reds[1:,:]))
    colors = np.concatenate([[colors[0,:]], colors, [colors[-1,:]]], 0)
    cmap, norm = from_levels_and_colors(levels, colors, extend='both')
    cmap_r, norm_r = from_levels_and_colors(levels, np.flipud(colors), extend='both')
    # Common parameters
    #hcell, wcell = 0.5, 0.6
    #hpad, wpad = 0, 0

    ## normal = plt.Normalize(-4.49, 4.49)
    ## cmap = plt.cm.get_cmap('seismic', 9) 
    #cmap = plt.cm.get_cmap('seismic', 15) 

    nrows, ncols = std_anom.index.size+1, std_anom.columns.size
    hcell, wcell = 0.5, 0.5
    hpad, wpad = 1, 1    
    fig=plt.figure(figsize=(ncols*wcell+wpad, nrows*hcell+hpad))
    ax = fig.add_subplot(111)
    ax.axis('off')
    #do the table
    header = ax.table(cellText=[['']],
                          colLabels=['-- NAFO division 2H --'],
                          loc='center'
                          )
    header.set_fontsize(13)
    #the_table=ax.table(cellText=valsbottom_clim_summer.png, rowLabels=std_anom.index, colLabels=std_anom.columns, 
    the_table=ax.table(cellText=vals, rowLabels=std_anom.index, colLabels=year_list,
                        loc='center', cellColours=cmap(normal(vals_color)), cellLoc='center',
                        bbox=[0, 0, 1, 0.5]
                        )
    # change font color to white where needed:
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(12.5)
    table_props = the_table.properties()
    #table_cells = table_props['child_artists']
    last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
    for key, cell in the_table.get_celld().items():
        cell_text = cell.get_text().get_text()
        if is_number(cell_text) == False:
            pass
        elif key[0] == 0: #year's row = no color
            pass
        elif key[1] in last_columns:
             cell._text.set_color('darkslategray')
        elif (float(cell_text) <= -1.5) | (float(cell_text) >= 1.5) :
            cell._text.set_color('white')
        elif (cell_text=='nan'):
            cell._set_facecolor('lightgray')
            cell._text.set_color('lightgray')

    plt.savefig("scorecards_fall_2H.png", dpi=300)
    os.system('convert -trim scorecards_fall_2H.png scorecards_fall_2H.png')

    # French table
    std_anom = std_anom.rename({r'$\rm T_{bot}$' : r'$\rm T_{fond}$', r'$\rm T_{bot_{<200m}}$' : r'$\rm T_{fond_{<200m}}$', r'$\rm Area_{>2^{\circ}C}$' : r'$\rm Aire_{>2^{\circ}C}$', r'$\rm Area_{<1^{\circ}C}$' : r'$\rm Aire_{<1^{\circ}C}$'})
    year_list[-1] = u'ET'

    header = ax.table(cellText=[['']],
                          colLabels=['-- Division 2H de l\'OPANO --'],
                          loc='center'
                          )
    header.set_fontsize(13)
    #the_table=ax.table(cellText=vals, rowLabels=std_anom.index, colLabels=std_anom.columns, 
    the_table=ax.table(cellText=vals, rowLabels=std_anom.index, colLabels=year_list,
                        loc='center', cellColours=cmap(normal(vals_color)), cellLoc='center',
                        bbox=[0, 0, 1, 0.5]
                        )
    # change font color to white where needed:
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(12.5)
    table_props = the_table.properties()
    #table_cells = table_props['child_artists']
    last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
    for key, cell in the_table.get_celld().items():
        cell_text = cell.get_text().get_text()
        if is_number(cell_text) == False:
            pass
        elif key[0] == 0: #year's row = no color
            pass
        elif key[1] in last_columns:
             cell._text.set_color('darkslategray')
        elif (float(cell_text) <= -1.5) | (float(cell_text) >= 1.5) :
            cell._text.set_color('white')
        elif (cell_text=='nan'):
            cell._set_facecolor('lightgray')
            cell._text.set_color('lightgray')

    plt.savefig("scorecards_fall_2H_FR.png", dpi=300)
    os.system('convert -trim scorecards_fall_2H_FR.png scorecards_fall_2H_FR.png')

 # 1.
    infile = 'stats_2J_fall.pkl'
    df = pd.read_pickle(infile)
    df.index = pd.to_datetime(df.index) # update index to datetime
    df = df[(df.index.year>=years[0]) & (df.index.year<=years[-1])]
    # Flag bad years (no or weak sampling):
    bad_years = np.array([1995, 2022])
    for i in bad_years:
        df[df.index.year==i]=np.nan
    df['area_colder0'] = df['area_colder0']/1000 # In 1000km
    df['area_colder1'] = df['area_colder1']/1000 # In 1000km
    df['area_warmer2'] = df['area_warmer2']/1000
    df_clim = df[(df.index.year>=clim_year[0]) & (df.index.year<=clim_year[1])]
    std_anom = (df-df_clim.mean(axis=0))/df_clim.std(axis=0)
    # Save in .csv for future use
    std_anom.to_csv('bottomT_stn_anom_2J_fall.csv', sep=',', float_format='%0.3f')
    # Transpose for table
    std_anom = std_anom.T
    std_anom['MEAN'] = df_clim.mean(axis=0)
    std_anom['SD'] = df_clim.std(axis=0)
    std_anom = std_anom.reindex(['Tmean', 'Tmean_sha200', 'area_warmer2', 'area_colder1'])
    std_anom = std_anom.rename({'Tmean': r'$\rm T_{bot}$', 'Tmean_sha200': r'$\rm T_{bot_{<200m}}$', 'area_warmer2': r'$\rm Area_{>2^{\circ}C}$', 'area_colder1': r'$\rm Area_{<1^{\circ}C}$'})
    std_anom.rename(columns={'MEAN': r'$\rm \overline{x}$', 'SD': r'sd'}, inplace=True)

    vals = np.around(std_anom.values,1)
    vals[vals==-0.] = 0.
    vals_color = vals.copy()
    vals_color[-1,] = vals_color[-1,]*-1
    vals_color[:,-1] = 0 # No color to last two columns (mean and STD)
    vals_color[:,-2] = 0 
    #normal = plt.Normalize(-4.49, 4.49)
    #cmap = plt.cm.get_cmap('seismic', 9) 
    nrows, ncols = std_anom.index.size, std_anom.columns.size
    fig=plt.figure(figsize=(ncols*wcell+wpad, nrows*hcell+hpad))
    ax = fig.add_subplot(111)
    ax.axis('off')
    #do the table
    header = ax.table(cellText=[['']],
                          colLabels=['-- NAFO division 2J --'],
                          loc='center'
                          )
    header.set_fontsize(12.5)
    the_table=ax.table(cellText=vals, rowLabels=std_anom.index, colLabels=None, 
                        loc='center', cellColours=cmap(normal(vals_color)), cellLoc='center',
                        bbox=[0, 0, 1.0, 0.50]
                        )
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(12.5)
    # change font color to white where needed:
    table_props = the_table.properties()
    #table_cells = table_props['child_artists']
    last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
    for key, cell in the_table.get_celld().items():
        cell_text = cell.get_text().get_text()
        if is_number(cell_text) == False:
            pass
        #elif key[0] == 0:# <--- remove when no years
        #    pass
        elif key[1] in last_columns:
             cell._text.set_color('darkslategray')
        elif (float(cell_text) <= -1.5) | (float(cell_text) >= 1.5) :
            cell._text.set_color('white')
        elif (cell_text=='nan'):
            cell._set_facecolor('lightgray')
            cell._text.set_color('lightgray')

    plt.savefig("scorecards_fall_2J.png", dpi=300)
    os.system('convert -trim scorecards_fall_2J.png scorecards_fall_2J.png')

    # French table
    std_anom = std_anom.rename({r'$\rm T_{bot}$' : r'$\rm T_{fond}$', r'$\rm T_{bot_{<200m}}$' : r'$\rm T_{fond_{<200m}}$', r'$\rm Area_{>2^{\circ}C}$' : r'$\rm Aire_{>2^{\circ}C}$', r'$\rm Area_{<1^{\circ}C}$' : r'$\rm Aire_{<1^{\circ}C}$'})
    header = ax.table(cellText=[['']],
                          colLabels=['-- Division 2J de l\'OPANO --'],
                          loc='center'
                          )
    #the_table=ax.table(cellText=vals, rowLabels=std_anom.index, colLabels=std_anom.columns, 
    header.set_fontsize(12.5)
    the_table=ax.table(cellText=vals, rowLabels=std_anom.index, colLabels=None, 
                        loc='center', cellColours=cmap(normal(vals_color)), cellLoc='center',
                        bbox=[0, 0, 1.0, 0.50]
                        )
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(12.5)
    # change font color to white where needed:
    table_props = the_table.properties()
    #table_cells = table_props['child_artists']
    last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
    for key, cell in the_table.get_celld().items():
        cell_text = cell.get_text().get_text()
        if is_number(cell_text) == False:
            pass
        #elif key[0] == 0:# <--- remove when no years
        #    pass
        elif key[1] in last_columns:
             cell._text.set_color('darkslategray')
        elif (float(cell_text) <= -1.5) | (float(cell_text) >= 1.5) :
            cell._text.set_color('white')
        elif (cell_text=='nan'):
            cell._set_facecolor('lightgray')
            cell._text.set_color('lightgray')

    plt.savefig("scorecards_fall_2J_FR.png", dpi=300)
    os.system('convert -trim scorecards_fall_2J_FR.png scorecards_fall_2J_FR.png')


    
    # 2.
    infile = 'stats_3K_fall.pkl'
    df = pd.read_pickle(infile)
    df.index = pd.to_datetime(df.index) # update index to datetime
    df = df[(df.index.year>=years[0]) & (df.index.year<=years[-1])]
    # Flag bad years (no or weak sampling):
    ## bad_years = np.array([2003])
    ## for i in bad_years:
    ##     df[df.index.year==i]=np.nan
    df['area_colder0'] = df['area_colder0']/1000 # In 1000km
    df['area_colder1'] = df['area_colder1']/1000 # In 1000km
    df['area_warmer2'] = df['area_warmer2']/1000
    df_clim = df[(df.index.year>=clim_year[0]) & (df.index.year<=clim_year[1])]
    std_anom = (df-df_clim.mean(axis=0))/df_clim.std(axis=0)
    # Save in .csv for future use
    std_anom.to_csv('bottomT_stn_anom_3K_fall.csv', sep=',', float_format='%0.3f')
    # Transpose for table
    std_anom = std_anom.T
    std_anom['MEAN'] = df_clim.mean(axis=0)
    std_anom['SD'] = df_clim.std(axis=0)
    std_anom = std_anom.reindex(['Tmean', 'Tmean_sha200', 'area_warmer2', 'area_colder1'])
    std_anom = std_anom.rename({'Tmean': r'$\rm T_{bot}$', 'Tmean_sha200': r'$\rm T_{bot_{<200m}}$', 'area_warmer2': r'$\rm Area_{>2^{\circ}C}$', 'area_colder1': r'$\rm Area_{<1^{\circ}C}$'})
    std_anom.rename(columns={'MEAN': r'$\rm \overline{x}$', 'SD': r'sd'}, inplace=True)

    vals = np.around(std_anom.values,1)
    vals[vals==-0.] = 0.
    vals_color = vals.copy()
    vals_color[-1,] = vals_color[-1,]*-1
    vals_color[:,-1] = 0 # No color to last two columns (mean and STD)
    vals_color[:,-2] = 0 
    #normal = plt.Normalize(-4.49, 4.49)
    #cmap = plt.cm.get_cmap('seismic', 9) 
    nrows, ncols = std_anom.index.size, std_anom.columns.size
    fig=plt.figure(figsize=(ncols*wcell+wpad, nrows*hcell+hpad))
    ax = fig.add_subplot(111)
    ax.axis('off')
    #do the table
    header = ax.table(cellText=[['']],
                          colLabels=['-- NAFO division 3K --'],
                          loc='center'
                          )
    #the_table=ax.table(cellText=vals, rowLabels=std_anom.index, colLabels=std_anom.columns, 
    header.set_fontsize(12.5)
    the_table=ax.table(cellText=vals, rowLabels=std_anom.index, colLabels=None, 
                        loc='center', cellColours=cmap(normal(vals_color)), cellLoc='center',
                        bbox=[0, 0, 1.0, 0.50]
                        )
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(12.5)
    # change font color to white where needed:
    table_props = the_table.properties()
    #table_cells = table_props['child_artists']
    last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
    for key, cell in the_table.get_celld().items():
        cell_text = cell.get_text().get_text()
        if is_number(cell_text) == False:
            pass
        #elif key[0] == 0:# <--- remove when no years
        #    pass
        elif key[1] in last_columns:
             cell._text.set_color('darkslategray')
        elif (float(cell_text) <= -1.5) | (float(cell_text) >= 1.5) :
            cell._text.set_color('white')
        elif (cell_text=='nan'):
            cell._set_facecolor('lightgray')
            cell._text.set_color('lightgray')

    plt.savefig("scorecards_fall_3K.png", dpi=300)
    os.system('convert -trim scorecards_fall_3K.png scorecards_fall_3K.png')

    # French table
    std_anom = std_anom.rename({r'$\rm T_{bot}$' : r'$\rm T_{fond}$', r'$\rm T_{bot_{<200m}}$' : r'$\rm T_{fond_{<200m}}$', r'$\rm Area_{>2^{\circ}C}$' : r'$\rm Aire_{>2^{\circ}C}$', r'$\rm Area_{<1^{\circ}C}$' : r'$\rm Aire_{<1^{\circ}C}$'})
    header = ax.table(cellText=[['']],
                          colLabels=['-- Division 3K de l\'OPANO --'],
                          loc='center'
                          )
    #the_table=ax.table(cellText=vals, rowLabels=std_anom.index, colLabels=std_anom.columns, 
    header.set_fontsize(12.5)
    the_table=ax.table(cellText=vals, rowLabels=std_anom.index, colLabels=None, 
                        loc='center', cellColours=cmap(normal(vals_color)), cellLoc='center',
                        bbox=[0, 0, 1.0, 0.50]
                        )
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(12.5)
    # change font color to white where needed:
    table_props = the_table.properties()
    #table_cells = table_props['child_artists']
    last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
    for key, cell in the_table.get_celld().items():
        cell_text = cell.get_text().get_text()
        if is_number(cell_text) == False:
            pass
        #elif key[0] == 0:# <--- remove when no years
        #    pass
        elif key[1] in last_columns:
             cell._text.set_color('darkslategray')
        elif (float(cell_text) <= -1.5) | (float(cell_text) >= 1.5) :
            cell._text.set_color('white')
        elif (cell_text=='nan'):
            cell._set_facecolor('lightgray')
            cell._text.set_color('lightgray')

    plt.savefig("scorecards_fall_3K_FR.png", dpi=300)
    os.system('convert -trim scorecards_fall_3K_FR.png scorecards_fall_3K_FR.png')

    # 3.
    infile = 'stats_3LNO_fall.pkl'
    df = pd.read_pickle(infile)
    df.index = pd.to_datetime(df.index) # update index to datetime
    df = df[(df.index.year>=years[0]) & (df.index.year<=years[-1])]
    # Flag bad years (no or weak sampling):
    bad_years = np.array([2021])
    for i in bad_years:
        df[df.index.year==i]=np.nan
    df['area_colder0'] = df['area_colder0']/1000 # In 1000km
    df['area_colder1'] = df['area_colder1']/1000 # In 1000km
    df['area_warmer2'] = df['area_warmer2']/1000
    df_clim = df[(df.index.year>=clim_year[0]) & (df.index.year<=clim_year[1])]
    std_anom = (df-df_clim.mean(axis=0))/df_clim.std(axis=0)
    # Save in .csv for future use
    std_anom.to_csv('bottomT_stn_anom_3LNO_fall.csv', sep=',', float_format='%0.3f')
    # Transpose for table
    std_anom = std_anom.T
    std_anom['MEAN'] = df_clim.mean(axis=0)
    std_anom['SD'] = df_clim.std(axis=0)
    std_anom = std_anom.reindex(['Tmean', 'Tmean_sha200', 'area_warmer2', 'area_colder0'])
    std_anom = std_anom.rename({'Tmean': r'$\rm T_{bot}$', 'Tmean_sha200': r'$\rm T_{bot_{<200m}}$', 'area_warmer2': r'$\rm Area_{>2^{\circ}C}$', 'area_colder0': r'$\rm Area_{<0^{\circ}C}$'})
    std_anom.rename(columns={'MEAN': r'$\rm \overline{x}$', 'SD': r'sd'}, inplace=True)

    vals = np.around(std_anom.values,1)
    vals[vals==-0.] = 0.
    vals_color = vals.copy()
    vals_color[-1,] = vals_color[-1,]*-1
    vals_color[:,-1] = 0 # No color to last two columns (mean and STD)
    vals_color[:,-2] = 0
    #normal = plt.Normalize(-4.49, 4.49)
    #cmap = plt.cm.get_cmap('seismic', 9) 
    fig=plt.figure(figsize=(ncols*wcell+wpad, nrows*hcell+hpad))
    ax = fig.add_subplot(111)
    ax.axis('off')
    #do the table
    header = ax.table(cellText=[['']],
                          colLabels=['-- NAFO division 3LNO --'],
                          loc='center'
                          )

    header.set_fontsize(12.5)

    the_table=ax.table(cellText=vals, rowLabels=std_anom.index, colLabels=None, 
                        loc='center', cellColours=cmap(normal(vals_color)), cellLoc='center',
                        bbox=[0, 0, 1.0, 0.50]
                        )
    # change font color to white where needed:
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(12.5)
    table_props = the_table.properties()
    #table_cells = table_props['child_artists']
    last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
    for key, cell in the_table.get_celld().items():
        cell_text = cell.get_text().get_text()
        if is_number(cell_text) == False:
            pass
        #elif key[0] == 0:# <--- remove when no years
        #    pass
        elif key[1] in last_columns:
             cell._text.set_color('darkslategray')
        elif (float(cell_text) <= -1.5) | (float(cell_text) >= 1.5) :
            cell._text.set_color('white')
        elif (cell_text=='nan'):
            cell._set_facecolor('lightgray')
            cell._text.set_color('lightgray')

    plt.savefig("scorecards_fall_3LNO.png", dpi=300)
    os.system('convert -trim scorecards_fall_3LNO.png scorecards_fall_3LNO.png')

    # French table
    std_anom = std_anom.rename({r'$\rm T_{bot}$' : r'$\rm T_{fond}$', r'$\rm T_{bot_{<200m}}$' : r'$\rm T_{fond_{<200m}}$', r'$\rm Area_{>2^{\circ}C}$' : r'$\rm Aire_{>2^{\circ}C}$', r'$\rm Area_{<1^{\circ}C}$' : r'$\rm Aire_{<1^{\circ}C}$'})
    header = ax.table(cellText=[['']],
                          colLabels=['-- Divisions 3LNO de l\'OPANO --'],
                          loc='center'
                          )

    header.set_fontsize(12.5)

    the_table=ax.table(cellText=vals, rowLabels=std_anom.index, colLabels=None, 
                        loc='center', cellColours=cmap(normal(vals_color)), cellLoc='center',
                        bbox=[0, 0, 1.0, 0.50]
                        )
    # change font color to white where needed:
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(12.5)
    table_props = the_table.properties()
    #table_cells = table_props['child_artists']
    last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
    for key, cell in the_table.get_celld().items():
        cell_text = cell.get_text().get_text()
        if is_number(cell_text) == False:
            pass
        #elif key[0] == 0:# <--- remove when no years
        #    pass
        elif key[1] in last_columns:
             cell._text.set_color('darkslategray')
        elif (float(cell_text) <= -1.5) | (float(cell_text) >= 1.5) :
            cell._text.set_color('white')
        elif (cell_text=='nan'):
            cell._set_facecolor('lightgray')
            cell._text.set_color('lightgray')

    plt.savefig("scorecards_fall_3LNO_FR.png", dpi=300)
    os.system('convert -trim scorecards_fall_3LNO_FR.png scorecards_fall_3LNO_FR.png')

    plt.close('all')
    # English
    os.system('montage  scorecards_fall_2H.png scorecards_fall_2J.png scorecards_fall_3K.png scorecards_fall_3LNO.png -tile 1x4 -geometry +1+1  -background white  scorecards_botT_fall.png') 
    # French
    os.system('montage  scorecards_fall_2H_FR.png scorecards_fall_2J.png scorecards_fall_3K_FR.png scorecards_fall_3LNO_FR.png -tile 1x4 -geometry +1+1  -background white  scorecards_botT_fall_FR.png') 




    #### ------------- For Spring ---------------- ####
    # 1.
    infile = 'stats_3LNO_spring.pkl'
    df = pd.read_pickle(infile)
    df.index = pd.to_datetime(df.index) # update index to datetime
    df = df[(df.index.year>=years[0]) & (df.index.year<=years[-1])]
    # Flag bad years (no or weak sampling):
    bad_years = np.array([2020, 2021])
    for i in bad_years:
        df[df.index.year==i]=np.nan
    year_list = df.index.year.astype('str')
    year_list = [i[2:4] for i in year_list] # 2-digit year
    df.index = pd.to_datetime(df.index) # update index to datetime
    df['area_colder0'] = df['area_colder0']/1000 # In 1000km
    df['area_colder1'] = df['area_colder1']/1000 # In 1000km
    df['area_warmer2'] = df['area_warmer2']/1000
    df_clim = df[(df.index.year>=clim_year[0]) & (df.index.year<=clim_year[1])]
    std_anom = (df-df_clim.mean(axis=0))/df_clim.std(axis=0)
    std_anom = (df-df_clim.mean(axis=0))/df_clim.std(axis=0)
    # Save in .csv for future use
    std_anom.to_csv('bottomT_stn_anom_3LNO_spring.csv', sep=',', float_format='%0.3f')
    # Transpose for table
    std_anom = std_anom.T
    std_anom['MEAN'] = df_clim.mean(axis=0)
    std_anom['SD'] = df_clim.std(axis=0)
    std_anom = std_anom.reindex(['Tmean', 'Tmean_sha200', 'area_warmer2', 'area_colder0'])
    std_anom = std_anom.rename({'Tmean': r'$\rm T_{bot}$', 'Tmean_sha200': r'$\rm T_{bot_{<200m}}$', 'area_warmer2': r'$\rm Area_{>2^{\circ}C}$', 'area_colder0': r'$\rm Area_{<0^{\circ}C}$'})
    std_anom.rename(columns={'MEAN': r'$\rm \overline{x}$', 'SD': r'sd'}, inplace=True)

    year_list.append(r'$\rm \overline{x}$') # add 2 extra columns
    year_list.append(r'sd')
    vals = np.around(std_anom.values,1)
    vals[vals==-0.] = 0.
    vals_color = vals.copy()
    vals_color[-1,] = vals_color[-1,]*-1
    vals_color[:,-1] = 0 # No color to last two columns (mean and STD)
    vals_color[:,-2] = 0
    #normal = plt.Normalize(-4.49, 4.49)
    #cmap = plt.cm.get_cmap('seismic', 9) 
    nrows, ncols = std_anom.index.size, std_anom.columns.size
    fig=plt.figure(figsize=(ncols*wcell+wpad, nrows*hcell+hpad))
    ax = fig.add_subplot(111)
    ax.axis('off')
    #do the table
    header = ax.table(cellText=[['']],
                          colLabels=['-- NAFO division 3LNO --'],
                          loc='center'
                          )

    header.set_fontsize(12.5)
    #the_table=ax.table(cellText=vals, rowLabels=std_anom.index, colLabels=std_anom.columns, 
    the_table=ax.table(cellText=vals, rowLabels=std_anom.index, colLabels=year_list, 
                        loc='center', cellColours=cmap(normal(vals_color)), cellLoc='center',
                        bbox=[0, 0, 1.0, 0.50]
                        )
    # change font color to white where needed:
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(12.5)
    table_props = the_table.properties()
    #table_cells = table_props['child_artists']
    last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
    for key, cell in the_table.get_celld().items():
        cell_text = cell.get_text().get_text()
        if is_number(cell_text) == False:
            pass
        elif key[0] == 0:
            pass
        elif key[1] in last_columns:
             cell._text.set_color('darkslategray')
        elif (float(cell_text) <= -1.5) | (float(cell_text) >= 1.5) :
            cell._text.set_color('white')
        elif (cell_text=='nan'):
            cell._set_facecolor('lightgray')
            cell._text.set_color('lightgray')

    plt.savefig("scorecards_spring_3LNO.png", dpi=300)
    os.system('convert -trim scorecards_spring_3LNO.png scorecards_spring_3LNO.png')

    # French table
    std_anom = std_anom.rename({r'$\rm T_{bot}$' : r'$\rm T_{fond}$', r'$\rm T_{bot_{<200m}}$' : r'$\rm T_{fond_{<200m}}$', r'$\rm Area_{>2^{\circ}C}$' : r'$\rm Aire_{>2^{\circ}C}$', r'$\rm Area_{<1^{\circ}C}$' : r'$\rm Aire_{<1^{\circ}C}$'})
    year_list[-1] = u'ET'

    header = ax.table(cellText=[['']],
                          colLabels=['-- Divisions 3LNO de l\'OPANO --'],
                          loc='center'
                          )
    header.set_fontsize(13)
    #the_table=ax.table(cellText=vals, rowLabels=std_anom.index, colLabels=std_anom.columns, 
    the_table=ax.table(cellText=vals, rowLabels=std_anom.index, colLabels=year_list,
                        loc='center', cellColours=cmap(normal(vals_color)), cellLoc='center',
                        bbox=[0, 0, 1, 0.5]
                        )
    # change font color to white where needed:
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(12.5)
    table_props = the_table.properties()
    #table_cells = table_props['child_artists']
    last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
    for key, cell in the_table.get_celld().items():
        cell_text = cell.get_text().get_text()
        if is_number(cell_text) == False:
            pass
        elif key[0] == 0: #year's row = no color
            pass
        elif key[1] in last_columns:
             cell._text.set_color('darkslategray')
        elif (float(cell_text) <= -1.5) | (float(cell_text) >= 1.5) :
            cell._text.set_color('white')
        elif (cell_text=='nan'):
            cell._set_facecolor('lightgray')
            cell._text.set_color('lightgray')

    plt.savefig("scorecards_spring_3LNO_FR.png", dpi=300)
    os.system('convert -trim scorecards_spring_3LNO_FR.png scorecards_spring_3LNO_FR.png')


    # 2.
    infile = 'stats_3Ps_spring.pkl'
    df = pd.read_pickle(infile)
    df.index = pd.to_datetime(df.index) # update index to datetime
    df = df[(df.index.year>=years[0]) & (df.index.year<=years[-1])]
    # Flag bad years (no or weak sampling):
    bad_years = np.array([1980, 1981, 1985, 1986, 1987, 1988, 1989, 1990, 1991, 1992, 2006, 2020])
    for i in bad_years:
        df[df.index.year==i]=np.nan
    df['area_colder0'] = df['area_colder0']/1000 # In 1000km
    df['area_colder1'] = df['area_colder1']/1000 # In 1000km
    df['area_warmer2'] = df['area_warmer2']/1000
    df_clim = df[(df.index.year>=clim_year[0]) & (df.index.year<=clim_year[1])]
    std_anom = (df-df_clim.mean(axis=0))/df_clim.std(axis=0)
    # Save in .csv for future use
    std_anom.to_csv('bottomT_stn_anom_3Ps_spring.csv', sep=',', float_format='%0.3f')
    # Transpose for table
    std_anom = std_anom.T
    std_anom['MEAN'] = df_clim.mean(axis=0)
    std_anom['SD'] = df_clim.std(axis=0)
    std_anom = std_anom.reindex(['Tmean', 'Tmean_sha200', 'area_warmer2', 'area_colder0'])
    std_anom = std_anom.rename({'Tmean': r'$\rm T_{bot}$', 'Tmean_sha200': r'$\rm T_{bot_{<200m}}$', 'area_warmer2': r'$\rm Area_{>2^{\circ}C}$', 'area_colder0': r'$\rm Area_{<0^{\circ}C}$'})
    std_anom.rename(columns={'MEAN': r'$\rm \overline{x}$', 'SD': r'sd'}, inplace=True)

    vals = np.around(std_anom.values,1)
    vals[vals==-0.] = 0.
    vals_color = vals.copy()
    vals_color[-1,] = vals_color[-1,]*-1
    vals_color[:,-1] = 0 # No color to last two columns (mean and STD)
    vals_color[:,-2] = 0 
    #normal = plt.Normalize(-4.49, 4.49)
    #cmap = plt.cm.get_cmap('seismic', 9) 
    nrows, ncols = std_anom.index.size, std_anom.columns.size
    fig=plt.figure(figsize=(ncols*wcell+wpad, nrows*hcell+hpad))
    ax = fig.add_subplot(111)
    ax.axis('off')
    #do the table
    header = ax.table(cellText=[['']],
                          colLabels=['-- NAFO division 3Ps --'],
                          loc='center'
                          )
    header.set_fontsize(12.5)
    #the_table=ax.table(cellText=vals, rowLabels=std_anom.index, colLabels=std_anom.columns, 
    the_table=ax.table(cellText=vals, rowLabels=std_anom.index, colLabels=None, 
                        loc='center', cellColours=cmap(normal(vals_color)), cellLoc='center',
                        bbox=[0, 0, 1.0, 0.50]
                        )
    # change font color to white where needed:
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(12.5)
    table_props = the_table.properties()
    #table_cells = table_props['child_artists']
    last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
    for key, cell in the_table.get_celld().items():
        cell_text = cell.get_text().get_text()
        if is_number(cell_text) == False:
            pass
        #elif key[0] == 0: # <--- remove when no years
        #    pass
        elif key[1] in last_columns:
             cell._text.set_color('darkslategray')
        elif (float(cell_text) <= -1.5) | (float(cell_text) >= 1.5) :
            cell._text.set_color('white')
        elif (cell_text=='nan'):
            cell._set_facecolor('lightgray')
            cell._text.set_color('lightgray')
            
    plt.savefig("scorecards_spring_3Ps.png", dpi=300)
    os.system('convert -trim scorecards_spring_3Ps.png scorecards_spring_3Ps.png')

    # French table
    std_anom = std_anom.rename({r'$\rm T_{bot}$' : r'$\rm T_{fond}$', r'$\rm T_{bot_{<200m}}$' : r'$\rm T_{fond_{<200m}}$', r'$\rm Area_{>2^{\circ}C}$' : r'$\rm Aire_{>2^{\circ}C}$', r'$\rm Area_{<1^{\circ}C}$' : r'$\rm Aire_{<1^{\circ}C}$'})
    header = ax.table(cellText=[['']],
                          colLabels=['-- Division 3Ps de l\'OPANO --'],
                          loc='center'
                          )

    header.set_fontsize(12.5)

    the_table=ax.table(cellText=vals, rowLabels=std_anom.index, colLabels=None, 
                        loc='center', cellColours=cmap(normal(vals_color)), cellLoc='center',
                        bbox=[0, 0, 1.0, 0.50]
                        )
    # change font color to white where needed:
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(12.5)
    table_props = the_table.properties()
    #table_cells = table_props['child_artists']
    last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
    for key, cell in the_table.get_celld().items():
        cell_text = cell.get_text().get_text()
        if is_number(cell_text) == False:
            pass
        #elif key[0] == 0:# <--- remove when no years
        #    pass
        elif key[1] in last_columns:
             cell._text.set_color('darkslategray')
        elif (float(cell_text) <= -1.5) | (float(cell_text) >= 1.5) :
            cell._text.set_color('white')
        elif (cell_text=='nan'):
            cell._set_facecolor('lightgray')
            cell._text.set_color('lightgray')
            
    plt.savefig("scorecards_spring_3Ps_FR.png", dpi=300)
    os.system('convert -trim scorecards_spring_3Ps_FR.png scorecards_spring_3Ps_FR.png')

    plt.close('all')
    # English montage
    os.system('montage  scorecards_spring_3LNO.png scorecards_spring_3Ps.png -tile 1x3 -geometry +1+1  -background white  scorecards_botT_spring.png')
    # French montage
    os.system('montage  scorecards_spring_3LNO_FR.png scorecards_spring_3Ps_FR.png -tile 1x3 -geometry +1+1  -background white  scorecards_botT_spring_FR.png')

    
def sfa_bottom_scorecards(years, clim_year=[2006, 2020]):


    '''
    Similar to bottom_scorcards, but for SAFs

    usage example (see azmp_genreport.py):
    azrt.sfa_bottom_scorecards(years=np.arange(2006, 2022), clim_year=[2006, 2020])

    '''

    #### ------------- For summer ---------------- ####
    # 0.
    infile = 'stats_sfa2_summer.pkl'
    df = pd.read_pickle(infile)
    df.index = pd.to_datetime(df.index) # update index to datetime
    df = df[(df.index.year>=years[0]) & (df.index.year<=years[-1])]
    # Flag bad years (no or weak sampling):
    bad_years = np.array([])
    for i in bad_years:
        df[df.index.year==i]=np.nan
    year_list = df.index.year.astype('str')
    year_list = [i[2:4] for i in year_list] # 2-digit year
    df['area_colder0'] = df['area_colder0']/1000 # In 1000km
    df['area_colder1'] = df['area_colder1']/1000 # In 1000km
    df['area_warmer2'] = df['area_warmer2']/1000
    df_clim = df[(df.index.year>=clim_year[0]) & (df.index.year<=clim_year[1])]
    std_anom = (df-df_clim.mean(axis=0))/df_clim.std(axis=0)
    std_anom = std_anom.T
    std_anom['MEAN'] = df_clim.mean(axis=0)
    std_anom['SD'] = df_clim.std(axis=0)
    std_anom = std_anom.reindex(['Tmean', 'Tmean_sha200', 'area_warmer2', 'area_colder1'])
    std_anom = std_anom.rename({'Tmean': r'$\rm T_{bot}$', 'Tmean_sha200': r'$\rm T_{bot_{<200m}}$', 'area_warmer2': r'$\rm Area_{>2^{\circ}C}$', 'area_colder1': r'$\rm Area_{<1^{\circ}C}$'})
    std_anom.rename(columns={'MEAN': r'$\rm \overline{x}$', 'SD': r'sd'}, inplace=True)
    # Save in .csv for future use
    std_anom.to_csv('bottomT_stn_anom_sfa2_summer.csv', sep=',', float_format='%0.3f')
        
    # Get text values +  cell color
    year_list.append(r'$\rm \overline{x}$') # add 2 extra columns
    year_list.append(r'sd')   
    vals = np.around(std_anom.values,1)
    vals[vals==-0.] = 0.
    vals_color = vals.copy()
    vals_color[-1,] = vals_color[-1,]*-1 # Reverse last row colorscale
    vals_color[:,-1] = 0 # No color to last two columns (mean and STD)
    vals_color[:,-2] = 0
    #vals_color[(vals_color<0.5) & (vals_color>-.5)] = 0.

    # Build the colormap
    vmin = -3.49
    vmax = 3.49
    midpoint = 0
    levels = np.linspace(vmin, vmax, 15)
    midp = np.mean(np.c_[levels[:-1], levels[1:]], axis=1)
    colvals = np.interp(midp, [vmin, midpoint, vmax], [-1, 0., 1])
    normal = plt.Normalize(-3.49, 3.49)
    reds = plt.cm.Reds(np.linspace(0,1, num=7))
    blues = plt.cm.Blues_r(np.linspace(0,1, num=7))
    whites = [(1,1,1,1)]*2
    colors = np.vstack((blues[0:-1,:], whites, reds[1:,:]))
    colors = np.concatenate([[colors[0,:]], colors, [colors[-1,:]]], 0)
    cmap, norm = from_levels_and_colors(levels, colors, extend='both')
    cmap_r, norm_r = from_levels_and_colors(levels, np.flipud(colors), extend='both')
    # Common parameters
    #hcell, wcell = 0.5, 0.6
    #hpad, wpad = 0, 0

    ## normal = plt.Normalize(-4.49, 4.49)
    ## cmap = plt.cm.get_cmap('seismic', 9) 
    #cmap = plt.cm.get_cmap('seismic', 15) 

    nrows, ncols = std_anom.index.size+1, std_anom.columns.size
    hcell, wcell = 0.5, 0.5
    hpad, wpad = 1, 1    
    fig=plt.figure(figsize=(ncols*wcell+wpad, nrows*hcell+hpad))
    ax = fig.add_subplot(111)
    ax.axis('off')
    #do the table
    header = ax.table(cellText=[['']],
                          colLabels=['-- SFA2 / EAZ --'],
                          loc='center'
                          )
    header.set_fontsize(13)
    #the_table=ax.table(cellText=vals, rowLabels=std_anom.index, colLabels=std_anom.columns, 
    the_table=ax.table(cellText=vals, rowLabels=std_anom.index, colLabels=year_list,
                        loc='center', cellColours=cmap(normal(vals_color)), cellLoc='center',
                        bbox=[0, 0, 1, 0.5]
                        )
    # change font color to white where needed:
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(12.5)
    table_props = the_table.properties()
    last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
    for key, cell in the_table.get_celld().items():
        cell_text = cell.get_text().get_text()
        if is_number(cell_text) == False:
            pass
        elif key[0] == 0: #year's row = no color
            pass
        elif key[1] in last_columns:
             cell._text.set_color('darkslategray')
        elif (np.float(cell_text) <= -1.5) | (np.float(cell_text) >= 1.5) :
            cell._text.set_color('white')
        elif (cell_text=='nan'):
            cell._set_facecolor('lightgray')
            cell._text.set_color('lightgray')

    plt.savefig("scorecards_summer_SFA2.png", dpi=300)
    os.system('convert -trim scorecards_summer_SFA2.png scorecards_summer_SFA2.png')

    # French table
    std_anom = std_anom.rename({r'$\rm T_{bot}$' : r'$\rm T_{fond}$', r'$\rm T_{bot_{<200m}}$' : r'$\rm T_{fond_{<200m}}$', r'$\rm Area_{>2^{\circ}C}$' : r'$\rm Aire_{>2^{\circ}C}~$ ', r'$\rm Area_{<1^{\circ}C}$' : r'$\rm Aire_{<1^{\circ}C}~$ '})
    year_list[-1] = u'ET'

    header = ax.table(cellText=[['']],
                          colLabels=['-- SFA2 / EAZ --'],
                          loc='center'
                          )
    header.set_fontsize(13)
    the_table=ax.table(cellText=vals, rowLabels=std_anom.index, colLabels=year_list,
                        loc='center', cellColours=cmap(normal(vals_color)), cellLoc='center',
                        bbox=[0, 0, 1, 0.5]
                        )
    # change font color to white where needed:
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(12.5)
    table_props = the_table.properties()
    last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
    for key, cell in the_table.get_celld().items():
        cell_text = cell.get_text().get_text()
        if is_number(cell_text) == False:
            pass
        elif key[0] == 0: #year's row = no color
            pass
        elif key[1] in last_columns:
             cell._text.set_color('darkslategray')
        elif (np.float(cell_text) <= -1.5) | (np.float(cell_text) >= 1.5) :
            cell._text.set_color('white')
        elif (cell_text=='nan'):
            cell._set_facecolor('lightgray')
            cell._text.set_color('lightgray')

    plt.savefig("scorecards_summer_SFA2_FR.png", dpi=300)
    os.system('convert -trim scorecards_summer_SFA2_FR.png scorecards_summer_SFA2_FR.png')

 # 1.
    infile = 'stats_sfa3_summer.pkl'
    df = pd.read_pickle(infile)
    df.index = pd.to_datetime(df.index) # update index to datetime
    df = df[(df.index.year>=years[0]) & (df.index.year<=years[-1])]
    # Flag bad years (no or weak sampling):
    bad_years = np.array([2006, 2008, 2010, 2012])
    for i in bad_years:
        df[df.index.year==i]=np.nan
    df['area_colder0'] = df['area_colder0']/1000 # In 1000km
    df['area_colder1'] = df['area_colder1']/1000 # In 1000km
    df['area_warmer2'] = df['area_warmer2']/1000
    df_clim = df[(df.index.year>=clim_year[0]) & (df.index.year<=clim_year[1])]
    std_anom = (df-df_clim.mean(axis=0))/df_clim.std(axis=0)
    std_anom = std_anom.T
    std_anom['MEAN'] = df_clim.mean(axis=0)
    std_anom['SD'] = df_clim.std(axis=0)
    std_anom = std_anom.reindex(['Tmean', 'Tmean_sha200', 'area_warmer2', 'area_colder1'])
    std_anom = std_anom.rename({'Tmean': r'$\rm T_{bot}$', 'Tmean_sha200': r'$\rm T_{bot_{<200m}}$', 'area_warmer2': r'$\rm Area_{>2^{\circ}C}$', 'area_colder1': r'$\rm Area_{<1^{\circ}C}$'})
    std_anom.rename(columns={'MEAN': r'$\rm \overline{x}$', 'SD': r'sd'}, inplace=True)
    # Save in .csv for future use
    std_anom.to_csv('bottomT_stn_anom_sfa3_summer.csv', sep=',', float_format='%0.3f')

    vals = np.around(std_anom.values,1)
    vals[vals==-0.] = 0.
    vals_color = vals.copy()
    vals_color[-1,] = vals_color[-1,]*-1
    vals_color[:,-1] = 0 # No color to last two columns (mean and STD)
    vals_color[:,-2] = 0 
    #normal = plt.Normalize(-4.49, 4.49)
    #cmap = plt.cm.get_cmap('seismic', 9) 
    nrows, ncols = std_anom.index.size, std_anom.columns.size
    fig=plt.figure(figsize=(ncols*wcell+wpad, nrows*hcell+hpad))
    ax = fig.add_subplot(111)
    ax.axis('off')
    #do the table
    header = ax.table(cellText=[['']],
                          colLabels=['-- SFA3 / WAZ --'],
                          loc='center'
                          )
    header.set_fontsize(12.5)
    the_table=ax.table(cellText=vals, rowLabels=std_anom.index, colLabels=None, 
                        loc='center', cellColours=cmap(normal(vals_color)), cellLoc='center',
                        bbox=[0, 0, 1.0, 0.50]
                        )
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(12.5)
    # change font color to white where needed:
    table_props = the_table.properties()
    #table_cells = table_props['child_artists']
    last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
    for key, cell in the_table.get_celld().items():
        cell_text = cell.get_text().get_text()
        if is_number(cell_text) == False:
            pass
        #elif key[0] == 0:# <--- remove when no years
        #    pass
        elif key[1] in last_columns:
             cell._text.set_color('darkslategray')
        elif (np.float(cell_text) <= -1.5) | (np.float(cell_text) >= 1.5) :
            cell._text.set_color('white')
        elif (cell_text=='nan'):
            cell._set_facecolor('lightgray')
            cell._text.set_color('lightgray')

    plt.savefig("scorecards_summer_SFA3.png", dpi=300)
    os.system('convert -trim scorecards_summer_SFA3.png scorecards_summer_SFA3.png')

    # French table
    std_anom = std_anom.rename({r'$\rm T_{bot}$' : r'$\rm T_{fond}$', r'$\rm T_{bot_{<200m}}$' : r'$\rm T_{fond_{<200m}}$', r'$\rm Area_{>2^{\circ}C}$' : r'$\rm Aire_{>2^{\circ}C}~$ ', r'$\rm Area_{<1^{\circ}C}$' : r'$\rm Aire_{<1^{\circ}C}~$ '})
    header = ax.table(cellText=[['']],
                          colLabels=['-- SFA3 / WAZ --'],
                          loc='center'
                          )
    #the_table=ax.table(cellText=vals, rowLabels=std_anom.index, colLabels=std_anom.columns, 
    header.set_fontsize(12.5)
    the_table=ax.table(cellText=vals, rowLabels=std_anom.index, colLabels=None, 
                        loc='center', cellColours=cmap(normal(vals_color)), cellLoc='center',
                        bbox=[0, 0, 1.0, 0.50]
                        )
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(12.5)
    # change font color to white where needed:
    table_props = the_table.properties()
    last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
    for key, cell in the_table.get_celld().items():
        cell_text = cell.get_text().get_text()
        if is_number(cell_text) == False:
            pass
        #elif key[0] == 0:# <--- remove when no years
        #    pass
        elif key[1] in last_columns:
             cell._text.set_color('darkslategray')
        elif (np.float(cell_text) <= -1.5) | (np.float(cell_text) >= 1.5) :
            cell._text.set_color('white')
        elif (cell_text=='nan'):
            cell._set_facecolor('lightgray')
            cell._text.set_color('lightgray')

    plt.savefig("scorecards_summer_SFA3_FR.png", dpi=300)
    os.system('convert -trim scorecards_summer_SFA3_FR.png scorecards_summer_SFA3_FR.png')


    # 2.
    infile = 'stats_sfa4_summer.pkl'
    df = pd.read_pickle(infile)
    df.index = pd.to_datetime(df.index) # update index to datetime
    df = df[(df.index.year>=years[0]) & (df.index.year<=years[-1])]
    # Flag bad years (no or weak sampling):
    bad_years = np.array([])
    for i in bad_years:
        df[df.index.year==i]=np.nan
    df['area_colder0'] = df['area_colder0']/1000 # In 1000km
    df['area_colder1'] = df['area_colder1']/1000 # In 1000km
    df['area_warmer2'] = df['area_warmer2']/1000
    df_clim = df[(df.index.year>=clim_year[0]) & (df.index.year<=clim_year[1])]
    std_anom = (df-df_clim.mean(axis=0))/df_clim.std(axis=0)
    std_anom = std_anom.T
    std_anom['MEAN'] = df_clim.mean(axis=0)
    std_anom['SD'] = df_clim.std(axis=0)
    std_anom = std_anom.reindex(['Tmean', 'Tmean_sha200', 'area_warmer2', 'area_colder1'])
    std_anom = std_anom.rename({'Tmean': r'$\rm T_{bot}$', 'Tmean_sha200': r'$\rm T_{bot_{<200m}}$', 'area_warmer2': r'$\rm Area_{>2^{\circ}C}$', 'area_colder1': r'$\rm Area_{<1^{\circ}C}$'})
    std_anom.rename(columns={'MEAN': r'$\rm \overline{x}$', 'SD': r'sd'}, inplace=True)
    # Save in .csv for future use
    std_anom.to_csv('bottomT_stn_anom_sfa4_summer.csv', sep=',', float_format='%0.3f')
        
    vals = np.around(std_anom.values,1)
    vals[vals==-0.] = 0.
    vals_color = vals.copy()
    vals_color[-1,] = vals_color[-1,]*-1
    vals_color[:,-1] = 0 # No color to last two columns (mean and STD)
    vals_color[:,-2] = 0 
    #normal = plt.Normalize(-4.49, 4.49)
    #cmap = plt.cm.get_cmap('seismic', 9) 
    nrows, ncols = std_anom.index.size, std_anom.columns.size
    fig=plt.figure(figsize=(ncols*wcell+wpad, nrows*hcell+hpad))
    ax = fig.add_subplot(111)
    ax.axis('off')
    #do the table
    header = ax.table(cellText=[['']],
                          colLabels=['-- SFA4 --'],
                          loc='center'
                          )
    header.set_fontsize(12.5)
    the_table=ax.table(cellText=vals, rowLabels=std_anom.index, colLabels=None, 
                        loc='center', cellColours=cmap(normal(vals_color)), cellLoc='center',
                        bbox=[0, 0, 1.0, 0.50]
                        )
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(12.5)
    # change font color to white where needed:
    table_props = the_table.properties()
    #table_cells = table_props['child_artists']
    last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
    for key, cell in the_table.get_celld().items():
        cell_text = cell.get_text().get_text()
        if is_number(cell_text) == False:
            pass
        #elif key[0] == 0:# <--- remove when no years
        #    pass
        elif key[1] in last_columns:
             cell._text.set_color('darkslategray')
        elif (np.float(cell_text) <= -1.5) | (np.float(cell_text) >= 1.5) :
            cell._text.set_color('white')
        elif (cell_text=='nan'):
            cell._set_facecolor('lightgray')
            cell._text.set_color('lightgray')

    plt.savefig("scorecards_summer_SFA4.png", dpi=300)
    os.system('convert -trim scorecards_summer_SFA4.png scorecards_summer_SFA4.png')

    # French table
    std_anom = std_anom.rename({r'$\rm T_{bot}$' : r'$\rm T_{fond}$', r'$\rm T_{bot_{<200m}}$' : r'$\rm T_{fond_{<200m}}$', r'$\rm Area_{>2^{\circ}C}$' : r'$\rm Aire_{>2^{\circ}C}~$ ', r'$\rm Area_{<1^{\circ}C}$' : r'$\rm Aire_{<1^{\circ}C}~$ '})
    header = ax.table(cellText=[['']],
                          colLabels=['-- SFA4 --'],
                          loc='center'
                          )
    #the_table=ax.table(cellText=vals, rowLabels=std_anom.index, colLabels=std_anom.columns, 
    header.set_fontsize(12.5)
    the_table=ax.table(cellText=vals, rowLabels=std_anom.index, colLabels=None, 
                        loc='center', cellColours=cmap(normal(vals_color)), cellLoc='center',
                        bbox=[0, 0, 1.0, 0.50]
                        )
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(12.5)
    # change font color to white where needed:
    table_props = the_table.properties()
    #table_cells = table_props['child_artists']
    last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
    for key, cell in the_table.get_celld().items():
        cell_text = cell.get_text().get_text()
        if is_number(cell_text) == False:
            pass
        #elif key[0] == 0:# <--- remove when no years
        #    pass
        elif key[1] in last_columns:
             cell._text.set_color('darkslategray')
        elif (np.float(cell_text) <= -1.5) | (np.float(cell_text) >= 1.5) :
            cell._text.set_color('white')
        elif (cell_text=='nan'):
            cell._set_facecolor('lightgray')
            cell._text.set_color('lightgray')

    plt.savefig("scorecards_summer_SFA4_FR.png", dpi=300)
    os.system('convert -trim scorecards_summer_SFA4_FR.png scorecards_summer_SFA4_FR.png')
    plt.close('all')

    
    ## **SFA4 in stand-alone (with years row)
    # rename back to English
    std_anom = std_anom.rename({r'$\rm T_{fond}$' : r'$\rm T_{bot}$', r'$\rm T_{fond_{<200m}}$' : r'$\rm T_{bot_{<200m}}$', r'$\rm Aire_{>2^{\circ}C}~$ ' : r'$\rm Area_{>2^{\circ}C}$', r'$\rm Aire_{<1^{\circ}C}~$ ' : r'$\rm Area_{<1^{\circ}C}$'})
    year_list[-1] = u'SD'
    #do the table
    fig=plt.figure(figsize=(ncols*wcell+wpad, nrows*hcell+hpad))
    ax = fig.add_subplot(111)
    ax.axis('off')
    header = ax.table(cellText=[['']],
                          colLabels=['-- SFA4 --'],
                          loc='center'
                          )
    header.set_fontsize(13)
    the_table=ax.table(cellText=vals, rowLabels=std_anom.index, colLabels=year_list, 
                        loc='center', cellColours=cmap(normal(vals_color)), cellLoc='center',
                        bbox=[0, 0, 1.0, 0.5]
                        )
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(12.5)
    # change font color to white where needed:
    table_props = the_table.properties()
    last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
    for key, cell in the_table.get_celld().items():
        cell_text = cell.get_text().get_text()
        if is_number(cell_text) == False:
            pass
        elif key[0] == 0: #year's row = no color
            pass
        elif key[1] in last_columns:
             cell._text.set_color('darkslategray')
        elif (np.float(cell_text) <= -1.5) | (np.float(cell_text) >= 1.5) :
            cell._text.set_color('white')
        elif (cell_text=='nan'):
            cell._set_facecolor('lightgray')
            cell._text.set_color('lightgray')

    plt.savefig("scorecards_botT_SFA4.png", dpi=300)
    os.system('convert -trim scorecards_botT_SFA4.png scorecards_botT_SFA4.png')

    # French table (SFA4 with years)
    #std_anom = std_anom.rename({r'$\rm T_{bot}$' : r'$\rm T_{fond}$', r'$\rm T_{bot_{<200m}}$' : r'$\rm T_{fond_{<200m}}$', r'$\rm Area_{>2^{\circ}C}$' : r'$\rm Aire_{>2^{\circ}C}$', r'$\rm Area_{<1^{\circ}C}$' : r'$\rm Aire_{<1^{\circ}C}$'})
    std_anom = std_anom.rename({r'$\rm T_{bot}$' : r'$\rm T_{fond}$', r'$\rm T_{bot_{<200m}}$' : r'$\rm T_{fond_{<200m}}$', r'$\rm Area_{>2^{\circ}C}$' : r'$\rm Aire_{>2^{\circ}C}~$ ', r'$\rm Area_{<1^{\circ}C}$' : r'$\rm Aire_{<1^{\circ}C}~$ '})
    year_list[-1] = u'ET'
    header = ax.table(cellText=[['']],
                          colLabels=['-- SFA4 --'],
                          loc='center'
                          )
    header.set_fontsize(12.5)
    the_table=ax.table(cellText=vals, rowLabels=std_anom.index, colLabels=year_list, 
                        loc='center', cellColours=cmap(normal(vals_color)), cellLoc='center',
                        bbox=[0, 0, 1.0, 0.50]
                        )
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(12.5)
    # change font color to white where needed:
    table_props = the_table.properties()
    last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
    for key, cell in the_table.get_celld().items():
        cell_text = cell.get_text().get_text()
        if is_number(cell_text) == False:
            pass
        elif key[0] == 0: #year's row = no color
            pass
        elif key[1] in last_columns:
             cell._text.set_color('darkslategray')
        elif (np.float(cell_text) <= -1.5) | (np.float(cell_text) >= 1.5) :
            cell._text.set_color('white')
        elif (cell_text=='nan'):
            cell._set_facecolor('lightgray')
            cell._text.set_color('lightgray')

    plt.savefig("scorecards_botT_SFA4_FR.png", dpi=300)
    os.system('convert -trim scorecards_botT_SFA4_FR.png scorecards_botT_SFA4_FR.png')
    plt.close('all')

    
    ## Montage all SFAs in subplot - English
    os.system('montage  scorecards_summer_SFA2.png scorecards_summer_SFA3.png scorecards_summer_SFA4.png -tile 1x3 -geometry +1+1  -background white  scorecards_botT_SFA2-4_summer.png') 
    # French
    os.system('montage  scorecards_summer_SFA2_FR.png scorecards_summer_SFA3_FR.png scorecards_summer_SFA4_FR.png -tile 1x3 -geometry +1+1  -background white  scorecards_botT_SFA2-4_summer_FR.png')
    # Same, but EAZ/WAZ only
    os.system('montage  scorecards_summer_SFA2.png scorecards_summer_SFA3.png -tile 1x2 -geometry +1+1  -background white  scorecards_botT_SFA2-3_summer.png')
    os.system('montage  scorecards_summer_SFA2_FR.png scorecards_summer_SFA3_FR.png -tile 1x2 -geometry +1+1  -background white  scorecards_botT_SFA2-3_summer_FR.png')    
    # remove leftowvers
    os.system('rm scorecards_summer_*.png')

def sfa_bottomS_scorecards(years, clim_year=[2006, 2021]):


    '''
    Similar to sfa_bottom_scorecards, but for salinity in SFAs

    usage example (see azmp_genreport.py):
    azrt.sfa_bottomS_scorecards(years=np.arange(2006, 2022), clim_year=[2006, 2020])

    ** Note that varuables inside .pkl files are name "Tmean" isntead of "Smean"
    '''

    #### ------------- For summer ---------------- ####
    # 0.
    infile = 'stats_sfa2_summer_salinity.pkl'
    df = pd.read_pickle(infile)
    df.index = pd.to_datetime(df.index) # update index to datetime
    df = df[(df.index.year>=years[0]) & (df.index.year<=years[-1])]
    # Flag bad years (no or weak sampling):
    bad_years = np.array([])
    for i in bad_years:
        df[df.index.year==i]=np.nan
    year_list = df.index.year.astype('str')
    year_list = [i[2:4] for i in year_list] # 2-digit year
    df_clim = df[(df.index.year>=clim_year[0]) & (df.index.year<=clim_year[1])]
    std_anom = (df-df_clim.mean(axis=0))/df_clim.std(axis=0)
    std_anom = std_anom.T
    std_anom['MEAN'] = df_clim.mean(axis=0)
    std_anom['SD'] = df_clim.std(axis=0)
    std_anom = std_anom.reindex(['Tmean', 'Tmean_sha200'])
    std_anom = std_anom.rename({'Tmean': r'$\rm S_{bot}$', 'Tmean_sha200': r'$\rm S_{bot_{<200m}}$'})
    std_anom.rename(columns={'MEAN': r'$\rm \overline{x}$', 'SD': r'sd'}, inplace=True)
    # Save in .csv for future use
    std_anom.to_csv('bottomS_stn_anom_sfa2_summer.csv', sep=',', float_format='%0.3f')
        
    # Get text values +  cell color
    year_list.append(r'$\rm \overline{x}$') # add 2 extra columns
    year_list.append(r'sd')   
    vals = np.around(std_anom.values,1)
    vals[vals==-0.] = 0.
    vals_color = vals.copy()
    vals_color[:,-1] = 0 # No color to last two columns (mean and STD)
    vals_color[:,-2] = 0
    #vals_color[(vals_color<0.5) & (vals_color>-.5)] = 0.

    # Build the colormap
    vmin = -3.49
    vmax = 3.49
    midpoint = 0
    levels = np.linspace(vmin, vmax, 15)
    midp = np.mean(np.c_[levels[:-1], levels[1:]], axis=1)
    colvals = np.interp(midp, [vmin, midpoint, vmax], [-1, 0., 1])
    normal = plt.Normalize(-3.49, 3.49)
    reds = plt.cm.Reds(np.linspace(0,1, num=7))
    blues = plt.cm.Blues_r(np.linspace(0,1, num=7))
    whites = [(1,1,1,1)]*2
    colors = np.vstack((blues[0:-1,:], whites, reds[1:,:]))
    colors = np.concatenate([[colors[0,:]], colors, [colors[-1,:]]], 0)
    cmap, norm = from_levels_and_colors(levels, colors, extend='both')
    cmap_r, norm_r = from_levels_and_colors(levels, np.flipud(colors), extend='both')
    # Common parameters
    nrows, ncols = std_anom.index.size+1, std_anom.columns.size
    hcell, wcell = 0.5, 0.5
    hpad, wpad = 1, 1    
    fig=plt.figure(figsize=(ncols*wcell+wpad, nrows*hcell+hpad))
    ax = fig.add_subplot(111)
    ax.axis('off')
    #do the table
    header = ax.table(cellText=[['']],
                          colLabels=['-- SFA2 / EAZ --'],
                          loc='center'
                          )
    header.set_fontsize(13)
    #the_table=ax.table(cellText=vals, rowLabels=std_anom.index, colLabels=std_anom.columns, 
    the_table=ax.table(cellText=vals, rowLabels=std_anom.index, colLabels=year_list,
                        loc='center', cellColours=cmap(normal(vals_color)), cellLoc='center',
                        bbox=[0, 0, 1, 0.5]
                        )
    # change font color to white where needed:
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(12.5)
    table_props = the_table.properties()
    #table_cells = table_props['child_artists']
    last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
    for key, cell in the_table.get_celld().items():
        cell_text = cell.get_text().get_text()
        if is_number(cell_text) == False:
            pass
        elif key[0] == 0: #year's row = no color
            pass
        elif key[1] in last_columns:
             cell._text.set_color('darkslategray')
        elif (np.float(cell_text) <= -1.5) | (np.float(cell_text) >= 1.5) :
            cell._text.set_color('white')
        elif (cell_text=='nan'):
            cell._set_facecolor('lightgray')
            cell._text.set_color('lightgray')

    plt.savefig("scorecards_summerS_SFA2.png", dpi=300)
    os.system('convert -trim scorecards_summerS_SFA2.png scorecards_summerS_SFA2.png')

    # French table
    std_anom = std_anom.rename({r'$\rm S_{bot}$' : r'$\rm S_{fond}$', r'$\rm S_{bot_{<200m}}$' : r'$\rm S_{fond_{<200m}}$'})
    year_list[-1] = u'ET'

    header = ax.table(cellText=[['']],
                          colLabels=['-- SFA2 / EAZ --'],
                          loc='center'
                          )
    header.set_fontsize(13)
    #the_table=ax.table(cellText=vals, rowLabels=std_anom.index, colLabels=std_anom.columns, 
    the_table=ax.table(cellText=vals, rowLabels=std_anom.index, colLabels=year_list,
                        loc='center', cellColours=cmap(normal(vals_color)), cellLoc='center',
                        bbox=[0, 0, 1, 0.5]
                        )
    # change font color to white where needed:
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(12.5)
    table_props = the_table.properties()
    #table_cells = table_props['child_artists']
    last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
    for key, cell in the_table.get_celld().items():
        cell_text = cell.get_text().get_text()
        if is_number(cell_text) == False:
            pass
        elif key[0] == 0: #year's row = no color
            pass
        elif key[1] in last_columns:
             cell._text.set_color('darkslategray')
        elif (np.float(cell_text) <= -1.5) | (np.float(cell_text) >= 1.5) :
            cell._text.set_color('white')
        elif (cell_text=='nan'):
            cell._set_facecolor('lightgray')
            cell._text.set_color('lightgray')

    plt.savefig("scorecards_summerS_SFA2_FR.png", dpi=300)
    os.system('convert -trim scorecards_summerS_SFA2_FR.png scorecards_summerS_SFA2_FR.png')

 # 1.
    infile = 'stats_sfa3_summer_salinity.pkl'
    df = pd.read_pickle(infile)
    df.index = pd.to_datetime(df.index) # update index to datetime
    df = df[(df.index.year>=years[0]) & (df.index.year<=years[-1])]
    # Flag bad years (no or weak sampling):
    bad_years = np.array([2006, 2008, 2010, 2012])
    for i in bad_years:
        df[df.index.year==i]=np.nan
    df_clim = df[(df.index.year>=clim_year[0]) & (df.index.year<=clim_year[1])]
    std_anom = (df-df_clim.mean(axis=0))/df_clim.std(axis=0)
    std_anom = std_anom.T
    std_anom['MEAN'] = df_clim.mean(axis=0)
    std_anom['SD'] = df_clim.std(axis=0)
    std_anom = std_anom.reindex(['Tmean', 'Tmean_sha200'])
    std_anom = std_anom.rename({'Tmean': r'$\rm S_{bot}$', 'Tmean_sha200': r'$\rm S_{bot_{<200m}}$'})
    std_anom.rename(columns={'MEAN': r'$\rm \overline{x}$', 'SD': r'sd'}, inplace=True)
    # Save in .csv for future use
    std_anom.to_csv('bottomS_stn_anom_sfa3_summer.csv', sep=',', float_format='%0.3f')

    vals = np.around(std_anom.values,1)
    vals[vals==-0.] = 0.
    vals_color = vals.copy()
    vals_color[:,-1] = 0 # No color to last two columns (mean and STD)
    vals_color[:,-2] = 0 
    nrows, ncols = std_anom.index.size, std_anom.columns.size
    fig=plt.figure(figsize=(ncols*wcell+wpad, nrows*hcell+hpad))
    ax = fig.add_subplot(111)
    ax.axis('off')
    #do the table
    header = ax.table(cellText=[['']],
                          colLabels=['-- SFA3 / WAZ --'],
                          loc='center'
                          )
    header.set_fontsize(12.5)
    the_table=ax.table(cellText=vals, rowLabels=std_anom.index, colLabels=None, 
                        loc='center', cellColours=cmap(normal(vals_color)), cellLoc='center',
                        bbox=[0, 0, 1.0, 0.50]
                        )
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(12.5)
    # change font color to white where needed:
    table_props = the_table.properties()
    #table_cells = table_props['child_artists']
    last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
    for key, cell in the_table.get_celld().items():
        cell_text = cell.get_text().get_text()
        if is_number(cell_text) == False:
            pass
        #elif key[0] == 0:# <--- remove when no years
        #    pass
        elif key[1] in last_columns:
             cell._text.set_color('darkslategray')
        elif (np.float(cell_text) <= -1.5) | (np.float(cell_text) >= 1.5) :
            cell._text.set_color('white')
        elif (cell_text=='nan'):
            cell._set_facecolor('lightgray')
            cell._text.set_color('lightgray')

    plt.savefig("scorecards_summerS_SFA3.png", dpi=300)
    os.system('convert -trim scorecards_summerS_SFA3.png scorecards_summerS_SFA3.png')

    # French table
    std_anom = std_anom.rename({r'$\rm S_{bot}$' : r'$\rm S_{fond}$', r'$\rm S_{bot_{<200m}}$' : r'$\rm S_{fond_{<200m}}$'})
    header = ax.table(cellText=[['']],
                          colLabels=['-- SFA3 / WAZ --'],
                          loc='center'
                          )
    header.set_fontsize(12.5)
    the_table=ax.table(cellText=vals, rowLabels=std_anom.index, colLabels=None, 
                        loc='center', cellColours=cmap(normal(vals_color)), cellLoc='center',
                        bbox=[0, 0, 1.0, 0.50]
                        )
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(12.5)
    # change font color to white where needed:
    table_props = the_table.properties()
    #table_cells = table_props['child_artists']
    last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
    for key, cell in the_table.get_celld().items():
        cell_text = cell.get_text().get_text()
        if is_number(cell_text) == False:
            pass
        #elif key[0] == 0:# <--- remove when no years
        #    pass
        elif key[1] in last_columns:
             cell._text.set_color('darkslategray')
        elif (np.float(cell_text) <= -1.5) | (np.float(cell_text) >= 1.5) :
            cell._text.set_color('white')
        elif (cell_text=='nan'):
            cell._set_facecolor('lightgray')
            cell._text.set_color('lightgray')

    plt.savefig("scorecards_summerS_SFA3_FR.png", dpi=300)
    os.system('convert -trim scorecards_summerS_SFA3_FR.png scorecards_summerS_SFA3_FR.png')


    # 2.
    infile = 'stats_sfa4_summer_salinity.pkl'
    df = pd.read_pickle(infile)
    df.index = pd.to_datetime(df.index) # update index to datetime
    df = df[(df.index.year>=years[0]) & (df.index.year<=years[-1])]
    # Flag bad years (no or weak sampling):
    bad_years = np.array([])
    for i in bad_years:
        df[df.index.year==i]=np.nan
    df_clim = df[(df.index.year>=clim_year[0]) & (df.index.year<=clim_year[1])]
    std_anom = (df-df_clim.mean(axis=0))/df_clim.std(axis=0)
    std_anom = std_anom.T
    std_anom['MEAN'] = df_clim.mean(axis=0)
    std_anom['SD'] = df_clim.std(axis=0)
    std_anom = std_anom.reindex(['Tmean', 'Tmean_sha200'])
    std_anom = std_anom.rename({'Tmean': r'$\rm S_{bot}$', 'Tmean_sha200': r'$\rm S_{bot_{<200m}}$'})
    std_anom.rename(columns={'MEAN': r'$\rm \overline{x}$', 'SD': r'sd'}, inplace=True)
    # Save in .csv for future use
    std_anom.to_csv('bottomS_stn_anom_sfa4_summer.csv', sep=',', float_format='%0.3f')
        
    vals = np.around(std_anom.values,1)
    vals[vals==-0.] = 0.
    vals_color = vals.copy()
    vals_color[:,-1] = 0 # No color to last two columns (mean and STD)
    vals_color[:,-2] = 0 
    #normal = plt.Normalize(-4.49, 4.49)
    #cmap = plt.cm.get_cmap('seismic', 9) 
    nrows, ncols = std_anom.index.size, std_anom.columns.size
    fig=plt.figure(figsize=(ncols*wcell+wpad, nrows*hcell+hpad))
    ax = fig.add_subplot(111)
    ax.axis('off')
    #do the table
    header = ax.table(cellText=[['']],
                          colLabels=['-- SFA4 --'],
                          loc='center'
                          )
    header.set_fontsize(12.5)
    the_table=ax.table(cellText=vals, rowLabels=std_anom.index, colLabels=None, 
                        loc='center', cellColours=cmap(normal(vals_color)), cellLoc='center',
                        bbox=[0, 0, 1.0, 0.50]
                        )
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(12.5)
    # change font color to white where needed:
    table_props = the_table.properties()
    #table_cells = table_props['child_artists']
    last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
    for key, cell in the_table.get_celld().items():
        cell_text = cell.get_text().get_text()
        if is_number(cell_text) == False:
            pass
        #elif key[0] == 0:# <--- remove when no years
        #    pass
        elif key[1] in last_columns:
             cell._text.set_color('darkslategray')
        elif (np.float(cell_text) <= -1.5) | (np.float(cell_text) >= 1.5) :
            cell._text.set_color('white')
        elif (cell_text=='nan'):
            cell._set_facecolor('lightgray')
            cell._text.set_color('lightgray')

    plt.savefig("scorecards_summerS_SFA4.png", dpi=300)
    os.system('convert -trim scorecards_summerS_SFA4.png scorecards_summerS_SFA4.png')

    # French table
    std_anom = std_anom.rename({r'$\rm S_{bot}$' : r'$\rm S_{fond}$', r'$\rm S_{bot_{<200m}}$' : r'$\rm S_{fond_{<200m}}$'})
    header = ax.table(cellText=[['']],
                          colLabels=['-- SFA4 --'],
                          loc='center'
                          )
    header.set_fontsize(12.5)
    the_table=ax.table(cellText=vals, rowLabels=std_anom.index, colLabels=None, 
                        loc='center', cellColours=cmap(normal(vals_color)), cellLoc='center',
                        bbox=[0, 0, 1.0, 0.50]
                        )
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(12.5)
    # change font color to white where needed:
    table_props = the_table.properties()
    #table_cells = table_props['child_artists']
    last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
    for key, cell in the_table.get_celld().items():
        cell_text = cell.get_text().get_text()
        if is_number(cell_text) == False:
            pass
        #elif key[0] == 0:# <--- remove when no years
        #    pass
        elif key[1] in last_columns:
             cell._text.set_color('darkslategray')
        elif (np.float(cell_text) <= -1.5) | (np.float(cell_text) >= 1.5) :
            cell._text.set_color('white')
        elif (cell_text=='nan'):
            cell._set_facecolor('lightgray')
            cell._text.set_color('lightgray')

    plt.savefig("scorecards_summerS_SFA4_FR.png", dpi=300)
    os.system('convert -trim scorecards_summerS_SFA4_FR.png scorecards_summerS_SFA4_FR.png')
    plt.close('all')


    ## **SFA4 in stand-alone (with years row)
    # rename back to English
    std_anom = std_anom.rename({r'$\rm S_{fond}$' : r'$\rm S_{bot}$', r'$\rm S_{fond_{<200m}}$' : r'$\rm S_{bot_{<200m}}$'})
    year_list[-1] = u'SD'
    #do the table
    fig=plt.figure(figsize=(ncols*wcell+wpad, nrows*hcell+hpad))
    ax = fig.add_subplot(111)
    ax.axis('off')
    header = ax.table(cellText=[['']],
                          colLabels=['-- SFA4 --'],
                          loc='center'
                          )
    header.set_fontsize(13)
    the_table=ax.table(cellText=vals, rowLabels=std_anom.index, colLabels=year_list, 
                        loc='center', cellColours=cmap(normal(vals_color)), cellLoc='center',
                        bbox=[0, 0, 1.0, 0.5]
                        )
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(12.5)
    # change font color to white where needed:
    table_props = the_table.properties()
    last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
    for key, cell in the_table.get_celld().items():
        cell_text = cell.get_text().get_text()
        if is_number(cell_text) == False:
            pass
        elif key[0] == 0: #year's row = no color
            pass
        elif key[1] in last_columns:
             cell._text.set_color('darkslategray')
        elif (np.float(cell_text) <= -1.5) | (np.float(cell_text) >= 1.5) :
            cell._text.set_color('white')
        elif (cell_text=='nan'):
            cell._set_facecolor('lightgray')
            cell._text.set_color('lightgray')

    plt.savefig("scorecards_botS_SFA4.png", dpi=300)
    os.system('convert -trim scorecards_botS_SFA4.png scorecards_botS_SFA4.png')

    # French table (SFA4 with years)
    std_anom = std_anom.rename({r'$\rm S_{bot}$' : r'$\rm S_{fond}$', r'$\rm S_{bot_{<200m}}$' : r'$\rm S_{fond_{<200m}}$'})
    year_list[-1] = u'ET'
    header = ax.table(cellText=[['']],
                          colLabels=['-- SFA4 --'],
                          loc='center'
                          )
    header.set_fontsize(12.5)
    the_table=ax.table(cellText=vals, rowLabels=std_anom.index, colLabels=year_list, 
                        loc='center', cellColours=cmap(normal(vals_color)), cellLoc='center',
                        bbox=[0, 0, 1.0, 0.50]
                        )
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(12.5)
    # change font color to white where needed:
    table_props = the_table.properties()
    last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
    for key, cell in the_table.get_celld().items():
        cell_text = cell.get_text().get_text()
        if is_number(cell_text) == False:
            pass
        elif key[0] == 0: #year's row = no color
            pass
        elif key[1] in last_columns:
             cell._text.set_color('darkslategray')
        elif (np.float(cell_text) <= -1.5) | (np.float(cell_text) >= 1.5) :
            cell._text.set_color('white')
        elif (cell_text=='nan'):
            cell._set_facecolor('lightgray')
            cell._text.set_color('lightgray')

    plt.savefig("scorecards_botS_SFA4_FR.png", dpi=300)
    os.system('convert -trim scorecards_botS_SFA4_FR.png scorecards_botS_SFA4_FR.png')
    plt.close('all')

    
    ## Montage in subplot - English
    os.system('montage  scorecards_summerS_SFA2.png scorecards_summerS_SFA3.png scorecards_summerS_SFA4.png -tile 1x3 -geometry +1+1  -background white  scorecards_botS_SFA2-4_summer.png') 
    # French
    os.system('montage  scorecards_summerS_SFA2_FR.png scorecards_summerS_SFA3_FR.png scorecards_summerS_SFA4_FR.png -tile 1x3 -geometry +1+1  -background white  scorecards_botS_SFA2-4_summer_FR.png') 
    # Same, but EAZ/WAZ only
    os.system('montage  scorecards_summerS_SFA2.png scorecards_summerS_SFA3.png -tile 1x2 -geometry +1+1  -background white  scorecards_botS_SFA2-3_summer.png')
    os.system('montage  scorecards_summerS_SFA2_FR.png scorecards_summerS_SFA3_FR.png -tile 1x2 -geometry +1+1  -background white  scorecards_botS_SFA2-3_summer_FR.png')
    # remove leftovers
    os.system('rm scorecards_summerS_*.png')


def SS_bottom_scorecards(years, clim_year=[1991, 2020]):


    '''
    To generate AZMP score cards for bottom temperature
    on the Scotian shelf (NAFO Div. 4).
    
    Data are from D. Hebert

    usage example:

    bottom_scorecards(years=[1980, 2021], clim_year=[1991, 2020]):

    df_4v = pd.read_csv('/home/cyrf0006/data/Hebert_timeseries/summerGroundfishBottomTemperature_4V.dat', delimiter=r",", index_col='year', header=10)
    df_4vn = pd.read_csv('/home/cyrf0006/data/Hebert_timeseries/summerGroundfishBottomTemperature_4Vn.dat', delimiter=r",", index_col='year', header=10)
    df_4vs = pd.read_csv('/home/cyrf0006/data/Hebert_timeseries/summerGroundfishBottomTemperature_4Vs.dat', delimiter=r",", index_col='year', header=10)
    df_4w = pd.read_csv('/home/cyrf0006/data/Hebert_timeseries/summerGroundfishBottomTemperature_4W.dat', delimiter=r",", index_col='year', header=10)
    df_4x = pd.read_csv('/home/cyrf0006/data/Hebert_timeseries/summerGroundfishBottomTemperature_4X.dat', delimiter=r",", index_col='year', header=10)


    UNFINISHED
    '''

    #### ------------- For summer ---------------- ####

    # 1. 4Vn
    df_4vn = pd.read_csv('/home/cyrf0006/data/Hebert_timeseries/summerGroundfishBottomTemperature_4Vn.dat', delimiter=r",", index_col='year', header=10)
    df = pd.read_pickle(infile)
    df.index = pd.to_datetime(df.index) # update index to datetime
    df = df[(df.index.year>=years[0]) & (df.index.year<=years[-1])]
    # Flag bad years (no or weak sampling):
    bad_years = np.array([1980, 1982, 1984, 1985, 1986, 1987, 1988, 1989, 1990, 1992, 1993, 1994, 1995, 1996, 2000, 2002, 2003, 2005, 2007, 2009])
    for i in bad_years:
        df[df.index.year==i]=np.nan
    year_list = df.index.year.astype('str')
    year_list = [i[2:4] for i in year_list] # 2-digit year
    df['area_colder0'] = df['area_colder0']/1000 # In 1000km
    df['area_colder1'] = df['area_colder1']/1000 # In 1000km
    df['area_warmer2'] = df['area_warmer2']/1000
    df_clim = df[(df.index.year>=clim_year[0]) & (df.index.year<=clim_year[1])]
    std_anom = (df-df_clim.mean(axis=0))/df_clim.std(axis=0)
    # Save in .csv for future use
    std_anom.to_csv('bottomT_stn_anom_2H_fall.csv', sep=',', float_format='%0.3f')
    # Transpose for table
    std_anom = std_anom.T
    std_anom['MEAN'] = df_clim.mean(axis=0)
    std_anom['SD'] = df_clim.std(axis=0)
    std_anom = std_anom.reindex(['Tmean', 'Tmean_sha200', 'area_warmer2', 'area_colder1'])
    std_anom = std_anom.rename({'Tmean': r'$\rm T_{bot}$', 'Tmean_sha200': r'$\rm T_{bot_{<200m}}$', 'area_warmer2': r'$\rm Area_{>2^{\circ}C}$', 'area_colder1': r'$\rm Area_{<1^{\circ}C}$'})
    std_anom.rename(columns={'MEAN': r'$\rm \overline{x}$', 'SD': r'sd'}, inplace=True)
    
    # Get text values +  cell color
    year_list.append(r'$\rm \overline{x}$') # add 2 extra columns
    year_list.append(r'sd')   
    vals = np.around(std_anom.values,1)
    vals[vals==-0.] = 0.
    vals_color = vals.copy()
    vals_color[-1,] = vals_color[-1,]*-1 # Reverse last row colorscale
    vals_color[:,-1] = 0 # No color to last two columns (mean and STD)
    vals_color[:,-2] = 0
    #vals_color[(vals_color<0.5) & (vals_color>-.5)] = 0.

    # Build the colormap
    vmin = -3.49
    vmax = 3.49
    midpoint = 0
    levels = np.linspace(vmin, vmax, 15)
    midp = np.mean(np.c_[levels[:-1], levels[1:]], axis=1)
    colvals = np.interp(midp, [vmin, midpoint, vmax], [-1, 0., 1])
    normal = plt.Normalize(-3.49, 3.49)
    reds = plt.cm.Reds(np.linspace(0,1, num=7))
    blues = plt.cm.Blues_r(np.linspace(0,1, num=7))
    whites = [(1,1,1,1)]*2
    colors = np.vstack((blues[0:-1,:], whites, reds[1:,:]))
    colors = np.concatenate([[colors[0,:]], colors, [colors[-1,:]]], 0)
    cmap, norm = from_levels_and_colors(levels, colors, extend='both')
    cmap_r, norm_r = from_levels_and_colors(levels, np.flipud(colors), extend='both')
    # Common parameters
    #hcell, wcell = 0.5, 0.6
    #hpad, wpad = 0, 0

    ## normal = plt.Normalize(-4.49, 4.49)
    ## cmap = plt.cm.get_cmap('seismic', 9) 
    #cmap = plt.cm.get_cmap('seismic', 15) 

    nrows, ncols = std_anom.index.size+1, std_anom.columns.size
    hcell, wcell = 0.5, 0.5
    hpad, wpad = 1, 1    
    fig=plt.figure(figsize=(ncols*wcell+wpad, nrows*hcell+hpad))
    ax = fig.add_subplot(111)
    ax.axis('off')
    #do the table
    header = ax.table(cellText=[['']],
                          colLabels=['-- NAFO division 2H --'],
                          loc='center'
                          )
    header.set_fontsize(13)
    #the_table=ax.table(cellText=vals, rowLabels=std_anom.index, colLabels=std_anom.columns, 
    the_table=ax.table(cellText=vals, rowLabels=std_anom.index, colLabels=year_list,
                        loc='center', cellColours=cmap(normal(vals_color)), cellLoc='center',
                        bbox=[0, 0, 1, 0.5]
                        )
    # change font color to white where needed:
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(12.5)
    table_props = the_table.properties()
    #table_cells = table_props['child_artists']
    last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
    for key, cell in the_table.get_celld().items():
        cell_text = cell.get_text().get_text()
        if is_number(cell_text) == False:
            pass
        elif key[0] == 0: #year's row = no color
            pass
        elif key[1] in last_columns:
             cell._text.set_color('darkslategray')
        elif (np.float(cell_text) <= -1.5) | (np.float(cell_text) >= 1.5) :
            cell._text.set_color('white')
        elif (cell_text=='nan'):
            cell._set_facecolor('lightgray')
            cell._text.set_color('lightgray')

    plt.savefig("scorecards_fall_2H.png", dpi=300)
    os.system('convert -trim scorecards_fall_2H.png scorecards_fall_2H.png')

    # French table
    std_anom = std_anom.rename({r'$\rm T_{bot}$' : r'$\rm T_{fond}$', r'$\rm T_{bot_{<200m}}$' : r'$\rm T_{fond_{<200m}}$', r'$\rm Area_{>2^{\circ}C}$' : r'$\rm Aire_{>2^{\circ}C}$', r'$\rm Area_{<1^{\circ}C}$' : r'$\rm Aire_{<1^{\circ}C}$'})
    year_list[-1] = u'ET'

    header = ax.table(cellText=[['']],
                          colLabels=['-- Division 2H de l\'OPANO --'],
                          loc='center'
                          )
    header.set_fontsize(13)
    #the_table=ax.table(cellText=vals, rowLabels=std_anom.index, colLabels=std_anom.columns, 
    the_table=ax.table(cellText=vals, rowLabels=std_anom.index, colLabels=year_list,
                        loc='center', cellColours=cmap(normal(vals_color)), cellLoc='center',
                        bbox=[0, 0, 1, 0.5]
                        )
    # change font color to white where needed:
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(12.5)
    table_props = the_table.properties()
    #table_cells = table_props['child_artists']
    last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
    for key, cell in the_table.get_celld().items():
        cell_text = cell.get_text().get_text()
        if is_number(cell_text) == False:
            pass
        elif key[0] == 0: #year's row = no color
            pass
        elif key[1] in last_columns:
             cell._text.set_color('darkslategray')
        elif (np.float(cell_text) <= -1.5) | (np.float(cell_text) >= 1.5) :
            cell._text.set_color('white')
        elif (cell_text=='nan'):
            cell._set_facecolor('lightgray')
            cell._text.set_color('lightgray')

    plt.savefig("scorecards_fall_2H_FR.png", dpi=300)
    os.system('convert -trim scorecards_fall_2H_FR.png scorecards_fall_2H_FR.png')
