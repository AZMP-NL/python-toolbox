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
sys.path.append(os.path.expanduser('~/github/AZMP-NL/python-toolbox/azmp_modules'))
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


def bottom_temperature(
    season,
    lonLims,
    latLims,
    year,
    clim_fill=True,
    time_adjust=True,
    netcdf_path='~/data/CABOTS/',
    CASTS_path='~/data/CASTS/',
    climato_file=''):
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

    if netcdf_path.endswith('.nc') == False:
        if season=='spring':
            year_file=netcdf_path+'CABOTS_bottomstats_spring.nc'
        elif season=='summer':
            year_file=netcdf_path+'CABOTS_bottomstats_summer.nc'
        elif season=='fall':
            year_file=netcdf_path+'CABOTS_bottomstats_fall.nc'
    else:
        year_file=netcdf_path

    #Determine if we're making NSRF plot
    if 'NSRF' in climato_file:
        NSRF_plot = True
        print('NSRF data detected.')
    else:
        NSRF_plot = False

    ## ---- Load Climato data ---- ##
    print('Load ' + climato_file)
    h5f = h5py.File(climato_file, 'r')
    Tbot_climato = h5f['Tbot'][:]
    lon_reg = h5f['lon_reg'][:][0,:]
    lat_reg = h5f['lat_reg'][:][:,0]
    Zitp = h5f['Zitp'][:]
    h5f.close()

    #Load in the CASTS data for cast locations
    ds = xr.open_dataset(os.path.expanduser(CASTS_path)+str(year)+'.nc')
    # Selection of a subset region
    ds = ds.sel(time=((ds.longitude>lonLims[0])*(ds.longitude<=lonLims[1])))
    ds = ds.sel(time=((ds.latitude>latLims[0])*(ds.latitude<=latLims[1])))
    # Select time (save several options here)
    if season == 'summer':
        ds = ds.sel(time=((ds['time.month']>=6)) & ((ds['time.month']<=9)))
    elif season == 'spring':
        ds = ds.sel(time=((ds['time.month']>=4)) & ((ds['time.month']<=6)))
    elif season == 'fall':
        ds = ds.sel(time=((ds['time.month']>=9)) & ((ds['time.month']<=12)))
    else:
        print('!! no season specified, used them all! !!')

    #Determine which temperature casts have measurements
    ds_temp_filt = ds.temperature.values
    ds_temp_filt = np.isnan(ds_temp_filt).sum(axis=1) != ds.level.size
    ds = ds.sel(time=ds_temp_filt)
    lons = ds.longitude.values
    lats = ds.latitude.values

    ## ---- NAFO divisions ---- ##
    nafo_div = azu.get_nafo_divisions()

    ## ---- SFA divisions ---- ##
    myshp = open(os.path.expanduser('~/github/AZMP-NL/utils/SFAs/SFAs_PANOMICS_Fall2020_shp/SFAs_PANOMICS_Fall2020.shp'), 'rb')
    mydbf = open(os.path.expanduser('~/github/AZMP-NL/utils/SFAs/SFAs_PANOMICS_Fall2020_shp/SFAs_PANOMICS_Fall2020.dbf'), 'rb')
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

    ## ---- Get bottom_temperature data --- ##
    print('Get ' + year_file)
    ds = xr.open_dataset(os.path.expanduser(year_file))
    #Isolate for the year of interest
    ds = ds.sel(TIME=ds['TIME.year']==int(year))
    ds = ds.mean('TIME')

    # Selection of a subset region
    ds = ds.sel(X=((ds.LONGITUDE[0,:]>=lonLims[0])*(ds.LONGITUDE[0,:]<=lonLims[1])).values)
    ds = ds.sel(Y=((ds.LATITUDE[:,0]>=latLims[0])*(ds.LATITUDE[:,0]<=latLims[1])).values)
    if time_adjust:
        Tbot = ds.BOTTOM_TEMPERATURE_ADJUSTED.values
    else:
        Tbot = ds.BOTTOM_TEMPERATURE.values

    # Use climatology to fill missing pixels
    if clim_fill:
        print('Fill NaNs with climatology')
        Tbot_orig = Tbot.copy()
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
    sfa2 = Polygon(shrimp_area['2'])
    sfa3 = Polygon(shrimp_area['3'])
    sfa4 = Polygon(shrimp_area['4'])
    sfa5 = Polygon(shrimp_area['5'])
    sfa6 = Polygon(shrimp_area['6'])
    sfa7 = Polygon(shrimp_area['7'])

    if NSRF_plot:
        if season == 'summer':
            for i, xx in enumerate(lon_reg):
                for j,yy in enumerate(lat_reg):
                    point = Point(lon_reg[i], lat_reg[j])
                    if sfa2.contains(point) | sfa3.contains(point) | sfa4.contains(point):
                        pass #nothing to do but cannot implement negative statement "if not" above
                    else:
                        Tbot[j,i] = np.nan
                        Tbot_orig[j,i] = np.nan
        if season == 'fall':
            for i, xx in enumerate(lon_reg):
                for j,yy in enumerate(lat_reg):
                    point = Point(lon_reg[i], lat_reg[j])
                    if sfa4.contains(point) | sfa5.contains(point) | sfa6.contains(point) | sfa7.contains(point):
                        pass #nothing to do but cannot implement negative statement "if not" above
                    else:
                        Tbot[j,i] = np.nan
                        Tbot_orig[j,i] = np.nan
    else:
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

        elif season == 'summer':
            for i, xx in enumerate(lon_reg):
                for j,yy in enumerate(lat_reg):
                    point = Point(lon_reg[i], lat_reg[j])


        else:
            print('no division mask, all data taken')

    print(' -> Done!')

    # Temperature anomaly:
    anom = Tbot-Tbot_climato
    if NSRF_plot:
        div_toplot = ['sfa2','sfa3','sfa4']
    else:
        div_toplot = ['2H', '2J', '3K', '3L', '3N', '3O', '3Ps', '4R']

    #Set up the map
    land_10m = cfeature.NaturalEarthFeature('physical','land','10m',
    facecolor='tan')

    ## ---- Plot Anomaly ---- ##
    fig = plt.figure(figsize=(6,9))
    ax = plt.subplot(projection=ccrs.Mercator())
    ax._autoscaleXon = False
    ax._autoscaleYon = False

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
    if NSRF_plot:
        gl = ax.gridlines(draw_labels=['bottom'], xlocs=np.arange(-68,-54,2), ylocs=np.arange(58,68,2),
            dms=True, x_inline=False, y_inline=False, linestyle='--')
    else:
        gl = ax.gridlines(draw_labels=['bottom'], xlocs=[-50,-55,-60], ylocs=[40, 45, 50, 55, 60],
            dms=True, x_inline=False, y_inline=False, linestyle='--')
    gl.xlabel_style = {'size': 8}

    #Add a colourbar
    cax = fig.add_axes([0.15, 0.07, 0.73, 0.0125])
    cb = plt.colorbar(c, cax=cax, orientation='horizontal')
    cb.set_label(r'$\rm T(^{\circ}C)$', fontsize=8, fontweight='normal')
    cb.ax.tick_params(labelsize=8)

    #Plot the divisions
    for div in div_toplot:
        if NSRF_plot:
            div_lon, div_lat = shrimp_area[div[-1]][:,0], shrimp_area[div[-1]][:,1]
            ax.plot(div_lon, div_lat, 'k', linewidth=1, transform=ccrs.PlateCarree())
            ax.text(np.mean(div_lon), np.mean(div_lat), div.upper(), fontsize=8, color='black', fontweight='bold', transform=ccrs.PlateCarree()) 
        else:
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
    ax._autoscaleXon = False
    ax._autoscaleYon = False

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
    if NSRF_plot:
        gl = ax.gridlines(draw_labels=['bottom'], xlocs=np.arange(-68,-54,2), ylocs=np.arange(58,68,2),
            dms=True, x_inline=False, y_inline=False, linestyle='--')
    else:
        gl = ax.gridlines(draw_labels=['bottom'], xlocs=[-50,-55,-60], ylocs=[40, 45, 50, 55, 60],
            dms=True, x_inline=False, y_inline=False, linestyle='--')
    gl.xlabel_style = {'size': 8}

    #Add a colourbar
    cax = fig.add_axes([0.15, 0.07, 0.73, 0.0125])
    cb = plt.colorbar(c, cax=cax, orientation='horizontal')
    cb.set_label(r'$\rm T(^{\circ}C)$', fontsize=8, fontweight='normal')
    cb.ax.tick_params(labelsize=8)

    #Plot the divisions
    for div in div_toplot:
        if NSRF_plot:
            div_lon, div_lat = shrimp_area[div[-1]][:,0], shrimp_area[div[-1]][:,1]
            ax.plot(div_lon, div_lat, 'k', linewidth=1, transform=ccrs.PlateCarree())
            ax.text(np.mean(div_lon), np.mean(div_lat), div.upper(), fontsize=8, color='black', fontweight='bold', transform=ccrs.PlateCarree()) 
        else:
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
    ax._autoscaleXon = False
    ax._autoscaleYon = False

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
    if NSRF_plot:
        gl = ax.gridlines(draw_labels=['bottom','left'], xlocs=np.arange(-68,-54,2), ylocs=np.arange(58,68,2),
            dms=True, x_inline=False, y_inline=False, linestyle='--')
    else:
        gl = ax.gridlines(draw_labels=['bottom','left'], xlocs=[-50,-55,-60], ylocs=[40, 45, 50, 55, 60],
            dms=True, x_inline=False, y_inline=False, linestyle='--')
    gl.xlabel_style = {'size': 8}
    gl.ylabel_style = {'size': 8}

    #Add a colourbar
    cax = fig.add_axes([0.15, 0.07, 0.73, 0.0125])
    cb = plt.colorbar(c, cax=cax, orientation='horizontal')
    cb.set_label(r'$\rm T(^{\circ}C)$', fontsize=8, fontweight='normal')
    cb.ax.tick_params(labelsize=8)


    #Plot the divisions
    for div in div_toplot:
        if NSRF_plot:
            div_lon, div_lat = shrimp_area[div[-1]][:,0], shrimp_area[div[-1]][:,1]
            ax.plot(div_lon, div_lat, 'k', linewidth=1, transform=ccrs.PlateCarree())
            ax.text(np.mean(div_lon), np.mean(div_lat), div.upper(), fontsize=8, color='black', fontweight='bold', transform=ccrs.PlateCarree()) 
        else:
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
    if NSRF_plot:
        save_end = 'sfa_'
    else:
        save_end = ''
    os.system('montage bottom_temp_climato_' + season + '_' + year + '.png bottom_temp_' + season + '_' + year + '.png bottom_temp_anomaly_' + season + '_' + year + '.png  -tile 3x1 -geometry +10+10  -background white  '+save_end+'bottomT_' + season + year + '.png')
    # in French
    os.system('montage bottom_temp_climato_' + season + '_' + year + '_FR.png bottom_temp_' + season + '_' + year + '_FR.png bottom_temp_anomaly_' + season + '_' + year + '_FR.png  -tile 3x1 -geometry +10+10  -background white  '+save_end+'bottomT_' + season + year + '_FR.png')
    # Move to year folder
    os.system('cp '+save_end+'bottomT_' + season + year + '.png '+save_end+'bottomT_' + season + year + '_FR.png ' + year)


#### bottom_salinity
def bottom_salinity(
    season,
    year,
    lonLims,
    latLims,
    clim_fill=True,
    time_adjust=True,
    netcdf_path='~/data/CABOTS/',
    CASTS_path='~/data/CASTS/',
    climato_file=''):
    '''
    Bottom salinity maps for AZMP ResDoc

    This function is adapted from azmp_bottomS.py

    * process in : /home/cyrf0006/AZMP/state_reports/bottomT

    usage ex:
    >> import azmp_report_tools as azrt
    >> azrt.bottom_salinity(season='spring', year='2019'):

    Frederic.Cyr@dfo-mpo.gc.ca - January 2020

    '''

    if len(climato_file) == 0:
        if season=='spring':
            climato_file='operation_files/Sbot_climato_spring_0.10.h5'
        elif season=='fall':
            climato_file='operation_files/Sbot_climato_fall_0.10.h5'
        elif season=='summer':
            climato_file='operation_files/Sbot_climato_summer_0.10.h6'

    if netcdf_path.endswith('.nc') == False:
        if season=='spring':
            year_file=netcdf_path+'CABOTS_bottomstats_spring.nc'
        elif season=='summer':
            year_file=netcdf_path+'CABOTS_bottomstats_summer.nc'
        elif season=='fall':
            year_file=netcdf_path+'CABOTS_bottomstats_fall.nc'
    else:
        year_file=netcdf_path

    #Determine if we're making NSRF plot
    if 'NSRF' in climato_file:
        NSRF_plot = True
        print('NSRF data detected.')
    else:
        NSRF_plot = False

    ## ---- Load Climato data ---- ##    
    print('Load ' + climato_file)
    h5f = h5py.File(climato_file, 'r')
    Sbot_climato = h5f['Sbot'][:]
    lon_reg = h5f['lon_reg'][:][0,:]
    lat_reg = h5f['lat_reg'][:][:,0]
    Zitp = h5f['Zitp'][:]
    h5f.close()

    #Load in the CASTS data for cast locations
    ds = xr.open_dataset(os.path.expanduser(CASTS_path)+str(year)+'.nc')
    # Selection of a subset region
    ds = ds.sel(time=((ds.longitude>lonLims[0])*(ds.longitude<=lonLims[1])))
    ds = ds.sel(time=((ds.latitude>latLims[0])*(ds.latitude<=latLims[1])))
    # Select time (save several options here)
    if season == 'summer':
        ds = ds.sel(time=((ds['time.month']>=6)) & ((ds['time.month']<=9)))
    elif season == 'spring':
        ds = ds.sel(time=((ds['time.month']>=4)) & ((ds['time.month']<=6)))
    elif season == 'fall':
        ds = ds.sel(time=((ds['time.month']>=9)) & ((ds['time.month']<=12)))
    else:
        print('!! no season specified, used them all! !!')

    #Determine which salinity casts have measurements
    ds_saln_filt = ds.salinity.values
    ds_saln_filt = np.isnan(ds_saln_filt).sum(axis=1) != ds.level.size
    ds = ds.sel(time=ds_saln_filt)
    lons = ds.longitude.values
    lats = ds.latitude.values

    ## ---- NAFO divisions ---- ##
    nafo_div = azu.get_nafo_divisions()

    ## ---- SFA divisions ---- ##
    myshp = open(os.path.expanduser('~/github/AZMP-NL/utils/SFAs/SFAs_PANOMICS_Fall2020_shp/SFAs_PANOMICS_Fall2020.shp'), 'rb')
    mydbf = open(os.path.expanduser('~/github/AZMP-NL/utils/SFAs/SFAs_PANOMICS_Fall2020_shp/SFAs_PANOMICS_Fall2020.dbf'), 'rb')
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

    ## ---- Get bottom_salinity data --- ##
    print('Get ' + year_file)
    ds = xr.open_dataset(os.path.expanduser(year_file))
    #Isolate for the year of interest
    ds = ds.sel(TIME=ds['TIME.year']==int(year))
    ds = ds.mean('TIME')

    # Selection of a subset region
    ds = ds.sel(X=((ds.LONGITUDE[0,:]>=lonLims[0])*(ds.LONGITUDE[0,:]<=lonLims[1])).values)
    ds = ds.sel(Y=((ds.LATITUDE[:,0]>=latLims[0])*(ds.LATITUDE[:,0]<=latLims[1])).values)
    if time_adjust:
        Sbot = ds.BOTTOM_SALINITY_ADJUSTED.values
    else:
        Sbot = ds.BOTTOM_SALINITY.values

    # Use climatology to fill missing pixels
    if clim_fill:
        print('Fill NaNs with climatology')
        Sbot_orig = Sbot.copy()
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
    sfa2 = Polygon(shrimp_area['2'])
    sfa3 = Polygon(shrimp_area['3'])
    sfa4 = Polygon(shrimp_area['4'])
    sfa5 = Polygon(shrimp_area['5'])
    sfa6 = Polygon(shrimp_area['6'])
    sfa7 = Polygon(shrimp_area['7'])

    if NSRF_plot:
        if season == 'summer':
            for i, xx in enumerate(lon_reg):
                for j,yy in enumerate(lat_reg):
                    point = Point(lon_reg[i], lat_reg[j])
                    if sfa2.contains(point) | sfa3.contains(point) | sfa4.contains(point):
                        pass #nothing to do but cannot implement negative statement "if not" above
                    else:
                        Sbot[j,i] = np.nan
                        Sbot_orig[j,i] = np.nan
        if season == 'fall':
            for i, xx in enumerate(lon_reg):
                for j,yy in enumerate(lat_reg):
                    point = Point(lon_reg[i], lat_reg[j])
                    if sfa4.contains(point) | sfa5.contains(point) | sfa6.contains(point) | sfa7.contains(point):
                        pass #nothing to do but cannot implement negative statement "if not" above
                    else:
                        Sbot[j,i] = np.nan
                        Sbot_orig[j,i] = np.nan
    else:
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
                        Sbot[j,i] = np.nan ### <--------------------- Do mask the fall / OR / 
                        Sbot_orig[j,i] = np.nan
                        #Tbot[j,i] = np.nan ### <--------------------- Do not mask the fall!!!!!

        elif season == 'summer':
            for i, xx in enumerate(lon_reg):
                for j,yy in enumerate(lat_reg):
                    point = Point(lon_reg[i], lat_reg[j])


        else:
            print('no division mask, all data taken')

    print(' -> Done!')

    # Salinity anomaly:
    anom = Sbot-Sbot_climato
    if NSRF_plot:
        div_toplot = ['sfa2','sfa3','sfa4']
    else:
        div_toplot = ['2H', '2J', '3K', '3L', '3N', '3O', '3Ps', '4R']

    #Set up the map
    land_10m = cfeature.NaturalEarthFeature('physical','land','10m',
    facecolor='tan')

    ## ---- Plot Anomaly ---- ##
    fig = plt.figure(figsize=(6,9))
    ax = plt.subplot(projection=ccrs.Mercator())
    ax._autoscaleXon = False
    ax._autoscaleYon = False

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
    if NSRF_plot:
        gl = ax.gridlines(draw_labels=['bottom'], xlocs=np.arange(-68,-54,2), ylocs=np.arange(58,68,2),
            dms=True, x_inline=False, y_inline=False, linestyle='--')
    else:
        gl = ax.gridlines(draw_labels=['bottom'], xlocs=[-50,-55,-60], ylocs=[40, 45, 50, 55, 60],
            dms=True, x_inline=False, y_inline=False, linestyle='--')
    gl.xlabel_style = {'size': 8}

    #Add a colourbar
    cax = fig.add_axes([0.15, 0.07, 0.73, 0.0125])
    cb = plt.colorbar(c, cax=cax, orientation='horizontal')
    cb.set_label(r'$\rm S$', fontsize=8, fontweight='normal')
    cb.ax.tick_params(labelsize=8)
    
    #Plot the divisions
    for div in div_toplot:
        if NSRF_plot:
            div_lon, div_lat = shrimp_area[div[-1]][:,0], shrimp_area[div[-1]][:,1]
            ax.plot(div_lon, div_lat, 'k', linewidth=1, transform=ccrs.PlateCarree())
            ax.text(np.mean(div_lon), np.mean(div_lat), div.upper(), fontsize=8, color='black', fontweight='bold', transform=ccrs.PlateCarree()) 
        else:
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
    ax._autoscaleXon = False
    ax._autoscaleYon = False

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
    if NSRF_plot:
        gl = ax.gridlines(draw_labels=['bottom'], xlocs=np.arange(-68,-54,2), ylocs=np.arange(58,68,2),
            dms=True, x_inline=False, y_inline=False, linestyle='--')
    else:
        gl = ax.gridlines(draw_labels=['bottom'], xlocs=[-50,-55,-60], ylocs=[40, 45, 50, 55, 60],
            dms=True, x_inline=False, y_inline=False, linestyle='--')
    gl.xlabel_style = {'size': 8}

    #Add a colourbar
    cax = fig.add_axes([0.15, 0.07, 0.73, 0.0125])
    cb = plt.colorbar(c, cax=cax, orientation='horizontal')
    cb.set_label(r'$\rm S$', fontsize=8, fontweight='normal')
    cb.ax.tick_params(labelsize=8)
    
    #Plot the divisions
    for div in div_toplot:
        if NSRF_plot:
            div_lon, div_lat = shrimp_area[div[-1]][:,0], shrimp_area[div[-1]][:,1]
            ax.plot(div_lon, div_lat, 'k', linewidth=1, transform=ccrs.PlateCarree())
            ax.text(np.mean(div_lon), np.mean(div_lat), div.upper(), fontsize=8, color='black', fontweight='bold', transform=ccrs.PlateCarree()) 
        else:
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
    ax._autoscaleXon = False
    ax._autoscaleYon = False

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

    if NSRF_plot:
        gl = ax.gridlines(draw_labels=['bottom','left'], xlocs=np.arange(-68,-54,2), ylocs=np.arange(58,68,2),
            dms=True, x_inline=False, y_inline=False, linestyle='--')
    else:
        gl = ax.gridlines(draw_labels=['bottom','left'], xlocs=[-50,-55,-60], ylocs=[40, 45, 50, 55, 60],
            dms=True, x_inline=False, y_inline=False, linestyle='--')
    gl.xlabel_style = {'size': 8}
    gl.ylabel_style = {'size': 8}

    #Add a colourbar
    cax = fig.add_axes([0.15, 0.07, 0.73, 0.0125])
    cb = plt.colorbar(c, cax=cax, orientation='horizontal')
    cb.set_label(r'$\rm S$', fontsize=8, fontweight='normal')
    cb.ax.tick_params(labelsize=8)
    
    #Plot the divisions
    for div in div_toplot:
        if NSRF_plot:
            div_lon, div_lat = shrimp_area[div[-1]][:,0], shrimp_area[div[-1]][:,1]
            ax.plot(div_lon, div_lat, 'k', linewidth=1, transform=ccrs.PlateCarree())
            ax.text(np.mean(div_lon), np.mean(div_lat), div.upper(), fontsize=8, color='black', fontweight='bold', transform=ccrs.PlateCarree()) 
        else:
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
    if NSRF_plot:
        save_end = 'sfa'
    else:
        save_end = ''
    os.system('montage bottom_sal_climato_' + season + '_' + year + '.png bottom_sal_' + season + '_' + year + '.png bottom_sal_anomaly_' + season + '_' + year + '.png  -tile 3x1 -geometry +10+10  -background white  '+save_end+'_bottomS_' + season + year + '.png') 
    # French
    os.system('montage bottom_sal_climato_' + season + '_' + year + '_FR.png bottom_sal_' + season + '_' + year + '_FR.png bottom_sal_anomaly_' + season + '_' + year + '_FR.png  -tile 3x1 -geometry +10+10  -background white  '+save_end+'_bottomS_' + season + year + '_FR.png') 
    # Move to year folder
    os.system('cp '+save_end+'_bottomS_' + season + year + '.png '+save_end+'_bottomS_' + season + year + '_FR.png ' + year)


#### bottom_stats
def bottom_stats(
    years,
    season,
    time_adjust=True,
    plot=False,
    netcdf_path='',
    var='',
    climato_file=''):

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

    #Finding the mean of an empty slice or double scalar errorys
    np.seterr(divide='ignore', invalid='ignore')
    warnings.simplefilter("ignore", category=RuntimeWarning)

    #Label either S or T
    if var == 'temperature':
        var_strt = 'T'
    elif var == 'salinity':
        var_strt = 'S'

    # load climato
    if len(climato_file) == 0:
        if season=='spring':
            climato_file='operation_files/'+var_strt+'bot_climato_spring_0.10.h5'
        elif season=='fall':
            climato_file='operation_files/'+var_strt+'bot_climato_fall_0.10.h5'
        elif season=='summer':
            climato_file='operation_files/'+var_strt+'bot_climato_summer_0.10.h5'
    else:
        print('Climato file provided')

    h5f = h5py.File(climato_file, 'r')
    if var == 'temperature':
        Tbot_climato = h5f['Tbot'][:]
    elif var == 'salinity':
        Tbot_climato = h5f['Sbot'][:]
    lon_reg = h5f['lon_reg'][:][0,:]
    lat_reg = h5f['lat_reg'][:][:,0]
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

    #Get the bottom temperature for each year
    year_file = netcdf_path
    if var == 'temperature':
        Tdict_together = azu.get_bottomT(year_file,years,season,climato_file,time_adjust=time_adjust,lab_mask=False)
    elif var == 'salinity':
        Tdict_together = azu.get_bottomS(year_file,years,season,climato_file,time_adjust=time_adjust,lab_mask=False)

    #Cycle through and fill with climatology
    for year in Tdict_together:
        Tdict_together[year][var_strt+'bot_orig'] = Tdict_together[year][var_strt+'bot']
        Tbot = Tdict_together[year][var_strt+'bot']
        place1 = Tbot_climato.copy()
        place1[~np.isnan(Tbot)] = Tbot[~np.isnan(Tbot)]
        Tbot = place1
        Tdict_together[year][var_strt+'bot'] = Tbot

    # NAFO division stats    
    dict_stats_2GH = azu.polygon_temperature_stats(Tdict_together, shape_2GH, var=var)
    dict_stats_2G = azu.polygon_temperature_stats(Tdict_together, shape_2G, var=var)
    dict_stats_2H = azu.polygon_temperature_stats(Tdict_together, shape_2H, var=var)
    dict_stats_2J = azu.polygon_temperature_stats(Tdict_together, shape_2J, var=var)
    dict_stats_2HJ = azu.polygon_temperature_stats(Tdict_together, shape_2HJ, var=var)
    dict_stats_3LNO = azu.polygon_temperature_stats(Tdict_together, shape_3LNO, var=var)
    dict_stats_3M = azu.polygon_temperature_stats(Tdict_together, shape_3M, var=var)
    dict_stats_3Ps = azu.polygon_temperature_stats(Tdict_together, shape_3Ps, var=var)
    dict_stats_3K = azu.polygon_temperature_stats(Tdict_together, shape_3K, var=var)
    dict_stats_3L = azu.polygon_temperature_stats(Tdict_together, shape_3L, var=var)
    dict_stats_3O = azu.polygon_temperature_stats(Tdict_together, shape_3O, var=var)

    # Append bottom temperature for multi-index export
    df_list = []
    for year in Tdict_together:
        df = pd.DataFrame(index=lat_reg, columns=lon_reg)
        df.index.name='latitude'
        df.columns.name='longitude'
        df[:] = Tdict_together[year][var_strt+'bot']
        df_list.append(df)

    #Convert to data frames
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

    if var == 'salinity':
        season = season + '_' + var
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
    df_mindex.to_pickle(season + '_bottom_'+var+'.pkl')

    #### bottom_stats
def sfa_bottom_stats(
    years,
    season,
    time_adjust=True,
    netcdf_path='',
    var='',
    sfas=[2,3,4],
    climato_file='',):

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
    lon_reg = h5f['lon_reg'][:][0,:]
    lat_reg = h5f['lat_reg'][:][:,0]
    Zitp = h5f['Zitp'][:]
    h5f.close()

    # Derive some map parameters
    lon_0 = np.round(np.mean(lon_reg))
    lat_0 = np.round(np.mean(lat_reg))
    lonLims = [lon_reg[0], lon_reg[-1]]
    latLims = [lat_reg[0], lat_reg[-1]]

    ## Get SFAs data
    myshp = open(os.path.expanduser('~/github/AZMP-NL/utils/SFAs/SFAs_PANOMICS_Fall2020_shp/SFAs_PANOMICS_Fall2020.shp'), 'rb')
    mydbf = open(os.path.expanduser('~/github/AZMP-NL/utils/SFAs/SFAs_PANOMICS_Fall2020_shp/SFAs_PANOMICS_Fall2020.dbf'), 'rb')
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

    #Get the bottom temp, saln for each year
    year_file = netcdf_path
    if var == 'temperature':
        Tdict_together = azu.get_bottomT(year_file,years,season,climato_file,time_adjust=time_adjust,lab_mask=False)
    elif var == 'salinity':
        Tdict_together = azu.get_bottomS(year_file,years,season,climato_file,time_adjust=time_adjust,lab_mask=False)

    #Cycle through and fill with climatology
    for year in Tdict_together:
        if var == 'temperature':
            var_l = 'T'
        elif var == 'salinity':
            var_l = 'S'
        Tdict_together[year][var_l+'bot_orig'] = Tdict_together[year][var_l+'bot']
        Tbot = Tdict_together[year][var_l+'bot']
        place1 = Tbot_climato.copy()
        place1[~np.isnan(Tbot)] = Tbot[~np.isnan(Tbot)]
        Tbot = place1
        Tdict_together[year][var_l+'bot'] = Tbot

    # NAFO division stats
    if season == 'summer':
        dict_stats_sfa2 = azu.polygon_temperature_stats(Tdict_together, sfa2, var=var)
        dict_stats_sfa3 = azu.polygon_temperature_stats(Tdict_together, sfa3, var=var)
        dict_stats_sfa4 = azu.polygon_temperature_stats(Tdict_together, sfa4, var=var)
    elif season=='fall':
        dict_stats_sfa4 = azu.polygon_temperature_stats(Tdict_together, sfa4, var=var)
        dict_stats_sfa5 = azu.polygon_temperature_stats(Tdict_together, sfa5, var=var)
        dict_stats_sfa6 = azu.polygon_temperature_stats(Tdict_together, sfa6, var=var)
        dict_stats_sfa7 = azu.polygon_temperature_stats(Tdict_together, sfa7, var=var)

    # Append bottom temperature for multi-index export
    df_list = []
    for year in Tdict_together:
        df = pd.DataFrame(index=lat_reg, columns=lon_reg)
        df.index.name='latitude'
        df.columns.name='longitude'
        df[:] = Tdict_together[year][var_l+'bot']
        df_list.append(df)

    #Convert to data frames
    df_sfa2 = pd.DataFrame.from_dict(dict_stats_sfa2, orient='index')
    df_sfa3 = pd.DataFrame.from_dict(dict_stats_sfa3, orient='index')
    df_sfa4 = pd.DataFrame.from_dict(dict_stats_sfa4, orient='index')
    df_sfa5 = pd.DataFrame.from_dict(dict_stats_sfa5, orient='index')
    df_sfa6 = pd.DataFrame.from_dict(dict_stats_sfa6, orient='index')
    df_sfa7 = pd.DataFrame.from_dict(dict_stats_sfa7, orient='index')

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
    year_index = pd.Series(years)
    year_index.name='year'
    df_mindex = pd.concat(df_list,keys=year_index)
    df_mindex.to_pickle(season + '_sfa_bottom_'+var+'.pkl')


def bottom_scorecards(years, clim_year=[1991, 2020]):


    '''
    To generate AZMP score cards for bottom temperature
    Uses pickled object generated by bottom_stats()

    usage example:


    bottom_scorecards(years=[1980, 2019], clim_year=[1981, 2010]):

    '''

    #### ------------- For fall ---------------- ####
    # 0.
    infile = 'operation_files/stats_2H_fall.pkl'
    df = pd.read_pickle(infile)
    df.index = pd.to_datetime(df.index) # update index to datetime
    df = df[(df.index.year>=years[0]) & (df.index.year<=years[-1])]
    percent_coverage = df.percent_coverage.values.copy().round(0)
    # Flag bad years (no or weak sampling):
    #bad_years = np.array([1980, 1982, 1984, 1985, 1986, 1987, 1988, 1989, 1990, 1992, 1993, 1994, 1995, 1996, 2000, 2002, 2003, 2005, 2007, 2009, 2022])
    bad_years = np.array([1982, 1984, 1985, 1986, 1988, 1989, 1990, 1992, 1993, 1994, 1995, 1996, 2000, 2002, 2003, 2005, 2007, 2009])
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
    std_anom = std_anom.reindex(['Tmean', 'Tmean_sha200', 'area_warmer2', 'area_colder1', 'percent_coverage'])
    std_anom = std_anom.rename({'Tmean': r'$\rm T_{bot}$', 'Tmean_sha200': r'$\rm T_{bot_{<200m}}$', 'area_warmer2': r'$\rm Area_{>2^{\circ}C}$', 'area_colder1': r'$\rm Area_{<1^{\circ}C}$', 'percent_coverage': '% Cov'})
    std_anom.rename(columns={'MEAN': r'$\rm \overline{x}$', 'SD': r'sd'}, inplace=True)

    #Re-fill the percent coverage
    std_anom.iloc[4][:-2] = percent_coverage
    std_anom.iloc[4][-2:] = np.nan

    # Get text values +  cell color
    year_list.append(r'$\rm \overline{x}$') # add 2 extra columns
    year_list.append(r'sd')   
    vals = np.around(std_anom.values,1)
    vals[vals==-0.] = 0.
    vals_color = vals.copy()
    vals_color[-2,] = vals_color[-2,]*-1 # Reverse last row colorscale
    vals_color[:,-1] = 0 # No color to last two columns (mean and STD)
    vals_color[:,-2] = 0
    vals_color[4,:] = 0
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
        elif key[0] == 5: #last row, percentages
            if cell_text == 'nan':
                cell._set_facecolor('lightgray')
                cell._text.set_color('lightgray')
            elif key[1] != -1:
                cell._text.set_text(int(float(cell_text)))
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
    std_anom = std_anom.rename({r'$\rm T_{bot}$' : r'$\rm T_{fond}$', r'$\rm T_{bot_{<200m}}$' : r'$\rm T_{fond_{<200m}}$', r'$\rm Area_{>2^{\circ}C}$' : r'$\rm Aire_{>2^{\circ}C}$', r'$\rm Area_{<1^{\circ}C}$' : r'$\rm Aire_{<1^{\circ}C}$', '% Coverage': '% Cov'})
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
        elif key[0] == 5:
            if cell_text == 'nan':
                cell._set_facecolor('lightgray')
                cell._text.set_color('lightgray')
            elif key[1] != -1:
                cell._text.set_text(int(float(cell_text)))
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
    infile = 'operation_files/stats_2J_fall.pkl'
    df = pd.read_pickle(infile)
    df.index = pd.to_datetime(df.index) # update index to datetime
    df = df[(df.index.year>=years[0]) & (df.index.year<=years[-1])]
    percent_coverage = df.percent_coverage.values.copy()
    # Flag bad years (no or weak sampling):
    #bad_years = np.array([1995, 2022])
    bad_years = np.array([1995])
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
    std_anom = std_anom.reindex(['Tmean', 'Tmean_sha200', 'area_warmer2', 'area_colder1', 'percent_coverage'])
    std_anom = std_anom.rename({'Tmean': r'$\rm T_{bot}$', 'Tmean_sha200': r'$\rm T_{bot_{<200m}}$', 'area_warmer2': r'$\rm Area_{>2^{\circ}C}$', 'area_colder1': r'$\rm Area_{<1^{\circ}C}$', 'percent_coverage': '% Cov'})
    std_anom.rename(columns={'MEAN': r'$\rm \overline{x}$', 'SD': r'sd'}, inplace=True)

    #Re-fill the percent coverage
    std_anom.iloc[4][:-2] = percent_coverage
    std_anom.iloc[4][-2:] = np.nan

    vals = np.around(std_anom.values,1)
    vals[vals==-0.] = 0.
    vals_color = vals.copy()
    vals_color[-2,] = vals_color[-2,]*-1
    vals_color[:,-1] = 0 # No color to last two columns (mean and STD)
    vals_color[:,-2] = 0 
    vals_color[4,:] = 0
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
        elif key[0] == 4:
            if cell_text == 'nan':
                cell._set_facecolor('lightgray')
                cell._text.set_color('lightgray')
            elif key[1] != -1:
                cell._text.set_text(int(float(cell_text)))
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
    std_anom = std_anom.rename({r'$\rm T_{bot}$' : r'$\rm T_{fond}$', r'$\rm T_{bot_{<200m}}$' : r'$\rm T_{fond_{<200m}}$', r'$\rm Area_{>2^{\circ}C}$' : r'$\rm Aire_{>2^{\circ}C}$', r'$\rm Area_{<1^{\circ}C}$' : r'$\rm Aire_{<1^{\circ}C}$', '% Coverage': '% Cov'})
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
        elif key[0] == 4:
            if cell_text == 'nan':
                cell._set_facecolor('lightgray')
                cell._text.set_color('lightgray')
            elif key[1] != -1:
                cell._text.set_text(int(float(cell_text)))
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
    infile = 'operation_files/stats_3K_fall.pkl'
    df = pd.read_pickle(infile)
    df.index = pd.to_datetime(df.index) # update index to datetime
    df = df[(df.index.year>=years[0]) & (df.index.year<=years[-1])]
    percent_coverage = df.percent_coverage.values.copy()
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
    std_anom = std_anom.reindex(['Tmean', 'Tmean_sha200', 'area_warmer2', 'area_colder1','percent_coverage'])
    std_anom = std_anom.rename({'Tmean': r'$\rm T_{bot}$', 'Tmean_sha200': r'$\rm T_{bot_{<200m}}$', 'area_warmer2': r'$\rm Area_{>2^{\circ}C}$', 'area_colder1': r'$\rm Area_{<1^{\circ}C}$', 'percent_coverage': '% Cov'})
    std_anom.rename(columns={'MEAN': r'$\rm \overline{x}$', 'SD': r'sd'}, inplace=True)

    #Re-fill the percent coverage
    std_anom.iloc[4][:-2] = percent_coverage
    std_anom.iloc[4][-2:] = np.nan

    vals = np.around(std_anom.values,1)
    vals[vals==-0.] = 0.
    vals_color = vals.copy()
    vals_color[-2,] = vals_color[-2,]*-1
    vals_color[:,-1] = 0 # No color to last two columns (mean and STD)
    vals_color[:,-2] = 0 
    vals_color[4,:] = 0
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
        elif key[0] == 4:
            if cell_text == 'nan':
                cell._set_facecolor('lightgray')
                cell._text.set_color('lightgray')
            elif key[1] != -1:
                cell._text.set_text(int(float(cell_text)))
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
    std_anom = std_anom.rename({r'$\rm T_{bot}$' : r'$\rm T_{fond}$', r'$\rm T_{bot_{<200m}}$' : r'$\rm T_{fond_{<200m}}$', r'$\rm Area_{>2^{\circ}C}$' : r'$\rm Aire_{>2^{\circ}C}$', r'$\rm Area_{<1^{\circ}C}$' : r'$\rm Aire_{<1^{\circ}C}$', '% Coverage': '% Cov'})
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
        elif key[0] == 4:
            if cell_text == 'nan':
                cell._set_facecolor('lightgray')
                cell._text.set_color('lightgray')
            elif key[1] != -1:
                cell._text.set_text(int(float(cell_text)))
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
    infile = 'operation_files/stats_3LNO_fall.pkl'
    df = pd.read_pickle(infile)
    df.index = pd.to_datetime(df.index) # update index to datetime
    df = df[(df.index.year>=years[0]) & (df.index.year<=years[-1])]
    percent_coverage = df.percent_coverage.values.copy()
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
    std_anom = std_anom.reindex(['Tmean', 'Tmean_sha200', 'area_warmer2', 'area_colder0','percent_coverage'])
    std_anom = std_anom.rename({'Tmean': r'$\rm T_{bot}$', 'Tmean_sha200': r'$\rm T_{bot_{<200m}}$', 'area_warmer2': r'$\rm Area_{>2^{\circ}C}$', 'area_colder0': r'$\rm Area_{<0^{\circ}C}$', 'percent_coverage': '% Cov'})
    std_anom.rename(columns={'MEAN': r'$\rm \overline{x}$', 'SD': r'sd'}, inplace=True)

    #Re-fill the percent coverage
    std_anom.iloc[4][:-2] = percent_coverage
    std_anom.iloc[4][-2:] = np.nan

    vals = np.around(std_anom.values,1)
    vals[vals==-0.] = 0.
    vals_color = vals.copy()
    vals_color[-2,] = vals_color[-2,]*-1
    vals_color[:,-1] = 0 # No color to last two columns (mean and STD)
    vals_color[:,-2] = 0
    vals_color[4,:] = 0
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
        elif key[0] == 4:
            if cell_text == 'nan':
                cell._set_facecolor('lightgray')
                cell._text.set_color('lightgray')
            elif key[1] != -1:
                cell._text.set_text(int(float(cell_text)))
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
    std_anom = std_anom.rename({r'$\rm T_{bot}$' : r'$\rm T_{fond}$', r'$\rm T_{bot_{<200m}}$' : r'$\rm T_{fond_{<200m}}$', r'$\rm Area_{>2^{\circ}C}$' : r'$\rm Aire_{>2^{\circ}C}$', r'$\rm Area_{<1^{\circ}C}$' : r'$\rm Aire_{<1^{\circ}C}$', '% Coverage': '% Cov'})
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
        elif key[0] == 4:
            if cell_text == 'nan':
                cell._set_facecolor('lightgray')
                cell._text.set_color('lightgray')
            elif key[1] != -1:
                cell._text.set_text(int(float(cell_text)))
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
    os.system('montage  scorecards_fall_2H.png scorecards_fall_2J.png scorecards_fall_3K.png scorecards_fall_3LNO.png -tile 1x4 -geometry +1+10  -background white  scorecards_botT_fall.png') 
    # French
    os.system('montage  scorecards_fall_2H_FR.png scorecards_fall_2J_FR.png scorecards_fall_3K_FR.png scorecards_fall_3LNO_FR.png -tile 1x4 -geometry +1+10  -background white  scorecards_botT_fall_FR.png') 




    #### ------------- For Spring ---------------- ####
    # 1.
    infile = 'operation_files/stats_3LNO_spring.pkl'
    df = pd.read_pickle(infile)
    df.index = pd.to_datetime(df.index) # update index to datetime
    df = df[(df.index.year>=years[0]) & (df.index.year<=years[-1])]
    percent_coverage = df.percent_coverage.values.copy()
    # Flag bad years (no or weak sampling):
    #bad_years = np.array([2020, 2021])
    bad_years = np.array([2020])
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
    std_anom = std_anom.reindex(['Tmean', 'Tmean_sha200', 'area_warmer2', 'area_colder0', 'percent_coverage'])
    std_anom = std_anom.rename({'Tmean': r'$\rm T_{bot}$', 'Tmean_sha200': r'$\rm T_{bot_{<200m}}$', 'area_warmer2': r'$\rm Area_{>2^{\circ}C}$', 'area_colder0': r'$\rm Area_{<0^{\circ}C}$', 'percent_coverage': '% Cov'})
    std_anom.rename(columns={'MEAN': r'$\rm \overline{x}$', 'SD': r'sd'}, inplace=True)

    #Re-fill the percent coverage
    std_anom.iloc[4][:-2] = percent_coverage
    std_anom.iloc[4][-2:] = np.nan

    year_list.append(r'$\rm \overline{x}$') # add 2 extra columns
    year_list.append(r'sd')
    vals = np.around(std_anom.values,1)
    vals[vals==-0.] = 0.
    vals_color = vals.copy()
    vals_color[-2,] = vals_color[-2,]*-1
    vals_color[:,-1] = 0 # No color to last two columns (mean and STD)
    vals_color[:,-2] = 0
    vals_color[4,:] = 0
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
        elif key[0] == 5:
            if cell_text == 'nan':
                cell._set_facecolor('lightgray')
                cell._text.set_color('lightgray')
            elif key[1] != -1:
                cell._text.set_text(int(float(cell_text)))
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
    std_anom = std_anom.rename({r'$\rm T_{bot}$' : r'$\rm T_{fond}$', r'$\rm T_{bot_{<200m}}$' : r'$\rm T_{fond_{<200m}}$', r'$\rm Area_{>2^{\circ}C}$' : r'$\rm Aire_{>2^{\circ}C}$', r'$\rm Area_{<1^{\circ}C}$' : r'$\rm Aire_{<1^{\circ}C}$', '% Coverage': '% Cov'})
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
        elif key[0] == 5:
            if cell_text == 'nan':
                cell._set_facecolor('lightgray')
                cell._text.set_color('lightgray')
            elif key[1] != -1:
                cell._text.set_text(int(float(cell_text)))
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
    infile = 'operation_files/stats_3Ps_spring.pkl'
    df = pd.read_pickle(infile)
    df.index = pd.to_datetime(df.index) # update index to datetime
    df = df[(df.index.year>=years[0]) & (df.index.year<=years[-1])]
    percent_coverage = df.percent_coverage.values.copy()
    # Flag bad years (no or weak sampling):
    #bad_years = np.array([1980, 1981, 1985, 1986, 1987, 1988, 1989, 1990, 1991, 1992, 2006, 2020])
    bad_years = np.array([1981, 1985, 1986, 1987, 1988, 1989, 1990, 1991, 1992, 2020])
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
    std_anom = std_anom.reindex(['Tmean', 'Tmean_sha200', 'area_warmer2', 'area_colder0', 'percent_coverage'])
    std_anom = std_anom.rename({'Tmean': r'$\rm T_{bot}$', 'Tmean_sha200': r'$\rm T_{bot_{<200m}}$', 'area_warmer2': r'$\rm Area_{>2^{\circ}C}$', 'area_colder0': r'$\rm Area_{<0^{\circ}C}$', 'percent_coverage': '% Cov'})
    std_anom.rename(columns={'MEAN': r'$\rm \overline{x}$', 'SD': r'sd'}, inplace=True)

    #Re-fill the percent coverage
    std_anom.iloc[4][:-2] = percent_coverage
    std_anom.iloc[4][-2:] = np.nan

    vals = np.around(std_anom.values,1)
    vals[vals==-0.] = 0.
    vals_color = vals.copy()
    vals_color[-2,] = vals_color[-2,]*-1
    vals_color[:,-1] = 0 # No color to last two columns (mean and STD)
    vals_color[:,-2] = 0 
    vals_color[4,:] = 0
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
        elif key[0] == 4:
            if cell_text == 'nan':
                cell._set_facecolor('lightgray')
                cell._text.set_color('lightgray')
            elif key[1] != -1:
                cell._text.set_text(int(float(cell_text)))
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
    std_anom = std_anom.rename({r'$\rm T_{bot}$' : r'$\rm T_{fond}$', r'$\rm T_{bot_{<200m}}$' : r'$\rm T_{fond_{<200m}}$', r'$\rm Area_{>2^{\circ}C}$' : r'$\rm Aire_{>2^{\circ}C}$', r'$\rm Area_{<1^{\circ}C}$' : r'$\rm Aire_{<1^{\circ}C}$', '% Coverage': '% Cov'})
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
        elif key[0] == 4:
            if cell_text == 'nan':
                cell._set_facecolor('lightgray')
                cell._text.set_color('lightgray')
            elif key[1] != -1:
                cell._text.set_text(int(float(cell_text)))
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
    os.system('montage  scorecards_spring_3LNO.png scorecards_spring_3Ps.png -tile 1x3 -geometry +1+10  -background white  scorecards_botT_spring.png')
    # French montage
    os.system('montage  scorecards_spring_3LNO_FR.png scorecards_spring_3Ps_FR.png -tile 1x3 -geometry +1+10  -background white  scorecards_botT_spring_FR.png')

    
def sfa_bottom_scorecards(years, clim_year=[2006, 2020]):


    '''
    Similar to bottom_scorcards, but for SAFs

    usage example (see azmp_genreport.py):
    azrt.sfa_bottom_scorecards(years=np.arange(2006, 2022), clim_year=[2006, 2020])

    '''

    #### ------------- For summer ---------------- ####
    # 0.
    infile = 'operation_files/stats_sfa2_summer.pkl'
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
        elif (float(cell_text) <= -1.5) | (float(cell_text) >= 1.5) :
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
    infile = 'operation_files/stats_sfa3_summer.pkl'
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
        elif (float(cell_text) <= -1.5) | (float(cell_text) >= 1.5) :
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
        elif (float(cell_text) <= -1.5) | (float(cell_text) >= 1.5) :
            cell._text.set_color('white')
        elif (cell_text=='nan'):
            cell._set_facecolor('lightgray')
            cell._text.set_color('lightgray')

    plt.savefig("scorecards_summer_SFA3_FR.png", dpi=300)
    os.system('convert -trim scorecards_summer_SFA3_FR.png scorecards_summer_SFA3_FR.png')


    # 2.
    infile = 'operation_files/stats_sfa4_summer.pkl'
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
        elif (float(cell_text) <= -1.5) | (float(cell_text) >= 1.5) :
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
        elif (float(cell_text) <= -1.5) | (float(cell_text) >= 1.5) :
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
        elif (float(cell_text) <= -1.5) | (float(cell_text) >= 1.5) :
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
        elif (float(cell_text) <= -1.5) | (float(cell_text) >= 1.5) :
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


def bottomS_scorecards(years, clim_year=[2006, 2021]):


    '''
    Similar to bottom_scorecards, but for salinity

    usage example (see azmp_genreport.py):
    azrt.sfa_bottomS_scorecards(years=np.arange(2006, 2022), clim_year=[2006, 2020])

    ** Note that varuables inside .pkl files are name "Tmean" isntead of "Smean"
    '''

    #### ------------- For summer ---------------- ####
    # 0.
    infile = 'operation_files/stats_2H_fall_salinity.pkl'
    df = pd.read_pickle(infile)
    df.index = pd.to_datetime(df.index) # update index to datetime
    df = df[(df.index.year>=years[0]) & (df.index.year<=years[-1])]
    percent_coverage = df.percent_coverage.values.copy().round(0)
    # Flag bad years (no or weak sampling):
    bad_years = df.index.year.values[percent_coverage < 80]
    for i in bad_years:
        df[df.index.year==i]=np.nan
    year_list = df.index.year.astype('str')
    year_list = [i[2:4] for i in year_list] # 2-digit year
    df_clim = df[(df.index.year>=clim_year[0]) & (df.index.year<=clim_year[1])]
    std_anom = (df-df_clim.mean(axis=0))/df_clim.std(axis=0)
    std_anom = std_anom.T
    std_anom['MEAN'] = df_clim.mean(axis=0)
    std_anom['SD'] = df_clim.std(axis=0)
    std_anom = std_anom.reindex(['Tmean', 'Tmean_sha200', 'percent_coverage'])
    std_anom = std_anom.rename({'Tmean': r'$\rm S_{bot}$', 'Tmean_sha200': r'$\rm S_{bot_{<200m}}$', 'percent_coverage': '% Cov'})
    std_anom.rename(columns={'MEAN': r'$\rm \overline{x}$', 'SD': r'sd'}, inplace=True)
    # Save in .csv for future use
    std_anom.to_csv('bottomS_stn_anom_2H_fall.csv', sep=',', float_format='%0.3f')

    #Re-fill the percent coverage
    std_anom.iloc[-1][:-2] = percent_coverage
    std_anom.iloc[-1][-2:] = np.nan

    # Get text values +  cell color
    year_list.append(r'$\rm \overline{x}$') # add 2 extra columns
    year_list.append(r'sd')   
    vals = np.around(std_anom.values,1)
    vals[vals==-0.] = 0.
    vals_color = vals.copy()
    vals_color[:,-1] = 0 # No color to last two columns (mean and STD)
    vals_color[:,-2] = 0
    vals_color[-1,:] = 0
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
        elif key[0] == 3:
            if cell_text == 'nan':
                cell._set_facecolor('lightgray')
                cell._text.set_color('lightgray')
            elif key[1] != -1:
                cell._text.set_text(int(float(cell_text)))
        elif key[1] in last_columns:
             cell._text.set_color('darkslategray')
        elif (float(cell_text) <= -1.5) | (float(cell_text) >= 1.5) :
            cell._text.set_color('white')
        elif (cell_text=='nan'):
            cell._set_facecolor('lightgray')
            cell._text.set_color('lightgray')

    plt.savefig("scorecards_fallS_2H.png", dpi=300)
    os.system('convert -trim scorecards_fallS_2H.png scorecards_fallS_2H.png')

    # French table
    std_anom = std_anom.rename({r'$\rm S_{bot}$' : r'$\rm S_{fond}$', r'$\rm S_{bot_{<200m}}$' : r'$\rm S_{fond_{<200m}}$', '% Coverage': '% Cov'})
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
        elif key[0] == 3:
            if cell_text == 'nan':
                cell._set_facecolor('lightgray')
                cell._text.set_color('lightgray')
            elif key[1] != -1:
                cell._text.set_text(int(float(cell_text)))
        elif key[1] in last_columns:
             cell._text.set_color('darkslategray')
        elif (float(cell_text) <= -1.5) | (float(cell_text) >= 1.5) :
            cell._text.set_color('white')
        elif (cell_text=='nan'):
            cell._set_facecolor('lightgray')
            cell._text.set_color('lightgray')

    plt.savefig("scorecards_fallS_2H_FR.png", dpi=300)
    os.system('convert -trim scorecards_fallS_2H_FR.png scorecards_fallS_2H_FR.png')

 # 1.
    infile = 'operation_files/stats_2J_fall_salinity.pkl'
    df = pd.read_pickle(infile)
    df.index = pd.to_datetime(df.index) # update index to datetime
    df = df[(df.index.year>=years[0]) & (df.index.year<=years[-1])]
    percent_coverage = df.percent_coverage.values.copy().round(0)
    # Flag bad years (no or weak sampling):
    bad_years = df.index.year.values[percent_coverage < 80]
    for i in bad_years:
        df[df.index.year==i]=np.nan
    df_clim = df[(df.index.year>=clim_year[0]) & (df.index.year<=clim_year[1])]
    std_anom = (df-df_clim.mean(axis=0))/df_clim.std(axis=0)
    std_anom = std_anom.T
    std_anom['MEAN'] = df_clim.mean(axis=0)
    std_anom['SD'] = df_clim.std(axis=0)
    std_anom = std_anom.reindex(['Tmean', 'Tmean_sha200', 'percent_coverage'])
    std_anom = std_anom.rename({'Tmean': r'$\rm S_{bot}$', 'Tmean_sha200': r'$\rm S_{bot_{<200m}}$', 'percent_coverage': '% Cov'})
    std_anom.rename(columns={'MEAN': r'$\rm \overline{x}$', 'SD': r'sd'}, inplace=True)
    # Save in .csv for future use
    std_anom.to_csv('bottomS_stn_anom_2J_fall.csv', sep=',', float_format='%0.3f')

    #Re-fill the percent coverage
    std_anom.iloc[-1][:-2] = percent_coverage
    std_anom.iloc[-1][-2:] = np.nan

    vals = np.around(std_anom.values,1)
    vals[vals==-0.] = 0.
    vals_color = vals.copy()
    vals_color[:,-1] = 0 # No color to last two columns (mean and STD)
    vals_color[:,-2] = 0
    vals_color[-1,:] = 0
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
        elif key[0] == 2:
            if cell_text == 'nan':
                cell._set_facecolor('lightgray')
                cell._text.set_color('lightgray')
            elif key[1] != -1:
                cell._text.set_text(int(float(cell_text)))
        elif key[1] in last_columns:
             cell._text.set_color('darkslategray')
        elif (float(cell_text) <= -1.5) | (float(cell_text) >= 1.5) :
            cell._text.set_color('white')
        elif (cell_text=='nan'):
            cell._set_facecolor('lightgray')
            cell._text.set_color('lightgray')

    plt.savefig("scorecards_fallS_2J.png", dpi=300)
    os.system('convert -trim scorecards_fallS_2J.png scorecards_fallS_2J.png')

    # French table
    std_anom = std_anom.rename({r'$\rm S_{bot}$' : r'$\rm S_{fond}$', r'$\rm S_{bot_{<200m}}$' : r'$\rm S_{fond_{<200m}}$', '% Coverage': '% Cov'})
    header = ax.table(cellText=[['']],
                          colLabels=['-- Division 2J de l\'OPANO --'],
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
        elif key[0] == 2:
            if cell_text == 'nan':
                cell._set_facecolor('lightgray')
                cell._text.set_color('lightgray')
            elif key[1] != -1:
                cell._text.set_text(int(float(cell_text)))
        elif key[1] in last_columns:
             cell._text.set_color('darkslategray')
        elif (float(cell_text) <= -1.5) | (float(cell_text) >= 1.5) :
            cell._text.set_color('white')
        elif (cell_text=='nan'):
            cell._set_facecolor('lightgray')
            cell._text.set_color('lightgray')

    plt.savefig("scorecards_fallS_2J_FR.png", dpi=300)
    os.system('convert -trim scorecards_fallS_2J_FR.png scorecards_fallS_2J_FR.png')


    # 2.
    infile = 'operation_files/stats_3K_fall_salinity.pkl'
    df = pd.read_pickle(infile)
    df.index = pd.to_datetime(df.index) # update index to datetime
    df = df[(df.index.year>=years[0]) & (df.index.year<=years[-1])]
    percent_coverage = df.percent_coverage.values.copy().round(0)
    # Flag bad years (no or weak sampling):
    bad_years = df.index.year.values[percent_coverage < 80]
    for i in bad_years:
        df[df.index.year==i]=np.nan
    df_clim = df[(df.index.year>=clim_year[0]) & (df.index.year<=clim_year[1])]
    std_anom = (df-df_clim.mean(axis=0))/df_clim.std(axis=0)
    std_anom = std_anom.T
    std_anom['MEAN'] = df_clim.mean(axis=0)
    std_anom['SD'] = df_clim.std(axis=0)
    std_anom = std_anom.reindex(['Tmean', 'Tmean_sha200', 'percent_coverage'])
    std_anom = std_anom.rename({'Tmean': r'$\rm S_{bot}$', 'Tmean_sha200': r'$\rm S_{bot_{<200m}}$', 'percent_coverage': '% Cov'})
    std_anom.rename(columns={'MEAN': r'$\rm \overline{x}$', 'SD': r'sd'}, inplace=True)
    # Save in .csv for future use
    std_anom.to_csv('bottomS_stn_anom_3K_fall.csv', sep=',', float_format='%0.3f')

    #Re-fill the percent coverage
    std_anom.iloc[-1][:-2] = percent_coverage
    std_anom.iloc[-1][-2:] = np.nan

    vals = np.around(std_anom.values,1)
    vals[vals==-0.] = 0.
    vals_color = vals.copy()
    vals_color[:,-1] = 0 # No color to last two columns (mean and STD)
    vals_color[:,-2] = 0
    vals_color[-1,:] = 0
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
        elif key[0] == 2:
            if cell_text == 'nan':
                cell._set_facecolor('lightgray')
                cell._text.set_color('lightgray')
            elif key[1] != -1:
                cell._text.set_text(int(float(cell_text)))
        elif key[1] in last_columns:
             cell._text.set_color('darkslategray')
        elif (float(cell_text) <= -1.5) | (float(cell_text) >= 1.5) :
            cell._text.set_color('white')
        elif (cell_text=='nan'):
            cell._set_facecolor('lightgray')
            cell._text.set_color('lightgray')

    plt.savefig("scorecards_fallS_3K.png", dpi=300)
    os.system('convert -trim scorecards_fallS_3K.png scorecards_fallS_3K.png')

    # French table
    std_anom = std_anom.rename({r'$\rm S_{bot}$' : r'$\rm S_{fond}$', r'$\rm S_{bot_{<200m}}$' : r'$\rm S_{fond_{<200m}}$', '% Coverage': '% Cov'})
    header = ax.table(cellText=[['']],
                          colLabels=['-- Division 3K de l\'OPANO --'],
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
        elif key[0] == 2:
            if cell_text == 'nan':
                cell._set_facecolor('lightgray')
                cell._text.set_color('lightgray')
            elif key[1] != -1:
                cell._text.set_text(int(float(cell_text)))
        elif key[1] in last_columns:
             cell._text.set_color('darkslategray')
        elif (float(cell_text) <= -1.5) | (float(cell_text) >= 1.5) :
            cell._text.set_color('white')
        elif (cell_text=='nan'):
            cell._set_facecolor('lightgray')
            cell._text.set_color('lightgray')

    plt.savefig("scorecards_fallS_3K_FR.png", dpi=300)
    os.system('convert -trim scorecards_fallS_3K_FR.png scorecards_fallS_3K_FR.png')
    plt.close('all')

    # 3.
    infile = 'operation_files/stats_3LNO_fall_salinity.pkl'
    df = pd.read_pickle(infile)
    df.index = pd.to_datetime(df.index) # update index to datetime
    df = df[(df.index.year>=years[0]) & (df.index.year<=years[-1])]
    percent_coverage = df.percent_coverage.values.copy().round(0)
    # Flag bad years (no or weak sampling):
    bad_years = df.index.year.values[percent_coverage < 80]
    for i in bad_years:
        df[df.index.year==i]=np.nan
    df_clim = df[(df.index.year>=clim_year[0]) & (df.index.year<=clim_year[1])]
    std_anom = (df-df_clim.mean(axis=0))/df_clim.std(axis=0)
    std_anom = std_anom.T
    std_anom['MEAN'] = df_clim.mean(axis=0)
    std_anom['SD'] = df_clim.std(axis=0)
    std_anom = std_anom.reindex(['Tmean', 'Tmean_sha200', 'percent_coverage'])
    std_anom = std_anom.rename({'Tmean': r'$\rm S_{bot}$', 'Tmean_sha200': r'$\rm S_{bot_{<200m}}$', 'percent_coverage': '% Cov'})
    std_anom.rename(columns={'MEAN': r'$\rm \overline{x}$', 'SD': r'sd'}, inplace=True)
    # Save in .csv for future use
    std_anom.to_csv('bottomS_stn_anom_3LNO_fall.csv', sep=',', float_format='%0.3f')

    #Re-fill the percent coverage
    std_anom.iloc[-1][:-2] = percent_coverage
    std_anom.iloc[-1][-2:] = np.nan

    vals = np.around(std_anom.values,1)
    vals[vals==-0.] = 0.
    vals_color = vals.copy()
    vals_color[:,-1] = 0 # No color to last two columns (mean and STD)
    vals_color[:,-2] = 0
    vals_color[-1,:] = 0
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
        elif key[0] == 2:
            if cell_text == 'nan':
                cell._set_facecolor('lightgray')
                cell._text.set_color('lightgray')
            elif key[1] != -1:
                cell._text.set_text(int(float(cell_text)))
        elif key[1] in last_columns:
             cell._text.set_color('darkslategray')
        elif (float(cell_text) <= -1.5) | (float(cell_text) >= 1.5) :
            cell._text.set_color('white')
        elif (cell_text=='nan'):
            cell._set_facecolor('lightgray')
            cell._text.set_color('lightgray')

    plt.savefig("scorecards_fallS_3LNO.png", dpi=300)
    os.system('convert -trim scorecards_fallS_3LNO.png scorecards_fallS_3LNO.png')

    # French table
    std_anom = std_anom.rename({r'$\rm S_{bot}$' : r'$\rm S_{fond}$', r'$\rm S_{bot_{<200m}}$' : r'$\rm S_{fond_{<200m}}$', '% Coverage': '% Cov'})
    header = ax.table(cellText=[['']],
                          colLabels=['-- Division 3LNO de l\'OPANO --'],
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
        elif key[0] == 2:
            if cell_text == 'nan':
                cell._set_facecolor('lightgray')
                cell._text.set_color('lightgray')
            elif key[1] != -1:
                cell._text.set_text(int(float(cell_text)))
        elif key[1] in last_columns:
             cell._text.set_color('darkslategray')
        elif (float(cell_text) <= -1.5) | (float(cell_text) >= 1.5) :
            cell._text.set_color('white')
        elif (cell_text=='nan'):
            cell._set_facecolor('lightgray')
            cell._text.set_color('lightgray')

    plt.savefig("scorecards_fallS_3LNO_FR.png", dpi=300)
    os.system('convert -trim scorecards_fallS_3LNO_FR.png scorecards_fallS_3LNO_FR.png')
    plt.close('all')



    ## Montage in subplot - English
    os.system('montage  scorecards_fallS_2H.png scorecards_fallS_2J.png scorecards_fallS_3K.png scorecards_fallS_3LNO.png -tile 1x4 -geometry +1+1  -background white  scorecards_botS_2HJ3KLNO_fall.png') 
    # French
    os.system('montage  scorecards_fallS_2H_FR.png scorecards_fallS_2J_FR.png scorecards_fallS_3K_FR.png scorecards_fallS_3LNO_FR.png -tile 1x4 -geometry +1+1  -background white  scorecards_botS_2HJ3KLNO_fall_FR.png') 
    # remove leftovers
    os.system('rm scorecards_fallS_*.png')






def sfa_bottomS_scorecards(years, clim_year=[2006, 2021]):


    '''
    Similar to sfa_bottom_scorecards, but for salinity in SFAs

    usage example (see azmp_genreport.py):
    azrt.sfa_bottomS_scorecards(years=np.arange(2006, 2022), clim_year=[2006, 2020])

    ** Note that varuables inside .pkl files are name "Tmean" isntead of "Smean"
    '''

    #### ------------- For summer ---------------- ####
    # 0.
    infile = 'operation_files/stats_sfa2_summer_salinity.pkl'
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
        elif (float(cell_text) <= -1.5) | (float(cell_text) >= 1.5) :
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
        elif (float(cell_text) <= -1.5) | (float(cell_text) >= 1.5) :
            cell._text.set_color('white')
        elif (cell_text=='nan'):
            cell._set_facecolor('lightgray')
            cell._text.set_color('lightgray')

    plt.savefig("scorecards_summerS_SFA2_FR.png", dpi=300)
    os.system('convert -trim scorecards_summerS_SFA2_FR.png scorecards_summerS_SFA2_FR.png')

 # 1.
    infile = 'operation_files/stats_sfa3_summer_salinity.pkl'
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
        elif (float(cell_text) <= -1.5) | (float(cell_text) >= 1.5) :
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
        elif (float(cell_text) <= -1.5) | (float(cell_text) >= 1.5) :
            cell._text.set_color('white')
        elif (cell_text=='nan'):
            cell._set_facecolor('lightgray')
            cell._text.set_color('lightgray')

    plt.savefig("scorecards_summerS_SFA3_FR.png", dpi=300)
    os.system('convert -trim scorecards_summerS_SFA3_FR.png scorecards_summerS_SFA3_FR.png')


    # 2.
    infile = 'operation_files/stats_sfa4_summer_salinity.pkl'
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
        elif (float(cell_text) <= -1.5) | (float(cell_text) >= 1.5) :
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
        elif (float(cell_text) <= -1.5) | (float(cell_text) >= 1.5) :
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
        elif (float(cell_text) <= -1.5) | (float(cell_text) >= 1.5) :
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
        elif (float(cell_text) <= -1.5) | (float(cell_text) >= 1.5) :
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
