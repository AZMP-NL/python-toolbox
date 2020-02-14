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
#from sys import version_info
# read/write tools
import netCDF4
import h5py
import xarray as xr
import pandas as pd
# maps
os.environ['PROJ_LIB'] = '/home/cyrf0006/anaconda3/share/proj'
from mpl_toolkits.basemap import Basemap
# interpolation tools
from scipy.interpolate import griddata
from scipy.interpolate import interp1d
# Shaping tools
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from shapely.ops import cascaded_union
## AZMP custom imports
import azmp_utils as azu
## for scorecards
import unicodedata
from matplotlib.colors import from_levels_and_colors


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
    
def bottom_temperature(season, year, zmin=0, zmax=1000, dz=5, proj='merc', netcdf_path='/home/cyrf0006/data/dev_database/netCDF/'):

    
    '''
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
    lon_reg = np.arange(lonLims[0]+dc/2, lonLims[1]-dc/2, dc)
    lat_reg = np.arange(latLims[0]+dc/2, latLims[1]-dc/2, dc)
    azu.get_bottomT_climato('/home/cyrf0006/data/dev_database/netCDF_1m_clim/*.nc', lon_reg, lat_reg, season='fall', h5_outputfile='Tbot_climato_fall_0.10.h5') 

    * see: /home/cyrf0006/AZMP/state_reports/bottomT

    usage ex:
    >> import azmp_report_tools as azrt
    >> azrt.bottom_temperature(season='spring', year='2019'):

    Frederic.Cyr@dfo-mpo.gc.ca - January 2020

    '''

    if season=='spring':
        climato_file = 'Tbot_climato_spring_0.10.h5'
    elif season=='fall':
        climato_file = 'Tbot_climato_fall_0.10.h5'
    elif season=='summer':
        climato_file = 'Tbot_climato_summer_0.10.h5'
    year_file = netcdf_path + year + '.nc'


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
    ds = xr.open_mfdataset(year_file)

    # Remome problematic datasets
    print('!!Remove MEDBA data!!')
    print('  ---> I Should be improme because I remove good data!!!!')
    ds = ds.where(ds.instrument_ID!='MEDBA', drop=True)
    ds = ds.where(ds.instrument_ID!='MEDTE', drop=True)

    # Selection of a subset region
    ds = ds.where((ds.longitude>lonLims[0]) & (ds.longitude<lonLims[1]), drop=True)
    ds = ds.where((ds.latitude>latLims[0]) & (ds.latitude<latLims[1]), drop=True)
    # Select time (save several options here)
    if season == 'summer':
        #ds = ds.sel(time=ds['time.season']=='JJA')
        ds = ds.sel(time=((ds['time.month']>=7)) & ((ds['time.month']<=9)))
    elif season == 'spring':
        #ds = ds.sel(time=ds['time.season']=='MAM')
        ds = ds.sel(time=((ds['time.month']>=4)) & ((ds['time.month']<=6)))
    elif season == 'fall':
        #ds = ds.sel(time=ds['time.season']=='SON')
        ds = ds.sel(time=((ds['time.month']>=9)) & ((ds['time.month']<=12)))
    else:
        print('!! no season specified, used them all! !!')

    # Vertical binning (on dataArray; more appropriate here)
    da_temp = ds['temperature']
    lons = np.array(ds.longitude)
    lats = np.array(ds.latitude)
    #bins = np.arange(dz/2.0, ds.level.max(), dz)
    bins = np.arange(dz/2.0, 1000, dz)
    da_temp = da_temp.groupby_bins('level', bins).mean(dim='level')
    #To Pandas Dataframe
    df_temp = da_temp.to_pandas()
    df_temp.columns = bins[0:-1] #rename columns with 'bins'
    # Remove empty columns
    idx_empty_rows = df_temp.isnull().all(1).nonzero()[0]
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
                idx_interp = np.arange(np.int(idx_good[0]),np.int(idx_good[-1]+1))
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

    # Mask data outside Nafo div.
    print('Mask according to NAFO division for ' + season)
    # Polygons
    polygon4R = Polygon(zip(nafo_div['4R']['lon'], nafo_div['4R']['lat']))
    polygon3K = Polygon(zip(nafo_div['3K']['lon'], nafo_div['3K']['lat']))
    polygon3L = Polygon(zip(nafo_div['3L']['lon'], nafo_div['3L']['lat']))
    polygon3N = Polygon(zip(nafo_div['3N']['lon'], nafo_div['3N']['lat']))
    polygon3O = Polygon(zip(nafo_div['3O']['lon'], nafo_div['3O']['lat']))
    polygon3Ps = Polygon(zip(nafo_div['3Ps']['lon'], nafo_div['3Ps']['lat']))
    polygon2J = Polygon(zip(nafo_div['2J']['lon'], nafo_div['2J']['lat']))

    # Contour of data to mask
    contour_mask = np.load('100m_contour_labrador.npy')
    polygon_mask = Polygon(contour_mask)

    if season == 'spring':
        for i, xx in enumerate(lon_reg):
            for j,yy in enumerate(lat_reg):
                point = Point(lon_reg[i], lat_reg[j])
                #if (~polygon3L.contains(point)) & (~polygon3N.contains(point)) & (~polygon3O.contains(point)) & (~polygon3Ps.contains(point)):
                if polygon3L.contains(point) | polygon3N.contains(point) | polygon3O.contains(point) | polygon3Ps.contains(point):
                    pass #nothing to do but cannot implement negative statement "if not" above
                else:
                    Tbot[j,i] = np.nan

    elif season == 'fall':
        for i, xx in enumerate(lon_reg):
            for j,yy in enumerate(lat_reg):
                point = Point(lon_reg[i], lat_reg[j])
                if polygon2J.contains(point) | polygon3K.contains(point) | polygon3L.contains(point) | polygon3N.contains(point) | polygon3O.contains(point):
                    pass #nothing to do but cannot implement negative statement "if not" above
                else:
                    Tbot[j,i] = np.nan ### <--------------------- Do mask the fall / OR / 
                    #Tbot[j,i] = np.nan ### <--------------------- Do not mask the fall!!!!!

                if polygon_mask.contains(point): # mask data near Labrador in fall
                    Tbot[j,i] = np.nan 

    elif season == 'summer':
        for i, xx in enumerate(lon_reg):
            for j,yy in enumerate(lat_reg):
                point = Point(lon_reg[i], lat_reg[j])
                # Just mask labrador
                if polygon_mask.contains(point): # mask data near Labrador in fall
                    Tbot[j,i] = np.nan 

    else:
        print('no division mask, all data taken')

    print(' -> Done!')    

    # Temperature anomaly:
    anom = Tbot-Tbot_climato

    ## ---- Plot Anomaly ---- ##
    fig, ax = plt.subplots(nrows=1, ncols=1)
    m = Basemap(ax=ax, projection='merc',lon_0=lon_0,lat_0=lat_0, llcrnrlon=lonLims[0],llcrnrlat=latLims[0],urcrnrlon=lonLims[1],urcrnrlat=latLims[1], resolution= 'i')
    levels = np.linspace(-3.5, 3.5, 8)
    #levels = np.linspace(-3.5, 3.5, 16)
    xi, yi = m(*np.meshgrid(lon_reg, lat_reg))
    c = m.contourf(xi, yi, anom, levels, cmap=plt.cm.RdBu_r, extend='both')
    cc = m.contour(xi, yi, -Zitp, [100, 500, 1000, 4000], colors='grey');
    plt.clabel(cc, inline=1, fontsize=10, fmt='%d')
    if season=='fall':
        plt.title('Fall Bottom Temperature ' + year + ' Anomaly')
    elif season=='spring':
        plt.title('Spring Bottom Temperature ' + year + ' Anomaly')
    else:
        plt.title('Bottom Temperature ' + year + '  Anomaly')
    m.fillcontinents(color='tan');
    m.drawparallels([40, 45, 50, 55, 60], labels=[0,0,0,0], fontsize=12, fontweight='normal');
    m.drawmeridians([-60, -55, -50, -45], labels=[0,0,0,1], fontsize=12, fontweight='normal');
    cax = fig.add_axes([0.16, 0.05, 0.7, 0.025])
    cb = plt.colorbar(c, cax=cax, orientation='horizontal')
    cb.set_label(r'$\rm T(^{\circ}C)$', fontsize=12, fontweight='normal')
    div_toplot = ['2J', '3K', '3L', '3N', '3O', '3Ps', '4R']
    for div in div_toplot:
        div_lon, div_lat = m(nafo_div[div]['lon'], nafo_div[div]['lat'])
        m.plot(div_lon, div_lat, 'k', linewidth=2)
        ax.text(np.mean(div_lon), np.mean(div_lat), div, fontsize=12, color='black', fontweight='bold')    
    # Save Figure
    fig.set_size_inches(w=6, h=9)
    fig.set_dpi(300)
    outfile = 'bottom_temp_anomaly_' + season + '_' + year + '.png'
    fig.savefig(outfile)
    os.system('convert -trim ' + outfile + ' ' + outfile)
    # Save French Figure
    plt.sca(ax)
    if season=='fall':
        plt.title(u'Anomalie de température au fond - Automne ' + year )
    elif season=='spring':
        plt.title(u'Anomalie de température au fond - Printemp ' + year )
    else:
        plt.title(u'Anomalie de température au fond ' + year )
    fig.set_size_inches(w=6, h=9)
    fig.set_dpi(300)
    outfile = 'bottom_temp_anomaly_' + season + '_' + year + '_FR.png'
    fig.savefig(outfile)
    os.system('convert -trim ' + outfile + ' ' + outfile)

    ## ---- Plot Temperature ---- ##
    fig, ax = plt.subplots(nrows=1, ncols=1)
    m = Basemap(ax=ax, projection='merc',lon_0=lon_0,lat_0=lat_0, llcrnrlon=lonLims[0],llcrnrlat=latLims[0],urcrnrlon=lonLims[1],urcrnrlat=latLims[1], resolution= 'i')
    #levels = np.linspace(-2, 6, 9)
    levels = np.linspace(-2, 6, 17)
    xi, yi = m(*np.meshgrid(lon_reg, lat_reg))
    c = m.contourf(xi, yi, Tbot, levels, cmap=plt.cm.RdBu_r, extend='both')
    cc = m.contour(xi, yi, -Zitp, [100, 500, 1000, 4000], colors='grey');
    plt.clabel(cc, inline=1, fontsize=10, fmt='%d')
    if season=='fall':
        plt.title('Fall Bottom Temperature ' + year)
    elif season=='spring':
        plt.title('Spring Bottom Temperature ' + year)
    else:
        plt.title('Bottom Temperature ' + year)
    m.fillcontinents(color='tan');
    m.drawparallels([40, 45, 50, 55, 60], labels=[0,0,0,0], fontsize=12, fontweight='normal');
    m.drawmeridians([-60, -55, -50, -45], labels=[0,0,0,1], fontsize=12, fontweight='normal');
    x, y = m(lons, lats)
    m.scatter(x,y, s=50, marker='.',color='k')
    cax = fig.add_axes([0.16, 0.05, 0.7, 0.025])
    #cax = plt.axes([0.85,0.15,0.04,0.7], facecolor='grey')
    cb = plt.colorbar(c, cax=cax, orientation='horizontal')
    cb.set_label(r'$\rm T(^{\circ}C)$', fontsize=12, fontweight='normal')
    div_toplot = ['2J', '3K', '3L', '3N', '3O', '3Ps', '4R']
    for div in div_toplot:
        div_lon, div_lat = m(nafo_div[div]['lon'], nafo_div[div]['lat'])
        m.plot(div_lon, div_lat, 'k', linewidth=2)
        ax.text(np.mean(div_lon), np.mean(div_lat), div, fontsize=12, color='black', fontweight='bold')
    # Save Figure
    fig.set_size_inches(w=6, h=9)
    fig.set_dpi(200)
    outfile = 'bottom_temp_' + season + '_' + year + '.png'
    fig.savefig(outfile)
    os.system('convert -trim ' + outfile + ' ' + outfile)
    # Save French Figure
    plt.sca(ax)
    if season=='fall':
        plt.title(u'Température au fond - Automne ' + year )
    elif season=='spring':
        plt.title(u'Température au fond - Printemp ' + year )
    else:
        plt.title(u'Température au fond ' + year )
    fig.set_size_inches(w=6, h=9)
    fig.set_dpi(300)
    outfile = 'bottom_temp_' + season + '_' + year + '_FR.png'
    fig.savefig(outfile)
    os.system('convert -trim ' + outfile + ' ' + outfile)


    ## ---- Plot Climato ---- ##
    fig, ax = plt.subplots(nrows=1, ncols=1)
    m = Basemap(ax=ax, projection='merc',lon_0=lon_0,lat_0=lat_0, llcrnrlon=lonLims[0],llcrnrlat=latLims[0],urcrnrlon=lonLims[1],urcrnrlat=latLims[1], resolution= 'i')
    #levels = np.linspace(-2, 6, 9)
    levels = np.linspace(-2, 6, 17)
    xi, yi = m(*np.meshgrid(lon_reg, lat_reg))
    c = m.contourf(xi, yi, Tbot_climato, levels, cmap=plt.cm.RdBu_r, extend='both')
    cc = m.contour(xi, yi, -Zitp, [100, 500, 1000, 4000], colors='grey');
    plt.clabel(cc, inline=1, fontsize=10, fmt='%d')
    if season=='fall':
        plt.title('Fall Bottom Temperature Climatology')
    elif season=='spring':
        plt.title('Spring Bottom Temperature Climatology')
    else:
        plt.title('Bottom Temperature Climatology')
    m.fillcontinents(color='tan');
    m.drawparallels([40, 45, 50, 55, 60], labels=[1,0,0,0], fontsize=12, fontweight='normal');
    m.drawmeridians([-60, -55, -50, -45], labels=[0,0,0,1], fontsize=12, fontweight='normal');
    cax = fig.add_axes([0.16, 0.05, 0.7, 0.025])
    cb = plt.colorbar(c, cax=cax, orientation='horizontal')
    cb.set_label(r'$\rm T(^{\circ}C)$', fontsize=12, fontweight='normal')
    div_toplot = ['2J', '3K', '3L', '3N', '3O', '3Ps', '4R']
    for div in div_toplot:
        div_lon, div_lat = m(nafo_div[div]['lon'], nafo_div[div]['lat'])
        m.plot(div_lon, div_lat, 'k', linewidth=2)
        ax.text(np.mean(div_lon), np.mean(div_lat), div, fontsize=12, color='black', fontweight='bold')    
    # Save Figure
    fig.set_size_inches(w=6, h=9)
    fig.set_dpi(300)
    outfile = 'bottom_temp_climato_' + season + '_' + year + '.png'
    fig.savefig(outfile)
    os.system('convert -trim ' + outfile + ' ' + outfile)
    # Save French Figure
    plt.sca(ax)
    if season=='fall':
        plt.title(u'Climatologie de température au fond - Automne ' + year )
    elif season=='spring':
        plt.title(u'Climatologie de température au fond - Printemp ' + year )
    else:
        plt.title(u'Climatologie de température au fond ' + year )
    fig.set_size_inches(w=6, h=9)
    fig.set_dpi(300)
    outfile = 'bottom_temp_climato_' + season + '_' + year + '_FR.png'
    fig.savefig(outfile)
    os.system('convert -trim ' + outfile + ' ' + outfile)


    # Convert to a subplot
    os.system('montage bottom_temp_climato_' + season + '_' + year + '.png bottom_temp_' + season + '_' + year + '.png bottom_temp_anomaly_' + season + '_' + year + '.png  -tile 3x1 -geometry +10+10  -background white  bottomT_' + season + year + '.png') 
    # in French
    os.system('montage bottom_temp_climato_' + season + '_' + year + '_FR.png bottom_temp_' + season + '_' + year + '_FR.png bottom_temp_anomaly_' + season + '_' + year + '_FR.png  -tile 3x1 -geometry +10+10  -background white  bottomT_' + season + year + '_FR.png') 
    # Move to year folder
    os.system('mv bottomT_' + season + year + '.png bottomT_' + season + year + '_FR.png ../' + year)


#### bottom_salinity
def bottom_salinity(season, year, zmin=0, zmax=1000, dz=5, proj='merc', netcdf_path='/home/cyrf0006/data/dev_database/netCDF/'):

    
    '''
    Bottom salinity maps for AZMP ResDoc

    This function is adapted from azmp_bottomS.py

    * process in : /home/cyrf0006/AZMP/state_reports/bottomT

    usage ex:
    >> import azmp_report_tools as azrt
    >> azrt.bottom_salinity(season='spring', year='2019'):

    Frederic.Cyr@dfo-mpo.gc.ca - January 2020

    '''

    if season=='spring':
        climato_file = 'Sbot_climato_spring_0.10.h5'
    elif season=='fall':
        climato_file = 'Sbot_climato_fall_0.10.h5'

    year_file = netcdf_path + year + '.nc'

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
    # Remome problematic datasets
    print('!!Remove MEDBA data!!')
    print('  ---> I Should be improme because I remove good data!!!!')
    ds = ds.where(ds.instrument_ID!='MEDBA', drop=True)
    ds = ds.where(ds.instrument_ID!='MEDTE', drop=True)
    # Selection of a subset region
    ds = ds.where((ds.longitude>lonLims[0]) & (ds.longitude<lonLims[1]), drop=True)
    ds = ds.where((ds.latitude>latLims[0]) & (ds.latitude<latLims[1]), drop=True)
    # Select time (save several options here)
    if season == 'summer':
        ds = ds.sel(time=((ds['time.month']>=7)) & ((ds['time.month']<=9)))
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
    # Vertical binning (on dataArray; more appropriate here
    da_sal = ds['salinity']
    bins = np.arange(dz/2.0, ds.level.max(), dz)
    da_sal = da_sal.groupby_bins('level', bins).mean(dim='level')
    #To Pandas Dataframe
    df_sal = da_sal.to_pandas()
    df_sal.columns = bins[0:-1] #rename columns with 'bins'
    # Remove empty columns & drop coordinates (for cast identification on map)
    idx_empty_rows = df_sal.isnull().all(1).nonzero()[0]
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
                idx_interp = np.arange(np.int(idx_good[0]),np.int(idx_good[-1]+1))
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

    # Mask data outside Nafo div.
    print('Mask according to NAFO division for ' + season)
    # Polygons
    polygon3K = Polygon(zip(nafo_div['3K']['lon'], nafo_div['3K']['lat']))
    polygon3L = Polygon(zip(nafo_div['3L']['lon'], nafo_div['3L']['lat']))
    polygon3N = Polygon(zip(nafo_div['3N']['lon'], nafo_div['3N']['lat']))
    polygon3O = Polygon(zip(nafo_div['3O']['lon'], nafo_div['3O']['lat']))
    polygon3Ps = Polygon(zip(nafo_div['3Ps']['lon'], nafo_div['3Ps']['lat']))
    polygon2J = Polygon(zip(nafo_div['2J']['lon'], nafo_div['2J']['lat']))

    # Contour of data to mask
    contour_mask = np.load('100m_contour_labrador.npy')
    polygon_mask = Polygon(contour_mask)

    if season == 'spring':
        for i, xx in enumerate(lon_reg):
            for j,yy in enumerate(lat_reg):
                point = Point(lon_reg[i], lat_reg[j])
                if polygon3L.contains(point) | polygon3N.contains(point) | polygon3O.contains(point) | polygon3Ps.contains(point):
                    pass #nothing to do but cannot implement negative statement "if not" above
                else:
                    Sbot[j,i] = np.nan

    elif season == 'fall':
        for i, xx in enumerate(lon_reg):
            for j,yy in enumerate(lat_reg):
                point = Point(lon_reg[i], lat_reg[j])
                if polygon2J.contains(point) | polygon3K.contains(point) | polygon3L.contains(point) | polygon3N.contains(point) | polygon3O.contains(point):
                    pass #nothing to do but cannot implement negative statement "if not" above
                else:
                    Sbot[j,i] = np.nan ### <--------------------- Do not mask the fall!!!!!

                if polygon_mask.contains(point): # mask data near Labrador in fall
                    Sbot[j,i] = np.nan 
    else:
        print('no division mask, all data taken')

    print(' -> Done!')    

    # Salinity anomaly:
    anom = Sbot-Sbot_climato

    ## ---- Plot Anomaly ---- ##
    fig, ax = plt.subplots(nrows=1, ncols=1)
    m = Basemap(ax=ax, projection='merc',lon_0=lon_0,lat_0=lat_0, llcrnrlon=lonLims[0],llcrnrlat=latLims[0],urcrnrlon=lonLims[1],urcrnrlat=latLims[1], resolution= 'i')
    #levels = np.linspace(-1, 1, 6)
    levels = np.array([-1, -.8, -.6, -.4, -.2, .2, .4, .6, .8, 1])
    xi, yi = m(*np.meshgrid(lon_reg, lat_reg))
    c = m.contourf(xi, yi, anom, levels, cmap=plt.cm.RdBu_r, extend='both')
    cc = m.contour(xi, yi, -Zitp, [100, 500, 1000, 4000], colors='grey');
    plt.clabel(cc, inline=1, fontsize=10, fmt='%d')
    if season=='fall':
        plt.title('Fall Bottom Salinity ' + year + ' Anomaly')
    elif season=='spring':
        plt.title('Spring Bottom Salinity ' + year + ' Anomaly')
    else:
        plt.title('Bottom Salinity ' + year + '  Anomaly')
    m.fillcontinents(color='tan');
    m.drawparallels([40, 45, 50, 55, 60], labels=[0,0,0,0], fontsize=12, fontweight='normal');
    m.drawmeridians([-60, -55, -50, -45], labels=[0,0,0,1], fontsize=12, fontweight='normal');
    cax = fig.add_axes([0.16, 0.05, 0.7, 0.025])
    cb = plt.colorbar(c, cax=cax, orientation='horizontal')
    cb.set_label(r'$\rm S$', fontsize=12, fontweight='normal')
    div_toplot = ['2J', '3K', '3L', '3N', '3O', '3Ps']
    for div in div_toplot:
        div_lon, div_lat = m(nafo_div[div]['lon'], nafo_div[div]['lat'])
        m.plot(div_lon, div_lat, 'k', linewidth=2)
        ax.text(np.mean(div_lon), np.mean(div_lat), div, fontsize=12, color='black', fontweight='bold')    
    # Save Figure
    fig.set_size_inches(w=6, h=9)
    fig.set_dpi(200)
    outfile = 'bottom_sal_anomaly_' + season + '_' + year + '.png'
    fig.savefig(outfile)
    os.system('convert -trim ' + outfile + ' ' + outfile)
    # Save French Figure
    plt.sca(ax)
    if season=='fall':
        plt.title(u'Anomalie de salinité au fond - Automne ' + year )
    elif season=='spring':
        plt.title(u'Anomalie de salinité au fond - Printemp ' + year )
    else:
        plt.title(u'Anomalie de salinité au fond ' + year )
    fig.set_size_inches(w=6, h=9)
    fig.set_dpi(300)
    outfile = 'bottom_sal_anomaly_' + season + '_' + year + '_FR.png'
    fig.savefig(outfile)
    os.system('convert -trim ' + outfile + ' ' + outfile)

    ## ---- Plot Salinity ---- ##
    fig, ax = plt.subplots(nrows=1, ncols=1)
    m = Basemap(ax=ax, projection='merc',lon_0=lon_0,lat_0=lat_0, llcrnrlon=lonLims[0],llcrnrlat=latLims[0],urcrnrlon=lonLims[1],urcrnrlat=latLims[1], resolution= 'i')
    levels = np.linspace(30, 36, 13)
    xi, yi = m(*np.meshgrid(lon_reg, lat_reg))
    c = m.contourf(xi, yi, Sbot, levels, cmap=plt.cm.RdBu_r, extend='both')
    cc = m.contour(xi, yi, -Zitp, [100, 500, 1000, 4000], colors='grey');
    plt.clabel(cc, inline=1, fontsize=10, fmt='%d')
    if season=='fall':
        plt.title('Fall Bottom Salinity ' + year)
    elif season=='spring':
        plt.title('Spring Bottom Salinity ' + year)
    else:
        plt.title('Bottom Salinity ' + year)
    m.fillcontinents(color='tan');
    m.drawparallels([40, 45, 50, 55, 60], labels=[0,0,0,0], fontsize=12, fontweight='normal');
    m.drawmeridians([-60, -55, -50, -45], labels=[0,0,0,1], fontsize=12, fontweight='normal');
    x, y = m(lons, lats)
    m.scatter(x,y, s=50, marker='.',color='k')
    cax = fig.add_axes([0.16, 0.05, 0.7, 0.025])
    cb = plt.colorbar(c, cax=cax, orientation='horizontal')
    cb.set_label(r'$\rm S$', fontsize=12, fontweight='normal')
    div_toplot = ['2J', '3K', '3L', '3N', '3O', '3Ps']
    for div in div_toplot:
        div_lon, div_lat = m(nafo_div[div]['lon'], nafo_div[div]['lat'])
        m.plot(div_lon, div_lat, 'k', linewidth=2)
        ax.text(np.mean(div_lon), np.mean(div_lat), div, fontsize=12, color='black', fontweight='bold')
    # Save Figure
    fig.set_size_inches(w=6, h=9)
    fig.set_dpi(200)
    outfile = 'bottom_sal_' + season + '_' + year + '.png'
    fig.savefig(outfile)
    os.system('convert -trim ' + outfile + ' ' + outfile)
    # Save French Figure
    plt.sca(ax)
    if season=='fall':
        plt.title(u'Salinité au fond - Automne ' + year )
    elif season=='spring':
        plt.title(u'Salinité au fond - Printemp ' + year )
    else:
        plt.title(u'Salinité au fond ' + year )
    fig.set_size_inches(w=6, h=9)
    fig.set_dpi(300)
    outfile = 'bottom_sal_' + season + '_' + year + '_FR.png'
    fig.savefig(outfile)
    os.system('convert -trim ' + outfile + ' ' + outfile)

    ## ---- Plot Climato ---- ##
    fig, ax = plt.subplots(nrows=1, ncols=1)
    m = Basemap(ax=ax, projection='merc',lon_0=lon_0,lat_0=lat_0, llcrnrlon=lonLims[0],llcrnrlat=latLims[0],urcrnrlon=lonLims[1],urcrnrlat=latLims[1], resolution= 'i')
    levels = np.linspace(30, 36, 13)
    #levels = np.linspace(30, 36, 7)
    xi, yi = m(*np.meshgrid(lon_reg, lat_reg))
    c = m.contourf(xi, yi, Sbot_climato, levels, cmap=plt.cm.RdBu_r, extend='both')
    cc = m.contour(xi, yi, -Zitp, [100, 500, 1000, 4000], colors='grey');
    plt.clabel(cc, inline=1, fontsize=10, fmt='%d')
    if season=='fall':
        plt.title('Fall Bottom Salinity Climatology')
    elif season=='spring':
        plt.title('Spring Bottom Salinity Climatology')
    else:
        plt.title('Bottom Salinity Climatology')
    m.fillcontinents(color='tan');
    m.drawparallels([40, 45, 50, 55, 60], labels=[1,0,0,0], fontsize=12, fontweight='normal');
    m.drawmeridians([-60, -55, -50, -45], labels=[0,0,0,1], fontsize=12, fontweight='normal');
    cax = fig.add_axes([0.16, 0.05, 0.7, 0.025])
    cb = plt.colorbar(c, cax=cax, orientation='horizontal')
    cb.set_label(r'$\rm S$', fontsize=12, fontweight='normal')
    div_toplot = ['2J', '3K', '3L', '3N', '3O', '3Ps']
    for div in div_toplot:
        div_lon, div_lat = m(nafo_div[div]['lon'], nafo_div[div]['lat'])
        m.plot(div_lon, div_lat, 'k', linewidth=2)
        ax.text(np.mean(div_lon), np.mean(div_lat), div, fontsize=12, color='black', fontweight='bold')    
    # Save Figure
    fig.set_size_inches(w=6, h=9)
    fig.set_dpi(200)
    outfile = 'bottom_sal_climato_' + season + '_' + year + '.png'
    fig.savefig(outfile)
    os.system('convert -trim ' + outfile + ' ' + outfile)
    # Save French Figure
    plt.sca(ax)
    if season=='fall':
        plt.title(u'Climatoligie de salinité au fond - Automne ' + year )
    elif season=='spring':
        plt.title(u'Climatologie de salinité au fond - Printemp ' + year )
    else:
        plt.title(u'Climatologie de salinité au fond ' + year )
    fig.set_size_inches(w=6, h=9)
    fig.set_dpi(300)
    outfile = 'bottom_sal_climato_' + season + '_' + year + '_FR.png'
    fig.savefig(outfile)
    os.system('convert -trim ' + outfile + ' ' + outfile)

    # Convert to a subplot
    os.system('montage bottom_sal_climato_' + season + '_' + year + '.png bottom_sal_' + season + '_' + year + '.png bottom_sal_anomaly_' + season + '_' + year + '.png  -tile 3x1 -geometry +10+10  -background white  bottomS_' + season + year + '.png') 
    # French
    os.system('montage bottom_sal_climato_' + season + '_' + year + '_FR.png bottom_sal_' + season + '_' + year + '_FR.png bottom_sal_anomaly_' + season + '_' + year + '_FR.png  -tile 3x1 -geometry +10+10  -background white  bottomS_' + season + year + '_FR.png') 
    # Move to year folder
    os.system('mv bottomS_' + season + year + '.png bottomS_' + season + year + '_FR.png ../' + year)


#### bottom_stats
def bottom_stats(years, season, proj='merc', plot=False, SPECIAL=None, netcdf_path='/home/cyrf0006/data/dev_database/netCDF/'):

    '''
        Function bottom_stats() based on script azmp_bottom_stats.py

        See the latter on how to generate the bottom climatology.
        See it also for specific usage such as Plaice - COSEWIC analysis.

        usage example:
        >> import azmp_report_tools as azrt
        >> import numpy as np
        >> azrt.bottom_stats(years=np.arange(1980, 2020), season='fall')

        Frederic.Cyr@dfo-mpo.gc.ca - January 2020
    '''

    # load climato
    if season == 'fall':
        climato_file = 'Tbot_climato_fall_0.10.h5'
    elif season == 'spring':
        climato_file = 'Tbot_climato_spring_0.10.h5'
    elif season == 'summer':
        if SPECIAL=='COSEWIC_north':
            climato_file = 'Tbot_climato_2GH_summer_0.10.h5'
        elif SPECIAL=='COSEWIC_south':
            climato_file = 'Tbot_climato_SA45_summer_0.10.h5'
        else:
            climato_file = 'Tbot_climato_summer_0.10.h5'
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
    shape_3LNO = cascaded_union(shape)
    shape_3M = Polygon(zip(nafo_div['3M']['lon'], nafo_div['3M']['lat']))
    shape_3Ps = Polygon(zip(nafo_div['3Ps']['lon'], nafo_div['3Ps']['lat']))
    shape_2G = Polygon(zip(nafo_div['2G']['lon'], nafo_div['2G']['lat']))
    shape_2H = Polygon(zip(nafo_div['2H']['lon'], nafo_div['2H']['lat']))
    shape_2J = Polygon(zip(nafo_div['2J']['lon'], nafo_div['2J']['lat']))
    shape_3K = Polygon(zip(nafo_div['3K']['lon'], nafo_div['3K']['lat']))
    shape_3L = Polygon(zip(nafo_div['3L']['lon'], nafo_div['3L']['lat']))
    shape_3O = Polygon(zip(nafo_div['3O']['lon'], nafo_div['3O']['lat']))
    shape = [shape_2J, shape_2H]
    shape_2HJ = cascaded_union(shape)
    shape = [shape_2G, shape_2H]
    shape_2GH = cascaded_union(shape)
    shape_4R = Polygon(zip(nafo_div['4R']['lon'], nafo_div['4R']['lat']))
    shape_4S = Polygon(zip(nafo_div['4S']['lon'], nafo_div['4S']['lat']))
    shape_4T = Polygon(zip(nafo_div['4T']['lon'], nafo_div['4T']['lat']))
    shape = [shape_4R, shape_4S]
    shape_4RS = cascaded_union(shape)
    shape = [shape_4R, shape_4S, shape_4T]
    shape_4RST = cascaded_union(shape)
    polygon4Vs = Polygon(zip(nafo_div['4Vs']['lon'], nafo_div['4Vs']['lat']))
    polygon4Vn = Polygon(zip(nafo_div['4Vn']['lon'], nafo_div['4Vn']['lat']))
    polygon4W = Polygon(zip(nafo_div['4W']['lon'], nafo_div['4W']['lat']))
    polygon4X = Polygon(zip(nafo_div['4X']['lon'], nafo_div['4X']['lat']))
    shape = [polygon4Vs.buffer(0), polygon4Vn.buffer(0), polygon4W.buffer(0), polygon4X.buffer(0)]
    shape_4VWX = cascaded_union(shape)
    shape_5Y = Polygon(zip(nafo_div['5Y']['lon'], nafo_div['5Y']['lat']))

    dict_stats_3LNO = {}
    dict_stats_3M = {}
    dict_stats_3Ps = {}
    dict_stats_3K = {}
    dict_stats_3L = {}
    dict_stats_3O = {}
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
        print(' ---- ' + np.str(year) + ' ---- ')
        year_file = netcdf_path + np.str(year) + '.nc'
        Tdict = azu.get_bottomT(year_file, season, climato_file)    
        Tbot = Tdict['Tbot']
        lons = Tdict['lons']
        lats = Tdict['lats']
        anom = Tbot-Tbot_climato


        # NAFO division stats    
        dict_stats_2GH[np.str(year)] = azu.polygon_temperature_stats(Tdict, shape_2GH)
        dict_stats_2J[np.str(year)] = azu.polygon_temperature_stats(Tdict, shape_2J)
        dict_stats_2HJ[np.str(year)] = azu.polygon_temperature_stats(Tdict, shape_2HJ)
        dict_stats_3LNO[np.str(year)] = azu.polygon_temperature_stats(Tdict, shape_3LNO)
        dict_stats_3M[np.str(year)] = azu.polygon_temperature_stats(Tdict, shape_3M)
        dict_stats_3Ps[np.str(year)] = azu.polygon_temperature_stats(Tdict, shape_3Ps)
        dict_stats_3K[np.str(year)] = azu.polygon_temperature_stats(Tdict, shape_3K)
        dict_stats_3L[np.str(year)] = azu.polygon_temperature_stats(Tdict, shape_3L)
        dict_stats_3O[np.str(year)] = azu.polygon_temperature_stats(Tdict, shape_3O)
        dict_stats_4R[np.str(year)] = azu.polygon_temperature_stats(Tdict, shape_4R)
        dict_stats_4S[np.str(year)] = azu.polygon_temperature_stats(Tdict, shape_4S)
        dict_stats_4RS[np.str(year)] = azu.polygon_temperature_stats(Tdict, shape_4RS)
        dict_stats_4RST[np.str(year)] = azu.polygon_temperature_stats(Tdict, shape_4RST)
        #dict_stats_4T[np.str(year)] = azu.polygon_temperature_stats(Tdict, shape_4T)
        dict_stats_4VWX[np.str(year)] = azu.polygon_temperature_stats(Tdict, shape_4VWX)
        #dict_stats_5Y[np.str(year)] = azu.polygon_temperature_stats(Tdict, shape_5Y)

        # Append bottom temperature for multi-index export
        df = pd.DataFrame(index=lat_reg, columns=lon_reg)
        df.index.name='latitude'
        df.columns.name='longitude'
        df[:] = Tbot
        df_list.append(df)


        if plot:
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
            cax = plt.axes([0.85,0.15,0.04,0.7], facecolor='grey')
            cb = plt.colorbar(c, cax=cax)
            cb.set_label(r'$\rm T(^{\circ}C)$', fontsize=12, fontweight='normal')
            div_toplot = ['2J', '3K', '3L', '3N', '3O', '3Ps', '4R', '4Vn', '4Vs', '4W', '4X']
            for div in div_toplot:
                div_lon, div_lat = m(nafo_div[div]['lon'], nafo_div[div]['lat'])
                m.plot(div_lon, div_lat, 'k', linewidth=2)
                ax.text(np.mean(div_lon), np.mean(div_lat), div, fontsize=12, color='black', fontweight='bold')    
            # Save Figure
            fig.set_size_inches(w=7, h=8)
            fig.set_dpi(200)
            outfile = 'bottom_temp_anomaly_' + season + '_' + np.str(year) + '.png'
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
            cax = plt.axes([0.85,0.15,0.04,0.7], facecolor='grey')
            cb = plt.colorbar(c, cax=cax)
            cb.set_label(r'$\rm T(^{\circ}C)$', fontsize=12, fontweight='normal')
            div_toplot = ['2J', '3K', '3L', '3N', '3O', '3Ps', '4R', '4Vn', '4Vs', '4W', '4X']
            for div in div_toplot:
                div_lon, div_lat = m(nafo_div[div]['lon'], nafo_div[div]['lat'])
                m.plot(div_lon, div_lat, 'k', linewidth=2)
                ax.text(np.mean(div_lon), np.mean(div_lat), div, fontsize=12, color='black', fontweight='bold')
            # Save Figure
            fig.set_size_inches(w=7, h=8)
            fig.set_dpi(200)
            outfile = 'bottom_temp_' + season + '_' + np.str(year) + '.png'
            fig.savefig(outfile)

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
    outname = 'stats_2J_' + season + '.pkl'
    df_2J.to_pickle(outname)
    outname = 'stats_2HJ_' + season + '.pkl'
    df_2HJ.to_pickle(outname)
    outname = 'stats_2GH_' + season + '.pkl'
    df_2GH.to_pickle(outname)
    outname = 'stats_4R_' + season + '.pkl'
    df_4R.to_pickle(outname)
    outname = 'stats_4S_' + season + '.pkl'
    df_4S.to_pickle(outname)
    outname = 'stats_4RS_' + season + '.pkl'
    df_4RS.to_pickle(outname)
    outname = 'stats_4RST_' + season + '.pkl'
    df_4RST.to_pickle(outname)
    outname = 'stats_4T_' + season + '.pkl'
    df_4T.to_pickle(outname)
    outname = 'stats_4VWX_' + season + '.pkl'
    df_4VWX.to_pickle(outname)
    #outname = 'stats_5Y_' + season + '.pkl'
    #df_5Y.to_pickle(outname)

    # Save in multi-index  dataFrame
    year_index = pd.Series(years)
    year_index.name='year'
    df_mindex = pd.concat(df_list,keys=year_index)
    df_mindex.to_pickle(season + '_bottom_temperature.pkl')


def bottom_scorecards(years, clim_year=[1981, 2010]):


    '''
    To generate AZMP score cards for bottom temperature
    Uses pickled object generated by bottom_stats()

    usage example:


    bottom_scorecards(years=[1980, 2019], clim_year=[1981, 2010]):

    '''

    #### ------------- For fall ---------------- ####
    # 1.
    infile = 'stats_2J_fall.pkl'
    df = pd.read_pickle(infile)
    df.index = pd.to_datetime(df.index) # update index to datetime
    df = df[(df.index.year>=years[0]) & (df.index.year<=years[-1])]
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
                          colLabels=['-- NAFO division 2J --'],
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
    table_cells = table_props['child_artists']
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
    plt.savefig("scorecards_fall_2J.png", dpi=300)
    os.system('convert -trim scorecards_fall_2J.png scorecards_fall_2J.png')

    # French table
    std_anom = std_anom.rename({r'$\rm T_{bot}$' : r'$\rm T_{fond}$', r'$\rm T_{bot_{<200m}}$' : r'$\rm T_{fond_{<200m}}$', r'$\rm Area_{>2^{\circ}C}$' : r'$\rm Aire_{>2^{\circ}C}$', r'$\rm Area_{<1^{\circ}C}$' : r'$\rm Aire_{<1^{\circ}C}$'})
    year_list[-1] = u'ET'

    header = ax.table(cellText=[['']],
                          colLabels=['-- Division 2J de l\'OPANO --'],
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
    table_cells = table_props['child_artists']
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
    plt.savefig("scorecards_fall_2J_FR.png", dpi=300)
    os.system('convert -trim scorecards_fall_2J_FR.png scorecards_fall_2J_FR.png')

    # 2.
    infile = 'stats_3K_fall.pkl'
    df = pd.read_pickle(infile)
    df.index = pd.to_datetime(df.index) # update index to datetime
    df = df[(df.index.year>=years[0]) & (df.index.year<=years[-1])]
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
    table_cells = table_props['child_artists']
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
    table_cells = table_props['child_artists']
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

    plt.savefig("scorecards_fall_3K_FR.png", dpi=300)
    os.system('convert -trim scorecards_fall_3K_FR.png scorecards_fall_3K_FR.png')

    # 3.
    infile = 'stats_3LNO_fall.pkl'
    df = pd.read_pickle(infile)
    df.index = pd.to_datetime(df.index) # update index to datetime
    df = df[(df.index.year>=years[0]) & (df.index.year<=years[-1])]
    df['area_colder0'] = df['area_colder0']/1000 # In 1000km
    df['area_colder1'] = df['area_colder1']/1000 # In 1000km
    df['area_warmer2'] = df['area_warmer2']/1000
    df_clim = df[(df.index.year>=clim_year[0]) & (df.index.year<=clim_year[1])]
    std_anom = (df-df_clim.mean(axis=0))/df_clim.std(axis=0)
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
    table_cells = table_props['child_artists']
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
    table_cells = table_props['child_artists']
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

    plt.savefig("scorecards_fall_3LNO_FR.png", dpi=300)
    os.system('convert -trim scorecards_fall_3LNO_FR.png scorecards_fall_3LNO_FR.png')


    # English
    os.system('montage  scorecards_fall_2J.png scorecards_fall_3K.png scorecards_fall_3LNO.png -tile 1x3 -geometry +1+1  -background white  scorecards_botT_fall.png') 
    # French
    os.system('montage  scorecards_fall_2J_FR.png scorecards_fall_3K_FR.png scorecards_fall_3LNO_FR.png -tile 1x3 -geometry +1+1  -background white  scorecards_botT_fall_FR.png') 




    #### ------------- For Spring ---------------- ####
    # 1.
    infile = 'stats_3LNO_spring.pkl'
    df = pd.read_pickle(infile)
    df.index = pd.to_datetime(df.index) # update index to datetime
    df = df[(df.index.year>=years[0]) & (df.index.year<=years[-1])]
    year_list = df.index.year.astype('str')
    year_list = [i[2:4] for i in year_list] # 2-digit year
    df.index = pd.to_datetime(df.index) # update index to datetime
    df['area_colder0'] = df['area_colder0']/1000 # In 1000km
    df['area_colder1'] = df['area_colder1']/1000 # In 1000km
    df['area_warmer2'] = df['area_warmer2']/1000
    df_clim = df[(df.index.year>=clim_year[0]) & (df.index.year<=clim_year[1])]
    std_anom = (df-df_clim.mean(axis=0))/df_clim.std(axis=0)
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
    table_cells = table_props['child_artists']
    last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
    for key, cell in the_table.get_celld().items():
        cell_text = cell.get_text().get_text()
        if is_number(cell_text) == False:
            pass
        elif key[0] == 0:
            pass
        elif key[1] in last_columns:
             cell._text.set_color('darkslategray')
        elif (np.float(cell_text) <= -1.5) | (np.float(cell_text) >= 1.5) :
            cell._text.set_color('white')
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
    table_cells = table_props['child_artists']
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
    plt.savefig("scorecards_spring_3LNO_FR.png", dpi=300)
    os.system('convert -trim scorecards_spring_3LNO_FR.png scorecards_spring_3LNO_FR.png')


    # 2.
    infile = 'stats_3Ps_spring.pkl'
    df = pd.read_pickle(infile)
    df.index = pd.to_datetime(df.index) # update index to datetime
    df = df[(df.index.year>=years[0]) & (df.index.year<=years[-1])]
    df['area_colder0'] = df['area_colder0']/1000 # In 1000km
    df['area_colder1'] = df['area_colder1']/1000 # In 1000km
    df['area_warmer2'] = df['area_warmer2']/1000
    df_clim = df[(df.index.year>=clim_year[0]) & (df.index.year<=clim_year[1])]
    std_anom = (df-df_clim.mean(axis=0))/df_clim.std(axis=0)
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
    table_cells = table_props['child_artists']
    last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
    for key, cell in the_table.get_celld().items():
        cell_text = cell.get_text().get_text()
        if is_number(cell_text) == False:
            pass
        #elif key[0] == 0: # <--- remove when no years
        #    pass
        elif key[1] in last_columns:
             cell._text.set_color('darkslategray')
        elif (np.float(cell_text) <= -1.5) | (np.float(cell_text) >= 1.5) :
            cell._text.set_color('white')

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
    table_cells = table_props['child_artists']
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

    plt.savefig("scorecards_spring_3Ps_FR.png", dpi=300)
    os.system('convert -trim scorecards_spring_3Ps_FR.png scorecards_spring_3Ps_FR.png')

    # English montage
    os.system('montage  scorecards_spring_3LNO.png scorecards_spring_3Ps.png -tile 1x3 -geometry +1+1  -background white  scorecards_botT_spring.png')
    # French montage
    os.system('montage  scorecards_spring_3LNO_FR.png scorecards_spring_3Ps_FR.png -tile 1x3 -geometry +1+1  -background white  scorecards_botT_spring_FR.png')
