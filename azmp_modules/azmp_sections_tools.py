"""Some tools to visualize NAFC custom pfiles (e.g., AZMP sections)

References
----------

Atlantic Zone Monitoring Program @NAFC:
https://azmp-nl.github.io/

"""

__author__ = 'Frederic.Cyr@dfo-mpo.gc.ca'
__version__ = '0.1'

## to do:
# - 
# - 
# - 

import re
import pandas as pd
import matplotlib.pyplot as plt
#import pfiles_basics
import numpy as np
#import time as tt
import xarray as xr
import netCDF4
#import os
#from sys import version_info
from math import radians, cos, sin, asin, sqrt
from math import radians, cos, sin, asin, sqrt
from scipy.interpolate import interp1d  # to remove NaNs in profiles
from scipy.interpolate import griddata

def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    km = 6367 * c
    return km


def section_bathymetry(section_name):
    """
    Retrieve high-resultion bathymetry (depth vs along-distance from 1st station on section) along AZMP-NL sections.
    Bathymetry taken from high-resolution data (https://github.com/AZMP-NL/bathymetry)


    !!!TO DO: - add an error if specified section name doesn't exist
              - add optional 'path' input. Default can wget from AZMP git repo...
    """

    standard_sections = ["3PS", "FC", "SEGB", "SWSP", "BB", "FI", "S27", "SESP", "WB", "BI", "MB", "SA", "SI"]
    bathy_files = ["3psline.txt", "fcline.txt", "segbline.txt", "swspbline.txt", "bbline.txt", \
                   "filine.txt", "s27line.txt", "sespbline.txt", "wbline.txt", "biline.txt", \
                   "mbline.txt", "saline.txt", "siline.txt"]

    bathy_dict = dict(zip(standard_sections, bathy_files))

    # This is my personal path, maybe find a way to generelize this.
    bathy_file = '/home/cyrf0006/github/AZMP-NL/bathymetry/bottom_profiles/' + bathy_dict[section_name]
    
    bathy = np.loadtxt(bathy_file, delimiter=",", unpack=False)
    bathy_x = bathy[:,0]/1000.0
    bathy_y = np.abs(bathy[:,1])
    
    bathy_x_close = np.append(bathy_x, np.array([bathy_x[-1], bathy_x[0], bathy_x[0]]))
    bathy_y_close = np.append(bathy_y ,np.array([bathy_y.max(), bathy_y.max(), bathy_y[0]]))
    bathymetry = zip(bathy_x_close, bathy_y_close)
    return bathymetry


def get_section(section_name, year, season, var_name, dlat=2, dlon=2, dc=.2, dz=1):
    """
    To get a transect plot of a certain variable along a certain section and for a defined year and season.

    Input:
    section_name - Abbreviation of the section ('BB', 'FC, 'SEGB', etc.)
    year - the year (e.g., 2018)
    season - 'spring', 'summer', 'fall', 'all'
    var_name - Name of variable in netCDF files (e.g., 'temperature', 'salinity')
    dlat, dlon - How far from the station (in degrees) we search (default 2) 
    dc - Disretization of the regulat grid for field 3D interpolation (e.g., the temperature cube around section, default 2x2)
    dz - vertical discretization (z bin size; default 1)

    Output:
    df_stn - DataFrame of stations vs depth by looking at station names only (work not well for the period before AZMP)
    df_itp - DataFrame of stations vs depth obtained from the interpolation of all available data.
    
    usage example:
    import azmp_sections_tools as azst
    df_stn, df_itp = azst.get_section('BB', 2018, 'summer', 'temperature')
    
    """
    ## ---- Get Stations ---- ## 
    df_stn = pd.read_excel('/home/cyrf0006/github/AZMP-NL/data/STANDARD_SECTIONS.xlsx')
    df_stn = df_stn.drop(['SECTION', 'LONG'], axis=1)
    df_stn = df_stn.rename(columns={'LONG.1': 'LON'})
    df_stn = df_stn.dropna()
    df_stn = df_stn[df_stn.STATION.str.contains(section_name+'-')]
    df_stn = df_stn.reset_index(drop=True)
    stn_list = df_stn.STATION.values

    # Derive regular grid
    latLims = np.array([df_stn.LAT.min() - dlat, df_stn.LAT.max() + dlat])
    lonLims = np.array([df_stn.LON.min() - dlon, df_stn.LON.max() + dlon])
    lon_reg = np.arange(lonLims[0]+dc/2, lonLims[1]-dc/2, dc)
    lat_reg = np.arange(latLims[0]+dc/2, latLims[1]-dc/2, dc)

    ## --------- Get Bathymetry -------- ####
    print('Get bathy...')
    dataFile = '/home/cyrf0006/data/GEBCO/GEBCO_2014_1D.nc' # Maybe find a better way to handle this file
    lon_grid, lat_grid = np.meshgrid(lon_reg,lat_reg)
    # Load data
    dataset = netCDF4.Dataset(dataFile)
    x = [-179-59.75/60, 179+59.75/60] # to correct bug in 30'' dataset?
    y = [-89-59.75/60, 89+59.75/60]
    spacing = dataset.variables['spacing']
    # Compute Lat/Lon
    nx = int((x[-1]-x[0])/spacing[0]) + 1  # num pts in x-dir
    ny = int((y[-1]-y[0])/spacing[1]) + 1  # num pts in y-dir
    lon = np.linspace(x[0],x[-1],nx)
    lat = np.linspace(y[0],y[-1],ny)
    # interpolate data on regular grid 
    # Reshape data
    zz = dataset.variables['z']
    Z = zz[:].reshape(ny, nx)
    Z = np.flipud(Z) # <------------ important!!!
    # Reduce data according to Region params
    idx_lon = np.where((lon>=lonLims[0]) & (lon<=lonLims[1]))
    idx_lat = np.where((lat>=latLims[0]) & (lat<=latLims[1]))
    Z = Z[idx_lat[0][0]:idx_lat[0][-1]+1, idx_lon[0][0]:idx_lon[0][-1]+1]
    lon = lon[idx_lon[0]]
    lat = lat[idx_lat[0]]
    # interpolate data on regular grid 
    lon_grid_bathy, lat_grid_bathy = np.meshgrid(lon,lat)
    lon_vec_bathy = np.reshape(lon_grid_bathy, lon_grid_bathy.size)
    lat_vec_bathy = np.reshape(lat_grid_bathy, lat_grid_bathy.size)
    z_vec = np.reshape(Z, Z.size)
    Zitp = griddata((lon_vec_bathy, lat_vec_bathy), z_vec, (lon_grid, lat_grid), method='linear')
    print(' -> Done!')


    ## -------- Get CTD data -------- ##
    year_file = '/home/cyrf0006/data/dev_database/' + str(year) + '.nc'
    print('Get ' + year_file)
    ds = xr.open_mfdataset(year_file)

    # Remame problematic datasets
    print('!!Remove MEDBA & MEDTE data!!')
    print('  ---> I Should be improme because I remove good data!!!!')
    ds = ds.where(ds.instrument_ID!='MEDBA', drop=True)
    ds = ds.where(ds.instrument_ID!='MEDTE', drop=True)

    # Select Region
    ds = ds.where((ds.longitude>lonLims[0]) & (ds.longitude<lonLims[1]), drop=True)
    ds = ds.where((ds.latitude>latLims[0]) & (ds.latitude<latLims[1]), drop=True)

    # Select time (save several options here)
    if season == 'summer':
        ds = ds.sel(time=((ds['time.month']>=7)) & ((ds['time.month']<=9)))
    elif season == 'spring':
        ds = ds.sel(time=((ds['time.month']>=3)) & ((ds['time.month']<=5)))
    elif season == 'fall':
        ds = ds.sel(time=((ds['time.month']>=10)) & ((ds['time.month']<=12)))
    else:
        print('!! no season specified, used them all! !!')

    # Extract variable
    da = ds[var_name]
    lons = np.array(ds.longitude)
    lats = np.array(ds.latitude)
    bins = np.arange(dz/2.0, 500, dz)
    da = da.groupby_bins('level', bins).mean(dim='level')
    #To Pandas Dataframe
    df = da.to_pandas()
    df.columns = bins[0:-1] #rename columns with 'bins'
    # Remove empty columns
    idx_empty_rows = df.isnull().all(1).nonzero()[0]
    df = df.dropna(axis=0,how='all')
    lons = np.delete(lons,idx_empty_rows)
    lats = np.delete(lats,idx_empty_rows)
    print(' -> Done!')


    ## --- fill 3D cube --- ##  
    print('Fill regular cube')
    z = df.columns.values
    V = np.full((lat_reg.size, lon_reg.size, z.size), np.nan)

    # Aggregate on regular grid
    for i, xx in enumerate(lon_reg):
        for j, yy in enumerate(lat_reg):    
            idx = np.where((lons>=xx-dc/2) & (lons<xx+dc/2) & (lats>=yy-dc/2) & (lats<yy+dc/2))
            tmp = np.array(df.iloc[idx].mean(axis=0))
            #idx_good = np.argwhere((~np.isnan(tmp)) & (tmp<30))
            idx_good = np.argwhere((~np.isnan(tmp)))
            if np.size(idx_good)==1:
                V[j,i,:] = np.array(df.iloc[idx].mean(axis=0))
            elif np.size(idx_good)>1: # vertical interpolation between pts
                interp = interp1d(np.squeeze(z[idx_good]), np.squeeze(tmp[idx_good]))  # <---- Attention: strange syntax but works
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
        if idx_good.size>7: # will ignore depth where no data exist
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


    ## ---- Extract section info (test 2 options) ---- ##
    temp_coords = np.array([[x,y] for x in lat_reg for y in lon_reg])
    VV = V.reshape(temp_coords.shape[0],V.shape[2])
    ZZ = Zitp.reshape(temp_coords.shape[0],1)
    section_only = []
    df_section_itp = pd.DataFrame(index=stn_list, columns=z)
    for stn in stn_list:
        # 1. Section only (by station name)
        ds_tmp = ds.where(ds.comments == stn, drop=True)  
        section_only.append(ds_tmp)

        #2.  From interpolated field (closest to station)
        station = df_stn[df_stn.STATION==stn]
        idx_opti = np.argmin(np.sum(np.abs(temp_coords - np.array(zip(station.LAT,station.LON))), axis=1))
        Tprofile = VV[idx_opti,:]
        # remove data below bottom
        bottom_depth = -ZZ[idx_opti]
        Tprofile[z>=bottom_depth]=np.nan
        # store in dataframe
        df_section_itp.loc[stn] = Tprofile

    df_section_itp.index.name = 'station'    
   
    # convert option #1 to dataframe    
    ds_section = xr.concat(section_only, dim='time')
    da = ds_section['temperature']
    df_section_stn = da.to_pandas()
    df_section_stn.index = ds_section.comments.values
    df_section_stn.index.name = 'station'    
    
    # Compute distance vector for option #1 - exact station
    distance_stn = np.full((df_section_stn.index.shape), np.nan)
    for i, stn in enumerate(df_section_stn.index):
        distance_stn[i] = haversine(df_stn.LON[0], df_stn.LAT[0],
                                    df_stn[df_stn.STATION==df_section_stn.index[i]].LON,
                                    df_stn[df_stn.STATION==df_section_stn.index[i]].LAT)
    
    # Compute distance vector for option #2 - interp field
    distance_itp = np.full((df_section_itp.index.shape), np.nan)
    for i, stn in enumerate(df_section_itp.index):
        distance_itp[i] = haversine(df_stn.LON[0], df_stn.LAT[0],
                                    df_stn[df_stn.STATION==df_section_itp.index[i]].LON,
                                    df_stn[df_stn.STATION==df_section_itp.index[i]].LAT)

    # Add a second index with distance
    df_section_stn['distance'] = distance_stn  
    df_section_stn.set_index('distance', append=True, inplace=True)
    df_section_itp['distance'] = distance_itp  
    df_section_itp.set_index('distance', append=True, inplace=True)
    
    return df_section_stn, df_section_itp

def standard_section_plot(nc_file, survey_name, section_name, var_name):
    """
    Contour plot on standard AZMP-NL sections for a certain year (specified with nc_file), season (specified as survey), section and variable.

    !!! Note: Maybe it would be better to just return the dataframe and do the plot (ylim, colormap,  etc. on a separate script?)
    Maybe the best way would be to keep this function, but make the plot=False by default and always return dataframe
    -> Need to calculate CIL area also (maybe another function using df output.)
    -> I think survey_name is non-existant in nc files... (rather cast IDs). Need to modify pfile_tools...
    (actually may bne working now) - edit 2018-05-29

    NOT FINISHED!
    ex:
    import azmp_sections_tools as azs
    azs.standard_section_plot('/home/cyrf0006/data/dev_database/2018.nc', TEL)


    """

    # For plots
    font = {'family' : 'normal',
            'weight' : 'bold',
            'size'   : 14}
    plt.rc('font', **font)

    # Open netCDF file using xarray
    ds = xr.open_dataset(nc_file)
    zVec = ds.level.to_masked_array()

    # Select Survey and switch to DataFrame
    #ds = ds.sel(time=ds['trip_ID'].values==survey_name)
    df = ds.to_dataframe()

    ## ----  plot FC-line ---- ##
    df_section =  df[df['trip_ID'].str.contains(survey_name)] # Restrict for survey-name
    df_section =  df_section[df_section['comments'].str.contains(section_name + '-')] # restrict for section name
    sr_var = df_section[var_name]
    sr_lon = df_section['longitude']
    sr_lat = df_section['latitude']
    sr_stn = df_section['comments']

    df_stn = sr_stn.unstack()
    df_var = sr_var.unstack()
    df_lat = sr_lat.unstack()
    df_lon = sr_lon.unstack()

    # compute distance
    lon_array = df_lon.values[0,]
    lat_array = df_lat.values[0,]
    distance = np.zeros(lon_array.shape)
    for i in range(lon_array.size):
        distance[i] = haversine(lon_array[0], lat_array[0], lon_array[i], lat_array[i])

    # Check which direction we are going (approaching St.27 or not)
    St27 = [47.550, -52.590]
    if haversine(lon_array[0], lat_array[0], St27[1], St27[0]) > haversine(lon_array[-1], lat_array[-1], St27[1], St27[0]):
        distance = np.abs(distance-distance.max())
        
    # retrieve bathymetry using function
    bathymetry = section_bathymetry(section_name)

    # Check maximum depth (for cast positioning on figure)
    bathy = np.array(bathymetry) # list to array
    cast_depth = []
    for i in distance:
        min_idx = np.argmin(np.abs(i-bathy[:,0]))
        cast_depth = np.append(cast_depth, bathy[:,1][min_idx])

        
    # Do the plot
    v = 20
    #v1 = 
    #v = np.linspace(-2,12,36)
    #v1 = np.linspace(-2,12,8)
    fig, ax = plt.subplots()
#    c = plt.contourf(distance, df_var.index, df_var, 20, cmap=plt.cm.RdBu_r)
    c = plt.contourf(distance, df_var.index, df_var, v, cmap=plt.cm.RdBu_r)
    plt.contour(distance, df_var.index, df_var, [0,], colors='k')
    ax.set_ylim([0, 400])
    #ax.set_xlim([0, distance.max()])
    ax.set_xlim([0, 300])
    ax.set_ylabel('Depth (m)', fontWeight = 'bold')
    ax.set_xlabel('Distance (km)', fontWeight = 'bold')
    fig_title = 'AZMP_' + survey_name + '_' + section_name + '_' + var_name
    ax.set_title(fig_title)
    ax.invert_yaxis()
    #cb = plt.colorbar(c, ticks=v1)
#    plt.colorbar(c, cax=cax, ticks=v1)
    plt.colorbar(c)

    Bgon = plt.Polygon(bathymetry,color=np.multiply([1,.9333,.6667],.4), alpha=1)
    ax.add_patch(Bgon)
    for i in range(0,len(distance)):
        plt.plot(np.array([distance[i], distance[i]]), np.array([df_var.index[0], cast_depth[i]]), '--k', linewidth=0.1)    

    fig.set_size_inches(w=9,h=4.5)
    fig_name = fig_title + '.png'
    fig.set_dpi(300)
    fig.savefig(fig_name)
    return None


def extract_section_casts(nc_file, section_name, year_lims=[], survey_name=[], nc_outfile='out.nc'):
    """
    To extract hydrographic data from a certain section.
    [Menu to be finished]

    Eventually, I could add conditions on 'year_lims' and 'survey_name' that are not considered right now.
    
    Initially created after a request by F. Li on SI data.

    A full usage example:
    import azmp_sections_tools as azst
    nc_file = '/home/cyrf0006/data/dev_database/201*.nc'
    azst.extract_section_casts(nc_file, section_name='SI', nc_outfile='SI.nc')
    
    Frederic.Cyr@dfo-mpo.gc.ca
    December 2018    
    
    """
    # Open netCDF file using xarray
    ds = xr.open_mfdataset(nc_file)

    # Read unique station names from section file
    df_stn = pd.read_excel('/home/cyrf0006/github/AZMP-NL/data/STANDARD_SECTIONS.xlsx')
    df_stn = df_stn.drop(['SECTION', 'LONG'], axis=1)
    df_stn = df_stn.rename(columns={'LONG.1': 'LON'})
    df_stn = df_stn.dropna()
    df_stn = df_stn[df_stn.STATION.str.contains(section_name)]
    df_stn = df_stn.reset_index(drop=True)
    stn_list = df_stn.STATION.values

    # loop on stations and append datasets in a list
    datasets = []
    for stn in stn_list:
        ds_section = ds.where(ds.comments == stn, drop=True)  
        datasets.append(ds_section)

    # concatenate the list of datasets
    ds_combined = xr.concat(datasets, dim='time')

    # save combined dataset in NetCDF
    ds_combined.to_netcdf(nc_outfile)
    
    
    return None

