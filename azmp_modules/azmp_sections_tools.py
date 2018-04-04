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
#import netCDF4 as nc
#import os
#from sys import version_info
from math import radians, cos, sin, asin, sqrt
from math import radians, cos, sin, asin, sqrt


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


def standard_section_plot(nc_file, survey_name, section_name, var_name):
    """
    Contour plot on standard AZMP-NL sections for a certain year (specified with nc_file), season (specified as survey), section and variable.

    !!! Note: Maybe it would be better to just return the dataframe and do the plot (ylim, colormap,  etc. on a separate script?)
    Maybe the best way would be to keep this function, but make the plot=False by default and always return dataframe
    - Need to calculate CIL area also (maybe another function using df output.)
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
    ds = ds.sel(time=ds['survey_ID'].values==survey_name)
    df = ds.to_dataframe()

    ## ----  plot FC-line ---- ##
    df_section =  df[df['comments'].str.contains(section_name + '-')]
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
    v = np.linspace(-2,12,36)
    v1 = np.linspace(-2,12,8)
    fig, ax = plt.subplots()
#    c = plt.contourf(distance, df_var.index, df_var, 20, cmap=plt.cm.RdBu_r)
    c = plt.contourf(distance, df_var.index, df_var, v, cmap=plt.cm.RdBu_r)
    plt.contour(distance, df_var.index, df_var, [0,], colors='k')
    ax.set_ylim([0, 400])
    ax.set_xlim([0, distance.max()])
    ax.set_xlim([0, 450])
    ax.set_ylabel('Depth (m)', fontWeight = 'bold')
    ax.set_xlabel('Distance (km)', fontWeight = 'bold')
    fig_title = 'AZMP_' + survey_name + '_' + section_name + '_' + var_name
    ax.set_title(fig_title)
    ax.invert_yaxis()
    cb = plt.colorbar(c, ticks=v1)
#    plt.colorbar(c, cax=cax, ticks=v1)
    
    Bgon = plt.Polygon(bathymetry,color=np.multiply([1,.9333,.6667],.4), alpha=1)
    ax.add_patch(Bgon)
    for i in range(0,len(distance)):
        plt.plot(np.array([distance[i], distance[i]]), np.array([df_var.index[0], cast_depth[i]]), '--k', linewidth=0.1)    

    fig.set_size_inches(w=9,h=4.5)
    fig_name = fig_title + '.png'
    fig.set_dpi(300)
    fig.savefig(fig_name)
    return None
