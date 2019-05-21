'''
This is a test with xarray, see if I can manipulate a lot of netCDF files.
Try this script in /home/cyrf0006/research/AZMP_database/2017_data
'''

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr
import datetime
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from scipy.interpolate import griddata # here should remove nans or empty profiles
import netCDF4
#import matplotlib
#matplotlib.interactive(True)

# For plots
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 14}
plt.rc('font', **font)

# This is a dataset
# in /home/cyrf0006/research/AZMP_database/netCDF_first_set
ds = xr.open_dataset('/home/cyrf0006/data/dev_database/2018.nc')

# Selection of a subset region
ds = ds.where((ds.longitude>-61) & (ds.longitude<-56), drop=True)
ds = ds.where((ds.latitude>44) & (ds.latitude<48), drop=True)

# Select a depth range
ds = ds.sel(level=ds['level']<500)

# Save NetCDF
ds.to_netcdf('LS_casts_2018.nc')

# SAve CSV T,S
da_temp = ds.temperature
df_temp = da_temp.to_pandas()
df_temp.to_csv('LC_temp_2018.csv', sep=',')

da_sal = ds.salinity
df_sal = da_sal.to_pandas()
df_sal.to_csv('LC_sal_2018.csv', sep=',')

# SAve CSV lat.lon
lon = ds.longitude.values
np.savetxt("LC_lon_.csv", lon, delimiter=",", fmt='%1.5f')
lat = ds.latitude.values
np.savetxt("LC_lat_.csv", lat, delimiter=",", fmt='%1.5f')

# To get the list of pfiles for Galbraith:
ids = ds.trip_ID.values
np.savetxt("LC_ids_.csv", ids, delimiter=",", fmt='%s')

# Then I ran this in /home/cyrf0006/data/pfiles2018/ on orwell:
#while read i; do cp $i* LC_galbraith; done < LC_ids_.csv 
#zip -r LC_galbraith.zip LC_galbraith/
