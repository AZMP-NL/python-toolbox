'''
Try to plot Viking buoy stuff from netCDF files
'''

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr
import datetime
import water_masses as wm

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 14}
plt.rc('font', **font)

# This is a dataset
ds = xr.open_mfdataset('2017_viking.nc')


# Select a depth range
ds = ds.sel(level=ds['level']<200)
ds = ds.sel(level=ds['level']>0)

# Sort time dimension (this takes time to display!!)
#ds = ds.isel(time=np.argsort(ds.time))

# Selection of only summer data
#ds_summer = ds.sel(time=ds['time.season']=='JJA')

# Monthly average (This takes a lot of time given the sorting above)
ds_sub = ds.resample('1M', dim='time', how='mean')


# To Pandas Dataframe
da_temp = ds_sub['temperature']
df_temp = da_temp.to_pandas()
da_sal = ds_sub['salinity']
df_sal = da_sal.to_pandas()


## --- CIL core --- ## 
plt.figure(1)
plt.clf()
plt.contourf(df_temp.index, df_temp.columns, df_temp.T)
plt.grid('on')
plt.xlabel('Time', fontsize=15, fontweight='bold')
plt.ylabel('Depth (m)', fontsize=15, fontweight='bold')
plt.gca().invert_yaxis()
plt.show()
