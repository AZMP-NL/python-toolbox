'''
This is a test with xarray, see if I can manipulate a lot of netCDF files.
Try this script in /home/cyrf0006/research/AZMP_database
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
ds = xr.open_mfdataset('/home/cyrf0006/research/AZMP_database/netCDF_first_set/2017.nc')
#ds = xr.open_mfdataset('/home/cyrf0006/research/AZMP_database/2017_data/2017_viking.nc')


# Select a depth range
ds = ds.sel(level=ds['level']<500)
#ds = ds.sel(level=ds['level']>10)

# Selection of a subset region
ds = ds.where((ds.longitude>-55) & (ds.longitude<-50), drop=True)
#ds = ds.where((ds.latitude>50) & (ds.latitude<55), drop=True)
ds = ds.where((ds.latitude>45) & (ds.latitude<50), drop=True)

# Sort time dimension (this takes time to display!!)
ds = ds.isel(time=np.argsort(ds.time))

# Selection of only summer data
#ds_summer = ds.sel(time=ds['time.season']=='JJA')

# Monthly average (This takes a lot of time given the sorting above)
ds_monthly = ds.resample('M', dim='time', how='mean')


# To Pandas Dataframe
da_temp = ds_monthly['temperature']
df_temp = da_temp.to_pandas()


# Only one month
df_temp_may = df_temp.loc[df_temp.index.month==5]
df_temp_june = df_temp.loc[df_temp.index.month==6]
df_temp_july = df_temp.loc[df_temp.index.month==7]
df_temp_aug = df_temp.loc[df_temp.index.month==8]
df_temp_mar = df_temp.loc[df_temp.index.month==3]
df_temp_apr = df_temp.loc[df_temp.index.month==4]

# save in csv format (for Peter)
df_temp_aug.T.to_csv('2017_aug.csv')
df_temp_mar.T.to_csv('2017_mar.csv')




## --- plot --- ## 
plt.figure(1)
plt.clf()
plt.plot(df_temp_mar.values.squeeze(), df_temp_mar.columns, 'k--', linewidth=4)
#plt.plot(df_temp_apr.values.squeeze(), df_temp_apr.columns, 'k--', linewidth=4)
plt.plot(df_temp_aug.values.squeeze(), df_temp_aug.columns, 'k-', linewidth=4)
plt.grid('on')
plt.legend(['Mar 2017', 'Aug 2017'])
plt.xlabel(r'T ($^{\circ}$C)', fontsize=15, fontweight='bold')
plt.ylabel(r'Depth (m)', fontsize=15, fontweight='bold')
plt.gca().invert_yaxis()
plt.show()
