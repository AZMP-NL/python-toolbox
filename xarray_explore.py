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
# This is a dataset
ds = xr.open_mfdataset('*.nc')


# Select a depth range
ds = ds.sel(level=ds['level']<500)
ds = ds.sel(level=ds['level']>10)

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
da_sal = ds_monthly['salinity']
df_sal = da_sal.to_pandas()

# Only one month
df_temp_may = df_temp.loc[df_temp.index.month==5]
df_temp_june = df_temp.loc[df_temp.index.month==6]
df_temp_july = df_temp.loc[df_temp.index.month==7]
df_concat = pd.concat((df_temp_may, df_temp_june, df_temp_july))
df_all = df_concat.resample('A').mean()

## --- CIL core --- ## 
plt.figure(1)
plt.clf()
plt.plot(df_temp_may.index, df_temp_may.min(axis=1), '-')
plt.plot(df_temp_june.index, df_temp_june.min(axis=1), '-')
plt.plot(df_temp_july.index, df_temp_july.min(axis=1), '-')
plt.plot(df_all.index, df_all.min(axis=1), 'k-', linewidth=2)


wm_def = wm.water_masses_def()

plt.legend(['May', 'June', 'July', 'mean'])
plt.ylabel(r'$T_{min}$ in monthly mean profile ($^{\circ}$C)')
plt.xlabel('Year')
plt.title('CIL core temperature')
plt.grid('on')
plt.show()


## --- No. of cast per year --- ##
years = np.arange(1950, 2016)
time_series = ds.time.to_pandas()
time_series.loc[df_temp.index.month==5]
year_count = np.zeros(years.shape)
for idx, year in enumerate(years):
    year_count[idx] = np.size(np.where(time_series.index.year==year))
    
plt.figure(2)
plt.clf()
plt.plot(years, year_count)
plt.show()


## --- basic T-S properties --- #
T = np.array(df_temp)
TT = np.reshape(T, T.size)
S = np.array(df_sal)
SS = np.reshape(S, S.size)

plt.figure(1)
plt.clf()
plt.plot(SS,TT, '.k')

for w in wm_def.keys():
    shape = np.array(wm_def.get(w))
    p = plt.plot(shape[:,0], shape[:,1])
    co = p[0].get_color() 
    plt.text(np.mean(shape[:,0]), np.mean(shape[:,1]), w, color=co)
    
plt.show()


keyboard







# 





plt.contourf(df_temp.index, df_temp.columns, df_temp.T, 20)


# HERE!!


# this is now a dataArray
da_temp = ds['temperature']
da_lat = ds['latitude']
da_lon = ds['longitude']


# Averrage profile of the whole timeseries
A = da_temp.mean(dim='time')

# Selection of only summer data
summer_temp = temp.sel(time=temp['time.season']=='JJA')
summer_lat = lat.sel(time=lat['time.season']=='JJA')
summer_lon = lon.sel(time=lon['time.season']=='JJA')




ds = xr.Dataset()
ds['data1'] = xr.DataArray(np.arange(100), coords={'t1': np.linspace(0, 1, 100)})
ds['data1b'] = xr.DataArray(np.arange(100, 200), coords={'t1': np.linspace(0, 1, 100)})
ds['data2'] = xr.DataArray(np.arange(200, 500), coords={'t2': np.linspace(0, 1, 300)})
ds['data2b'] = xr.DataArray(np.arange(600, 900), coords={'t2': np.linspace(0, 1, 300)})
ds.where(ds.data1 < 50, drop=True)



