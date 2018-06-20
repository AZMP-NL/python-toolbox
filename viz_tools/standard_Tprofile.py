'''
This is a test with xarray, see if I can manipulate a lot of netCDF files.
Try this script in /home/cyrf0006/research/AZMP_database
'''

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr
import datetime

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 14}
plt.rc('font', **font)

# This is a dataset
#ds = xr.open_mfdataset('/home/cyrf0006/AZMP/database/netCDF_first_set/2017.nc')
#ds = xr.open_mfdataset('/home/cyrf0006/research/AZMP_database/2017_data/2017_viking.nc')
ds = xr.open_mfdataset('/home/cyrf0006/data/dev_database/*.nc')


# Select a depth range
ds = ds.sel(level=ds['level']<500)
#ds = ds.sel(level=ds['level']>10)

# Selection of a subset region
ds = ds.where((ds.longitude>-55) & (ds.longitude<-50), drop=True)
ds = ds.where((ds.latitude>45) & (ds.latitude<50), drop=True)

#ds = ds.where((ds.longitude>-52) & (ds.longitude<-50), drop=True)
#ds = ds.where((ds.latitude>48) & (ds.latitude<50), drop=True)



# only 2017
ds = ds.sel(time=ds['time.year']>=2015)

# Sort time dimension (this takes time to display!!)
ds = ds.sortby('time')


# Selection of only summer data
#ds_summer = ds.sel(time=ds['time.season']=='JJA')

# Monthly average (This takes a lot of time given the sorting above)
ds_monthly = ds.resample(time="M").mean('time') 

# To Pandas Dataframe
da_temp = ds_monthly['temperature']
df_temp = da_temp.to_pandas()
da_sal = ds_monthly['salinity']
df_sal = da_sal.to_pandas()


# --- Climatology --- (useful only if multi-year are used (otehrwise should not change the result)) 
df_temp_clim = df_temp.groupby(df_temp.index.month).mean()  
df_temp_clim.index=pd.date_range('2000-01-31', periods=12, freq='M')
df_temp_clim.index.name = 'time'

df_sal_clim = df_sal.groupby(df_sal.index.month).mean()  
df_sal_clim.index=pd.date_range('2000-01-31', periods=12, freq='M')
df_sal_clim.index.name = 'time'

df_temp = df_temp_clim # By-pass df_temp
df_sal = df_sal_clim # By-pass df_sal




# Only one month
df_temp_may = df_temp.loc[df_temp.index.month==5]
df_temp_june = df_temp.loc[df_temp.index.month==6]
df_temp_july = df_temp.loc[df_temp.index.month==7]
df_temp_aug = df_temp.loc[df_temp.index.month==8]
df_temp_mar = df_temp.loc[df_temp.index.month==3]
df_temp_apr = df_temp.loc[df_temp.index.month==4]
df_temp_nov = df_temp.loc[df_temp.index.month==11]

df_sal_may = df_sal.loc[df_sal.index.month==5]
df_sal_june = df_sal.loc[df_sal.index.month==6]
df_sal_july = df_sal.loc[df_sal.index.month==7]
df_sal_aug = df_sal.loc[df_sal.index.month==8]
df_sal_mar = df_sal.loc[df_sal.index.month==3]
df_sal_apr = df_sal.loc[df_sal.index.month==4]
df_sal_nov = df_sal.loc[df_sal.index.month==11]


# save in csv format (for Peter)
## df_temp_aug.T.to_csv('2017_temp_aug.csv')
## df_temp_mar.T.to_csv('2017_temp_mar.csv')
## df_temp_apr.T.to_csv('2017_temp_apr.csv')

## df_sal_aug.T.to_csv('2017_sal_aug.csv')
## df_sal_mar.T.to_csv('2017_sal_mar.csv')
## df_sal_apr.T.to_csv('2017_sal_apr.csv')

df_temp_aug.T.to_csv('clim_temp_aug.csv')
df_temp_mar.T.to_csv('clim_temp_mar.csv')
df_temp_apr.T.to_csv('clim_temp_apr.csv')

df_sal_aug.T.to_csv('clim_sal_aug.csv')
df_sal_mar.T.to_csv('clim_sal_mar.csv')
df_sal_apr.T.to_csv('clim_sal_apr.csv')


## --- plot --- ## 
fig = plt.figure(1)
plt.clf()
plt.plot(df_temp_mar.values.squeeze(), df_temp_mar.columns, linewidth=4)
plt.plot(df_temp_june.values.squeeze(), df_temp_june.columns, linewidth=4)
#plt.plot(df_temp_may.values.squeeze(), df_temp_may.columns, 'k--', linewidth=4)
plt.plot(df_temp_aug.values.squeeze(), df_temp_aug.columns, linewidth=4)
plt.plot(df_temp_nov.values.squeeze(), df_temp_nov.columns, linewidth=4)
plt.grid('on')
plt.legend(['Mar', 'June', 'Aug', 'Nov'])
#plt.legend(['Apr 2017', 'Aug 2017'])
plt.xlabel(r'T ($^{\circ}$C)', fontsize=15, fontweight='bold')
plt.ylabel(r'Depth (m)', fontsize=15, fontweight='bold')
plt.gca().invert_yaxis()
fig.set_size_inches(w=7,h=8)
fig_name = 'clim_profile_temp.png'
fig.set_dpi(300)
fig.savefig(fig_name)



## --- plot --- ## 
fig = plt.figure(2)
plt.clf()
plt.plot(df_sal_mar.values.squeeze(), df_sal_mar.columns, 'k--', linewidth=4)
#plt.plot(df_sal_may.values.squeeze(), df_sal_may.columns, 'k--', linewidth=4)
#plt.plot(df_sal_apr.values.squeeze(), df_sal_apr.columns, 'k--', linewidth=4)
plt.plot(df_sal_aug.values.squeeze(), df_sal_aug.columns, 'k-', linewidth=4)
plt.grid('on')
plt.legend(['March', 'August'])
plt.xlabel(r'S', fontsize=15, fontweight='bold')
plt.ylabel(r'Depth (m)', fontsize=15, fontweight='bold')
plt.gca().invert_yaxis()
fig.set_size_inches(w=7,h=8)
fig_name = 'clim_profile_sal.png'
fig.set_dpi(300)
fig.savefig(fig_name)

