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
df_all.to_pickle('temp_summer_1948-2017.pkl')

## --- Monthly T profile ---- #




## --- CIL core --- ## 
plt.figure(1)
plt.clf()
plt.plot(df_temp_may.index, df_temp_may.min(axis=1), '-')
plt.plot(df_temp_june.index, df_temp_june.min(axis=1), '-')
plt.plot(df_temp_july.index, df_temp_july.min(axis=1), '-')
plt.plot(df_all.index, df_all.min(axis=1), 'k-', linewidth=2)


wm_def = wm.water_masses_def_petrie()

plt.legend(['May', 'June', 'July', 'mean'], fontsize=15)
plt.ylabel(r'$T_{min}$ in monthly mean profile ($^{\circ}$C)', fontsize=15, fontweight='bold')
plt.xlabel('Year', fontsize=15, fontweight='bold')
plt.title('CIL core temperature', fontsize=15, fontweight='bold')
plt.grid('on')
plt.show()


## --- No. of cast per year (monthly) --- ##
years = np.arange(1912, 2017)
time_series = ds.time.to_pandas()
month05 = time_series.loc[time_series.index.month==5]
month06 = time_series.loc[time_series.index.month==6]
month07 = time_series.loc[time_series.index.month==7]
year_count_05 = np.zeros(years.shape)
year_count_06 = np.zeros(years.shape)
year_count_07 = np.zeros(years.shape)
year_count = np.zeros(years.shape)

for idx, year in enumerate(years):
    year_count_05[idx] = np.size(np.where(month05.index.year==year))
    year_count_06[idx] = np.size(np.where(month06.index.year==year))
    year_count_07[idx] = np.size(np.where(month07.index.year==year))
    year_count[idx] = np.size(np.where(time_series.index.year==year))

plt.figure(2)
plt.clf()
plt.plot(years, year_count_05)
plt.plot(years, year_count_06)
plt.plot(years, year_count_07)
plt.plot(years, np.mean(np.stack((year_count_05, year_count_06, year_count_07)), axis=0), 'k-', linewidth=2)
plt.legend(['May', 'June', 'July', 'mean'])
plt.ylabel('No. profiles in monthly average')
plt.show()

# Histogram, yearly
plt.figure(2)
plt.clf()
plt.bar(years, year_count, align='center')
plt.ylabel('Number of stations', fontsize=14, fontweight='bold')
plt.xlabel('Year', fontsize=14, fontweight='bold')
plt.annotate('1992;\nground fish\nmoratorium', xy=(1991, 2000), xytext=(1940, 2000),
            arrowprops=dict(facecolor='black', width=1, shrink=0.05), fontsize=14
            )
plt.annotate('1999;\nAZMP', xy=(1999, 750), xytext=(1999, 1800),
            arrowprops=dict(facecolor='black', width=1, shrink=0.05), fontsize=14
            )

plt.show()



## --- basic T-S properties --- #
T = np.array(df_temp)
TT = np.reshape(T, T.size)
S = np.array(df_sal)
SS = np.reshape(S, S.size)
TT = TT[~np.isnan(SS)]
SS = SS[~np.isnan(SS)]
SS = SS[~np.isnan(TT)]
TT = TT[~np.isnan(TT)]
TT = TT[(SS>20) & (SS<36)]
SS = SS[(SS>20) & (SS<36)]

plt.figure(1)
plt.clf()
plt.plot(SS,TT, '.k')

for w in wm_def.keys():
    shape = np.array(wm_def.get(w))
    p = plt.plot(shape[:,0], shape[:,1])
    co = p[0].get_color() 
    plt.text(np.mean(shape[:,0]), np.mean(shape[:,1]), w, color=co)
    
plt.show()


# T-S seasonal timeseries
df_temp_season = df_temp.resample('3M').mean()
df_sal_season = df_sal.resample('3M').mean()
df_sal_year = df_sal.resample('A').mean()

T = np.array(df_temp_season)
TT = np.reshape(T, T.size)
S = np.array(df_sal_season)

Tbin = np.arange(-2,10)
Slist = []
for idx in np.arange(0, T.shape[0]):
    TVec = T[idx,:]
    SVec = T[idx,:]  
    digitized = np.digitize(TVec, Tbin) #<- this is awesome! 
    Slist.append([SVec[digitized == i].mean() for i in range(0, len(Tbin))])


plt.figure(2)
plt.clf()
#plt.contourf(df_temp.index, Tbin, np.array(Slist).T, 20)  
plt.pcolormesh(df_temp_season.index, Tbin, np.array(Slist).T, vmin=-2, vmax=10)  
plt.show()

## Tcontour-seasonal
fig = plt.figure(3)
plt.clf()
v = np.arange(-2,10)
plt.contourf(df_temp_season.index, df_temp_season.columns, df_temp_season.T, v, cmap=plt.cm.RdBu_r)  
plt.xlim([pd.Timestamp('1950-01-01'), pd.Timestamp('2017-12-01')])
plt.gca().invert_yaxis()
plt.ylabel('Depth (m)', fontsize=14, fontweight='bold')
plt.xlabel('years', fontsize=14, fontweight='bold')
cb = plt.colorbar()
plt.title(r'$\rm T(^{\circ}C)$', fontsize=14, fontweight='bold')
fig.set_size_inches(w=12,h=9)
fig_name = 'Tcontours_1950-2016.png'
fig.savefig(fig_name)


## Scontour-seasonal
fig = plt.figure(4)
plt.clf()
v = np.linspace(30, 35, 20)
plt.contourf(df_sal_season.index, df_sal_season.columns, df_sal_season.T, v, cmap=plt.cm.RdBu_r)  
plt.xlim([pd.Timestamp('1950-01-01'), pd.Timestamp('2017-12-01')])
plt.gca().invert_yaxis()
plt.ylabel('Depth (m)', fontsize=14, fontweight='bold')
plt.xlabel('years', fontsize=14, fontweight='bold')
cb = plt.colorbar()
plt.title(r'$\rm S$', fontsize=14, fontweight='bold')
fig.set_size_inches(w=12,h=9)
fig_name = 'Scontours_1950-2016.png'
fig.savefig(fig_name)

# Salinity timeseries
plt.figure(5)
plt.clf()
plt.plot(df_sal_season.index, df_sal_season.iloc[:,df_sal_season.columns<300].mean(axis=1), '-')
#plt.plot(df_sal_year.index, df_sal_year.iloc[:,df_sal_year.columns<200].mean(axis=1), '-r')
A = df_sal_season.iloc[:,df_sal_season.columns<300].mean(axis=1)
plt.plot(df_sal_season.index, A.rolling(16, center=True).mean(), 'r')
df_sal_season.to_pickle('salinity_1948-2017.pkl')

## plt.plot(df_sal_season.index, df_sal_season.iloc[:,((df_sal_season.columns>150)&(df_sal_season.columns<300))].mean(axis=1), '-')
## #plt.plot(df_sal_year.index, df_sal_year.iloc[:,((df_sal_season.columns>200)&(df_sal_season.columns<300))].mean(axis=1), '-r')
## B =  df_sal_season.iloc[:,((df_sal_season.columns>150)&(df_sal_season.columns<300))].mean(axis=1)
## plt.plot(df_sal_season.index, B.rolling(24).mean(), 'r')


plt.ylabel(r'$S', fontsize=15, fontweight='bold')
plt.xlabel('Year', fontsize=15, fontweight='bold')
plt.title(r'$<S>_{0-300m}$', fontsize=15, fontweight='bold')
plt.grid('on')
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



