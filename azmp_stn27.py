# -*- coding: utf-8 -*-
'''
Station 27 analysis for CSAS/NAFO ResDocs. This scripts compute:

- vertically average T-S
- CIL T, coreT and thickness
- 

** see azmp_stn27_explore.py for more options and ways to explore the dataset

Frederic.Cyr@dfo-mpo.gc.ca - June 2019

'''

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr
import datetime
import os

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 14}
plt.rc('font', **font)

## ---- Some custom parameters ---- ##
s27 = [47.55,-52.59]
dc = .1
year_clim = [1981, 2010]
current_year = 2018
variable = 'temperature'
use_viking = True
XLIM = [datetime.date(current_year, 1, 1), datetime.date(current_year, 12, 31)]
if variable == 'temperature':
    V = np.arange(-2, 20, .5)
elif variable == 'salinity':
    V = np.arange(29.5, 34, .25)


## ---- Open data and select ---- ##
ds = xr.open_mfdataset('/home/cyrf0006/data/dev_database/*.nc')

# Remome GTS datasets
ds = ds.where(ds.instrument_ID!='MEDBA', drop=True) # BATHY GTS message 
ds = ds.where(ds.instrument_ID!='MEDTE', drop=True) # TESAC GTS message 

# Select a depth range
ds = ds.sel(level=ds['level']<180)
ds = ds.sel(level=ds['level']>0)

# Select stn27 data according to lat-lon in a box [47.55,-52.59]
ds = ds.where((ds.longitude>s27[1]-dc/2) & (ds.longitude<s27[1]+dc/2), drop=True) # original one
ds = ds.where((ds.latitude>s27[0]-dc/2) & (ds.latitude<s27[0]+dc/2), drop=True)
ds = ds.sortby('time')

da = ds[variable]
df = da.to_pandas()

## ---- Open Viking data and concatenate ---- ##
if use_viking:
    # open dataset
    ds_vik = xr.open_mfdataset('/home/cyrf0006/data/dev_database/viking_nc/*_viking.nc')
    # Select a depth range
    ds_vik = ds_vik.sel(level=ds_vik['level']<180)
    ds_vik = ds_vik.sel(level=ds_vik['level']>0)
    # get variable
    da_vik = ds_vik[variable]
    df_vik = da_vik.to_pandas()
    df_vik = df_vik.dropna(how='all')
    # concatenate (better if they sahre same vertical resolution)
    df = pd.concat([df, df_vik]) 


## ---- 1. Climatologies ---- ##

# Monthly average (15th of the month)
df_monthly = df.resample('MS', loffset=pd.Timedelta(14, 'd')).mean()
# Select years for climato
df_clim_period = df_monthly[(df_monthly.index.year>=year_clim[0]) & (df_monthly.index.year<=year_clim[1])]
# Monthly clim
monthly_clim = df_clim_period.groupby(df_clim_period.index.month).mean()
# set index (year 1900)
monthly_clim.index = pd.to_datetime(monthly_clim.index.values, format='%m')
monthly_clim.to_pickle('s27_' + variable + '_monthly_clim.pkl')

# Weekly average (previous version of climatology)
#df_weekly = df.resample('W').mean()
# Monthly clim
#weekly_clim_raw = df_weekly.groupby(df_weekly.index.week).mean()
# set index (year 1900)
#weekly_clim.index = pd.to_datetime(weekly_clim.index.values, format='%W')

# Weekly clim (upsample monthly clim to weekly)
weekly_clim = monthly_clim.resample('W').mean().interpolate(method='linear') 
weekly_clim.to_pickle('s27_' + variable + '_weekly_clim.pkl')

# Update climatology index to current year
weekly_clim.index = pd.to_datetime('2018-' +  weekly_clim.index.month.astype(np.str) + '-' + weekly_clim.index.day.astype(np.str))
weekly_clim = weekly_clim.dropna(how='all', axis=1)

# plot
fig, ax = plt.subplots(nrows=1, ncols=1)
c = plt.contourf(weekly_clim.index, weekly_clim.columns, weekly_clim.values.T, V, extend='both', cmap=plt.cm.RdBu_r)
#cc = plt.contour(weekly_clim.index, weekly_clim.columns, weekly_clim.values.T, V, colors='k')
#plt.plot(weekly_clim.index, np.repeat(0, df_temp.index.size), '|k', markersize=20)    
cax = fig.add_axes([0.91, .15, 0.01, 0.7])
cb = plt.colorbar(c, cax=cax, orientation='vertical')
if variable == 'temperature':
    ccc = plt.contour(weekly_clim.index, weekly_clim.columns, weekly_clim.values.T, [0,], colors='k', linewidths=3)
    cb.set_label(r'$\rm T(^{\circ}C)$', fontsize=12, fontweight='normal')
elif variables== 'salinity':
    cb.set_label(r'S', fontsize=12, fontweight='normal')
plt.ylim([0, 175])
plt.xlim([XLIM[0], XLIM[1]])
#plt.clabel(cc, inline=1, fontsize=10, colors='k', fmt='%1.1f')
plt.grid('on')
plt.xlabel('Time', fontsize=15, fontweight='bold')
plt.ylabel('Depth (m)', fontsize=15, fontweight='bold')
plt.ylim([0, 175])
plt.gca().invert_yaxis()
# Save Figure
fig.set_size_inches(w=12, h=6)
fig.set_dpi(200)
outfile = 's27_' + variable + '_clim.png'
fig.savefig(outfile)
os.system('convert -trim ' + outfile + ' ' + outfile)



## ---- 2. Year average and anomaly ---- ##
df_year = df[df.index.year==current_year]
df_weekly = df_year.resample('W').mean().interpolate(method='linear') 
df_weekly.dropna(how='all')


anom = df_weekly - weekly_clim





# HERE!!!





keyboard


da = ds['time']
df = da.to_pandas()
df = df[df.index.year>=1999]
plt.plot(df.index.year, df.index.weekofyear, '.k')
plt.ylabel('Week of year')
#plt.ylabel('Month')
plt.xlabel('Year')
plt.xlabel('Station 27 occupation')
plt.show()



# To Pandas Dataframe
#a_temp = ds_monthly['temperature']
#a_temp = ds_season['temperature']
df_temp = da_temp.to_pandas()
#da_sal = ds_season['salinity']
da_sal = ds_monthly['salinity']
df_sal = da_sal.to_pandas()




#here
fig = plt.figure(3)
plt.clf()
v = np.arange(-2,10)
plt.contourf(df_temp.index, df_temp.columns, df_temp.T, v, cmap=plt.cm.RdBu_r, extend='both')  
plt.gca().invert_yaxis()
plt.ylabel('Depth (m)', fontsize=14, fontweight='bold')
cb = plt.colorbar()
plt.title(r'$\rm T(^{\circ}C)$', fontsize=14, fontweight='bold')
plt.show()



keyboard

## # Pickle data
## df_temp.to_pickle('historical_temp_1948-2017.pkl') # 
## df_sal.to_pickle('historical_sal_1948-2017.pkl') # 
## print(' -> Done!')


# Compute climatology
df_temp_may = df_temp.loc[df_temp.index.month==5]
df_temp_june = df_temp.loc[df_temp.index.month==6]
df_temp_july = df_temp.loc[df_temp.index.month==7]
df_concat = pd.concat((df_temp_may, df_temp_june, df_temp_july))
df_all = df_concat.resample('As').mean()
#df_all.to_pickle('temp_summer_1948-2017.pkl') # 



## --- CIL core --- ## 
fig = plt.figure(1)
plt.clf()
plt.plot(df_temp_may.index, df_temp_may.min(axis=1), '.')
plt.plot(df_temp_june.index, df_temp_june.min(axis=1), '.')
plt.plot(df_temp_july.index, df_temp_july.min(axis=1), '.')
plt.plot(df_all.index, df_all.min(axis=1).rolling(5,center=True).mean(), 'k-', linewidth=3)

plt.legend(['May', 'June', 'July', '5y moving ave'], fontsize=15)
plt.ylabel(r'$T_{min}$ in monthly mean profile ($^{\circ}$C)', fontsize=15, fontweight='bold')
plt.xlabel('Year', fontsize=15, fontweight='bold')
plt.title('CIL core temperature', fontsize=15, fontweight='bold')
plt.xlim([pd.Timestamp('1948-01-01'), pd.Timestamp('2017-01-01')])
plt.ylim([-2, 1])
plt.xticks(pd.date_range('1950-01-01', periods=7, freq='10Y'))
plt.grid('on')


fig.set_size_inches(w=9,h=6)
fig_name = 'CIL_core_1948-2018.png'
fig.set_dpi(300)
fig.savefig(fig_name)



keyboard

## --- No. of cast per year --- ##
years = np.arange(1912, 2019)
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

fig = plt.figure(2)
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

fig.set_size_inches(w=9,h=6)
fig_name = 'no_profile_histo.png'
fig.set_dpi(300)
fig.savefig(fig_name)


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

## Tcontour-seasonal
fig = plt.figure(3)
plt.clf()
v = np.arange(-2,10)
plt.contourf(df_temp_season.index, df_temp_season.columns, df_temp_season.T, v, cmap=plt.cm.RdBu_r)  
plt.xlim([pd.Timestamp('1948-01-01'), pd.Timestamp('2017-01-01')])
plt.xticks(pd.date_range('1950-01-01', periods=7, freq='10Y'))
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



## Salinity timeseries
fig = plt.figure(5)
plt.clf()
plt.plot(df_sal_season.index, df_sal_season.iloc[:,df_sal_season.columns<200].mean(axis=1), '-')
#plt.plot(df_sal_year.index, df_sal_year.iloc[:,df_sal_year.columns<200].mean(axis=1), '-r')
A = df_sal_season.iloc[:,df_sal_season.columns<200].mean(axis=1)
#A = df_sal_season.iloc[:,(df_sal_season.columns>75) & (df_sal_season.columns<150)].mean(axis=1)
plt.plot(df_sal_season.index, A.rolling(20, center=True).mean(), 'r')
df_sal_season.to_pickle('salinity_1948-2017.pkl')

## plt.plot(df_sal_season.index, df_sal_season.iloc[:,((df_sal_season.columns>150)&(df_sal_season.columns<300))].mean(axis=1), '-')
## #plt.plot(df_sal_year.index, df_sal_year.iloc[:,((df_sal_season.columns>200)&(df_sal_season.columns<300))].mean(axis=1), '-r')
## B =  df_sal_season.iloc[:,((df_sal_season.columns>150)&(df_sal_season.columns<300))].mean(axis=1)
## plt.plot(df_sal_season.index, B.rolling(24).mean(), 'r')


plt.ylabel(r'$\rm <S>_{0-200m}$', fontsize=15, fontweight='bold')
plt.xlabel('Year', fontsize=15, fontweight='bold')
#plt.title(r'$<S>_{0-300m}$', fontsize=15, fontweight='bold')
plt.xlim([pd.Timestamp('1948-01-01'), pd.Timestamp('2017-01-01')])
plt.xticks(pd.date_range('1950-01-01', periods=7, freq='10Y'))
plt.ylim([31.75, 33.50])
plt.grid('on')

fig.set_size_inches(w=9,h=6)
fig_name = 'sal_0-200m_1948-2017.png'
fig.set_dpi(300)
fig.savefig(fig_name)



