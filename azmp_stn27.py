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
import matplotlib.dates as mdates
from matplotlib.ticker import NullFormatter
from matplotlib.dates import MonthLocator, DateFormatter
import cmocean

## font = {'family' : 'normal',
##         'weight' : 'bold',
##         'size'   : 14}
## plt.rc('font', **font)

## ---- Some custom parameters ---- ##
s27 = [47.55,-52.59]
dc = .1
year_clim = [1981, 2010]
current_year = 2018
variable = 'sigma-t'
use_viking = True
XLIM = [datetime.date(current_year, 1, 1), datetime.date(current_year, 12, 31)]

# Derived parameter
if variable == 'temperature':
    V = np.arange(-1, 16, 1)
    #Vanom = np.linspace(-5.5, 5.5, 12)
    #Vanom = np.array([-6, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 6])
    Vanom = np.linspace(-5.5, 5.5, 23)
    Vanom_ticks = np.array([-5, -3, -1, 1, 3, 5])
    Vanom = np.delete(Vanom, np.where(Vanom==0))
    CMAP = cmocean.tools.lighten(cmocean.cm.thermal, .9)
    #CMAP = cmocean.cm.thermal
elif variable == 'salinity':
    V = np.arange(30, 33.5, .25)
    Vanom = np.linspace(-1, 1, 21) 
    Vanom = np.delete(Vanom, np.where(Vanom==0))
    Vanom_ticks = np.linspace(-1, 1, 11) 
    CMAP = cmocean.cm.haline
else:
    V = 20
    Vanom = 20
    Vanom_ticks = np.linspace(-1, 1, 11) 
    CMAP = cmocean.cm.thermal

## months = mdates.MonthLocator()  # every month
## month_fmt = mdates.DateFormatter('%b')

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
df_hydro = da.to_pandas()

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
    df = pd.concat([df_hydro, df_vik]) 
    #df = df.interpolate(axis=1, method='linear', limit_area='inside') # doesnt seems to work as expected
    df = df.interpolate(axis=1).where(df.bfill(axis=1).notnull()) # here interpolate and leave NaNs.
else:
    df = df_hydro.copy()
    
## ---- 1. Climatologies ---- ##
# Monthly average (15th of the month) +  pickle for further analysis
df_monthly = df.resample('MS', loffset=pd.Timedelta(14, 'd')).mean()
df_monthly.to_pickle('S27_' + variable + '_monthly.pkl')
# Select years for climato
df_clim_period = df_monthly[(df_monthly.index.year>=year_clim[0]) & (df_monthly.index.year<=year_clim[1])]
# Monthly clim
monthly_clim = df_clim_period.groupby(df_clim_period.index.month).mean()
# set index (year 1900)
monthly_clim.index = pd.to_datetime(monthly_clim.index.values, format='%m')
monthly_clim.to_pickle('S27_' + variable + '_monthly_clim.pkl')

# Weekly average (previous version of climatology)
#df_weekly = df.resample('W').mean()
# Monthly clim
#weekly_clim_raw = df_weekly.groupby(df_weekly.index.week).mean()
# set index (year 1900)
#weekly_clim.index = pd.to_datetime(weekly_clim.index.values, format='%W')

# Weekly clim (upsample monthly clim to weekly)
weekly_clim = monthly_clim.resample('W').mean().interpolate(method='linear') 
weekly_clim.to_pickle('S27_' + variable + '_weekly_clim.pkl')

# Update climatology index to current year
weekly_clim.index = pd.to_datetime('2018-' +  weekly_clim.index.month.astype(np.str) + '-' + weekly_clim.index.day.astype(np.str))
#weekly_clim.dropna(how='all', axis=1, inplace=True)

# plot
fig, ax = plt.subplots(nrows=1, ncols=1)
c = plt.contourf(weekly_clim.index, weekly_clim.columns, weekly_clim.values.T, V, extend='both', cmap=CMAP)
#cc = plt.contour(weekly_clim.index, weekly_clim.columns, weekly_clim.values.T, V, colors='k')
#plt.plot(weekly_clim.index, np.repeat(0, df_temp.index.size), '|k', markersize=20)    
plt.ylim([0, 175])
plt.xlim([XLIM[0], XLIM[1]])
#plt.clabel(cc, inline=1, fontsize=10, colors='k', fmt='%1.1f')
#plt.grid('on')
#plt.xlabel('Time', fontsize=15, fontweight='bold')
plt.ylabel('Depth (m)', fontsize=15, fontweight='bold')
plt.ylim([0, 175])
plt.title('Climatology')
ax.invert_yaxis()
# format date ticks
ax.xaxis.set_major_locator(MonthLocator())
ax.xaxis.set_minor_locator(MonthLocator(bymonthday=15))
ax.xaxis.set_major_formatter(NullFormatter())
ax.xaxis.set_minor_formatter(NullFormatter())
# format the ticks
#ax.xaxis.set_major_locator(months)
#ax.xaxis.set_major_formatter(month_fmt)
cax = fig.add_axes([0.91, .15, 0.01, 0.7])
cb = plt.colorbar(c, cax=cax, orientation='vertical')
if variable == 'temperature':
    ccc = ax.contour(weekly_clim.index, weekly_clim.columns, weekly_clim.values.T, [0,], colors='k', linewidths=3)
    cb.set_label(r'$\rm T(^{\circ}C)$', fontsize=12, fontweight='normal')
elif variable == 'salinity':
    cb.set_label(r'S', fontsize=12, fontweight='normal')
ax.xaxis.label.set_visible(False)
# Save Figure
fig.set_size_inches(w=12, h=6)
outfile_clim = 's27_' + variable + '_clim.png'
fig.savefig(outfile_clim, dpi=200)
os.system('convert -trim ' + outfile_clim + ' ' + outfile_clim)

# Save French Figure
ax.set_ylabel('Profondeur (m)', fontsize=15, fontweight='bold')
fig.set_size_inches(w=12, h=6)
outfile_climFR = 's27_' + variable + '_clim_FR.png'
fig.savefig(outfile_climFR, dpi=200)
os.system('convert -trim ' + outfile_climFR + ' ' + outfile_climFR)


## ---- 2. Year average and anomaly ---- ##
df_year = df[df.index.year==current_year]
df_weekly = df_year.resample('W').mean().interpolate(method='linear') 
#df_weekly.dropna(how='all', axis=1, inplace=True)
anom = df_weekly - weekly_clim

# plot current year
fig, ax = plt.subplots(nrows=1, ncols=1)
c = plt.contourf(df_weekly.index, df_weekly.columns, df_weekly.values.T, V, extend='both', cmap=CMAP)
plt.plot(df_year.index, np.repeat(0, df_year.index.size), '|k', markersize=20)    
plt.ylim([0, 175])
plt.xlim([XLIM[0], XLIM[1]])
#plt.clabel(cc, inline=1, fontsize=10, colors='k', fmt='%1.1f')
#plt.grid('on')
#plt.xlabel('Time', fontsize=15, fontweight='bold')
plt.ylabel('Depth (m)', fontsize=15, fontweight='bold')
plt.title(str(current_year) + ' observations')
plt.ylim([0, 175])
ax.invert_yaxis()
# format date ticks
ax.xaxis.set_major_locator(MonthLocator())
ax.xaxis.set_minor_locator(MonthLocator(bymonthday=15))
ax.xaxis.set_major_formatter(NullFormatter())
ax.xaxis.set_minor_formatter(NullFormatter())
cax = fig.add_axes([0.91, .15, 0.01, 0.7])
cb = plt.colorbar(c, cax=cax, orientation='vertical')
if variable == 'temperature':
    ccc = ax.contour(df_weekly.index, df_weekly.columns, df_weekly.values.T, [0,], colors='k', linewidths=3)
    cb.set_label(r'$\rm T(^{\circ}C)$', fontsize=12, fontweight='normal')
elif variable == 'salinity':
    cb.set_label(r'S', fontsize=12, fontweight='normal')
ax.xaxis.label.set_visible(False)
# Save Figure
fig.set_size_inches(w=12, h=6)
outfile_year = 's27_' + variable + '_' + str(current_year) + '.png'
fig.savefig(outfile_year, dpi=200)
os.system('convert -trim ' + outfile_year + ' ' + outfile_year)

# Save French Figure
ax.set_ylabel('Profondeur (m)', fontsize=15, fontweight='bold')
fig.set_size_inches(w=12, h=6)
outfile_yearFR = 's27_' + variable + '_' + str(current_year) + '_FR.png'
fig.savefig(outfile_yearFR, dpi=200)
os.system('convert -trim ' + outfile_yearFR + ' ' + outfile_yearFR)

# plot anomaly
fig, ax = plt.subplots(nrows=1, ncols=1)
c = plt.contourf(anom.index, anom.columns, anom.values.T, Vanom, extend='both', cmap=cmocean.cm.balance)
plt.ylim([0, 175])
plt.xlim([XLIM[0], XLIM[1]])
plt.ylabel('Depth (m)', fontsize=15, fontweight='bold')
plt.title(str(current_year) + ' anomaly')
plt.ylim([0, 175])
ax.invert_yaxis()
# format date ticks
ax.xaxis.set_major_locator(MonthLocator())
ax.xaxis.set_minor_locator(MonthLocator(bymonthday=15))
ax.xaxis.set_major_formatter(NullFormatter())
ax.xaxis.set_minor_formatter(DateFormatter('%b'))
cax = fig.add_axes([0.91, .15, 0.01, 0.7])
cb = plt.colorbar(c, cax=cax, ticks=Vanom_ticks, orientation='vertical')
if variable == 'temperature':
    cb.set_label(r'$\rm T(^{\circ}C)$', fontsize=12, fontweight='normal')
elif variable == 'salinity':
    cb.set_label(r'S', fontsize=12, fontweight='normal')
# Save Figure
fig.set_size_inches(w=12, h=6)
outfile_anom = 's27_' + variable + '_anom_' + str(current_year) + '.png'
fig.savefig(outfile_anom, dpi=200)
os.system('convert -trim ' + outfile_anom + ' ' + outfile_anom)

# Save French Figure
ax.xaxis.set_minor_formatter(DateFormatter('%m'))
ax.set_ylabel('Profondeur (m)', fontsize=15, fontweight='bold')
fig.set_size_inches(w=12, h=6)
outfile_anomFR = 's27_' + variable + '_anom_' + str(current_year) + '_FR.png'
fig.savefig(outfile_anomFR, dpi=200)
os.system('convert -trim ' + outfile_anomFR + ' ' + outfile_anomFR)

# Convert to a subplot
os.system('montage ' + outfile_clim + ' ' +  outfile_year + ' ' + outfile_anom  + ' -tile 1x3 -geometry +10+10  -background white  s27_' + variable + '_subplot_' + str(current_year) + '.png') 
os.system('montage ' + outfile_climFR + ' ' +  outfile_yearFR + ' ' + outfile_anomFR  + ' -tile 1x3 -geometry +10+10  -background white  s27_' + variable + '_subplot_' + str(current_year) + '_FR.png') 

## ---- Station occupation plot ---- ##
da_occu = ds['time']
df_occu = da_occu.to_pandas()
df_occu = df_occu[df_occu.index.year>=1945]

# plot 1 - weekly occupations only
fig, ax = plt.subplots(nrows=1, ncols=1)
plt.clf()
plt.plot(df_occu.index.year, df_occu.index.weekofyear, '.k')
plt.ylabel('Week of year')
plt.xlabel('Station 27 occupation')
# Save Figure
fig.set_size_inches(w=12, h=6)
outfile_occu = 's27_occupation.png'
fig.savefig(outfile_occu, dpi=200)
os.system('convert -trim ' + outfile_occu + ' ' + outfile_occu)

# plot 2 - weekly occupations + no of week per year
W = df_occu.resample('w').count()
W = W[W>0]  
W = W.resample('Y').count()
fig = plt.figure()
# ax1
ax1 = plt.subplot2grid((2, 1), (0, 0))
plt.bar(W.index.year, W.values)
ax1.xaxis.label.set_visible(False)
ax1.tick_params(labelbottom='off')
ax1.set_ylabel('No. of weekly occupations')
plt.title('Station 27 occupation')
# ax2
ax2 = plt.subplot2grid((2, 1), (1, 0))
plt.plot(df_occu.index.year, df_occu.index.weekofyear, '.k')
ax2.set_ylabel('week of year')
# Save Figure
fig.set_size_inches(w=6, h=7)
outfile_occu = 's27_occupation_stats.png'
fig.savefig(outfile_occu, dpi=200)
os.system('convert -trim ' + outfile_occu + ' ' + outfile_occu)

