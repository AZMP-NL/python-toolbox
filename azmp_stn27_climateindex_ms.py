# -*- coding: utf-8 -*-
'''
Station 27 analysis for manuscript on climate index

Frederic.Cyr@dfo-mpo.gc.ca - October 2020

Run in '/home/cyrf0006/research/ms_climate_index'

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
import gsw

# Adjust fontsize/weight
font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 14}
plt.rc('font', **font)

def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)

## ---- Some custom parameters ---- ##
#year_clim = [1981, 2010]
year_clim = [1991, 2020]
current_year = 2021
XLIM = [datetime.date(1945, 1, 1), datetime.date(2020, 12, 31)]
french_months = ['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D']
months = mdates.MonthLocator()  # every month
month_fmt = mdates.DateFormatter('%b')
Vsig = np.append(np.arange(21,27), 26.5)
Vsig = np.append(Vsig, 27)

# Load pickled data
df_temp = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/stn27/S27_temperature_monthly.pkl')
df_sal = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/stn27/S27_salinity_monthly.pkl')
# start in 1950
df_temp = df_temp[df_temp.index.year>=1951]
df_sal = df_sal[df_sal.index.year>=1951]
# Flag years with less than 8 months
df_temp[df_temp.index.year==1980] = np.nan
df_temp[df_temp.index.year==1981] = np.nan
df_temp[df_temp.index.year==2020] = np.nan
df_sal[df_sal.index.year==1980] = np.nan
df_sal[df_sal.index.year==1981] = np.nan
df_sal[df_sal.index.year==2020] = np.nan
# Climatology period
df_temp_clim_period = df_temp[(df_temp.index.year>=year_clim[0]) & (df_temp.index.year<=year_clim[1])]
df_sal_clim_period = df_sal[(df_sal.index.year>=year_clim[0]) & (df_sal.index.year<=year_clim[1])]



## 1. ---- Monthly clim ---- ##
monthly_climT = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/stn27/S27_temperature_monthly_clim.pkl')
monthly_climS = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/stn27/S27_salinity_monthly_clim.pkl')
# add 13th months
monthly_climT.loc[pd.to_datetime('1901-01-01')] = monthly_climT.loc[pd.to_datetime('1900-01-01')]
monthly_climS.loc[pd.to_datetime('1901-01-01')] = monthly_climS.loc[pd.to_datetime('1900-01-01')]
# Calculate density
SA = gsw.SA_from_SP(monthly_climS, monthly_climS.columns, -52, 47)
CT = gsw.CT_from_t(SA, monthly_climT, monthly_climT.columns)
SIG = gsw.sigma0(SA, CT)
monthly_climSIG = pd.DataFrame(SIG, index=monthly_climT.index, columns=monthly_climT.columns)

# plot T
fig, ax = plt.subplots(nrows=1, ncols=1)
CMAP = cmocean.tools.lighten(cmocean.cm.thermal, .9)
V = np.arange(-2, 14, 1)
c = plt.contourf(monthly_climT.index, monthly_climT.columns, monthly_climT.values.T, V, extend='max', cmap=CMAP)
cc = plt.contour(monthly_climSIG.index, monthly_climSIG.columns, monthly_climSIG.values.T, Vsig, colors='dimgray')
plt.clabel(cc, inline=1, fontsize=10, colors='dimgray', fmt='%1.1f')
plt.ylim([0, 175])
plt.ylabel('Depth (m)', fontsize=14)
plt.ylim([0, 175])
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
ccc = ax.contour(monthly_climT.index, monthly_climT.columns, monthly_climT.T, [0,], colors='k', linewidths=3)
cb.set_label(r'$\rm T(^{\circ}C)$', fontsize=14, fontweight='normal')
ax.xaxis.label.set_visible(False)
# Save Figure
fig.set_size_inches(w=13, h=6)
outfile_clim = 's27_temp_clim.png'
fig.savefig(outfile_clim, dpi=200)
os.system('convert -trim ' + outfile_clim + ' ' + outfile_clim)


# plot S
fig, ax = plt.subplots(nrows=1, ncols=1)
V = np.arange(30, 33.5, .25)
CMAP = cmocean.cm.haline
c = plt.contourf(monthly_climS.index, monthly_climS.columns, monthly_climS.values.T, V, extend='both', cmap=CMAP)
cc = plt.contour(monthly_climSIG.index, monthly_climSIG.columns, monthly_climSIG.values.T, Vsig, colors='dimgray')
plt.clabel(cc, inline=1, fontsize=10, colors='dimgray', fmt='%1.1f')
plt.ylim([0, 175])
plt.ylabel('Depth (m)', fontsize=14)
plt.ylim([0, 175])
ax.invert_yaxis()
# format date ticks
ax.xaxis.set_major_locator(MonthLocator())
ax.xaxis.set_minor_locator(MonthLocator(bymonthday=15))
ax.xaxis.set_major_formatter(NullFormatter())
ax.xaxis.set_minor_formatter(NullFormatter())
# format the ticks
ax.xaxis.set_major_locator(months)
ax.xaxis.set_major_formatter(month_fmt)
cax = fig.add_axes([0.91, .15, 0.01, 0.7])
cb = plt.colorbar(c, cax=cax, orientation='vertical')
ccc = ax.contour(monthly_climS.index, monthly_climS.columns, monthly_climS.T, [0,], colors='k', linewidths=3)
cb.set_label(r'S', fontsize=14, fontweight='normal')
# Save Figure
fig.set_size_inches(w=13, h=6)
outfile_clim = 's27_sal_clim.png'
fig.savefig(outfile_clim, dpi=200)
os.system('convert -trim ' + outfile_clim + ' ' + outfile_clim)

# Convert to subplot
os.system('montage  s27_temp_clim.png s27_sal_clim.png -tile 1x2 -geometry +10+10  -background white  s27_TSclim.png')

## ---- 2.  Vertically averaged conditions ---- ##
# Temperature
my_ts = df_temp.mean(axis=1)
# New from March 2020 - Annual anomaly = mean of monthly anomalies
ts_stack = my_ts.groupby([(my_ts.index.year),(my_ts.index.month)]).mean()
ts_unstack = ts_stack.unstack()
ts_clim_period = my_ts[(my_ts.index.year>=year_clim[0]) & (my_ts.index.year<=year_clim[1])]
ts_monthly_stack = ts_clim_period.groupby([(ts_clim_period.index.year),(ts_clim_period.index.month)]).mean()
ts_monthly_clim = ts_monthly_stack.mean(level=1)
ts_monthly_std = ts_monthly_stack.std(level=1)
monthly_anom = ts_unstack - ts_monthly_clim 
monthly_stdanom = (ts_unstack - ts_monthly_clim) /  ts_monthly_std
anom_std = monthly_stdanom.mean(axis=1)
anom_std.index = pd.to_datetime(anom_std.index, format='%Y')
anom = monthly_anom.mean(axis=1) 
anom.index = pd.to_datetime(anom.index, format='%Y')
# Plot anomaly
df1 = anom_std[anom_std>0]
df2 = anom_std[anom_std<0]
fig, ax = plt.subplots()
width = 200
p1 = plt.bar(df1.index, np.squeeze(df1.values), width, alpha=0.8, color='indianred', zorder=10)
p2 = plt.bar(df2.index, np.squeeze(df2.values), width, bottom=0, alpha=0.8, color='steelblue', zorder=10)
plt.fill_between([anom_std.index[0], anom_std.index[-1]], [-.5, -.5], [.5, .5], facecolor='gray', alpha=.2)
plt.ylabel('Normalized anomaly')
plt.title('Average temperature (0-176m)')
#plt.xlim(XLIM)
plt.ylim([-2.1, 2.1])
plt.grid()
ax.tick_params(labelbottom=False)
# Save Figure
fig.set_size_inches(w=8,h=4)
fig_name = 's27_vert_temp_anomaly.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)

# Salinity
my_ts = df_sal.mean(axis=1)
# New from March 2020 - Annual anomaly = mean of monthly anomalies
ts_stack = my_ts.groupby([(my_ts.index.year),(my_ts.index.month)]).mean()
ts_unstack = ts_stack.unstack()
ts_clim_period = my_ts[(my_ts.index.year>=year_clim[0]) & (my_ts.index.year<=year_clim[1])]
ts_monthly_stack = ts_clim_period.groupby([(ts_clim_period.index.year),(ts_clim_period.index.month)]).mean()
ts_monthly_clim = ts_monthly_stack.mean(level=1)
ts_monthly_std = ts_monthly_stack.std(level=1)
monthly_anom = ts_unstack - ts_monthly_clim 
monthly_stdanom = (ts_unstack - ts_monthly_clim) /  ts_monthly_std
anom_std = monthly_stdanom.mean(axis=1) 
anom_std.index = pd.to_datetime(anom_std.index, format='%Y')
anom = monthly_anom.mean(axis=1) 
anom.index = pd.to_datetime(anom.index, format='%Y')
# Plot anomaly
df1 = anom_std[anom_std>0]
df2 = anom_std[anom_std<0]
fig, ax = plt.subplots()
width = 200
p1 = plt.bar(df1.index, np.squeeze(df1.values), width, alpha=0.8, color='indianred', zorder=10)
p2 = plt.bar(df2.index, np.squeeze(df2.values), width, bottom=0, alpha=0.8, color='steelblue', zorder=10)
plt.fill_between([anom_std.index[0], anom_std.index[-1]], [-.5, -.5], [.5, .5], facecolor='gray', alpha=.2)
plt.ylabel('Normalized anomaly')
plt.title('Average salinity (0-176m)')
#plt.xlim(XLIM)
plt.ylim([-2.5, 2.5])
plt.grid()
ax.tick_params(labelbottom=False)
# Save Figure
fig.set_size_inches(w=8,h=4)
fig_name = 's27_vert_sal_anomaly.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)


## ---- 2.  CIL summer statistics ---- ##
# Select summer (MJJ)
df_temp_summer = df_temp[(df_temp.index.month>=6) & (df_temp.index.month<=8)]
df_temp_summer = df_temp_summer.resample('As').mean()
df_temp_summer = df_temp_summer[df_temp_summer.columns[(df_temp_summer.columns>=10)]] # remove top 10m to be sure

cil_temp = np.full(df_temp_summer.index.shape, np.nan)
cil_core = np.full(df_temp_summer.index.shape, np.nan)
cil_coredepth = np.full(df_temp_summer.index.shape, np.nan)
cil_thickness = np.full(df_temp_summer.index.shape, np.nan)

zbin = df_temp_summer.columns[2] - df_temp_summer.columns[1]

for idx, YEAR in enumerate(df_temp_summer.index.year):
    # Get single year
    tmp =  df_temp_summer.iloc[idx]
    # CIL core
    cil_core[idx] = np.nanmin(tmp.values)
    # CIL core depth
    z_idx = np.where(tmp==np.nanmin(tmp.values))
    cil_coredepth[idx] = np.mean(tmp.index[z_idx])
    # CIL thickness
    CIL_idxs = np.squeeze(np.where(tmp<=0))
    cil_thickness[idx] =  np.size(CIL_idxs)*zbin
    # CIL temp
    cil_temp[idx] =  tmp[CIL_idxs].mean()
    
# Convert to pandas series    
cil_temp = pd.Series(cil_temp, index=df_temp_summer.index)
cil_core = pd.Series(cil_core, index=df_temp_summer.index)
cil_coredepth = pd.Series(cil_coredepth, index=df_temp_summer.index)
cil_thickness = pd.Series(cil_thickness, index=df_temp_summer.index)
cil_summer_stats = pd.concat([cil_temp, cil_core, cil_coredepth, cil_thickness], axis=1, keys=['CIL temp', 'CIL core T', 'CIL core depth', 'CIL thickness'])

# Calculate anomalies
my_clim = cil_summer_stats[(cil_summer_stats.index.year>=year_clim[0]) & (cil_summer_stats.index.year<=year_clim[1])]
clim_mean = my_clim.mean()
clim_std = my_clim.std()
cil_summer_anom = (cil_summer_stats - clim_mean) / clim_std

# plot CIL index (HERE TO CHOOSE WHICH WE USE)
#anom_std = cil_summer_anom[['CIL temp', 'CIL core T']].mean(axis=1)
anom_std = cil_summer_anom[['CIL core T']].mean(axis=1)
df1 = anom_std[anom_std>0]
df2 = anom_std[anom_std<0]
fig, ax = plt.subplots()
width = 200
p1 = plt.bar(df1.index, np.squeeze(df1.values), width, alpha=0.8, color='indianred', zorder=10)
p2 = plt.bar(df2.index, np.squeeze(df2.values), width, bottom=0, alpha=0.8, color='steelblue', zorder=10)
plt.fill_between([anom_std.index[0], anom_std.index[-1]], [-.5, -.5], [.5, .5], facecolor='gray', alpha=.2)
plt.ylabel('Normalized anomaly')
plt.title('CIL core temperature')
plt.ylim([-2, 4.5])
plt.grid()
# Save Figure
fig.set_size_inches(w=8,h=4)
fig_name = 's27_CILtemp_anomaly.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)    

# Convert to subplot
os.system('montage s27_vert_temp_anomaly.png s27_vert_sal_anomaly.png  s27_CILtemp_anomaly.png -tile 1x3 -geometry +10+30  -background white  s27_anom.png')
