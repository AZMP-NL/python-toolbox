# -*- coding: utf-8 -*-
'''
Station 27 analysis for CSAS/NAFO ResDocs. This scripts compute:

- vertically average T-S
- CIL T, coreT and thickness
- 

** see azmp_stn27.py and azmp_stn27_explore.py for more options and ways to explore the dataset


Frederic.Cyr@dfo-mpo.gc.ca - July 2019

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

## font = {'family' : 'normal',
##         'weight' : 'bold',
##         'size'   : 14}
## plt.rc('font', **font)


def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)

## ---- Some custom parameters ---- ##
year_clim = [1981, 2010]
current_year = 2018
XLIM = [datetime.date(1945, 1, 1), datetime.date(2020, 12, 31)]
french_months = ['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D']

# Load pickled data
df_temp = pd.read_pickle('S27_temperature_monthly.pkl')
df_sal = pd.read_pickle('S27_salinity_monthly.pkl')
df_temp = df_temp[df_temp.index.year>=1950]
df_sal = df_sal[df_sal.index.year>=1950]

df_temp_clim_period = df_temp[(df_temp.index.year>=year_clim[0]) & (df_temp.index.year<=year_clim[1])]
df_sal_clim_period = df_sal[(df_sal.index.year>=year_clim[0]) & (df_sal.index.year<=year_clim[1])]

# Compute density
Z = df_temp.columns
SP = df_sal.values
PT = df_temp.values
SA = gsw.SA_from_SP(SP, Z, -50, 47)
CT = gsw.CT_from_pt(SA, PT)
RHO = gsw.rho(SA, CT, Z)
SIG0 = gsw.sigma0(SA, CT)
df_rho = pd.DataFrame(RHO, index=df_temp.index, columns=df_temp.columns)
df_sig = pd.DataFrame(SIG0, index=df_temp.index, columns=df_temp.columns)

# Compute N2 and MLD
print('Compute N2 and MLD for every index (make take some time...)')
N2 = np.full((df_rho.index.size, df_rho.columns.size-1), np.nan)
MLD = np.full((df_rho.index.size), np.nan)
for i,idx in enumerate(df_rho.index):     
        N2_tmp, pmid = gsw.Nsquared(SA[i,:], CT[i,:], Z, 47)
        N2[i,:] =  N2_tmp
        N2_tmp[np.where((pmid<=10) | (pmid>=100))] = np.nan
        if ~np.isnan(N2_tmp).all():
                MLD[i] = pmid[np.nanargmax(N2_tmp)]
print('  Done!')
df_N2 = pd.DataFrame(N2, index=df_temp.index, columns=pmid)
MLD = pd.Series(MLD, index=df_temp.index)
MLD.to_pickle('S27_MLD_monthly.pkl')



## ---- 1.  Vertically averaged conditions ---- ##
# Temperature
my_ts = df_temp.mean(axis=1)
my_ts = my_ts.resample('As').mean()
# Temperature climatology
clim_vertical = df_temp_clim_period.mean(axis=1)
clim_vertical = clim_vertical.resample('As').mean()
clim_mean = clim_vertical.mean()
clim_std = clim_vertical.std()
# Anomaly
anom = (my_ts - clim_mean)
anom_std = anom / clim_std    
# Plot anomaly
df1 = anom_std[anom_std>0]
df2 = anom_std[anom_std<0]
fig = plt.figure(1)
fig.clf()
width = 200
p1 = plt.bar(df1.index, np.squeeze(df1.values), width, alpha=0.8, color='indianred', zorder=10)
p2 = plt.bar(df2.index, np.squeeze(df2.values), width, bottom=0, alpha=0.8, color='steelblue', zorder=10)
plt.fill_between([anom_std.index[0], anom_std.index[-1]], [-.5, -.5], [.5, .5], facecolor='gray', alpha=.2)
plt.ylabel('Normalized anomaly')
plt.title('Station 27 - Average temperature (0-176m)')
#plt.xlim(XLIM)
plt.ylim([-3, 3])
plt.grid()
# Save Figure
fig.set_size_inches(w=7,h=4)
fig_name = 's27_vert_temp_anomaly.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)
# Save French Figure
plt.ylabel(u'Anomalie standardizée')
plt.title(u'Station 27 - Température moyenne (0-176m)')
fig_name = 's27_vert_temp_anomalyFR.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)

# Salinity
my_ts = df_sal.mean(axis=1)
my_ts = my_ts.resample('As').mean()
# Salinity climatology
clim_vertical = df_sal_clim_period.mean(axis=1)
clim_vertical = clim_vertical.resample('As').mean()
clim_mean = clim_vertical.mean()
clim_std = clim_vertical.std()
# Anomaly
anom = (my_ts - clim_mean)
anom_std = anom / clim_std    
# Plot anomaly
df1 = anom_std[anom_std>0]
df2 = anom_std[anom_std<0]
fig = plt.figure(1)
fig.clf()
width = 200
p1 = plt.bar(df1.index, np.squeeze(df1.values), width, alpha=0.8, color='indianred', zorder=10)
p2 = plt.bar(df2.index, np.squeeze(df2.values), width, bottom=0, alpha=0.8, color='steelblue', zorder=10)
plt.fill_between([anom_std.index[0], anom_std.index[-1]], [-.5, -.5], [.5, .5], facecolor='gray', alpha=.2)
plt.ylabel('Normalized anomaly')
plt.title('Station 27 - Average salinity (0-176m)')
#plt.xlim(XLIM)
plt.ylim([-3, 3])
plt.grid()
# Save Figure
fig.set_size_inches(w=7,h=4)
fig_name = 's27_vert_sal_anomaly.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)
# Save French Figure
plt.ylabel(u'Anomalie standardizée')
plt.title(u'Station 27 - Salinité moyenne (0-176m)')
fig_name = 's27_vert_sal_anomalyFR.png'
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
cil_summer_stats.to_pickle('S27_CIL_summer_stats.pkl')

# **plot CIL T**
my_ts = cil_temp
my_clim = my_ts[(my_ts.index.year>=year_clim[0]) & (my_ts.index.year<=year_clim[1])]
# Temperature climatology
clim_mean = my_clim.mean()
clim_std = my_clim.std()
# Anomaly
anom = (my_ts - clim_mean)
anom_std = anom / clim_std    
# Plot anomaly
df1 = anom_std[anom_std>0]
df2 = anom_std[anom_std<0]
fig = plt.figure(1)
fig.clf()
width = 200
p1 = plt.bar(df1.index, np.squeeze(df1.values), width, alpha=0.8, color='indianred', zorder=10)
p2 = plt.bar(df2.index, np.squeeze(df2.values), width, bottom=0, alpha=0.8, color='steelblue', zorder=10)
plt.ylabel('Normalized anomaly')
plt.title('Station 27 - CIL mean temperature')
plt.ylim([-5, 5])
plt.grid()
# Save Figure
fig.set_size_inches(w=7,h=4)
fig_name = 's27_CILtemp_anomaly.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)    
# Save French Figure
plt.ylabel(u'Anomalie standardizée')
plt.title(u'Station 27 - Température moyenne de la CIF')
fig_name = 's27_CILtemp_anomalyFR.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)

# **plot core T**
my_ts = cil_core
my_clim = my_ts[(my_ts.index.year>=year_clim[0]) & (my_ts.index.year<=year_clim[1])]
# Temperature climatology
clim_mean = my_clim.mean()
clim_std = my_clim.std()
# Anomaly
anom = (my_ts - clim_mean)
anom_std = anom / clim_std    
# Plot anomaly
df1 = anom_std[anom_std>0]
df2 = anom_std[anom_std<0]
fig = plt.figure(1)
fig.clf()
width = 200
p1 = plt.bar(df1.index, np.squeeze(df1.values), width, alpha=0.8, color='indianred', zorder=10)
p2 = plt.bar(df2.index, np.squeeze(df2.values), width, bottom=0, alpha=0.8, color='steelblue', zorder=10)
plt.ylabel('Normalized anomaly')
plt.title('Station 27 - CIL core temperature')
plt.ylim([-5, 5])
plt.grid()
# Save Figure
fig.set_size_inches(w=7,h=4)
fig_name = 's27_CILcore_anomaly.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)    
# Save French Figure
plt.ylabel(u'Anomalie standardizée')
plt.title(u'Station 27 - Température du coeur de la CIF')
fig_name = 's27_CILcore_anomalyFR.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)

# **plot core depth**
my_ts = cil_coredepth
my_clim = my_ts[(my_ts.index.year>=year_clim[0]) & (my_ts.index.year<=year_clim[1])]
# Temperature climatology
clim_mean = my_clim.mean()
clim_std = my_clim.std()
# Anomaly
anom = (my_ts - clim_mean)
anom_std = anom / clim_std    
# Plot anomaly
df1 = anom_std[anom_std>0]
df2 = anom_std[anom_std<0]
fig = plt.figure(1)
fig.clf()
width = 200
p1 = plt.bar(df1.index, np.squeeze(df1.values), width, alpha=0.8, color='indianred', zorder=10)
p2 = plt.bar(df2.index, np.squeeze(df2.values), width, bottom=0, alpha=0.8, color='steelblue', zorder=10)
plt.ylabel('Normalized anomaly')
plt.title('Station 27 - CIL core depth')
plt.ylim([-4, 4])
plt.grid()
# Save Figure
fig.set_size_inches(w=7,h=4)
fig_name = 's27_CILcoredepth_anomaly.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)    
# Save French Figure
plt.ylabel(u'Anomalie standardizée')
plt.title(u'Station 27 - Profondeur du coeur de la CIF')
fig_name = 's27_CILcoredepth_anomalyFR.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)

# **plot thickness**
my_ts = cil_thickness
my_clim = my_ts[(my_ts.index.year>=year_clim[0]) & (my_ts.index.year<=year_clim[1])]
# Temperature climatology
clim_mean = my_clim.mean()
clim_std = my_clim.std()
# Anomaly
anom = (my_ts - clim_mean)
anom_std = anom / clim_std    
# Plot anomaly
df1 = anom_std[anom_std<0]
df2 = anom_std[anom_std>0]
fig = plt.figure(1)
fig.clf()
width = 200
p1 = plt.bar(df1.index, np.squeeze(df1.values), width, alpha=0.8, color='indianred', zorder=10)
p2 = plt.bar(df2.index, np.squeeze(df2.values), width, bottom=0, alpha=0.8, color='steelblue', zorder=10)
plt.ylabel('Normalized anomaly')
plt.title('Station 27 - CIL thickness')
#plt.xlim(XLIM)
plt.ylim([-4, 4])
plt.grid()
# Save Figure
fig.set_size_inches(w=7,h=4)
fig_name = 's27_CILthickness_anomaly.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)
# Save French Figure
plt.ylabel(u'Anomalie standardizée')
plt.title(u'Station 27 - Épaisseur de la CIF')
fig_name = 's27_CILthickness_anomalyFR.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)

# **plot CIL timeseries**
fig, host = plt.subplots()
fig.subplots_adjust(right=0.75)
par1 = host.twinx()
par2 = host.twinx()
# Offset the right spine of par2.  The ticks and label have already been
# placed on the right by twinx above.
par2.spines["right"].set_position(("axes", 1.1))
# Having been created by twinx, par2 has its frame off, so the line of its
# detached spine is invisible.  First, activate the frame but make the patch
# and spines invisible.
make_patch_spines_invisible(par2)
# Second, show the right spine.
par2.spines["right"].set_visible(True)
p1, = host.plot(cil_core.index.year, cil_core.rolling(5).mean().values, "steelblue", label="CIL core temperature", zorder=10)
p2, = par1.plot(cil_coredepth.index.year, cil_coredepth.rolling(5).mean().values, "darkorange", label="CIL core depth", zorder=10)
p3, = par2.plot(cil_thickness.index.year, cil_thickness.rolling(5).mean().values, "indianred", label="CIL thickness", zorder=10)
#host.set_xlim(-2, 1)
host.set_ylim(-2, 1)
par1.set_ylim(90, 150)
par2.set_ylim(20, 140)
host.set_xlabel("")
host.set_ylabel(r"$\rm ^{\circ}C$")
par1.set_ylabel("Z (m)")
par2.set_ylabel("H (m)")
host.yaxis.label.set_color(p1.get_color())
par1.yaxis.label.set_color(p2.get_color())
par2.yaxis.label.set_color(p3.get_color())
tkw = dict(size=4, width=1.5)
host.tick_params(axis='y', colors=p1.get_color(), **tkw)
par1.tick_params(axis='y', colors=p2.get_color(), **tkw)
par2.tick_params(axis='y', colors=p3.get_color(), **tkw)
host.tick_params(axis='x', **tkw)
lines = [p1, p2, p3]
host.legend(lines, [l.get_label() for l in lines], loc=2)
plt.grid()
fig.set_size_inches(w=14,h=8)
fig_name = 's27_CIL_stats.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)


## ---- 3. Stratification ---- ##
rho_5m = df_rho.iloc[:,int(np.squeeze(np.where(df_rho.columns==5)))]
rho_50m = df_rho.iloc[:,int(np.squeeze(np.where(df_rho.columns==50)))]
strat_monthly = (rho_50m-rho_5m)/45
strat_monthly.to_pickle('S27_stratif_monthly.pkl')
    
# A) Barplot
strat = strat_monthly.resample('As').mean()
strat_clim = strat[(strat.index.year>=year_clim[0]) & (strat.index.year<=year_clim[1])]
clim_mean = strat_clim.mean()
clim_std = strat_clim.std()
# Anomaly
anom = (strat - clim_mean)*1000
anom_std = anom / (clim_std*1000)    
# Plot anomaly
df1 = anom_std[anom_std>0]
df2 = anom_std[anom_std<0]
fig = plt.figure(1)
fig.clf()
width = 200
p1 = plt.bar(df1.index, np.squeeze(df1.values), width, alpha=0.8, color='indianred', zorder=10)
p2 = plt.bar(df2.index, np.squeeze(df2.values), width, bottom=0, alpha=0.8, color='steelblue', zorder=10)
plt.ylabel('Normalized anomaly')
plt.title('Station 27 - Stratification (5-50m)')
#plt.xlim(XLIM)
plt.ylim([-3.2, 3.2])
plt.grid()
# Save Figure
fig.set_size_inches(w=7,h=4)
fig_name = 's27_stratif_bar.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)
# Save French Figure
plt.ylabel(u'Anomalie standardizée')
fig_name = 's27_stratif_bar_FR.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)

# B) timeseries
fig = plt.figure(1)
plt.clf()
anom.plot(color='gray')
anom.rolling(5).mean().plot(color='k', linewidth=3, linestyle='--')
#ticks = ax.xaxis.get_ticklocs() 
#plt.fill_between([ticks[0]-1, ticks[-1]+1], [-clim_std*1000, -clim_std*1000], [clim_std*1000, clim_std*1000], facecolor='gray', alpha=.2)
plt.grid()
plt.ylabel(r'Stratification anomaly $\rm (g\,m^{-4})$')
plt.xlabel(' ')
# Save Figure
fig.set_size_inches(w=7,h=4)
fig_name = 's27_stratif_plot.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)
# Same in French
plt.ylabel(r'Anomalie de stratification $\rm (g\,m^{-4})$')
fig_name = 's27_stratif_plotFR.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)

# C) Monthly bar plot for current year
strat_clim_period = strat_monthly[(strat_monthly.index.year>=year_clim[0]) & (strat_monthly.index.year<=year_clim[1])]
strat_monthly_stack = strat_clim_period.groupby([(strat_clim_period.index.year),(strat_clim_period.index.month)]).mean()
strat_monthly_clim = strat_monthly_stack.mean(level=1)
strat_monthly_std = strat_monthly_stack.std(level=1)
strat_current_year = strat_monthly[strat_monthly.index.year==current_year]
year_index = strat_current_year.index # backup index
strat_current_year.index=strat_monthly_std.index # reset index
monthly_anom = strat_current_year - strat_monthly_clim
monthly_std_anom = monthly_anom/strat_monthly_std
monthly_std_anom.index = year_index.strftime('%b') # replace index (by text)
monthly_anom.index = year_index.strftime('%b') # replace index (by text)
# plot
ind = np.arange(len(monthly_anom.keys()))  # the x locations for the groups
width = 0.35  # the width of the bars
fig, ax = plt.subplots()
rects1 = ax.bar(ind - width/2, strat_monthly_clim.values*1000, width, yerr=strat_monthly_std.values*.5*1000,
                label='1981-2010', color='lightblue', zorder=1)
rects2 = ax.bar(ind + width/2, np.squeeze(strat_current_year.values*1000), width, yerr=None,
                label=str(current_year), color='steelblue', zorder=1)
# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel(r'$\rm \frac{\Delta \rho}{\Delta z} (g\,m^{-4})$')
ax.set_title('Station 27 - Stratification')
ax.set_xticks(ind)
ax.set_xticklabels(monthly_anom.index)
ax.legend()
ax.yaxis.grid() # horizontal lines
# Save Figure
fig.set_size_inches(w=6,h=3)
fig_name = 's27_stratif_monthly.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)
# Save French Figure
ax.set_xticklabels(french_months)
fig_name = 's27_stratif_monthly_FR.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)



## ---- 4. MLD ---- ##
# A) Barplot
mld_monthly = MLD
mld = mld_monthly.resample('As').mean()
mld_clim = mld[(mld.index.year>=year_clim[0]) & (mld.index.year<=year_clim[1])]
clim_mean = mld_clim.mean()
clim_std = mld_clim.std()
# Anomaly
anom = (mld - clim_mean)
anom_std = anom / clim_std    
# Plot anomaly
df1 = anom_std[anom_std>0]
df2 = anom_std[anom_std<0]
fig = plt.figure(1)
fig.clf()
width = 200
p1 = plt.bar(df1.index, np.squeeze(df1.values), width, alpha=0.8, color='indianred', zorder=10)
p2 = plt.bar(df2.index, np.squeeze(df2.values), width, bottom=0, alpha=0.8, color='steelblue', zorder=10)
plt.ylabel(r'Normalized anomaly')
plt.title('Station 27 - Mixed layer depth')
#plt.xlim(XLIM)
plt.ylim([-3, 3])
plt.grid()
fig.set_size_inches(w=7,h=4)
fig_name = 's27_mld_bar.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)
# French Figure
plt.ylabel(u'Anomalie standardizée')
plt.title(u'Station 27 - Profondeur couche de mélange')
fig_name = 's27_mld_barFR.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)


# B) timeseries
fig = plt.figure(1)
plt.clf()
anom.plot(color='gray')
anom.rolling(5).mean().plot(color='k', linewidth=3, linestyle='--')
plt.grid()
plt.ylabel(r'MLD anomaly (m)')
plt.xlabel(' ')
fig.set_size_inches(w=7,h=4)
fig_name = 's27_mld_plot.png'
fig.savefig(fig_name, dpi=200)
os.system('convert -trim ' + fig_name + ' ' + fig_name)

# French Figure
plt.ylabel(r'Anomalie de la PCM (m)')
fig_name = 's27_mld_plotFR.png'
fig.savefig(fig_name, dpi=200)
os.system('convert -trim ' + fig_name + ' ' + fig_name)


# C) Monthly bar plot for current year
mld_clim_period = mld_monthly[(mld_monthly.index.year>=year_clim[0]) & (mld_monthly.index.year<=year_clim[1])]
mld_monthly_stack = mld_clim_period.groupby([(mld_clim_period.index.year),(mld_clim_period.index.month)]).mean()
mld_monthly_clim = mld_monthly_stack.mean(level=1)
mld_monthly_std = mld_monthly_stack.std(level=1)
mld_current_year = mld_monthly[mld_monthly.index.year==current_year]
year_index = mld_current_year.index # backup index
mld_current_year.index=mld_monthly_std.index # reset index
monthly_anom = mld_current_year - mld_monthly_clim
monthly_std_anom = monthly_anom/mld_monthly_std
monthly_std_anom.index = year_index.strftime('%b') # replace index (by text)
monthly_anom.index = year_index.strftime('%b') # replace index (by text)
# plot
ind = np.arange(len(monthly_anom.keys()))  # the x locations for the groups
width = 0.35  # the width of the bars
fig, ax = plt.subplots()
rects1 = ax.bar(ind - width/2, mld_monthly_clim.values, width, yerr=mld_monthly_std.values*.5,
                label='1981-2010', color='lightblue', zorder=1)
rects2 = ax.bar(ind + width/2, np.squeeze(mld_current_year.values), width, yerr=None,
                label=str(current_year), color='steelblue', zorder=1)
# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel(r'MLD (m)')
ax.set_title('Station 27 - Mixed layer depth')
ax.set_xticks(ind)
ax.set_xticklabels(monthly_anom.index)
ax.legend()
ax.yaxis.grid() # horizontal lines
#plt.ylim([0, 330])
fig.set_size_inches(w=6,h=3)
fig_name = 's27_mld_monthly.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)
# French Figure
ax.set_ylabel(r'PCM (m)')
ax.set_title(u'Station 27 - Profondeur couche de mélange')
ax.set_xticklabels(french_months)
fig_name = 's27_mld_monthly_FR.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)


## ---- Montage figures ---- ##
# T-S English
os.system('montage  s27_vert_temp_anomaly.png s27_vert_sal_anomaly.png -tile 1x2 -geometry +1+25  -background white  s27_TS_subplots.png')
# T-S French
os.system('montage  s27_vert_temp_anomalyFR.png s27_vert_sal_anomalyFR.png -tile 1x2 -geometry +1+25  -background white  s27_TS_subplotsFR.png')

# CIL English
os.system('montage  s27_CILtemp_anomaly.png s27_CILthickness_anomaly.png s27_CILcore_anomaly.png s27_CILcoredepth_anomaly.png  -tile 2x2 -geometry +20+25  -background white  s27_CIL_subplots.png')
# CIL French
os.system('montage  s27_CILtemp_anomalyFR.png s27_CILthickness_anomalyFR.png s27_CILcore_anomalyFR.png s27_CILcoredepth_anomalyFR.png  -tile 2x2 -geometry +20+25  -background white  s27_CIL_subplotsFR.png')



