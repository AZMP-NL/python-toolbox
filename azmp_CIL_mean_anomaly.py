# -*- coding: utf-8 -*-
'''

df_CIL_*_summer.pkl are generated by azmp_CIL_stats.py

* run in /home/cyrf0006/AZMP/state_reports/sections_plots/CIL



'''

import numpy as  np
import matplotlib.pyplot as plt
import pandas as pd
import os

# Adjust fontsize/weight
font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 14}
plt.rc('font', **font)

clim_year = [1991, 2020]
years = [1950, 2023]

baditp_SI = [1998,2022]
baditp_BB = []
baditp_FC = []


#Mark years where stn will be replaced with stn_man
rplstn_SI = np.concatenate((np.arange(1950,1992),np.array([1992,2000,2003,2004])))
rplstn_BB = np.concatenate((np.arange(1950,1992),np.array([1992,1997,1998,1999,2000,2001,2002,2022])))
rplstn_FC = np.concatenate((np.arange(1950,1992),np.array([1992,1997,1998,2000,2001,2002,2006,2007,2022])))

#### ---- Load the data and compute anomalies ---- ####
df_SI = pd.read_pickle('operation_files/df_CIL_SI_summer.pkl').astype('float')
df_BB = pd.read_pickle('operation_files/df_CIL_BB_summer.pkl').astype('float')
df_FC = pd.read_pickle('operation_files/df_CIL_FC_summer.pkl').astype('float')

# Set problem years equal to nan
df_SI['vol_stn'].loc[baditp_SI] = np.nan
df_SI['core_stn'].loc[baditp_SI] = np.nan
df_BB['vol_stn'].loc[baditp_BB] = np.nan
df_BB['core_stn'].loc[baditp_BB] = np.nan
df_FC['vol_stn'].loc[baditp_FC] = np.nan
df_FC['core_stn'].loc[baditp_FC] = np.nan


# Set replace years to interpolate
df_SI['vol_stn'].loc[rplstn_SI] = df_SI['vol_itp'].loc[rplstn_SI]
#df_SI['core_stn'].loc[rplstn_SI] = df_SI['core_itp'].loc[rplstn_SI]
df_BB['vol_stn'].loc[rplstn_BB] = df_BB['vol_itp'].loc[rplstn_BB]
#df_BB['core_stn'].loc[rplstn_BB] = df_BB['core_itp'].loc[rplstn_BB]
df_FC['vol_stn'].loc[rplstn_FC] = df_FC['vol_itp'].loc[rplstn_FC]
#df_FC['core_stn'].loc[rplstn_FC] = df_FC['core_itp'].loc[rplstn_FC]


#
df_years = df_SI[(df_SI.index>=years[0]) & (df_SI.index<=years[1])]
df_clim = df_SI[(df_SI.index>=clim_year[0]) & (df_SI.index<=clim_year[1])]
std_anom_SI = (df_years - df_clim.mean()) / df_clim.std()

df_years = df_BB[(df_BB.index>=years[0]) & (df_BB.index<=years[1])]
df_clim = df_BB[(df_BB.index>=clim_year[0]) & (df_BB.index<=clim_year[1])]
std_anom_BB = (df_years - df_clim.mean()) / df_clim.std()

df_years = df_FC[(df_FC.index>=years[0]) & (df_FC.index<=years[1])]
df_clim = df_FC[(df_FC.index>=clim_year[0]) & (df_FC.index<=clim_year[1])]
std_anom_FC = (df_years - df_clim.mean()) / df_clim.std()

# average quatities
std_core = pd.concat([std_anom_SI.core_stn, std_anom_BB.core_stn, std_anom_FC.core_stn], axis=1).mean(axis=1)
std_vol = pd.concat([std_anom_SI.vol_stn, std_anom_BB.vol_stn, std_anom_FC.vol_stn], axis=1).mean(axis=1)
#std_core.index = pd.to_datetime(std_core.index, format='%Y')
#std_vol.index = pd.to_datetime(std_vol.index, format='%Y')


# Save for climate index
cil_section_index = pd.concat([std_core, std_vol], axis=1, keys=['core','volume'])
cil_section_index.to_pickle('section_cil_index.pkl')


## ---- plot index ---- ##
width=.7
fig = plt.figure(4)
fig.clf()
ax = plt.subplot2grid((2, 1), (0, 0))
sign=std_core>0
std_core.plot(kind='bar', color=sign.map({False: 'steelblue', True: 'indianred'}), width = width, zorder=10)
n = 5 # xtick every n years
ticks = ax.xaxis.get_ticklocs()
ticklabels = [l.get_text() for l in ax.xaxis.get_ticklabels()]
ax.xaxis.set_ticks(ticks[::n])
ax.xaxis.set_ticklabels(ticklabels[::n])
plt.fill_between([ticks[0], ticks[-1]], [-.5, -.5], [.5, .5], facecolor='gray', alpha=.2)
plt.ylabel('Standardized Anomaly', weight='bold', fontsize=14)
plt.xlabel(' ')
plt.grid()
plt.ylim([-2.5,2.5])
plt.title('CIL core - Sections SI, BB & FC', weight='bold', fontsize=14)
ax.tick_params(labelbottom=False)

ax2 = plt.subplot2grid((2, 1), (1, 0))
sign=std_vol>0
std_vol.plot(kind='bar', color=sign.map({True: 'steelblue', False: 'indianred'}), width = width, zorder=10)
n = 5 # xtick every n years
ticks = ax2.xaxis.get_ticklocs()
ticklabels = [l.get_text() for l in ax2.xaxis.get_ticklabels()]
ax2.xaxis.set_ticks(ticks[::n])
ax2.xaxis.set_ticklabels(ticklabels[::n])
plt.fill_between([ticks[0], ticks[-1]], [-.5, -.5], [.5, .5], facecolor='gray', alpha=.2)
plt.ylabel('Standardized Anomaly', weight='bold', fontsize=14)
plt.xlabel(' ')
plt.grid()
plt.ylim([-2.5,2.5])
plt.title('CIL area - Sections SI, BB & FC', weight='bold', fontsize=14)
ax2.tick_params(axis='x', rotation=0)

fig_name = 'section_CIL_anomaly.png'
fig.set_size_inches(w=15,h=7)
fig.savefig(fig_name, dpi=200)
os.system('convert -trim ' + fig_name + ' ' + fig_name)


## ---- Figure in French ---- ##
fig = plt.figure(4)
fig.clf()
ax = plt.subplot2grid((2, 1), (0, 0))
sign=std_core>0
std_core.plot(kind='bar', color=sign.map({False: 'steelblue', True: 'indianred'}), width = width, zorder=10)
n = 5 # xtick every n years
ticks = ax.xaxis.get_ticklocs()
ticklabels = [l.get_text() for l in ax.xaxis.get_ticklabels()]
ax.xaxis.set_ticks(ticks[::n])
ax.xaxis.set_ticklabels(ticklabels[::n])
plt.ylabel('Anomalie normalisée', weight='bold', fontsize=14)
plt.xlabel(' ')
plt.grid()
plt.ylim([-2.5,2.5])
plt.title('Coeur de la CIF - Sections SI, BB & FC', weight='bold', fontsize=14)
ax.tick_params(labelbottom=False)

ax2 = plt.subplot2grid((2, 1), (1, 0))
sign=std_vol>0
std_vol.plot(kind='bar', color=sign.map({True: 'steelblue', False: 'indianred'}), width = width, zorder=10)
n = 5 # xtick every n years
ticks = ax2.xaxis.get_ticklocs()
ticklabels = [l.get_text() for l in ax2.xaxis.get_ticklabels()]
ax2.xaxis.set_ticks(ticks[::n])
ax2.xaxis.set_ticklabels(ticklabels[::n])
plt.ylabel('Anomalie normalisée', weight='bold', fontsize=14)
plt.xlabel(' ')
plt.grid()
plt.ylim([-2.5,2.5])
plt.title('Aire de la CIF - Sections SI, BB & FC', weight='bold', fontsize=14)
ax2.tick_params(axis='x', rotation=0)

fig_name = 'section_CIL_anomaly_FR.png'
fig.set_size_inches(w=15,h=7)
fig.savefig(fig_name, dpi=200)
os.system('convert -trim ' + fig_name + ' ' + fig_name)

## ---- CIL area for ms_climate index ---- ##
# Build the colormap
from matplotlib.colors import from_levels_and_colors
YlGn = plt.cm.YlGn(np.linspace(0,1, num=12))
YlGn = YlGn[1:,]
cmap, norm = from_levels_and_colors(np.arange(0,10), YlGn, extend='both') 

std_vol_std = pd.concat([std_anom_SI.vol_stn, std_anom_BB.vol_stn, std_anom_FC.vol_stn], axis=1) / 3
std_vol_std.columns = ['SI', 'BB', 'FC']
n = 5 # xtick every n years
#ax = std_vol_std.plot(kind='bar', stacked=True, cmap='tab10')
ax = std_vol_std.plot(kind='bar', stacked=True, cmap=cmap)
ticks = ax.xaxis.get_ticklocs()
ticklabels = [l.get_text() for l in ax.xaxis.get_ticklabels()]
ax.xaxis.set_ticks(ticks[::n])
ax.xaxis.set_ticklabels(ticklabels[::n])
plt.fill_between([ticks[0], ticks[-1]], [-.5, -.5], [.5, .5], facecolor='gray', alpha=.2)
plt.grid('on')
ax.set_ylabel(r'Normalized anomaly')
ax.set_title('CIL volume')
plt.ylim([-3,3])
fig = ax.get_figure()
fig.set_size_inches(w=12,h=8)
fig_name = 'CIL_volume.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)

# Save in French
ax.set_ylabel(r'Anomalie normalisée')
ax.set_title('Anomalies du volume de la CIF')
fig_name = 'CIL_volume_FR.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)

## ---- CIL area for ms_climate index [with scorecards] ---- ##
# Build the colormap - Stack
from matplotlib.colors import from_levels_and_colors
YlGn = plt.cm.YlGn(np.linspace(0,1, num=12))
YlGn = YlGn[1:,]
cmap_stack, norm_stack = from_levels_and_colors(np.arange(0,10), YlGn, extend='both') 
# Build the colormap - Scorecard
vmin = -3.49
vmax = 3.49
midpoint = 0
levels = np.linspace(vmin, vmax, 15)
midp = np.mean(np.c_[levels[:-1], levels[1:]], axis=1)
colvals = np.interp(midp, [vmin, midpoint, vmax], [-1, 0., 1])
normal = plt.Normalize(-3.49, 3.49)
reds = plt.cm.Reds(np.linspace(0,1, num=7))
blues = plt.cm.Blues_r(np.linspace(0,1, num=7))
whites = [(1,1,1,1)]*2
colors = np.vstack((blues[0:-1,:], whites, reds[1:,:]))
colors = np.concatenate([[colors[0,:]], colors, [colors[-1,:]]], 0)
cmap, norm = from_levels_and_colors(levels, colors, extend='both')
cmap_r, norm_r = from_levels_and_colors(levels, np.flipud(colors), extend='both')

# Concat and normalize
std_vol = pd.concat([std_anom_SI.vol_stn, std_anom_BB.vol_stn, std_anom_FC.vol_stn], axis=1)
#std_vol_std = pd.concat([std_anom_SI.vol_itp, std_anom_BB.vol_itp, std_anom_FC.vol_itp], axis=1) / 3
std_vol_norm = std_vol.divide((3 - std_vol.isna().sum(axis=1)).values, axis=0)
std_vol_norm.columns = ['SI', 'BB', 'FC']

# Plot
fig, ax = plt.subplots(nrows=1, ncols=1)
n = 5 # xtick every n years
ax = std_vol_norm.plot(kind='bar', stacked=True, cmap=cmap_stack)
ticks = ax.xaxis.get_ticklocs()
ticklabels = [l.get_text() for l in ax.xaxis.get_ticklabels()]
ax.xaxis.set_ticks(ticks[::n])
ax.xaxis.set_ticklabels(ticklabels[::n])
plt.fill_between([ticks[0]-1, ticks[-1]+1], [-.5, -.5], [.5, .5], facecolor='gray', alpha=.2)
plt.grid('on')
ax.set_xlabel(r'')
ax.set_ylabel(r'Normalized anomaly')
ax.set_title('CIL area')
plt.ylim([-3.3,2.6])

colors = cmap(normal(np.nansum(std_vol_norm.values*-1, axis=1)))
cell_text = np.nansum(std_vol_norm.values, axis=1).round(1).astype('str')
the_table = ax.table(cellText=[np.nansum(std_vol_norm.values, axis=1).round(1)],
        rowLabels=['CIL area subindex'],
        colLabels=None,
        cellColours = [colors],
        cellLoc = 'center', rowLoc = 'center',
        loc='bottom', bbox=[0, -0.14, 1, 0.05])
the_table.auto_set_font_size (False)
the_table.set_fontsize(9)

for key, cell in the_table.get_celld().items():
    if key[1] == -1:
        cell.set_linewidth(0)
        cell.set_fontsize(10)
    else:
        cell._text.set_rotation(90)
        
fig = ax.get_figure()
fig.set_size_inches(w=13,h=9.5)
fig_name = 'CIL_volume_climateindex.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)

## In French ##
# Plot
fig, ax = plt.subplots(nrows=1, ncols=1)
n = 5 # xtick every n years
ax = std_vol_norm.plot(kind='bar', stacked=True, cmap=cmap_stack)
ticks = ax.xaxis.get_ticklocs()
ticklabels = [l.get_text() for l in ax.xaxis.get_ticklabels()]
ax.xaxis.set_ticks(ticks[::n])
ax.xaxis.set_ticklabels(ticklabels[::n])
plt.fill_between([ticks[0]-1, ticks[-1]+1], [-.5, -.5], [.5, .5], facecolor='gray', alpha=.2)
plt.grid('on')
ax.set_xlabel(r'')
ax.set_ylabel(r'Anomalie normalisée')
ax.set_title('Aire de la CIF')
plt.ylim([-3.3,2.6])

colors = cmap(normal(np.nansum(std_vol_norm.values*-1, axis=1)))
cell_text = np.nansum(std_vol_norm.values, axis=1).round(1).astype('str')
the_table = ax.table(cellText=[np.nansum(std_vol_norm.values, axis=1).round(1)],
        rowLabels=['sous-indice aire CIF'],
        colLabels=None,
        cellColours = [colors],
        cellLoc = 'center', rowLoc = 'center',
        loc='bottom', bbox=[0, -0.14, 1, 0.05])
the_table.auto_set_font_size (False)
the_table.set_fontsize(9)

for key, cell in the_table.get_celld().items():
    if key[1] == -1:
        cell.set_linewidth(0)
        cell.set_fontsize(10)
    else:
        cell._text.set_rotation(90)
        
fig = ax.get_figure()
fig.set_size_inches(w=13,h=9.5)
fig_name = 'CIL_volume_climateindex_FR.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)



## ## ---- CIL area for ms_climate index V2 (RED/BLUE) ---- ##
## vol_stack = pd.concat([std_anom_SI.vol_itp[std_anom_SI.vol_itp<0], std_anom_SI.vol_itp[std_anom_SI.vol_itp>0],
##                          std_anom_BB.vol_itp[std_anom_BB.vol_itp<0], std_anom_BB.vol_itp[std_anom_BB.vol_itp>0],
##                          std_anom_FC.vol_itp[std_anom_FC.vol_itp<0], std_anom_FC.vol_itp[std_anom_FC.vol_itp>0]], axis=1) / 3


## fig, ax = plt.subplots(nrows=1, ncols=1)
## n = 5 # xtick every n years
## sign=vol_stack>0
## vol_stack.plot(kind='bar', stacked=True, color=['lightcoral', 'skyblue', 'indianred', 'steelblue', 'firebrick', 'royalblue'], width = width, zorder=10, ax=ax, legend=False)
## ticks = ax.xaxis.get_ticklocs()
## ticklabels = [l.get_text() for l in ax.xaxis.get_ticklabels()]
## ax.xaxis.set_ticks(ticks[::n])
## ax.xaxis.set_ticklabels(ticklabels[::n])
## ax.set_ylabel(r'Normalized anomaly')
## ax.set_title('CIL volume')
## plt.fill_between([ticks[0]-1, ticks[-1]+1], [-.5, -.5], [.5, .5], facecolor='gray', alpha=.2)
## plt.grid()
## plt.ylim([-3.5,3.5])
## # Custom legend
## import matplotlib.lines as mlines
## legend_elements = [mlines.Line2D([],[], marker='s',linestyle='None', color='lightcoral', markersize=12, label='\nSI'),
##                     mlines.Line2D([],[], marker='s',linestyle='None', color='skyblue', markersize=12, label=' '),
##                     mlines.Line2D([],[], marker='s',linestyle='None', color='indianred', markersize=12, label='\nBB'),
##                     mlines.Line2D([],[], marker='s',linestyle='None', color='steelblue', markersize=12, label=' '),
##                     mlines.Line2D([],[], marker='s',linestyle='None', color='firebrick', markersize=12, label='\nFC'),
##                     mlines.Line2D([],[], marker='s',linestyle='None', color='royalblue', markersize=12, label=' ')]
## ax.legend(handles=legend_elements)
      
## fig.set_size_inches(w=15,h=7)
## fig_name = 'CIL_volume_stacked.png'
## fig.savefig(fig_name, dpi=300)
## os.system('convert -trim -bordercolor White -border 10x10 ' + fig_name + ' ' + fig_name)

