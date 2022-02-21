# -*- coding: utf-8 -*-
'''
To generate figure for Redfish assessment

Frederic.Cyr@dfo-mpo.gc.ca
Feb. 2022
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

clim_year=[1991, 2020]


#### ------------- Load data ---------------- ####
# Peter's SST
df_sst = pd.read_csv('/home/cyrf0006/data/AZMP/Galbraith_data/AZMP_SST_Seasonal.dat', delimiter=r"\s+", index_col='#Yr', header=27)
df_sst.index.name = 'year'
df_sst = df_sst.replace(-99.00, np.nan)
#df_sst = df_sst[df_sst.index<=yearf]
df_sst.loc[1981] = np.nan      
df_sst.loc[1980] = np.nan      
df_sst.sort_index(inplace=True)

# select subset
df_sst = df_sst[['3P', '4Vn', '4Vs']]
# climatology
df_clim = df_sst[(df_sst.index>=clim_year[0]) & (df_sst.index<=clim_year[1])]
std_anom = (df_sst-df_clim.mean(axis=0))/df_clim.std(axis=0)


#### ------------- Stacked bar plot with scorecard ---------------- ####
# Concat and normalize
stack = std_anom.copy()
stack_norm = stack.divide((3 - stack.isna().sum(axis=1)).values, axis=0)
stack_norm.columns = ['3P', '4Vn', '4Vs']

# Build the colormap - Stack
from matplotlib.colors import from_levels_and_colors
YlGn = plt.cm.YlGn(np.linspace(0,1, num=12))
YlGn = YlGn[4:,]
cmap_stack, norm_stack = from_levels_and_colors(np.arange(0,7), YlGn, extend='both') 
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

# Plot
fig, ax = plt.subplots(nrows=1, ncols=1)
n = 5 # xtick every n years
ax = stack_norm.plot(kind='bar', stacked=True, cmap=cmap_stack)
ticks = ax.xaxis.get_ticklocs()
ticklabels = [l.get_text() for l in ax.xaxis.get_ticklabels()]
ax.xaxis.set_ticks(ticks[::n])
ax.xaxis.set_ticklabels(ticklabels[::n])
plt.fill_between([ticks[0]-.5, ticks[-1]+.5], [-.5, -.5], [.5, .5], facecolor='gray', alpha=.2)
plt.grid('on')
ax.set_xlabel(r'')
ax.set_ylabel(r'Normalized anomaly')
ax.set_title('Sea Surface Temperature')
plt.ylim([-2,2.25])


colors = cmap(normal(np.nansum(stack_norm.values, axis=1)))
cell_text = np.nansum(stack_norm.values, axis=1).round(1).astype('str')
the_table = ax.table(cellText=[np.nansum(stack_norm.values, axis=1).round(1)],
        rowLabels=['norm. anomaly'],
        colLabels=None,
        cellColours = [colors],
        cellLoc = 'center', rowLoc = 'center',
        loc='bottom', bbox=[0, -0.14, 1, 0.05])
the_table.auto_set_font_size (False)
the_table.set_fontsize(11)

for key, cell in the_table.get_celld().items():
    if key[1] == -1:
        cell.set_linewidth(0)
        cell.set_fontsize(10)
    else:
        cell._text.set_rotation(90)
        
fig = ax.get_figure()
fig.set_size_inches(w=13,h=9.5)
fig_name = 'SST_anomalies_redfish.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim -bordercolor White -border 10x10 ' + fig_name + ' ' + fig_name)


# Plot French Figure
stack_normFR = stack_norm.copy()
fig, ax = plt.subplots(nrows=1, ncols=1)
n = 5 # xtick every n years
ax = stack_normFR.plot(kind='bar', stacked=True, cmap=cmap_stack)
ticks = ax.xaxis.get_ticklocs()
ticklabels = [l.get_text() for l in ax.xaxis.get_ticklabels()]
ax.xaxis.set_ticks(ticks[::n])
ax.xaxis.set_ticklabels(ticklabels[::n])
plt.fill_between([ticks[0]-.5, ticks[-1]+.5], [-.5, -.5], [.5, .5], facecolor='gray', alpha=.2)
plt.grid('on')
ax.set_xlabel(r'')
ax.set_ylabel(r'Anomalie Normalisée')
ax.set_title('Température de surface')
plt.ylim([-2,2.25])

colors = cmap(normal(np.nansum(stack_normFR.values, axis=1)))
cell_text = np.nansum(stack_normFR.values, axis=1).round(1).astype('str')
the_table = ax.table(cellText=[np.nansum(stack_normFR.values, axis=1).round(1)],
        rowLabels=['sub-indice temp fond'],
        colLabels=None,
        cellColours = [colors],
        cellLoc = 'center', rowLoc = 'center',
        loc='bottom', bbox=[0, -0.14, 1, 0.05])
the_table.auto_set_font_size (False)
the_table.set_fontsize(11)

for key, cell in the_table.get_celld().items():
    if key[1] == -1:
        cell.set_linewidth(0)
        cell.set_fontsize(10)
    else:
        cell._text.set_rotation(90)
        
fig = ax.get_figure()
fig.set_size_inches(w=13,h=9.5)
fig_name = 'SST_anomalies_redfish_FR.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim -bordercolor White -border 10x10 ' + fig_name + ' ' + fig_name)
