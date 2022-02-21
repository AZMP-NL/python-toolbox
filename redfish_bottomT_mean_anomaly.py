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


#### ------------- Load data NL ---------------- ####
# 3Ps
infile = 'stats_3Ps_spring.pkl'
df = pd.read_pickle(infile)
df.index = pd.to_datetime(df.index) # update index to datetime
# Flag bad years (no or weak sampling):
bad_years = np.array([1980, 1981, 1985, 1986, 1987, 1988, 1989, 1990, 1991, 1992, 2006])
for i in bad_years:
    df[df.index.year==i]=np.nan
df['area_colder0'] = df['area_colder0']/1000 # In 1000km
df_clim = df[(df.index.year>=clim_year[0]) & (df.index.year<=clim_year[1])]
std_anom = (df-df_clim.mean(axis=0))/df_clim.std(axis=0)
std_anom3Ps = std_anom[['Tmean', 'Tmean_sha100']]

# Take only mean T
df_3Ps = std_anom3Ps.Tmean
df_3Ps.index = df_3Ps.index.year



#### ------------- Load data SS ---------------- ####
# 4Vn
df_4vn = pd.read_csv('/home/cyrf0006/data/Hebert_timeseries/summerGroundfishBottomTemperature_4Vn.dat', delimiter=r",", index_col='year', header=10)
df_4vn = df_4vn.anomaly
df_4vn = df_4vn.replace('    NA', np.nan)


# 4Vs
df_4vs = pd.read_csv('/home/cyrf0006/data/Hebert_timeseries/summerGroundfishBottomTemperature_4Vs.dat', delimiter=r",", index_col='year', header=10)
df_4vs = df_4vs.anomaly
df_4vs = df_4vs.replace('   NA', np.nan)

#### ------------- Stacked bar plot with scorecard ---------------- ####
# Concat and normalize
bottomT_stack = pd.concat([df_3Ps, df_4vn.astype('float'), df_4vs.astype('float')], axis=1)
bottomT_stack_norm = bottomT_stack.divide((3 - bottomT_stack.isna().sum(axis=1)).values, axis=0)
bottomT_stack_norm.columns = ['3Ps - Spring', '4Vn - Summer', '4Vs - Summer']

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
ax = bottomT_stack_norm.plot(kind='bar', stacked=True, cmap=cmap_stack)
ticks = ax.xaxis.get_ticklocs()
ticklabels = [l.get_text() for l in ax.xaxis.get_ticklabels()]
ax.xaxis.set_ticks(ticks[::n])
ax.xaxis.set_ticklabels(ticklabels[::n])
plt.fill_between([ticks[0], ticks[-1]], [-.5, -.5], [.5, .5], facecolor='gray', alpha=.2)
plt.grid('on')
ax.set_xlabel(r'')
ax.set_ylabel(r'Normalized anomaly')
ax.set_title('Bottom Temperature')
plt.ylim([-2,2.25])

colors = cmap(normal(np.nansum(bottomT_stack_norm.values, axis=1)))
cell_text = np.nansum(bottomT_stack_norm.values, axis=1).round(1).astype('str')
the_table = ax.table(cellText=[np.nansum(bottomT_stack_norm.values, axis=1).round(1)],
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
fig_name = 'bottomT_anomalies_redfish.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim -bordercolor White -border 10x10 ' + fig_name + ' ' + fig_name)

bottomT_stack_norm.columns = ['3Ps - Spring', '4Vn - Summer', '4Vs - Summer']

# Plot French Figure
bottomT_stack_normFR = bottomT_stack_norm.rename(columns={'3Ps - Spring':'3Ps - Printemps','4Vn - Summer':'4Vn -  Été','4Vs - Summer':'4Vs -  Été'})
fig, ax = plt.subplots(nrows=1, ncols=1)
n = 5 # xtick every n years
ax = bottomT_stack_normFR.plot(kind='bar', stacked=True, cmap=cmap_stack)
ticks = ax.xaxis.get_ticklocs()
ticklabels = [l.get_text() for l in ax.xaxis.get_ticklabels()]
ax.xaxis.set_ticks(ticks[::n])
ax.xaxis.set_ticklabels(ticklabels[::n])
plt.fill_between([ticks[0], ticks[-1]], [-.5, -.5], [.5, .5], facecolor='gray', alpha=.2)
plt.grid('on')
ax.set_xlabel(r'')
ax.set_ylabel(r'Anomalie Normalisée')
ax.set_title('Température fond')
plt.ylim([-2,2.25])

colors = cmap(normal(np.nansum(bottomT_stack_normFR.values, axis=1)))
cell_text = np.nansum(bottomT_stack_normFR.values, axis=1).round(1).astype('str')
the_table = ax.table(cellText=[np.nansum(bottomT_stack_normFR.values, axis=1).round(1)],
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
fig_name = 'bottomT_anomalies_redfish_FR.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim -bordercolor White -border 10x10 ' + fig_name + ' ' + fig_name)
