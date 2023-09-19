# -*- coding: utf-8 -*-
'''
To generate AZMP score cards for bottom temperature

Uses pickled object generated by azmp_bottom_stats.py

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

width = 0.7
clim_year = [1991, 2020]
years = [1980, 2023]

#### ------------- For fall ---------------- ####
# 0.
infile = 'bottom_temp_stats/stats_2H_fall.pkl'
df = pd.read_pickle(infile)
df.index = pd.to_datetime(df.index) # update index to datetime
# Flag bad years (no or weak sampling):
bad_years = np.array([1980, 1982, 1984, 1985, 1986, 1987, 1988, 1989, 1990, 1992, 1993, 1994, 1995, 1996, 2000, 2002, 2003, 2005, 2007, 2009, 2022])
for i in bad_years:
    df[df.index==i]=np.nan
df['area_colder0'] = df['area_colder0']/1000 # In 1000km
df_clim = df[(df.index.year>=clim_year[0]) & (df.index.year<=clim_year[1])]
std_anom = (df-df_clim.mean(axis=0))/df_clim.std(axis=0)
std_anom2H = std_anom[['Tmean', 'Tmean_sha200']]

# 1.
infile = 'bottom_temp_stats/stats_2J_fall.pkl'
df = pd.read_pickle(infile)
df.index = pd.to_datetime(df.index) # update index to datetime
# Flag bad years (no or weak sampling):
bad_years = np.array([1995, 2022])
for i in bad_years:
    df[df.index==i]=np.nan
df['area_colder0'] = df['area_colder0']/1000 # In 1000km
df_clim = df[(df.index.year>=clim_year[0]) & (df.index.year<=clim_year[1])]
std_anom = (df-df_clim.mean(axis=0))/df_clim.std(axis=0)
std_anom2J = std_anom[['Tmean', 'Tmean_sha200']]

# 2.
infile = 'bottom_temp_stats/stats_3K_fall.pkl'
df = pd.read_pickle(infile)
df.index = pd.to_datetime(df.index) # update index to datetime
df['area_colder0'] = df['area_colder0']/1000 # In 1000km
df_clim = df[(df.index.year>=clim_year[0]) & (df.index.year<=clim_year[1])]
std_anom = (df-df_clim.mean(axis=0))/df_clim.std(axis=0)
std_anom3K = std_anom[['Tmean', 'Tmean_sha300']]

# 3.
infile = 'bottom_temp_stats/stats_3LNO_fall.pkl'
df = pd.read_pickle(infile)
df.index = pd.to_datetime(df.index) # update index to datetime
# Flag bad years (no or weak sampling):
bad_years = np.array([2021])
for i in bad_years:
    df[df.index==i]=np.nan
df['area_colder0'] = df['area_colder0']/1000 # In 1000km
df_clim = df[(df.index.year>=clim_year[0]) & (df.index.year<=clim_year[1])]
std_anom = (df-df_clim.mean(axis=0))/df_clim.std(axis=0)
std_anom3LNO = std_anom[['Tmean', 'Tmean_sha100']]

# 4. plot all merged together
#std_anom_all = pd.concat([std_anom2J, std_anom3K, std_anom3LNO], axis=1)
std_anom_all = pd.concat([std_anom2H, std_anom2J, std_anom3K, std_anom3LNO], axis=1)
df = std_anom_all.mean(axis=1)
df = df[(df.index.year>=years[0]) & (df.index.year<=years[-1])]
df.index = df.index.year
# Save for climate index
df.to_pickle('bottomT_index_fall.pkl')

fig, ax = plt.subplots(nrows=1, ncols=1)
n = 5 # xtick every n years
sign=df>0
ax = df.plot(kind='bar', color=sign.map({True: 'indianred', False: 'steelblue'}), width = width)
ticks = ax.xaxis.get_ticklocs()
ticklabels = [l.get_text() for l in ax.xaxis.get_ticklabels()]
ax.xaxis.set_ticks(ticks[::n])
ax.xaxis.set_ticklabels(ticklabels[::n])
plt.ylabel('Normalized Anomaly', weight='bold', fontsize=14)
plt.fill_between([ticks[0]-1, ticks[-1]+1], [-.5, -.5], [.5, .5], facecolor='gray', alpha=.2)
#plt.xlabel('Year')
plt.title('Mean bottom temperature anomaly for 2HJ3KLNO - Fall', weight='bold', fontsize=14)
plt.grid()
plt.ylim([-2.5,2.5])
ax.yaxis.set_ticks(np.arange(-2.5, 3, .5))
#plt.gca().set_xlim([pd.to_datetime('1979-01-01'), pd.to_datetime('2018-01-01')]) 
fig.set_size_inches(w=15,h=7)
fig_name = 'mean_anomalies_fall.png'
fig.savefig(fig_name, dpi=200)
os.system('convert -trim mean_anomalies_fall.png mean_anomalies_fall.png')

# Save French Figure
plt.sca(ax)
plt.ylabel(u'Anomalie normalisée', weight='bold', fontsize=14)
plt.title(u'Température de fond moyenne pour 2HJ3KLNO - Automne', weight='bold', fontsize=14)
#fig.set_size_inches(w=15,h=7)
fig_name = 'mean_anomalies_fall_FR.png'
fig.savefig(fig_name, dpi=200)
os.system('convert -trim ' + fig_name + ' ' + fig_name)


#### ------------- For Spring ---------------- ####
# 1.
infile = 'bottom_temp_stats/stats_3LNO_spring.pkl'
df = pd.read_pickle(infile)
df.index = pd.to_datetime(df.index) # update index to datetime
# Flag bad years (no or weak sampling):
bad_years = np.array([2020, 2021])
for i in bad_years:
    df[df.index.year==i]=np.nan
df['area_colder0'] = df['area_colder0']/1000 # In 1000km
df_clim = df[(df.index.year>=clim_year[0]) & (df.index.year<=clim_year[1])]
std_anom = (df-df_clim.mean(axis=0))/df_clim.std(axis=0)
std_anom3LNO = std_anom[['Tmean', 'Tmean_sha100']]

# 2.
infile = 'bottom_temp_stats/stats_3Ps_spring.pkl'
df = pd.read_pickle(infile)
df.index = pd.to_datetime(df.index) # update index to datetime
# Flag bad years (no or weak sampling):
bad_years = np.array([1980, 1981, 1985, 1986, 1987, 1988, 1989, 1990, 1991, 1992, 2006, 2020])
for i in bad_years:
    df[df.index.year==i]=np.nan
df['area_colder0'] = df['area_colder0']/1000 # In 1000km
df_clim = df[(df.index.year>=clim_year[0]) & (df.index.year<=clim_year[1])]
std_anom = (df-df_clim.mean(axis=0))/df_clim.std(axis=0)
std_anom3Ps = std_anom[['Tmean', 'Tmean_sha100']]

# 4. plot all merged together
std_anom_all = pd.concat([std_anom3LNO, std_anom3Ps], axis=1)
df = std_anom_all.mean(axis=1)
df = df[(df.index.year>=years[0]) & (df.index.year<=years[-1])]
df.index = df.index.year
df.to_pickle('bottomT_index_spring.pkl')


fig, ax = plt.subplots(nrows=1, ncols=1)
n = 5 # xtick every n years
sign=df>0
ax = df.plot(kind='bar', color=sign.map({True: 'indianred', False: 'steelblue'}), width = width, zorder=10)
ticks = ax.xaxis.get_ticklocs()
ticklabels = [l.get_text() for l in ax.xaxis.get_ticklabels()]
ax.xaxis.set_ticks(ticks[::n])
ax.xaxis.set_ticklabels(ticklabels[::n])
plt.ylabel('Normalized Anomaly', weight='bold', fontsize=14)
plt.fill_between([ticks[0]-1, ticks[-1]+1], [-.5, -.5], [.5, .5], facecolor='gray', alpha=.2)
#plt.xlabel('Year')
plt.title('Mean bottom temperature anomaly for 3LNOPs - Spring', weight='bold', fontsize=14)
plt.grid()
plt.ylim([-2.5,2.5])
ax.yaxis.set_ticks(np.arange(-2.5, 3, .5))
#plt.gca().set_xlim([pd.to_datetime('1979-01-01'), pd.to_datetime('2018-01-01')]) 
fig.set_size_inches(w=15,h=7)
fig_name = 'mean_anomalies_spring.png'
fig.savefig(fig_name, dpi=200)
os.system('convert -trim mean_anomalies_spring.png mean_anomalies_spring.png')

# Save French Figure
plt.sca(ax)
plt.ylabel(u'Anomalie normalisée', weight='bold', fontsize=14)
plt.title(u'Température de fond moyenne pour 2HJ3KLNO - Printemps', weight='bold', fontsize=14)
#fig.set_size_inches(w=15,h=7)
fig_name = 'mean_anomalies_spring_FR.png'
fig.savefig(fig_name, dpi=200)
os.system('convert -trim ' + fig_name + ' ' + fig_name)


#### ------------- Merge spring and bottom ---------------- ####
bottomT_spring = pd.read_pickle('bottomT_index_spring.pkl')
bottomT_fall = pd.read_pickle('bottomT_index_fall.pkl')
bottomT = pd.concat([bottomT_spring, bottomT_fall], axis=1).mean(axis=1)

fig = plt.figure()
plt.plot(bottomT_spring)
plt.plot(bottomT_fall)
plt.plot(bottomT_fall.rolling(window=3, center=True).mean(), 'k', linewidth=2.5)
plt.grid()
plt.ylabel('Normalized anomaly')
plt.legend(['Spring', 'Fall', '3Y running mean'])
plt.title('Bottom Temperature')
ticks = plt.gca().xaxis.get_ticklocs()
plt.fill_between([ticks[0]-1, ticks[-1]+1], [-.5, -.5], [.5, .5], facecolor='gray', alpha=.2)
plt.xlim([1979.5, 2022.5])
fig_name = 'bottom_temp_anomaly.png'
fig.set_size_inches(w=12,h=7)
fig.savefig(fig_name, dpi=200)
os.system('convert -trim ' + fig_name + ' ' + fig_name)

#### ------------- Stacked spring and bottom for NL climate index ---------------- ####
# Concat and normalize
bottomT_stack = pd.concat([bottomT_spring, bottomT_fall], axis=1)
bottomT_stack_norm = bottomT_stack.divide((2 - bottomT_stack.isna().sum(axis=1)).values, axis=0)
bottomT_stack_norm.columns = ['3LNOPs - Spring', '2HJ3KLNO - Fall']

fig, ax = plt.subplots(nrows=1, ncols=1)
n = 5 # xtick every n years
sign=bottomT_stack_norm>0
bottomT_stack_norm.plot(kind='bar', stacked=True, color=['lightcoral', 'skyblue', 'indianred', 'steelblue'], width = width, zorder=10, ax=ax, legend=False)
ticks = ax.xaxis.get_ticklocs()
ticklabels = [l.get_text() for l in ax.xaxis.get_ticklabels()]
ax.xaxis.set_ticks(ticks[::n])
ax.xaxis.set_ticklabels(ticklabels[::n])
plt.ylabel('Normalized Anomaly', weight='bold', fontsize=14)
plt.fill_between([ticks[0]-1, ticks[-1]+1], [-.5, -.5], [.5, .5], facecolor='gray', alpha=.2)
#plt.xlabel('Year')
plt.title('Mean bottom temperature anomaly', weight='bold', fontsize=14)
plt.grid()
plt.ylim([-2.25,2.25])
#ax.yaxis.set_ticks(np.arange(-2.5, 3, .5))
# Custom legend
import matplotlib.lines as mlines
legend_elements = [mlines.Line2D([],[], marker='s',linestyle='None', color='lightcoral', markersize=15, label='\n3LNOPs - Spring'),
                    mlines.Line2D([],[], marker='s',linestyle='None', color='skyblue', markersize=15, label=' '),
                    mlines.Line2D([],[], marker='s',linestyle='None', color='indianred', markersize=15, label='\n2HJ3KLNO - Fall'),
                    mlines.Line2D([],[], marker='s',linestyle='None', color='steelblue', markersize=15, label=' ')]
ax.legend(handles=legend_elements)
      
#plt.gca().set_xlim([pd.to_datetime('1979-01-01'), pd.to_datetime('2018-01-01')]) 
fig.set_size_inches(w=15,h=7)
fig_name = 'bottomT_anomalies_stacked.png'
fig.savefig(fig_name, dpi=200)
os.system('convert -trim bottomT_anomalies_stacked.png bottomT_anomalies_stacked.png')

## ---- Bottom T ms_climate index [with scorecards] ---- ##
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
plt.ylim([-2,3])

colors = cmap(normal(np.nansum(bottomT_stack_norm.values, axis=1)))
cell_text = np.nansum(bottomT_stack_norm.values, axis=1).round(1).astype('str')
the_table = ax.table(cellText=[np.nansum(bottomT_stack_norm.values, axis=1).round(1)],
        rowLabels=['Bot. temp. subindex'],
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
fig_name = 'bottomT_anomalies_climateindex.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim -bordercolor White -border 10x10 ' + fig_name + ' ' + fig_name)


# Plot French Figure
bottomT_stack_normFR = bottomT_stack_norm.rename(columns={'3LNOPs - Spring':'3LNOPs - Printemps','2HJ3KLNO - Fall':'2HJ3KLNO - Automne'})
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
plt.ylim([-2,3])

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
fig_name = 'bottomT_anomalies_climateindex_FR.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim -bordercolor White -border 10x10 ' + fig_name + ' ' + fig_name)
