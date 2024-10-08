# -*- coding: utf-8 -*-
'''
To generate Colbourne's and STACFIS composite anomalies


Uses this pickled DataFrame:

/home/cyrf0006/AZMP/state_reports/SSTs/SSTs_merged_monthly.pkl
generated by from azmp_sst_scorecards.py


'''

import numpy as  np
import matplotlib.pyplot as plt
import pandas as pd
import os
import unicodedata
from matplotlib.colors import from_levels_and_colors
import seaborn as sn
import cmocean as cmo

# Adjust fontsize/weight
font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 13}
plt.rc('font', **font)


YEAR_MIN = 1950
YEAR_MAX = 2021
#clim_year = [1981, 2010]
clim_year = [1991, 2020]
width = 0.7


#### ---- LOAD THE DATA (and prepare) ---- ####
# 1. NAO
nao = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/airTemp/NAO_winter.pkl')
nao.name = 'Wint. NAO'
nao = nao[nao.index<YEAR_MAX+1]
# normalize NAO with std.
#nao_clim = nao[(nao.index>=clim_year[0]) & (nao.index<=clim_year[1])]
#nao = nao / nao_clim.std()
nao_natural = nao.copy()
nao = nao*-1

# 2. Air temp
if clim_year[0] == 1981:
    air = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/airTemp/airT_std_81anom.pkl')
else:
    air = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/airTemp/airT_std_anom.pkl')
air.index.name='Year'
air = air.mean(axis=1)
air.name = 'Air Temp'

# 3. Sea Ice (And icebergs)
if clim_year[0] == 1981:
    ice = pd.read_csv('/home/cyrf0006/data/AZMP/Galbraith_data/ice-area-thick.NL.ClimateIndex.1981-2010.dat', header=None, sep=' ')
else:
    ice = pd.read_csv('/home/cyrf0006/data/AZMP/Galbraith_data/ice-area-thick.NL.ClimateIndex.dat', header=None, sep=' ')
ice.set_index(0, inplace=True)
ice.index.name='Year'
ice.rename(columns={ice.columns[0]: "Sea Ice" }, inplace = True)
ice = ice[['Sea Ice']]
ice_natural = ice.copy()
ice = ice*-1

# previous version:
#ice = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/ice/ice_index.pkl')
#ice_natural = ice.copy()
#ice = ice*-1

# 4. Icebergs
if clim_year[0] == 1981:
    bergs = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/bergs/bergs_std_anom_1981clim.pkl')
else:
    bergs = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/bergs/bergs_std_anom.pkl')
bergs.index.name='Year'
bergs.name = 'Icebergs'
bergs_natural = bergs.copy()
bergs = bergs*-1

# 5. SSTs
if clim_year[0] == 1981:
    sst = pd.read_csv('/home/cyrf0006/data/AZMP/Galbraith_data/NL_ClimateIndex_SST.1981-2010.dat', header=None, sep=' ')
else:
    sst = pd.read_csv('/home/cyrf0006/data/AZMP/Galbraith_data/NL_ClimateIndex_SST.1991-2020.dat', header=0, sep=' ')
sst.set_index('year', inplace=True)
sst.index.name='Year'
sst.rename(columns={sst.columns[0]: "SST" }, inplace = True)
sst = sst[['SST']]
sst.index = sst.index 

# previous version
#sst = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/SSTs/SST_anom.pkl')

# 6. Stn27 (0-176m, 0-50m, 150-176m)
if clim_year[0] == 1981:
    s27_temp = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/stn27/1981clim/s27_temp_std_anom.pkl')
    s27_sal = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/stn27/1981clim/s27_sal_std_anom.pkl')    
else:
    s27_temp = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/stn27/s27_temp_std_anom.pkl')
    s27_sal = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/stn27/s27_sal_std_anom.pkl')

# Flag years with less than 8 months (!!!Should be in pkl object!!!)
s27_temp[s27_temp.index.year==1950] = np.nan
s27_temp[s27_temp.index.year==1980] = np.nan
s27_temp[s27_temp.index.year==1981] = np.nan
s27_temp[s27_temp.index.year==2020] = np.nan
s27_sal[s27_sal.index.year==1950] = np.nan
s27_sal[s27_sal.index.year==1980] = np.nan
s27_sal[s27_sal.index.year==1981] = np.nan
s27_sal[s27_sal.index.year==2020] = np.nan

# Set index and rename
s27_temp.index = s27_temp.index.year
#s27_temp.to_csv('S27_stn_anom_3fields.csv', float_format='%.2f')
s27_sal.index = s27_sal.index.year
# Take only 0-176m average (correlation is 0.94, see below) 
#s27_temp = s27_temp.mean(axis=1)
#s27_sal = s27_sal.mean(axis=1)
s27_temp = s27_temp['Temp 0-176m']
s27_sal = s27_sal['Sal 0-176m']
s27_temp.name = 'S27 T'
s27_sal.name = 'S27 S'
s27_sal_natural = s27_sal.copy()
s27_sal = s27_sal*-1 # assume fresh = warm

#  compa previous and new T index
## s27_temp1 = s27_temp.mean(axis=1)
## s27_temp2 = s27_temp['Temp 0-176m']
## A = pd.concat([s27_temp1, s27_temp2], axis=1)  
## s27_sal1 = s27_sal.mean(axis=1)
## s27_sal2 = s27_sal['Sal 0-176m']
## B = pd.concat([s27_sal1, s27_sal2], axis=1)  
## A.plot()
## plt.grid()
## plt.legend(['prelim (all 3 depths)', 'new (0-176m)'])
## B.plot()
## plt.grid()
## plt.legend(['prelim (all 3 depths)', 'new (0-176m)'])

# s27_cil not std_anom yet
s27_cil = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/stn27/S27_CIL_summer_stats.pkl')
s27_cil.index = s27_cil.index.year
s27_cil_clim = s27_cil[(s27_cil.index>=clim_year[0]) & (s27_cil.index<=clim_year[1])]
s27_cil = (s27_cil-s27_cil_clim.mean(axis=0))/s27_cil_clim.std(axis=0)

# Flag years with less than 8 months
s27_cil[s27_cil.index==1950] = np.nan
s27_cil[s27_cil.index==1980] = np.nan
s27_cil[s27_cil.index==1981] = np.nan
s27_cil[s27_cil.index==2020] = np.nan

#  compa previous and new CIL index
## C = pd.concat([s27_cil[['CIL temp', 'CIL core T']].mean(axis=1), s27_cil[['CIL core T']].mean(axis=1)], axis=1)  
## C.plot()
## plt.grid()
## plt.legend(['prelim (coreT, meanT)', 'coreT only'])
#s27_cil = s27_cil[['CIL temp', 'CIL core T']].mean(axis=1) # previous version
s27_cil = s27_cil[['CIL core T']].mean(axis=1) # new version
s27_cil.name = 'S27 CIL'

# 7. Section CIL (only area) [.pkl files form azmp_CIL_mean_anomaly.py]
if clim_year[0] == 1981:
    section_cil = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/sections_plots/CIL/section_cil_index_1981clim.pkl')
else:
    section_cil = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/sections_plots/CIL/section_cil_index.pkl')
# Select volume and/or coreT
section_cil = section_cil['volume'] # volume and/or core
section_cil.name = 'CIL area' 
section_cil_natural = section_cil.copy()
section_cil = section_cil*-1

# 8. bottomT [.pkl files from azmp_bottomT_mean_anomaly.py]
if clim_year[0] == 1981:
    bottomT_spring = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/bottomT/bottomT_index_spring_1981clim.pkl')
    bottomT_fall = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/bottomT/bottomT_index_fall_1981clim.pkl')
else:
    bottomT_spring = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/bottomT/bottomT_index_spring.pkl')
    bottomT_fall = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/bottomT/bottomT_index_fall.pkl')
bottomT = pd.concat([bottomT_spring, bottomT_fall], axis=1).mean(axis=1)
bottomT.name = 'Bottom T'

#### ----- Merge the data ---- ####
climate_index = pd.concat([nao, air, ice, bergs, sst, s27_temp,  s27_sal, s27_cil, section_cil, bottomT], axis=1)
climate_index = climate_index[climate_index.index>=YEAR_MIN].sort_index()
climate_index.loc[1950] = climate_index.loc[1950]*np.nan
climate_index_sc = climate_index.copy() # for scorecards at top

# keep a copy with Natural signs
climate_index_natural = pd.concat([nao_natural, air, ice_natural, bergs_natural, sst, s27_temp,  s27_sal_natural, s27_cil, section_cil_natural, bottomT], axis=1)
climate_index_natural = climate_index_natural[climate_index_natural.index>=YEAR_MIN].sort_index()

## restrict time series and normalize for plots.
climate_index = climate_index[climate_index.index<YEAR_MAX+1]
climate_index_norm = climate_index.divide((10 - climate_index.isna().sum(axis=1)).values, axis=0)
climate_index_norm_ns = climate_index_sc.divide((10 - climate_index_sc.isna().sum(axis=1)).values, axis=0)

#### ----- Plot climate index [1] ---- ####
climate_index_norm.reset_index(inplace=True)

# reset index
year = climate_index_norm['Year']
climate_index_norm = climate_index_norm.drop('Year', 1)

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


n = 5 # xtick every n years
#ax = climate_index.plot(kind='bar', stacked=True, cmap='gist_rainbow')
fig, ax = plt.subplots(nrows=1, ncols=1)
climate_index_norm.plot(ax=ax, kind='bar', stacked=True, cmap='nipy_spectral', zorder=10, legend=False)
ticks = ax.xaxis.get_ticklocs()
ticklabels = [l.get_text() for l in ax.xaxis.get_ticklabels()]
ax.xaxis.set_ticks(ticks[::n])
ax.xaxis.set_ticklabels(ticklabels[::n])
plt.fill_between([ticks[0]-1, ticks[-1]+1], [-.5, -.5], [.5, .5], facecolor='gray', alpha=.2)
plt.grid('on')
ax.set_ylabel(r'Normalized anomaly')
ax.set_title('NL Climate Index')
#ax.legend(loc=2, bbox_to_anchor=(0.45, 1), fontsize=12) 
ax.set_xticklabels(year[::n], rotation=0)
ax.set_ylim([-1.55, 1.55])
ax.set_xlim([ticks[0]-1, ticks[-1]+1])


#Save 1
fig_name = 'NLCI_capelin_keynote1.png'
fig.set_size_inches(w=12,h=8)
fig.savefig(fig_name, dpi=400)
os.system('convert -trim -bordercolor White -border 10x10 ' + fig_name + ' ' + fig_name)

# add the smooth
ax2 = ax.twinx()
climate_index_norm.mean(axis=1).rolling(2).mean().cumsum().plot(ax=ax2, color='magenta', linewidth=5)
ax2.set_yticks([])
YLIM = ax.get_ylim() 
ax2.set_ylim(YLIM)


#Save 2
fig_name = 'NLCI_capelin_keynote2.png'
fig.set_size_inches(w=12,h=8)
fig.savefig(fig_name, dpi=400)
os.system('convert -trim -bordercolor White -border 10x10 ' + fig_name + ' ' + fig_name)

# Breaking points based on the NLCI:
years_list = [
[-1, 21],
[21, 26],
[26, 32],
[32, 48],
[48, 64],
[64, 67],
[67, 72]
]

years_colors = [
['green'],
['red'],
['green'],
['red'],
['green'],
['red'],
['green']
]
        
for years in years_list:
    plt.plot([years[0], years[0]], [-1.55, 1.55], '--k')

ax2.set_ylim(YLIM)

#Save 3
fig_name = 'NLCI_capelin_keynote3.png'
fig.set_size_inches(w=12,h=8)
fig.savefig(fig_name, dpi=400)
os.system('convert -trim -bordercolor White -border 10x10 ' + fig_name + ' ' + fig_name)

#Save 4
for idx, years in enumerate(years_list):
    c = plt.fill_between([years[0], years[1]], [-1.55, -1.55], [1.55, 1.55], facecolor=years_colors[idx], alpha=.1, zorder=-1)
fig_name = 'NLCI_capelin_keynote4.png'
fig.set_size_inches(w=12,h=8)
fig.savefig(fig_name, dpi=400)
os.system('convert -trim -bordercolor White -border 10x10 ' + fig_name + ' ' + fig_name)

## Manuscript Figure
fig, ax = plt.subplots(nrows=1, ncols=1)
climate_index_norm.plot(ax=ax, kind='bar', stacked=True, cmap='nipy_spectral', zorder=10, legend=False)
ticks = ax.xaxis.get_ticklocs()
ticklabels = [l.get_text() for l in ax.xaxis.get_ticklabels()]
ax.xaxis.set_ticks(ticks[::n])
ax.xaxis.set_ticklabels(ticklabels[::n])
plt.fill_between([ticks[0]-1, ticks[-1]+1], [-.5, -.5], [.5, .5], facecolor='gray', alpha=.2)
plt.grid('on')
ax.set_ylabel(r'Normalized anomaly')
ax.set_title('NL Climate Index')
ax.legend(loc=2, bbox_to_anchor=(0, 1), fontsize=7) 
ax.set_xticklabels(year[::n], rotation=0)
ax.set_ylim([-1.55, 1.55])


colors = cmap(normal(np.nansum(climate_index_norm.values, axis=1)))
nlci_text = np.nansum(climate_index_norm.values, axis=1).round(1).astype('str')
if nlci_text[0] == '0.0':
    nlci_text[0] = 'nan'
the_table = ax.table(cellText=[nlci_text],
        rowLabels=['NL climate index'],
        colLabels=None,
        cellColours = [colors],
        cellLoc = 'center', rowLoc = 'center',
        loc='bottom', bbox=[0, -0.14, 1, 0.05])
the_table.auto_set_font_size (False)
the_table.set_fontsize(11)

for key, cell in the_table.get_celld().items():
    cell_text = cell.get_text().get_text() 
    if key[1] == -1:
        cell.set_linewidth(0)
        cell.set_fontsize(10)
    elif cell_text=='nan':
        cell._set_facecolor('darkgray')
        cell._text.set_color('darkgray')
        cell.set_fontsize(0)
    else:
        cell._text.set_rotation(90)


# Add regimes
for idx, years in enumerate(years_list):
    c = plt.fill_between([years[0], years[1]], [-1.55, -1.55], [1.55, 1.55], facecolor=years_colors[idx], alpha=.1, zorder=-1)
for years in years_list:
    plt.plot([years[0], years[0]], [-1.55, 1.55], '--k')

# add the smooth
ax2 = ax.twinx()
climate_index_norm.sum(axis=1).rolling(2).mean().cumsum().plot(ax=ax2, color='magenta', linewidth=5)
ax2.set_ylabel('NLCI Cumsum', color='magenta')
ax2.set_yticks(np.linspace(0,10,6))
ax2.tick_params(axis='y', colors='magenta')
plt.xlim([ticks[0]-1, ticks[-1]+1])
fig_name = 'NLCI_bottom-up_ms_sc.png'
fig.set_size_inches(w=12,h=8)
fig.savefig(fig_name, dpi=400)
os.system('convert -trim -bordercolor White -border 10x10 ' + fig_name + ' ' + fig_name)

## plot the smooth only
fig, ax = plt.subplots(nrows=1, ncols=1)
climate_index_norm.mean(axis=1).rolling(2).mean().cumsum().plot(ax=ax, color='black', linewidth=5)
ticks = np.arange(0, 75, n)
ax.xaxis.set_ticks(ticks)
ax.set_xticklabels(year[::n], rotation=0)
plt.grid('on')
ax.set_ylabel(r'')
ax.set_title('Cummulative NLCI')
ax.set_ylim([0,1.5])
for years in years_list:
    plt.plot([years[0], years[0]], [0, 1.5], '--k')
fig_name = 'NLCI_smooth_only.png'
fig.set_size_inches(w=12,h=8)
fig.savefig(fig_name, dpi=200)
os.system('convert -trim -bordercolor White -border 10x10 ' + fig_name + ' ' + fig_name)

## Loop on year and add periods of growth and reduction
for idx, years in enumerate(years_list):
    c = plt.fill_between([years[0], years[1]], [0, 0], [1.5, 1.5], facecolor=years_colors[idx], alpha=.2)
    fig_name = 'NLCI_smooth_only_period_' + str(idx) + '.png'
    fig.set_size_inches(w=12,h=8)
    fig.savefig(fig_name, dpi=200)
    os.system('convert -trim -bordercolor White -border 10x10 ' + fig_name + ' ' + fig_name)
    c.remove()

    
for idx, years in enumerate(years_list):
    c = plt.fill_between([years[0], years[1]], [0, 0], [1.5, 1.5], facecolor=years_colors[idx], alpha=.2)

fig_name = 'NLCI_smooth_only_period_all.png'
fig.set_size_inches(w=12,h=8)
fig.savefig(fig_name, dpi=200)
os.system('convert -trim -bordercolor White -border 10x10 ' + fig_name + ' ' + fig_name)



## NLCI vs bloom initiation
# Load bloom param and compute metrics
df_bloom = pd.read_csv('/home/cyrf0006/research/Wu_revisited/bloom_max_timing_anom.csv')
# plot
fig, ax = plt.subplots(nrows=1, ncols=1)
climate_index_norm.plot(ax=ax, kind='bar', stacked=True, cmap='nipy_spectral', zorder=10, legend=False)
ticks = ax.xaxis.get_ticklocs()
ticklabels = [l.get_text() for l in ax.xaxis.get_ticklabels()]
ax.xaxis.set_ticks(ticks[::n])
ax.xaxis.set_ticklabels(ticklabels[::n])
plt.fill_between([ticks[0]-1, ticks[-1]+1], [-.5, -.5], [.5, .5], facecolor='gray', alpha=.2)
plt.grid('on')
ax.set_ylabel(r'Normalized anomaly')
ax.set_title('NL Climate Index')
ax.set_xticklabels(year[::n], rotation=0)
ax.set_ylim([-1.55, 1.55])
# add bloom params
ax2 = ax.twinx()
df_bloom.plot(ax=ax2, color='black', linewidth=5)
ax2.set_yticks([])
YLIM = ax.get_ylim() 
ax2.set_ylim(YLIM)
#Save
fig_name = 'NLCI_bloom_capelin.png'
fig.set_size_inches(w=12,h=8)
fig.savefig(fig_name, dpi=400)
os.system('convert -trim -bordercolor White -border 10x10 ' + fig_name + ' ' + fig_name)
