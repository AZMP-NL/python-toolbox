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
        'size'   : 14}
plt.rc('font', **font)


YEAR_MIN = 1951
YEAR_MAX = 2019
clim_year = [1981, 2010]
width = 0.7


#### ---- LOAD THE DATA (and prepare) ---- ####
# 1. NAO
nao = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/NAO/NAO_winter.pkl')
nao = nao.rename(columns={'Value':'Wint. NAO'})
nao = nao[nao.index<2020]
nao_natural = nao.copy()
nao = nao*-1

# 2. Air temp
air = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/airTemp/airT_std_anom.pkl')
air.index.name='Year'
air = air.mean(axis=1)
air.name = 'Air Temp'

# 3. Sea Ice (And icebergs)
#ice = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/ice/ice_index.pkl')
#ice_natural = ice.copy()
#ice = ice*-1

# Peter's version
ice = pd.read_csv('/home/cyrf0006/data/AZMP/Galbraith_data/ice-area-thick.NL.ClimateIndex.dat', header=None, sep=' ')
ice.set_index(0, inplace=True)
ice.index.name='Year'
ice.rename(columns={ice.columns[0]: "Sea Ice" }, inplace = True)
ice = ice[['Sea Ice']]
ice_natural = ice.copy()
ice = ice*-1

# 4. Icebergs
bergs = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/bergs/bergs_std_anom.pkl')
bergs.index.name='Year'
bergs.name = 'Icebergs'
bergs_natural = bergs.copy()
bergs = bergs*-1

# 5. SSTs
#sst = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/SSTs/SST_anom.pkl')
# Peter's version
sst = pd.read_csv('/home/cyrf0006/data/AZMP/Galbraith_data/NL_ClimateIndex_SST.dat', header=None, sep=' ')
sst.set_index(0, inplace=True)
sst.index.name='Year'
sst.rename(columns={sst.columns[0]: "SST" }, inplace = True)
sst = sst[['SST']]

# 6. Stn27 (0-176m, 0-50m, 150-176m)
s27_temp = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/stn27/s27_temp_std_anom.pkl')
s27_temp.index = s27_temp.index.year
#s27_temp.to_csv('S27_stn_anom_3fields.csv', float_format='%.2f')
s27_sal = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/stn27/s27_sal_std_anom.pkl')
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
#  compa previous and new CIL index
## C = pd.concat([s27_cil[['CIL temp', 'CIL core T']].mean(axis=1), s27_cil[['CIL core T']].mean(axis=1)], axis=1)  
## C.plot()
## plt.grid()
## plt.legend(['prelim (coreT, meanT)', 'coreT only'])
#s27_cil = s27_cil[['CIL temp', 'CIL core T']].mean(axis=1) # previous version
s27_cil = s27_cil[['CIL core T']].mean(axis=1) # new version
s27_cil.name = 'S27 CIL'

# 7. Section CIL (only area)
section_cil = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/CIL/section_cil_index.pkl')
section_cil = section_cil['volume'] # volume and/or core
section_cil.name = 'CIL area' 
section_cil_natural = section_cil.copy()
section_cil = section_cil*-1

# 8. bottomT
bottomT_spring = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/bottomT/bottomT_index_spring.pkl')
bottomT_fall = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/bottomT/bottomT_index_fall.pkl')
bottomT = pd.concat([bottomT_spring, bottomT_fall], axis=1).mean(axis=1)
bottomT.name = 'Bottom T'

#### ----- Merge the data ---- ####
climate_index = pd.concat([nao, air, ice, bergs, sst, s27_temp,  s27_sal, s27_cil, section_cil, bottomT], axis=1)
climate_index = climate_index[climate_index.index>=1950]
climate_index.loc[1950] = climate_index.loc[1950]*np.nan
climate_index = climate_index[climate_index.index>=1950]

# keep a copy with Natural signs
climate_index_natural = pd.concat([nao_natural, air, ice_natural, bergs_natural, sst, s27_temp,  s27_sal, s27_cil, section_cil_natural, bottomT], axis=1)
climate_index_natural.loc[1950] = climate_index_natural.loc[1950]*np.nan

#### ----- Plot climate index ---- ####
climate_index = climate_index[climate_index.index<2020]
#climate_index_norm = climate_index / 10
climate_index_norm = climate_index.divide((10 - climate_index.isna().sum(axis=1)).values, axis=0)

n = 5 # xtick every n years
#ax = climate_index.plot(kind='bar', stacked=True, cmap='gist_rainbow')
ax = climate_index_norm.plot(kind='bar', stacked=True, cmap='nipy_spectral', zorder=10)
ticks = ax.xaxis.get_ticklocs()
ticklabels = [l.get_text() for l in ax.xaxis.get_ticklabels()]
ax.xaxis.set_ticks(ticks[::n])
ax.xaxis.set_ticklabels(ticklabels[::n])
plt.fill_between([ticks[0]-1, ticks[-1]+1], [-.5, -.5], [.5, .5], facecolor='gray', alpha=.2)
plt.grid('on')
ax.set_ylabel(r'Normalized anomaly')
ax.set_title('NL Climate Index')
#ax.legend(loc=2, bbox_to_anchor=(0.45, 1), fontsize=12) 
ax.legend(bbox_to_anchor=(1.0, 1), loc='upper left', fontsize=12)
#ax.legend(fontsize=12) 
fig = ax.get_figure()
fig.set_size_inches(w=15.5,h=8)
fig_name = 'NL_climate_index.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim -bordercolor White -border 10x10 ' + fig_name + ' ' + fig_name)

## Save index (in English)
# 1. All fields.
climate_index = climate_index[(climate_index.index>=YEAR_MIN) & (climate_index.index<=YEAR_MAX)]
climate_index.index.name = 'Year'
climate_index.to_csv('NL_climate_index_all_fields.csv', float_format='%.2f')
# 2. All fields natural signs
climate_index_natural = climate_index_natural[(climate_index_natural.index>=YEAR_MIN) & (climate_index_natural.index<=YEAR_MAX)]
climate_index_natural.index.name = 'Year'
climate_index_natural.to_csv('NL_climate_index_all_fields_natural_signs.csv', float_format='%.2f')
# 3. Mean index.
climate_index_mean = climate_index.mean(axis=1)
climate_index_mean.index.name = 'Year'
climate_index_mean = climate_index_mean.rename('Climate index').to_frame()  
climate_index_mean.to_csv('NL_climate_index.csv', float_format='%.2f')

#### ----- Plot climate index (French) ---- ####
climate_index.rename(columns={
    'Winter NAO':'ONA hiver',
    'Air Temp':'Temp Air',
    'Sea Ice':'Glace',
    'SST':'SST',
    'S27 CIL':'S27 CIF',
    'CIL area':'Aire CIF',
    'Bottom T.': 'Temp Fond'
    }, inplace=True)

ax = climate_index_norm.plot(kind='bar', stacked=True, cmap='nipy_spectral', zorder=10)
ticks = ax.xaxis.get_ticklocs()
ticklabels = [l.get_text() for l in ax.xaxis.get_ticklabels()]
ax.xaxis.set_ticks(ticks[::n])
ax.xaxis.set_ticklabels(ticklabels[::n])
plt.fill_between([ticks[0]-1, ticks[-1]+1], [-.5, -.5], [.5, .5], facecolor='gray', alpha=.2)
plt.grid('on')
ax.set_ylabel(r'Anomalie normalisée cummulée')
ax.set_title('Indice climatique pour TNL')
fig = ax.get_figure()
fig.set_size_inches(w=12,h=8)
fig_name = 'NL_climate_index_FR.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim -bordercolor White -border 10x10 ' + fig_name + ' ' + fig_name)



#### ----- Comparison with Eugene's CEI ---- ####
    
df_cei = pd.read_excel('/home/cyrf0006/data/AZMP/ColbourneStuff/composite_climate_index_english.xlsx')
df_cei = df_cei.T
df_cei = df_cei.dropna(how='all')
df_cei.columns = df_cei.iloc[0]
df_cei.drop(df_cei.index[0], inplace=True)
df_cei.drop(df_cei.index[-1], inplace=True)

#df_cei = df_cei.COMPOSITE
plt.close('all')
fig = plt.figure()
plt.plot(df_cei.mean(axis=1), linewidth=2)
plt.plot(climate_index.mean(axis=1), linewidth=2)
plt.grid()
plt.text(1950,-1.5,'r=0.87', fontsize=20)
plt.legend(['CEI (Petrie et al., 2007)','NL climate index (this study)'])
fig.set_size_inches(w=12,h=8)
fig_name = 'CEI_climate_compa.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim -bordercolor White -border 10x10 ' + fig_name + ' ' + fig_name)

# double axis
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax2 = ax1.twinx()
ax1.plot(df_cei.mean(axis=1), linewidth=2, color='tab:blue')
ax2.plot(climate_index.mean(axis=1), linewidth=2, color='tab:orange')
ax1.set_ylabel('CEI (Petrie et al., 2007)', color='tab:blue')
ax2.set_ylabel('NL climate index (this study)', color='tab:orange')
ax1.grid()
ax1.text(1951, -1.5, 'r=0.87', fontsize=16, color='dimgray', weight='bold')
ax1.tick_params(axis='y', colors='tab:blue')
ax2.tick_params(axis='y', colors='tab:orange')
#plt.legend(['CEI (Petrie et al., 2007)','NL climate index (this study)'])
fig.set_size_inches(w=12,h=8)
fig_name = 'CEI_climate_compa_scaled.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim -bordercolor White -border 10x10 ' + fig_name + ' ' + fig_name)



# scatterplot
dfs = pd.concat([df_cei.mean(axis=1), climate_index.mean(axis=1)], axis=1)
fig = plt.figure()
plt.scatter(dfs[0].values, dfs[1].values, c=dfs.index)
plt.grid()
cb = plt.colorbar()
cb.ax.set_ylabel('Year', fontsize=14)
plt.xlabel('CEI (Petrie et al., 2007)')
plt.ylabel('NL climate index (this study)')
plt.text(-2,1,'r=0.87', fontsize=20)
fig.set_size_inches(w=12,h=8)
fig_name = 'CEI_climate_scatter.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim -bordercolor White -border 10x10 ' + fig_name + ' ' + fig_name)



# correlation
df_merged = pd.concat([df_cei.mean(axis=1), climate_index.mean(axis=1)], axis=1)
df_merged.corr('pearson')

#### ----- Comparison stn 27 temperatures ---- ####
s27_temp_all = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/stn27/s27_temp_std_anom.pkl')
s27_temp_all.index = s27_temp_all.index.year
df_stn27 = pd.concat([s27_temp_all.mean(axis=1), s27_temp_all['Temp 0-176m']], axis=1)   
plt.figure()
df_stn27.plot()
plt.legend(['S27 Temp index','S27 0-176m'])

# correlation
df_stn27.corr('pearson')
        

#### ----- Cross-correlation between sub-indices ---- ####
# Build the colormap
vmin = -1.0
vmax = 1.0
midpoint = 0.0
levels = np.linspace(vmin, vmax, 21)
midp = np.mean(np.c_[levels[:-1], levels[1:]], axis=1)
colvals = np.interp(midp, [vmin, midpoint, vmax], [-1, 0., 1])
normal = plt.Normalize(-1.0, 1.0)
reds = plt.cm.Reds(np.linspace(0,1, num=6))
blues = plt.cm.Blues_r(np.linspace(0,1, num=6))
whites = [(.95,.95,.95,.95)]*9
colors = np.vstack((blues[0:-1,:], whites, reds[1:,:]))
colors = np.concatenate([[colors[0,:]], colors, [colors[-1,:]]], 0)
cmap, norm = from_levels_and_colors(midp, colors, extend='both')

## keyboard
## # Test bloom
## bloom = pd.read_csv('/home/cyrf0006/research/ecosystem_considerations/SpringBloomAnomalies.csv')
## bloom = bloom[(bloom.box=='Avalon_Channel') | (bloom.box=='Northern_Grand_Bank') |(bloom.box=='Southeast_Shoal')  ]
## bloom.set_index('year', inplace=True)
## bloom = bloom['initiation']
## bloom = bloom.groupby('year').mean()
## climate_index_natural['bloom'] = bloom

        
from scipy.stats import pearsonr
def calculate_pvalues(df):
    df = df.dropna()._get_numeric_data()
    dfcols = pd.DataFrame(columns=df.columns)
    pvalues = dfcols.transpose().join(dfcols, how='outer')
    for r in df.columns:
        for c in df.columns:
            pvalues[r][c] = round(pearsonr(df[r], df[c])[1], 4)
    return pvalues

# add climate index.
climate_index_natural['Climate Index'] = climate_index.mean(axis=1)
# correlation matric and pvalues
corrMatrix = climate_index_natural.corr().round(2)
pvalues = calculate_pvalues(climate_index_natural)
annot_text  = corrMatrix.astype('str')
corrMatrix_text = corrMatrix.copy()

for i in np.arange(11):
    for j in np.arange(11):
        if pvalues.iloc[i,j]>=.05:
            #annot_text.iloc[i,j] = annot_text.iloc[i,j]+'*'
            corrMatrix.iloc[i,j] = 0            
            corrMatrix_text.iloc[i,j] = ' '
            
plt.close('all')
fig = plt.figure(3)
mask = np.zeros_like(corrMatrix)
mask[np.triu_indices_from(mask)] = True
np.fill_diagonal(mask, 0)
sn.heatmap(corrMatrix, annot=corrMatrix_text.astype('str'), fmt='s', mask=mask, linewidths=.2, cmap=cmap, cbar=None, vmin=-1.05, vmax=1.05)
#sn.heatmap(corrMatrix, annot=True, mask=mask, linewidths=.2, cmap=cmap, cbar=None, linecolor='k')
plt.title('Pearson correlation coefficients')
fig.set_size_inches(w=13,h=14)
fig_name = 'NL_climate_index_correlation.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim -bordercolor White -border 10x10 ' + fig_name + ' ' + fig_name)
