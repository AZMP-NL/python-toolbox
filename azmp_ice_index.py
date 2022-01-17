'''
A script to read and plot data from Ice EXtent Excel sheet south of 55N on NL shelves

Note from Nov. 2021: Now deprecated. PSG calculates the new sea-ice subindex.

'''
# A first test to read Excel nutrient file and export to Pandas.

# Check in:
#  /home/cyrf0006/research/AZMP_database/biochem

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import datetime
from scipy import stats
import os

# Adjust fontsize/weight
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 18}
plt.rc('font', **font)

clim_year = [1981, 2010]
years = [1969, 2018]
width = 0.7
bergs_in = False

## ----  Load data a make index ---- ##
# Peter duration
df = pd.read_csv('/home/cyrf0006/data/seaIce_IML/OccurrenceOfIce.GEC.dat', header=None, delimiter=r"\s+", error_bad_lines=False, names=['Year', 'Region', 'First', 'Last', 'Duration', 'Date_first', 'Date_last'])
df.replace(-99, np.nan, inplace=True)  
df = df.set_index('Year', drop=True)
df_lab = df[df.Region==3]
df_nfl = df[df.Region==4]


# Peter's volume
df_vol = pd.read_csv('/home/cyrf0006/data/seaIce_IML/IceVolumeRegionsMax.GEC.dat', header=1, delimiter=r"\s+", error_bad_lines=False)
df_vol.rename(columns={'##1':'Year'}, inplace=True)
df_vol = df_vol.set_index('Year', drop=True)
vol_lab = df_vol['11'] 
vol_nfl = df_vol['14'] 

# Icebergs
#bergs = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/bergs/bergs_std_anom.pkl')
bergs = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/bergs/bergs_annual.pkl')
bergs.index.name='Year'

# Average and climatology
if bergs_in:
    df_ice = pd.concat([df_lab.Duration, df_nfl.Duration, vol_lab, vol_nfl, bergs], axis=1)
else:
    df_ice = pd.concat([df_lab.Duration, df_nfl.Duration, vol_lab, vol_nfl], axis=1)

df_clim = df_ice[(df_ice.index>=clim_year[0]) & (df_ice.index<=clim_year[1])]
df_anom = df_ice-df_clim.mean()
df_std_anom = df_anom/df_clim.std()

# Save index
ice_index = df_std_anom.mean(axis=1)
if bergs_in:
    ice_index.to_pickle('ice_bergs_index.pkl')
else:
    ice_index.to_pickle('ice_index.pkl')
    ice_index.loc[1968] = np.nan
    ice_index.loc[1967] = np.nan
    ice_index.loc[1966] = np.nan
    ice_index.loc[1965] = np.nan
    ice_index = ice_index.sort_index()                                           

# Limit to later than 1950
ice_index = ice_index[ice_index.index>=1965]


## ---- plot index ---- ##
fig, ax = plt.subplots(nrows=1, ncols=1)
sign=ice_index>0
ice_index.plot(kind='bar', color=sign.map({True: 'steelblue', False: 'indianred'}), width = width, zorder=10)
n = 5 # xtick every n years
ticks = ax.xaxis.get_ticklocs()
ticklabels = [l.get_text() for l in ax.xaxis.get_ticklabels()]
ax.xaxis.set_ticks(ticks[::n])
ax.xaxis.set_ticklabels(ticklabels[::n])
plt.ylabel('Normalized Anomaly', weight='bold', fontsize=14)
plt.grid()
plt.ylim([-2.5,2.5])
plt.fill_between([ticks[0]-1, ticks[-1]+1], [-.5, -.5], [.5, .5], facecolor='gray', alpha=.2)
if bergs_in:
    plt.title('NL Sea Ice and Icebergs index', weight='bold', fontsize=14)
    fig_name = 'ice_bergs_index.png'
else:
    plt.title('NL Sea Ice index', weight='bold', fontsize=14)
    fig_name = 'ice_index.png'
fig.set_size_inches(w=15,h=7)
fig.savefig(fig_name, dpi=200)
os.system('convert -trim ' + fig_name + ' ' + fig_name)

# Save French Figure
plt.sca(ax)
plt.ylabel(u'Anomalie normalisée', weight='bold', fontsize=14)
plt.title(u'Indice glace de mer - T-N-L', weight='bold', fontsize=14)
#fig.set_size_inches(w=15,h=7)
fig_name = 'ice_index_FR.png'
fig.savefig(fig_name, dpi=200)
os.system('convert -trim ' + fig_name + ' ' + fig_name)

