'''
This script extract bottom temperature from pickled object (created with azmp_bottom_stats.py) and store them in a .csv file.
This scipt was created in preparation of 2019 snow crab stock assessment.

-> see email request from D. Mullowney on 14-01-2019

Frederic.Cyr@dfo-mpo.gc.ca

'''

import numpy as  np
import matplotlib.pyplot as plt
import pandas as pd
import os

width = 0.35
clim_year = [1981, 2010]



## ----  Prepare the data ---- ##
# load from Excel sheets
s27_sst = '/home/cyrf0006/AZMP/S27/S27_SST.xlsx'
s27_bot = '/home/cyrf0006/AZMP/S27/S27_botT.xlsx'
df_sst = pd.read_excel(s27_sst)
df_bot = pd.read_excel(s27_bot)

# Rename columns
col_names = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
df_sst.columns = col_names
df_bot.columns = col_names

# Stack months under Years (pretty cool!)
df_sst = df_sst.stack() 
df_bot = df_bot.stack() 

# Transform to a series with values based the 15th of each month (had to convert years to string)
df_sst.index = pd.to_datetime('15-' + df_sst.index.get_level_values(1) + '-' + df_sst.index.get_level_values(0).values.astype(np.str))
df_bot.index = pd.to_datetime('15-' + df_bot.index.get_level_values(1) + '-' + df_bot.index.get_level_values(0).values.astype(np.str))

# Concatenate all timeseries
df_monthly = pd.concat([df_sst, df_bot], axis=1)
df_monthly.columns = ['surface', 'bottom']

# Compute annual mean
df = df_monthly.resample('As').mean()
df_annual = df.copy()
df_annual.index = df_annual.index.year # update index

# Save to csv
df_monthly.to_csv('S27_temp_monthly.csv', float_format='%.4f')
df_annual.to_csv('S27_temp_annual.csv', float_format='%.4f')

# compute anomalies
df_clim = df[(df.index.year>=clim_year[0]) & (df.index.year<=clim_year[1])]
std_anom = (df - df_clim.mean())/df_clim.std()
std_anom.index = std_anom.index.year

## ----  plot annual timemsereis ---- ##
fig = plt.figure(4)
fig.clf()
sign=std_anom.surface<0
ax = std_anom.surface.plot(kind='bar', color=sign.map({True: 'cornflowerblue', False: 'lightcoral'}), width = width, position=0)
sign2=std_anom.bottom<0
ax = std_anom.bottom.plot(kind='bar', color=sign2.map({True: 'navy', False: 'darkred'}), width = width, ax=ax, position=1)
plt.grid('on')

n = 5 # xtick every n years
ticks = ax.xaxis.get_ticklocs()
ticklabels = [l.get_text() for l in ax.xaxis.get_ticklabels()]
ax.xaxis.set_ticks(ticks[::n])
ax.xaxis.set_ticklabels(ticklabels[::n])
plt.grid('on')
#plt.legend(['Surface', '', 'Bottom', ''])
ax.set_ylabel(r'Standardized anomaly')
ax.set_title('Station 27')

# custom legend
from matplotlib.lines import Line2D
custom_lines = [Line2D([0], [0], color='cornflowerblue', lw=4),
                Line2D([0], [0], color='lightcoral', lw=4),
                Line2D([0], [0], color='navy', lw=4),
                Line2D([0], [0], color='darkred', lw=4)]

ax.legend(custom_lines, ['Surface', '', 'Bottom', ''])

fig = ax.get_figure()
fig.set_size_inches(w=12,h=8)
fig_name = 'S27_anom.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)



## CHECK FOR INTERACTIVE PLOTS: https://stackoverflow.com/questions/43387013/bar-chart-pandas-dataframe-with-bokeh
