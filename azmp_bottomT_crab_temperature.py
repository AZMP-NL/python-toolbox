'''
To generate timeseries of 2-4 degC bottom habitat.

Uses pickled object generated by azmp_bottom_stats.py

'''

import numpy as  np
import matplotlib.pyplot as plt
import pandas as pd
import os

# plot tuning:
font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 13}
plt.rc('font', **font)
width = 0.7

clim_year = [1981, 2010]

#### ------------- For fall ---------------- ####
# 1.
infile = 'stats_2J_fall.pkl'
df = pd.read_pickle(infile)
df.index = pd.to_datetime(df.index) # update index to datetime
df = df['Tmean']
df = df/1000
df_clim = df[(df.index.year>=clim_year[0]) & (df.index.year<=clim_year[1])]
std_anom = (df-df_clim.mean(axis=0))/df_clim.std(axis=0)
std_anom2J = std_anom

# 2.
infile = 'stats_3K_fall.pkl'
df = pd.read_pickle(infile)
df.index = pd.to_datetime(df.index) # update index to datetime
df = df['Tmean']
df = df/1000
df_clim = df[(df.index.year>=clim_year[0]) & (df.index.year<=clim_year[1])]
std_anom = (df-df_clim.mean(axis=0))/df_clim.std(axis=0)
std_anom3K = std_anom


# 3.
infile = 'stats_3LNO_fall.pkl'
df = pd.read_pickle(infile)
df.index = pd.to_datetime(df.index) # update index to datetime
df = df['Tmean']
df = df/1000
df_clim = df[(df.index.year>=clim_year[0]) & (df.index.year<=clim_year[1])]
std_anom = (df-df_clim.mean(axis=0))/df_clim.std(axis=0)
std_anom3LNO = std_anom


# 4. plot all merged together
std_anom_all = pd.concat([std_anom2J, std_anom3K, std_anom3LNO], axis=1)
df_2J3KLNO = std_anom_all.mean(axis=1)
df_2J3KLNO.index = df_2J3KLNO.index.year


#### ------------- For Spring ---------------- ####
# 1.
infile = 'stats_3LNO_spring.pkl'
df = pd.read_pickle(infile)
df.index = pd.to_datetime(df.index) # update index to datetime
df = df['Tmean']
df = df/1000
df_clim = df[(df.index.year>=clim_year[0]) & (df.index.year<=clim_year[1])]
std_anom = (df-df_clim.mean(axis=0))/df_clim.std(axis=0)
std_anom3LNO = std_anom

# 2.
infile = 'stats_3Ps_spring.pkl'
df = pd.read_pickle(infile)
df.index = pd.to_datetime(df.index) # update index to datetime
df = df['Tmean']
df = df/1000
df_clim = df[(df.index.year>=clim_year[0]) & (df.index.year<=clim_year[1])]
std_anom = (df-df_clim.mean(axis=0))/df_clim.std(axis=0)
std_anom3Ps = std_anom

# 4. plot all merged together
std_anom_all = pd.concat([std_anom3LNO, std_anom3Ps], axis=1)
df_3LNOPs = std_anom_all.mean(axis=1)
df_3LNOPs.index = df_3LNOPs.index.year


#### ------------- For Summer ---------------- ####
# 1.
infile = 'stats_4R_summer.pkl'
df = pd.read_pickle(infile)
df.index = pd.to_datetime(df.index) # update index to datetime
df = df['Tmean']
df = df/1000
# !!! Flagging some years!
flags = np.array([1980, 1981, 1982, 1990, 1995, 1996, 1997, 1998, 2015, 2018])
df[df.index.year.isin(flags)]=np.nan

df_clim = df[(df.index.year>=clim_year[0]) & (df.index.year<=clim_year[1])]
std_anom = (df-df_clim.mean(axis=0))/df_clim.std(axis=0)
std_anom4R = std_anom

# 4. plot all merged together
df_4R = std_anom4R
df_4R.index = df.index.year

df_all = pd.concat([df_2J3KLNO, df_3LNOPs], axis=1)
df_all = df_all.mean(axis=1)


fig = plt.figure(6)
fig.clf()
#ax = df_4R.plot()
#df_3LNOPs.plot(ax=ax)
ax = df_3LNOPs.plot()
df_2J3KLNO.plot(ax=ax)
plt.plot(df_all.index, df_all.rolling(5, center=True).mean(), 'k--', linewidth=3)
#plt.legend(['4R (summer)', '3LNOPs (spring)', '2J3KLNO (fall)', '5yr running mean 2J3KNLOPs'])
plt.legend(['3LNOPs (spring)', '2J3KLNO (fall)', '5yr running mean'])
plt.ylabel('Mean Standardized Anomaly', weight='bold', fontsize=14)
#plt.xlabel('Year')
plt.title(r'Bottom temperature anomaly', weight='bold', fontsize=14)
plt.grid()
plt.ylim([-2.5,2.5])
#plt.gca().set_xlim([pd.to_datetime('1979-01-01'), pd.to_datetime('2018-01-01')]) 
fig.set_size_inches(w=15,h=7)
fig_name = 'bottomT_trends.png'
#fig_name = 'crab_bottomT_anonlaie_all.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)
