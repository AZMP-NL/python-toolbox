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

width = 0.7
clim_year = [1981, 2010]

#### ------------- For fall ---------------- ####
# 1.
infile = 'stats_2HJ_fall.pkl'
df = pd.read_pickle(infile)
year_list = df.index # save former index
year_list = [i[2:4] for i in year_list] # 2-digit year
df.index = pd.to_datetime(df.index) # update index to datetime
df['area_colder2'] = df['area_colder2']/1000 # In 1000km
df = df[['area_colder2', 'area_colder2_perc', 'Tmean']] # keep only these 3 columns

df_clim = df[(df.index.year>=clim_year[0]) & (df.index.year<=clim_year[1])]
std_anom = (df-df_clim.mean(axis=0))/df_clim.std(axis=0)
std_anom = std_anom.rename(index=str, columns={"area_colder2": "area_std_anom", "area_colder2_perc": "area_perc_std_anom", "Tmean": "Tmean_std_anom"})
std_anom = std_anom.drop(columns=['area_std_anom']) # drop (area and percentage area same anomaly)

frames = [df, std_anom]
df = pd.concat(frames, axis=1)
df.index = df.index.year

df_area = df.area_perc_std_anom
fig = plt.figure(4)
fig.clf()
sign=df_area>0
ax = df_area.plot(kind='bar', color=sign.map({True: 'steelblue', False: 'indianred'}), width = width)
n = 5 # xtick every n years
ticks = ax.xaxis.get_ticklocs()
ticklabels = [l.get_text() for l in ax.xaxis.get_ticklabels()]
ax.xaxis.set_ticks(ticks[::n])
ax.xaxis.set_ticklabels(ticklabels[::n])
plt.ylabel('Normalized Anomaly', weight='bold', fontsize=14)
plt.fill_between([ticks[0]-1, ticks[-1]+1], [-.5, -.5], [.5, .5], facecolor='gray', alpha=.2)
plt.title(u'Area colder than 2$^{\circ}$C (anomaly) for 2HJ - Fall', weight='bold', fontsize=14)
plt.grid()
plt.ylim([-2.5,2.5])
fig.set_size_inches(w=15,h=7)
fig_name = '2HJ_area_anomalies_fall.png'
fig.savefig(fig_name)
os.system('convert -trim 2HJ_area_anomalies_fall.png 2HJ_area_anomalies_fall.png')

# save in csv
df.to_csv('bottomT_2HJ_fall.csv', float_format='%.4f')



# 2.
infile = 'stats_3K_fall.pkl'
df = pd.read_pickle(infile)
year_list = df.index # save former index
year_list = [i[2:4] for i in year_list] # 2-digit year
df.index = pd.to_datetime(df.index) # update index to datetime
df['area_colder2'] = df['area_colder2']/1000 # In 1000km
df = df[['area_colder2', 'area_colder2_perc', 'Tmean']] # keep only these 3 columns

df_clim = df[(df.index.year>=clim_year[0]) & (df.index.year<=clim_year[1])]
std_anom = (df-df_clim.mean(axis=0))/df_clim.std(axis=0)
std_anom = std_anom.rename(index=str, columns={"area_colder2": "area_std_anom", "area_colder2_perc": "area_perc_std_anom", "Tmean": "Tmean_std_anom"})
std_anom = std_anom.drop(columns=['area_std_anom']) # drop (area and percentage area same anomaly)

frames = [df, std_anom]
df = pd.concat(frames, axis=1)
df.index = df.index.year

df_area = df.area_perc_std_anom
fig = plt.figure(4)
fig.clf()
sign=df_area>0
ax = df_area.plot(kind='bar', color=sign.map({True: 'steelblue', False: 'indianred'}), width = width)
n = 5 # xtick every n years
ticks = ax.xaxis.get_ticklocs()
ticklabels = [l.get_text() for l in ax.xaxis.get_ticklabels()]
ax.xaxis.set_ticks(ticks[::n])
ax.xaxis.set_ticklabels(ticklabels[::n])
plt.ylabel('Normalized Anomaly', weight='bold', fontsize=14)
plt.fill_between([ticks[0]-1, ticks[-1]+1], [-.5, -.5], [.5, .5], facecolor='gray', alpha=.2)
plt.title(u'Area colder than 2$^{\circ}$C (anomaly) for 3K - Fall', weight='bold', fontsize=14)
plt.grid()
plt.ylim([-2.5,2.5])
fig.set_size_inches(w=15,h=7)
fig_name = '3K_area_anomalies_fall.png'
fig.savefig(fig_name)
os.system('convert -trim 3K_area_anomalies_fall.png 3K_area_anomalies_fall.png')

# save in csv
df.to_csv('bottomT_3K_fall.csv', float_format='%.4f')


# 3.
infile = 'stats_3LNO_fall.pkl'
df = pd.read_pickle(infile)
year_list = df.index # save former index
year_list = [i[2:4] for i in year_list] # 2-digit year
df.index = pd.to_datetime(df.index) # update index to datetime
df['area_colder2'] = df['area_colder2']/1000 # In 1000km
df = df[['area_colder2', 'area_colder2_perc', 'Tmean']] # keep only these 3 columns

df_clim = df[(df.index.year>=clim_year[0]) & (df.index.year<=clim_year[1])]
std_anom = (df-df_clim.mean(axis=0))/df_clim.std(axis=0)
std_anom = std_anom.rename(index=str, columns={"area_colder2": "area_std_anom", "area_colder2_perc": "area_perc_std_anom", "Tmean": "Tmean_std_anom"})
std_anom = std_anom.drop(columns=['area_std_anom']) # drop (area and percentage area same anomaly)

frames = [df, std_anom]
df = pd.concat(frames, axis=1)
df.index = df.index.year

df_area = df.area_perc_std_anom
fig = plt.figure(4)
fig.clf()
sign=df_area>0
ax = df_area.plot(kind='bar', color=sign.map({True: 'steelblue', False: 'indianred'}), width = width)
n = 5 # xtick every n years
ticks = ax.xaxis.get_ticklocs()
ticklabels = [l.get_text() for l in ax.xaxis.get_ticklabels()]
ax.xaxis.set_ticks(ticks[::n])
ax.xaxis.set_ticklabels(ticklabels[::n])
plt.ylabel('Normalized Anomaly', weight='bold', fontsize=14)
plt.fill_between([ticks[0]-1, ticks[-1]+1], [-.5, -.5], [.5, .5], facecolor='gray', alpha=.2)
plt.title(u'Area colder than 2$^{\circ}$C (anomaly) for 3LNO - Fall', weight='bold', fontsize=14)
plt.grid()
plt.ylim([-2.5,2.5])
fig.set_size_inches(w=15,h=7)
fig_name = '3LNO_area_anomalies_fall.png'
fig.savefig(fig_name)
os.system('convert -trim 3LNO_area_anomalies_fall.png 3LNO_area_anomalies_fall.png')

# save in csv
df.to_csv('bottomT_3LNO_fall.csv', float_format='%.4f')



#### ------------- For Spring ---------------- ####
# 1.
infile = 'stats_3LNO_spring.pkl'
df = pd.read_pickle(infile)
year_list = df.index # save former index
year_list = [i[2:4] for i in year_list] # 2-digit year
df.index = pd.to_datetime(df.index) # update index to datetime
df['area_colder2'] = df['area_colder2']/1000 # In 1000km
df = df[['area_colder2', 'area_colder2_perc', 'Tmean']] # keep only these 3 columns

df_clim = df[(df.index.year>=clim_year[0]) & (df.index.year<=clim_year[1])]
std_anom = (df-df_clim.mean(axis=0))/df_clim.std(axis=0)
std_anom = std_anom.rename(index=str, columns={"area_colder2": "area_std_anom", "area_colder2_perc": "area_perc_std_anom", "Tmean": "Tmean_std_anom"})
std_anom = std_anom.drop(columns=['area_std_anom']) # drop (area and percentage area same anomaly)

frames = [df, std_anom]
df = pd.concat(frames, axis=1)
df.index = df.index.year

df_area = df.area_perc_std_anom
fig = plt.figure(4)
fig.clf()
sign=df_area>0
ax = df_area.plot(kind='bar', color=sign.map({True: 'steelblue', False: 'indianred'}), width = width)
n = 5 # xtick every n years
ticks = ax.xaxis.get_ticklocs()
ticklabels = [l.get_text() for l in ax.xaxis.get_ticklabels()]
ax.xaxis.set_ticks(ticks[::n])
ax.xaxis.set_ticklabels(ticklabels[::n])
plt.ylabel('Normalized Anomaly', weight='bold', fontsize=14)
plt.fill_between([ticks[0]-1, ticks[-1]+1], [-.5, -.5], [.5, .5], facecolor='gray', alpha=.2)
plt.title(u'Area colder than 2$^{\circ}$C (anomaly) for 3LNO - Spring', weight='bold', fontsize=14)
plt.grid()
plt.ylim([-2.5,2.5])
fig.set_size_inches(w=15,h=7)
fig_name = '3LNO_area_anomalies_spring.png'
fig.savefig(fig_name)
os.system('convert -trim 3LNO_area_anomalies_spring.png 3LNO_area_anomalies_spring.png')

# save in csv
df.to_csv('bottomT_3LNO_spring.csv', float_format='%.4f')

# 2.
infile = 'stats_3Ps_spring.pkl'
df = pd.read_pickle(infile)
year_list = df.index # save former index
year_list = [i[2:4] for i in year_list] # 2-digit year
df.index = pd.to_datetime(df.index) # update index to datetime
df['area_colder2'] = df['area_colder2']/1000 # In 1000km
df = df[['area_colder2', 'area_colder2_perc', 'Tmean']] # keep only these 3 columns

df_clim = df[(df.index.year>=clim_year[0]) & (df.index.year<=clim_year[1])]
std_anom = (df-df_clim.mean(axis=0))/df_clim.std(axis=0)
std_anom = std_anom.rename(index=str, columns={"area_colder2": "area_std_anom", "area_colder2_perc": "area_perc_std_anom", "Tmean": "Tmean_std_anom"})
std_anom = std_anom.drop(columns=['area_std_anom']) # drop (area and percentage area same anomaly)

frames = [df, std_anom]
df = pd.concat(frames, axis=1)
df.index = df.index.year

df_area = df.area_perc_std_anom
fig = plt.figure(4)
fig.clf()
sign=df_area>0
ax = df_area.plot(kind='bar', color=sign.map({True: 'steelblue', False: 'indianred'}), width = width)
n = 5 # xtick every n years
ticks = ax.xaxis.get_ticklocs()
ticklabels = [l.get_text() for l in ax.xaxis.get_ticklabels()]
ax.xaxis.set_ticks(ticks[::n])
ax.xaxis.set_ticklabels(ticklabels[::n])
plt.ylabel('Normalized Anomaly', weight='bold', fontsize=14)
plt.fill_between([ticks[0]-1, ticks[-1]+1], [-.5, -.5], [.5, .5], facecolor='gray', alpha=.2)
plt.title(u'Area colder than 2$^{\circ}$C (anomaly) for 3Ps - Spring', weight='bold', fontsize=14)
plt.grid()
plt.ylim([-2.5,2.5])
fig.set_size_inches(w=15,h=7)
fig_name = '3Ps_area_anomalies_spring.png'
fig.savefig(fig_name)
os.system('convert -trim 3Ps_area_anomalies_spring.png 3Ps_area_anomalies_spring.png')

# save in csv
df.to_csv('bottomT_3Ps_spring.csv', float_format='%.4f')

