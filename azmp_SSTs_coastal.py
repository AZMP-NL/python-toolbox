'''
AZMP script to extract SST data in standard reporting boxes from the BIO's remote sensing server:
ftp://ftp.dfo-mpo.gc.ca/bometrics
ftp://ftp.dfo-mpo.gc.ca/bometrics/noaa/stats/boxes/

For oofline use, you can download all data with:
wget -m ftp://ftp.dfo-mpo.gc.ca/bometrics/noaa/stats/boxes/*.stat

check in : /home/cyrf0006/AZMP/annual_meetings/2019
http://www.bio.gc.ca/science/data-donnees/base/data-donnees/sst-en.php
Frederic.Cyr@dfo-mpo.gc.ca - February 2019

'''

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import datetime
import matplotlib.dates as mdates
import pandas as pd
import os
from sys import version_info
import re

# Adjust fontsize/weight
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 18}
plt.rc('font', **font)

## ----  Parameters to edit (getting data infile or offline) ---- ##
#prefix = 'ftp://ftp.dfo-mpo.gc.ca/bometrics/noaa/stats/boxes/'
prefix = '/home/cyrf0006/data/BIO_remote/bometrics/noaa/stats/boxes/'
# Get region names
df_box = pd.read_excel('/home/cyrf0006/github/AZMP-NL/utils/SST_boxes.xslx')

## ---- Loop on NL regions and store in a dataFrame ---- ##
dfs = []
df_labels = []
for box in df_box[(df_box.region=='NL') | (df_box.region=='GSL')].box_name.values:
    df = pd.read_csv(prefix + box +'_sst.stat', delimiter='\s+')
    df = df.rename(columns={'date-id':'date'})
    # Set index (need first to swap (a,b) by (7,15))
    date_tmp = df.date.str[-1].apply(lambda x: re.sub('a','07',x))
    date_day = date_tmp.str[-1].apply(lambda x: re.sub('b','21',x)) 
    date = df.date.map(lambda x: str(x)[:-1]) # strip last character
    df['date'] = date.astype(str) + date_day
    df = df.set_index('date')
    df.index = pd.to_datetime(df.index, format='%Y%b%d')

    dfs.append(df)
    df_labels.append(box)

# convert the list of DataFrames into a single multiindex DataFrame
df_all = pd.concat(dfs, keys=df_labels, axis=0)    


## ---- Just mean SST now ---- ##
df_sst = df_all.mean_sst 
df_sst = df_sst.unstack(level=0)
df_sst = df_sst.replace(-999.00000, np.NaN)
df_sst = df_sst[df_sst.index.year<=2018]
df_sst = df_sst.resample('Qs').mean()

# SUMMER ONLY!!!!!!!!
df_sst = df_sst[df_sst.index.month==7]

fig = plt.figure(1)
plt.clf
plt.plot(df_sst.mean(axis=1), color='k', linewidth=2)
plt.plot(df_sst.mean(axis=1).rolling(12).mean(), color='r', linewidth=5)
plt.ylabel('SST')
plt.xlabel('Year')
plt.grid()
fig.set_size_inches(w=12,h=9)
fig_name = 'SST_boxes.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)

 
df_sst = df_all.mean_sst
df_sst = df_sst.unstack(level=0)
df_sst = df_sst.replace(-999.00000, np.NaN)
df_sst = df_sst[(df_sst.index.year>=1998) & (df_sst.index.year<=2018)]
df_sst = df_sst.resample('As').mean()

fig = plt.figure(2)
plt.clf
plt.plot(df_sst.mean(axis=1), color='k', linewidth=2)
plt.ylabel('SST')
plt.xlabel('Year')
plt.grid()
fig.set_size_inches(w=12,h=9)
fig_name = 'SST_boxes_annual.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)


## Select coastel boxes only
df_coastal = df_sst[['Avalon_Channel', 'Cabot_Strait', 'Green-St._Pierre_Bank','Northeast_GSL', 'St.Anthony_Basin', 'Northeast_Nfld_Shelf']]

## Climatology over Shorter timeseries!!!
df = df_coastal.mean(axis=1)

clim_year = [1998, 2018]
clim = df[(df.index.year>=clim_year[0]) & (df.index.year<=clim_year[1])].mean()
std = df[(df.index.year>=clim_year[0]) & (df.index.year<=clim_year[1])].std()
std_anom = (df - clim)/std
std_anom.index = std_anom.index.year

std_anom = std_anom.sort_index()

fig = plt.figure(4)
fig.clf()
sign=std_anom>0
width = .7
n = 5 # xtick every n years
ax = std_anom.plot(kind='bar', color=sign.map({True: 'indianred', False: 'steelblue'}), width = width)
ticks = ax.xaxis.get_ticklocs()
ticklabels = [l.get_text() for l in ax.xaxis.get_ticklabels()]
ax.xaxis.set_ticks(ticks[::n])
ax.xaxis.set_ticklabels(ticklabels[::n])
plt.ylabel('Mean Standardized Anomaly', weight='bold', fontsize=14)
plt.title(u'SSTs', weight='bold', fontsize=14)
plt.grid()
plt.ylim([-2.5,2.5])
fig.set_size_inches(w=15,h=7)
fig_name = 'SST_anomalies_coastal_summer.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)

# save in csv
# std_anom.to_csv('SST_anom_short.csv', float_format='%.4f')


## Stackbar plotdddddf

clim_year = [1998, 2018]
clim_coast = df_coastal[(df_coastal.index.year>=clim_year[0]) & (df_coastal.index.year<=clim_year[1])].mean()
std_coast = df_coastal[(df_coastal.index.year>=clim_year[0]) & (df_coastal.index.year<=clim_year[1])].std()
std_anom_coastal = (df_coastal - clim_coast)/std_coast
anom_coastal = (df_coastal - clim_coast)
std_anom_coastal.index = std_anom_coastal.index.year
anom_coastal.index = anom_coastal.index.year

fig = plt.figure(5)
fig.clf()
ax = anom_coastal.plot(kind='bar', stacked=True, cmap='YlGn')
plt.grid('on')
ax.set_ylabel(r'[$^{\circ}$C]')

n = 5 # xtick every n years
ticks = ax.xaxis.get_ticklocs()
ticklabels = [l.get_text() for l in ax.xaxis.get_ticklabels()]
ax.xaxis.set_ticks(ticks[::n])
ax.xaxis.set_ticklabels(ticklabels[::n])
plt.grid('on')
plt.legend(['AC', 'CS', 'SPB', 'NEGSL', 'SAB', 'NENL'])
fig = ax.get_figure()
fig.set_size_inches(w=12,h=8)
fig_name = 'SST_coastal_anom.png'
#plt.annotate('data source: www.ncdc.noaa.gov/teleconnections/', xy=(.58, .01), xycoords='figure fraction', annotation_clip=False, FontSize=12)
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)
