'''
AZMP script to extract SST data in standard reporting boxes from the BIO's remote sensing server:
ftp://ftp.dfo-mpo.gc.ca/bometrics
ftp://ftp.dfo-mpo.gc.ca/bometrics/noaa/stats/boxes/

For oofline use, you can download all data with:
wget -m ftp://ftp.dfo-mpo.gc.ca/bometrics/noaa/stats/boxes/*.stat

check in : /home/cyrf0006/AZMP/state_reports/SSTs
(formerly in /home/cyrf0006/AZMP/annual_meetings/2019)


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
prefix_excel = '/home/cyrf0006/data/SSTs/'
# Get region names
df_box = pd.read_excel('/home/cyrf0006/github/AZMP-NL/utils/SST_boxes.xslx')

## ---- Loop on NL regions and store in a dataFrame ---- ##
col_names = ['Year', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
dfs = []
df_labels = []
dfs_excel = []
df_labels = []
for box in df_box[df_box.region=='NL'].box_name.values:
    
    # From BOMETRICS
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

    # From Eugene Excel
    df_excel = pd.read_excel(prefix_excel + 'SST_' + box +'.xlsx')
    df_excel.columns = col_names
    df_excel = df_excel.set_index('Year', drop=True)

    df_excel = df_excel.stack() 
    df_excel.index = pd.to_datetime('15-' + df_excel.index.get_level_values(1) + '-' + df_excel.index.get_level_values(0).values.astype(np.str))
    dfs_excel.append(df_excel)

    df_labels.append(box)

# convert the list of DataFrames into a single multiindex DataFrame
df_all = pd.concat(dfs, keys=df_labels, axis=0)    
df_all_excel = pd.concat(dfs_excel, keys=df_labels, axis=0)    

# Remove data where coverage is less than 15%
df_all = df_all[df_all['%coverage']>=15]

## ---- Just mean SST now ---- ##
df_sst = df_all.mean_sst 
df_sst = df_sst.unstack(level=0)
df_sst = df_sst.replace(-999.00000, np.NaN)
df_sst = df_sst[df_sst.index.year<=2019]
df_sst = df_sst.resample('MS', loffset=pd.Timedelta(14, 'd')).mean() #re-averaged bi-weekly on 15th of the month

df_sst_excel = df_all_excel.unstack(level=0)

# merge 2 monthly timeseries and save for Scorecards
df_monthly = pd.concat([df_sst, df_sst_excel]) # merge
df_monthly = df_monthly.resample('MS', loffset=pd.Timedelta(14, 'd')).mean() # re-average
df_monthly.to_pickle('SSTs_merged_monthly.pkl')


# Quick plot to check the merge
fig = plt.figure(1)
plt.clf
plt.plot(df_sst.mean(axis=1), color='k', linewidth=2)
plt.plot(df_sst.mean(axis=1).rolling(12).mean(), color='r', linewidth=5)
plt.plot(df_sst_excel.mean(axis=1), color='gray', linewidth=2)
plt.plot(df_sst_excel.mean(axis=1).rolling(12).mean(), color='m', linewidth=5)
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
df_sst = df_sst[(df_sst.index.year>=1998) & (df_sst.index.year<=2019)]
df_sst = df_sst.resample('As').mean()

df_sst_excel = df_sst_excel.resample('As').mean()


fig = plt.figure(2)
plt.clf
plt.plot(df_sst_excel.mean(axis=1), color='k', linewidth=2)
plt.plot(df_sst.mean(axis=1), color='r', linewidth=2)
plt.ylabel('SST')
plt.xlabel('Year')
plt.grid()
fig.set_size_inches(w=12,h=9)
fig_name = 'SST_boxes_annual.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)



## ---- Merging two timeseries + anomaly calculation ---- ##
df_merged = pd.concat([df_sst, df_sst_excel]) # merge
df_merged = df_merged.resample('As').mean() # re-average
#df = pd.concat([df_sst.mean(axis=1), df_sst_excel.mean(axis=1)], axis=1)    
df = df_merged.mean(axis=1)

clim_year = [1981, 2010]
clim = df[(df.index.year>=clim_year[0]) & (df.index.year<=clim_year[1])].mean()
std = df[(df.index.year>=clim_year[0]) & (df.index.year<=clim_year[1])].std()
std_anom = (df - clim)/std
std_anom.index = std_anom.index.year

# a little tweak for plot:
std_anom[std_anom.index==1981]=np.nan   
std_anom = std_anom.append(pd.Series([np.nan], index=[1980]))
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
plt.ylabel('Mean Normalized Anomaly', weight='bold', fontsize=14)
plt.title(u'NL SST Index', weight='bold', fontsize=14)
plt.grid()
plt.ylim([-2,2])
plt.fill_between([ticks[0]-1, ticks[-1]+1], [-.5, -.5], [.5, .5], facecolor='gray', alpha=.2)
fig.set_size_inches(w=15,h=7)
fig_name = 'SST_index.png'
fig.savefig(fig_name, dpi=200)
os.system('convert -trim ' + fig_name + ' ' + fig_name)

# Save French Figure
plt.sca(ax)
plt.ylabel(u'Anomalie normalisée', weight='bold', fontsize=14)
plt.title(u'Indice de température de surface - T-N-L', weight='bold', fontsize=14)
#fig.set_size_inches(w=15,h=7)
fig_name = 'SST_index_FR.png'
fig.savefig(fig_name, dpi=200)
os.system('convert -trim ' + fig_name + ' ' + fig_name)



# save in csv & pkl
std_anom.to_csv('SST_anom.csv', float_format='%.4f')
std_anom.to_pickle('SST_anom.pkl')
