#wget https://www.ncdc.noaa.gov/teleconnections/nao/data.csv
# check in : /home/cyrf0006/AZMP/annual_meetings/2019

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import datetime
import matplotlib.dates as mdates
import pandas as pd
import os
from sys import version_info

# Adjust fontsize/weight
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 18}
plt.rc('font', **font)

# Download and save up-to-date  NAO index from NOAA (data.csv) if needed
#url = 'https://www.ncdc.noaa.gov/teleconnections/nao/data.csv' (until 2020...)
url = 'https://www.cpc.ncep.noaa.gov/products/precip/CWlink/pna/norm.nao.monthly.b5001.current.ascii.table'
nao_file = '/home/cyrf0006/data/AZMP/indices/data.csv'
if os.path.exists(nao_file):
    py3 = version_info[0] > 2 #creates boolean value for test that Python major version > 2        
    response_isnt_good = True
    while response_isnt_good:
        if py3:
            response = input('Do you what to update '  + nao_file + '? [y/n]')
        else:
            response = raw_input('Do you what to update '  + nao_file + '? [y/n]')
        
        if response == 'y':
            import urllib3
            http = urllib3.PoolManager()
            r = http.request('GET', url)
            open('/home/cyrf0006/data/AZMP/indices/data.csv', 'wb').write(r.data)
            response_isnt_good = False
        elif response == 'n':
            response_isnt_good = False
        else:
            print(' -> Please answer "y" or "n"')
            
# Reload using pandas
col_names = ["Year", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
df = pd.read_csv(nao_file, header=None, delimiter=r"\s+", error_bad_lines=False)
df.set_axis(col_names, axis=1, inplace=True)                              

# Set index
df = df.set_index('Year')
df = df.stack()
df.index = pd.to_datetime(df.index.get_level_values(1) + '-' + df.index.get_level_values(0).astype('str'))

## ----  plot Winter NAO ---- ####
# Select only DJF
df_winter_djf = df[(df.index.month==12) | (df.index.month==1) | (df.index.month==2)]
df_winter = df[(df.index.month==12) | (df.index.month==1) | (df.index.month==2) | (df.index.month==3)]
df_summer = df[(df.index.month>=6) & (df.index.month<=9)]

# Start Dec-1950
df_winter_djf = df_winter_djf[df_winter_djf.index>pd.to_datetime('1950-10-01')]
df_winter = df_winter[df_winter.index>pd.to_datetime('1950-10-01')]
df_summer = df_summer[df_summer.index>pd.to_datetime('1950-01-01')]


# Average 3 consecutive values (DJF average); We loose index.
df_winter_djf = df_winter_djf.groupby(np.arange(len(df_winter_djf))//3).mean()
df_winter = df_winter.groupby(np.arange(len(df_winter))//4).mean()
df_summer = df_summer.groupby(np.arange(len(df_summer))//3).mean()

# Reset index using years only
year_unique = pd.unique(df.index.year)[1:,]
df_winter = df_winter.iloc[np.arange(0, year_unique.size)] # reduce if last month is december (belongs to following year)
df_winter_djf = df_winter_djf.iloc[np.arange(0, year_unique.size)] 
df_winter.index = year_unique
df_winter_djf.index = year_unique

df_summer = df_summer.iloc[np.arange(0, year_unique.size)] 
df_summer.index = year_unique

# pickle DataFrame for scorecards:
df.to_pickle('NAO_monthly.pkl')
df.to_csv('NAO_monthly.csv')
df_annual = df.resample('As').mean()
df_annual.index = df_annual.index.year
df_annual.to_pickle('NAO_annual.pkl')
df_winter.to_pickle('NAO_winter.pkl')
df_summer.to_pickle('NAO_summer.pkl')

## ---- plot winter NAO bar plots ---- ##
#df_winter[df_winter.index==2021]=np.nan # Remove 2021 for 2020 ResDoc
df1 = df_winter[df_winter>0]
df2 = df_winter[df_winter<0]


## ---- plot winter NAO bar plots #2 ---- ##
fig = plt.figure(4)
fig.clf()
width = .9
p1 = plt.bar(df1.index, np.squeeze(df1.values), width, alpha=0.8, color='steelblue')
p2 = plt.bar(df2.index, np.squeeze(df2.values), width, bottom=0, alpha=0.8, color='indianred')
#p1 = plt.bar(df1.index[-1], np.squeeze(df1.values[-1]), width, alpha=.5, color='white')
plt.ylabel('NAO subindex')
plt.title('Winter NAO average (DJFM)')
ticks = plt.gca().xaxis.get_ticklocs()
plt.fill_between([ticks[0]-1, ticks[-1]+1], [-.5, -.5], [.5, .5], facecolor='gray', alpha=.2)
plt.xlim([1950, 2023])
plt.grid()
fig.set_size_inches(w=15,h=9)
fig_name = 'NAO_winter_1950-2022.png'
#plt.annotate('data source: www.ncdc.noaa.gov/teleconnections/', xy=(.58, .01), xycoords='figure fraction', annotation_clip=False, FontSize=12)
fig.savefig(fig_name, dpi=300)
os.system('convert -trim -bordercolor White -border 10x10 ' + fig_name + ' ' + fig_name)


# French Figure
fig = plt.figure(4)
fig.clf()
width = .9
p1 = plt.bar(df1.index, np.squeeze(df1.values), width, alpha=0.8, color='steelblue')
p2 = plt.bar(df2.index, np.squeeze(df2.values), width, bottom=0, alpha=0.8, color='indianred')
#p1 = plt.bar(df1.index[-1], np.squeeze(df1.values[-1]), width, alpha=.5, color='white')
plt.ylabel('indice ONA')
plt.title('Oscillation Nord-Atlantique hivernale (DJFM)')
plt.grid()
fig.set_size_inches(w=15,h=9)
fig_name = 'NAO_winter_1950-2022_FR.png'
plt.annotate('source données: www.ncdc.noaa.gov/teleconnections/', xy=(.58, .01), xycoords='figure fraction', annotation_clip=False, FontSize=12)
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)


## ## ---- plot summer NAO bar plots ---- ##
df1 = df_summer[df_summer>0]
df2 = df_summer[df_summer<0]
fig = plt.figure(4)
fig.clf()
width = .9
p1 = plt.bar(df1.index, np.squeeze(df1.values), width, alpha=0.8, color='steelblue')
p2 = plt.bar(df2.index, np.squeeze(df2.values), width, bottom=0, alpha=0.8, color='indianred')
#p1 = plt.bar(df1.index[-1], np.squeeze(df1.values[-1]), width, alpha=.3, color='black')
plt.ylabel('NAO index')
#plt.xlabel('Year')
plt.title('Summer NAO average (JJAS)')
plt.grid()
fig.set_size_inches(w=15,h=9)
fig_name = 'NAO_summer_bar2_1950-2022.png'
#plt.annotate('data source: NCDC/NOAA', xy=(.75, .01), xycoords='figure fraction', annotation_clip=False, FontSize=12)
plt.annotate('data source: www.ncdc.noaa.gov/teleconnections/', xy=(.58, .01), xycoords='figure fraction', annotation_clip=False, FontSize=12)
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)
