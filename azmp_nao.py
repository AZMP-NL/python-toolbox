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
url = 'https://www.ncdc.noaa.gov/teleconnections/nao/data.csv'
nao_file = '/home/cyrf0006/data/AZMP/indices/data.csv'
if os.path.exists(nao_file):

    py3 = version_info[0] > 2 #creates boolean value for test that Python major version > 2        
    response_isnt_good = True
    while response_isnt_good:
        if py3:
            response = input("Do you want to remove it? (yes/no?): ")
        else:
            response = raw_input('Do you what to update'  + nao_file + '? [y/n]')
        
        if response == 'y':
            import urllib2
            open('/home/cyrf0006/data/AZMP/indices/data.csv', 'wb').write(urllib2.urlopen(url).read())
            response_isnt_good = False
        elif response == 'n':
            response_isnt_good = False
        else:
            print ' -> Please answer "yes" or "no"'
            
# Reload using pandas
df = pd.read_csv(nao_file, header=1)

# Set index
df = df.set_index('Date')
df.index = pd.to_datetime(df.index, format='%Y%m')


## ## ----  plot Monthly NAO + 5year running mean ---- ##
## fig = plt.figure(1)
## plt.plot(df)
## plt.plot(df.rolling(window=36, center=True).mean())
## plt.ylabel('NAO')
## plt.xlabel('Year')
## plt.grid()
## fig.set_size_inches(w=12,h=9)
## fig_name = 'NAO_monthly_1950-2019.png'
## fig.savefig(fig_name, dpi=300)



## ----  plot Winter NAO ---- ####
# Select only DJF
df_winter = df[(df.index.month==12) | (df.index.month==1) | (df.index.month==2)]

# Start Dec-1950
df_winter = df_winter[df_winter.index>pd.to_datetime('1950-10-01')]
#df_winter.resample('3M').mean()


# Average 3 consecutive values (DJF average); We loose index.
df_winter = df_winter.groupby(np.arange(len(df_winter))//3).mean()

# Reset index using years only
year_unique = pd.unique(df.index.year)[1:,]
df_winter = df_winter.iloc[np.arange(0, year_unique.size)] # reduce if last month is december (belongs to following year)
df_winter.index = year_unique

## fig = plt.figure(2)
## plt.plot(df_winter)
## plt.plot(df_winter.rolling(window=5, center=True).mean())
## plt.ylabel('Winter NAO')
## plt.xlabel('Year')
## plt.grid()
## fig.set_size_inches(w=12,h=9)
## fig_name = 'NAO_winter_1950-2019.png'
## fig.savefig(fig_name, dpi=300)


## ## ---- plot winter NAO bar plots ---- ##
df1 = df_winter[df_winter>0]
df2 = df_winter[df_winter<0]

## fig = plt.figure(3)
## plt.plot(df_winter, 'k')
## plt.fill_between(df1.index, np.squeeze(df1.values), color='b')
## plt.fill_between(df2.index, np.squeeze(df2.values), color='r')
## plt.ylabel('Winter NAO')
## plt.xlabel('Year')
## plt.grid()
## fig.set_size_inches(w=12,h=9)
## fig_name = 'NAO_winter_bar_1950-2019.png'
## fig.savefig(fig_name, dpi=300)

## ---- plot winter NAO bar plots #2 ---- ##
fig = plt.figure(4)
fig.clf()
width = .9
p1 = plt.bar(df1.index, np.squeeze(df1.values), width, alpha=0.8, color='steelblue')
p2 = plt.bar(df2.index, np.squeeze(df2.values), width, bottom=0, alpha=0.8, color='indianred')
p1 = plt.bar(df1.index[-1], np.squeeze(df1.values[-1]), width, alpha=.3, color='black')
plt.ylabel('NAO index')
#plt.xlabel('Year')
plt.title('Winter NAO average (DJF)')
plt.grid()
fig.set_size_inches(w=15,h=9)
fig_name = 'NAO_winter_bar2_1950-2019.png'
#plt.annotate('data source: NCDC/NOAA', xy=(.75, .01), xycoords='figure fraction', annotation_clip=False, FontSize=12)
plt.annotate('data source: www.ncdc.noaa.gov/teleconnections/', xy=(.58, .01), xycoords='figure fraction', annotation_clip=False, FontSize=12)
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)
