#wget https://www.esrl.noaa.gov/psd/data/correlation/amon.us.data
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

# Download and save up-to-date  AMO index from NOAA (data.csv) if needed
url = 'https://www.esrl.noaa.gov/psd/data/correlation/amon.us.data'
amo_file = '/home/cyrf0006/data/AZMP/indices/amon.us.data'
if os.path.exists(amo_file):
    py3 = version_info[0] > 2 #creates boolean value for test that Python major version > 2        
    response_isnt_good = True
    while response_isnt_good:
        if py3:
            response = input('Do you what to update '  + amo_file + '? [y/n]')
        else:
            response = raw_input('Do you what to update '  + amo_file + '? [y/n]')
        
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
df = pd.read_csv(amo_file, header=1, delimiter=r"\s+", error_bad_lines=False, names=col_names)
# Remove last rows
df = df.iloc[0:-4]

# Set index
df = df.set_index('Year')
df = df.stack()      
df.index = pd.to_datetime(df.index.get_level_values(1) + '-' + df.index.get_level_values(0))

# Remove wrong values
df = df.astype(float)
df = df[df>-99]


## ----  plot Monthly AMO + 5year running mean ---- ##
fig = plt.figure(1)
fig.clf()
plt.plot(df)
plt.plot(df.rolling(window=60, center=True).mean(), linewidth=3)
plt.ylabel('AMO')
plt.xlabel('Year')
plt.grid()
plt.title('AMO unsmoothed, detrended from the Kaplan SST V2')
plt.legend(['monthly', '5-year smooth'])
fig_name = 'AMO_monthly_1950-2020.png'
fig.set_size_inches(w=12,h=9)
fig.savefig(fig_name, dpi=300)


# pickle DataFrame for scorecards:
df.to_pickle('AMO_monthly.pkl')
df_annual = df.resample('As').mean()
df_annual.index = df_annual.index.year
df_annual.to_pickle('AMO_annual.pkl')


## ## ---- plot winter AMO bar plots ---- ##
df1 = df_annual[df_annual>0]
df2 = df_annual[df_annual<0]

## ---- plot winter AMO bar plots #2 ---- ##
fig = plt.figure(4)
fig.clf()
width = .9
p1 = plt.bar(df1.index, np.squeeze(df1.values), width, alpha=0.8, color='steelblue')
p2 = plt.bar(df2.index, np.squeeze(df2.values), width, bottom=0, alpha=0.8, color='indianred')
p1 = plt.bar(df1.index[-1], np.squeeze(df1.values[-1]), width, alpha=.3, color='black')
plt.ylabel('AMO index')
#plt.xlabel('Year')
plt.title('AMO unsmoothed, detrended from the Kaplan SST V2')
plt.grid()
fig.set_size_inches(w=15,h=9)
fig_name = 'AMO_bar_1950-2020.png'
plt.annotate('data source: http://www.esrl.noaa.gov/psd/data/timeseries/AMO/', xy=(.58, .01), xycoords='figure fraction', annotation_clip=False, FontSize=12)
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)
