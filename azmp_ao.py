#wget -O /home/cyrf0006/data/AZMP/indices/ao_data.csv https://www.ncdc.noaa.gov/teleconnections/ao/data.csv
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

# Download and save up-to-date  AO index from NOAA (data.csv) if needed
url = 'https://www.ncdc.noaa.gov/teleconnections/ao/data.csv'
ao_file = '/home/cyrf0006/data/AZMP/indices/ao_data.csv'
if os.path.exists(ao_file):
    py3 = version_info[0] > 2 #creates boolean value for test that Python major version > 2        
    response_isnt_good = True
    while response_isnt_good:
        if py3:
            response = input('Do you what to update '  + ao_file + '? [y/n]')
        else:
            response = raw_input('Do you what to update '  + ao_file + '? [y/n]')
        
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
df = pd.read_csv(ao_file, header=1)

# Set index
df = df.set_index('Date')
df.index = pd.to_datetime(df.index, format='%Y%m')

# Resample
# pickle DataFrame for scorecards:
#df.to_pickle('AO_monthly.pkl')
df_annual = df.resample('As').mean()
df_annual.index = df_annual.index.year
df_annual.to_pickle('AO_annual.pkl')

## ## ---- plot winter AO bar plots ---- ##
df1 = df_annual[df_annual>0]
df2 = df_annual[df_annual<0]

## ---- plot winter AO bar plots #2 ---- ##
fig = plt.figure(4)
fig.clf()
width = .9
p1 = plt.bar(df1.index, np.squeeze(df1.values), width, alpha=0.8, color='steelblue')
p2 = plt.bar(df2.index, np.squeeze(df2.values), width, bottom=0, alpha=0.8, color='indianred')
p1 = plt.bar(df1.index[-1], np.squeeze(df1.values[-1]), width, alpha=.3, color='black')
plt.ylabel('AO index')
#plt.xlabel('Year')
plt.title('Annual AO average')
plt.grid()
fig.set_size_inches(w=15,h=9)
fig_name = 'AO_1950-2019.png'
plt.annotate('data source: www.ncdc.noaa.gov/teleconnections/', xy=(.58, .01), xycoords='figure fraction', annotation_clip=False, FontSize=12)
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)
