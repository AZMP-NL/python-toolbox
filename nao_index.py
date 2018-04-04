'''
NAO index computed from this index:
https://www.ncdc.noaa.gov/teleconnections/nao/data.csv

The North Atlantic Oscillation (NAO) index is based on the surface sea-level pressure difference between the Subtropical (Azores) High and the Subpolar Low. The positive phase of the NAO reflects below-normal heights and pressure across the high latitudes of the North Atlantic and above-normal heights and pressure over the central North Atlantic, the eastern United States and western Europe. The negative phase reflects an opposite pattern of height and pressure anomalies over these regions. Both phases of the NAO are associated with basin-wide changes in the intensity and location of the North Atlantic jet stream and storm track, and in large-scale modulations of the normal patterns of zonal and meridional heat and moisture transport, which in turn results in changes in temperature and precipitation patterns often extending from eastern North America to western and central Europe.

Strong positive phases of the NAO tend to be associated with above-normal temperatures in the eastern United States and across northern Europe and below-normal temperatures in Greenland and oftentimes across southern Europe and the Middle East. They are also associated with above-normal precipitation over northern Europe and Scandinavia and below-normal precipitation over southern and central Europe. Opposite patterns of temperature and precipitation anomalies are typically observed during strong negative phases of the NAO. During particularly prolonged periods dominated by one particular phase of the NAO, abnormal height and temperature patterns are also often seen extending well into central Russia and north-central Siberia. The NAO exhibits considerable interseasonal and interannual variability, and prolonged periods (several months) of both positive and negative phases of the pattern are common.

The NAO index is obtained by projecting the NAO loading pattern to the daily anomaly 500 millibar height field over 0-90N. The NAO loading pattern has been chosen as the first mode of a Rotated Empirical Orthogonal Function (EOF) analysis using monthly mean 500 millibar height anomaly data from 1950 to 2000 over 0-90N latitude.


!!! TO DO !!!
For now the file is found here: /home/cyrf0006/research/AZMP_stateReports/data/indices/data.csv

But I should write a shell script to update it using a crontab and but it in ~/Data


Frederic.Cyr@dfo-mpo.gc.ca - March 2018

'''


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import datetime
import matplotlib.dates as mdates

# Adjust fontsize/weight
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 18}
plt.rc('font', **font)

# Load data
nao_file = '/home/cyrf0006/research/AZMP_stateReports/data/indices/data.csv'
df = pd.read_csv(nao_file, header=1)

# Set index
df = df.set_index('Date')
df.index = pd.to_datetime(df.index, format='%Y%m')




## ----  plot Monthly NAO + 5year running mean ---- ##
fig = plt.figure(1)
plt.plot(df)
plt.plot(df.rolling(window=36, center=True).mean())
plt.ylabel('NAO')
plt.xlabel('Year')
plt.grid()
fig.set_size_inches(w=12,h=9)
fig_name = 'NAO_monthly_1950-2018.png'
fig.savefig(fig_name)



## ----  plot Winter NAO ---- ####
# Select only DJF
df_winter = df[(df.index.month==12) | (df.index.month==1) | (df.index.month==2)]

# Start Dec-1950
df_winter = df_winter[df_winter.index>pd.to_datetime('1950-10-01')]
df_winter.resample('3M').mean()


# Average 3 consecutive values (DJF average); We loose index.
df_winter = df_winter.groupby(np.arange(len(df_winter))//3).mean()

# Reset index using years only
df_winter.index = pd.unique(df.index.year)[1:,]

fig = plt.figure(2)
plt.plot(df_winter)
plt.plot(df_winter.rolling(window=5, center=True).mean())
plt.ylabel('Winter NAO')
plt.xlabel('Year')
plt.grid()
fig.set_size_inches(w=12,h=9)
fig_name = 'NAO_winter_1950-2018.png'
fig.savefig(fig_name)


## ---- plot winter NAO bar plots ---- ##
df1 = df_winter[df_winter>0]
df2 = df_winter[df_winter<0]

fig = plt.figure(3)
plt.plot(df_winter, 'k')
plt.fill_between(df1.index, np.squeeze(df1.values), color='b')
plt.fill_between(df2.index, np.squeeze(df2.values), color='r')
plt.ylabel('Winter NAO')
plt.xlabel('Year')
plt.grid()
fig.set_size_inches(w=12,h=9)
fig_name = 'NAO_winter_bar_1950-2018.png'
fig.savefig(fig_name)

## ---- plot winter NAO bar plots #2 ---- ##
fig = plt.figure(4)
fig.clf()
width = .9
p1 = plt.bar(df1.index, np.squeeze(df1.values), width, alpha=0.8)
p2 = plt.bar(df2.index, np.squeeze(df2.values), width, bottom=0, alpha=0.8)
p1 = plt.bar(df1.index[-1], np.squeeze(df1.values[-1]), width, alpha=.3, color='black')
plt.ylabel('NAO index')
#plt.xlabel('Year')
plt.title('Winter NAO average (DJF)')
plt.grid()
fig.set_size_inches(w=15,h=9)
fig_name = 'NAO_winter_bar2_1950-2018.png'
#plt.annotate('data source: NCDC/NOAA', xy=(.75, .01), xycoords='figure fraction', annotation_clip=False, FontSize=12)
plt.annotate('data source: www.ncdc.noaa.gov/teleconnections/', xy=(.58, .01), xycoords='figure fraction', annotation_clip=False, FontSize=12)
fig.savefig(fig_name)

