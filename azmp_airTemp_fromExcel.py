'''
AZMP reporting - Air temperature from Colbourne's Excel sheets
(script ran in /home/cyrf0006/AZMP/annual_meetings/2019)

Using data from (see ~/research/PeopleStuff/ColbourneStuff):
AZMP_AIR_TEMP_COMPOSITE_2018.xlsx
Fred built:
AZMP_AIR_TEMP_COMPOSITE_BONAVISTA.xlsx
AZMP_AIR_TEMP_COMPOSITE_CARTWRIGHT.xlsx
AZMP_AIR_TEMP_COMPOSITE_IQALUIT.xlsx
AZMP_AIR_TEMP_COMPOSITE_NUUK.xlsx
AZMP_AIR_TEMP_COMPOSITE_STJOHNS.xlsx
that are loaded and plotted here. 

Ideally, I would use directly data from EC  Homogenized Temperature: ftp://ccrp.tor.ec.gc.ca/pub/AHCCD/Homog_monthly_mean_temp.zip

but since some data are delayed or unavailable for NL stations (NUUK is in Greenland, Bonavista N/A and Cartwright stops in 2015), Eugne used to got them from :
http://climate.weather.gc.ca/prods_servs/cdn_climate_summary_e.html
http://climate.weather.gc.ca/prods_servs/cdn_climate_summary_report_e.html?intYear=2018&intMonth=2&prov=NL&dataFormat=csv&btnSubmit=Download+data

and update the Excel files.

Eventually, I could find a way to update directly from server (see azmp_airTemp.py).

I took NUUK temperature here:
https://www.dmi.dk/publikationer/
https://www.dmi.dk/vejrarkiv/
using file 4250_2014_2018.csv
( this one ends in 2013: https://crudata.uea.ac.uk/cru/data/greenland/nuuk.dat)

I generated historical data from here (see azmp_dmi_nuukAirT.py):
JOURNAL OF GEOPHYSICAL RESEARCH, VOL. 111, D11105, doi:10.1029/2005JD006810, 2006
https://www.dmi.dk/fileadmin/user_upload/Rapporter/TR/2018/DMIRep18-05.zip
https://www.dmi.dk/publikationer/

Frederic.Cyr@dfo-mpo.gc.ca - February 2019

'''


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import datetime
import matplotlib.dates as mdates
from scipy.interpolate import griddata
import os

# Adjust fontsize/weight
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 14}
plt.rc('font', **font)

clim_year = [1981, 2010]
current_year = 2019

## ----  Prepare the data ---- ##
# load from Excel sheets
df_BB = pd.read_excel('/home/cyrf0006/research/PeopleStuff/ColbourneStuff/AZMP_AIR_TEMP_COMPOSITE_BONAVISTA.xlsx', header=0)
df_CA = pd.read_excel('/home/cyrf0006/research/PeopleStuff/ColbourneStuff/AZMP_AIR_TEMP_COMPOSITE_CARTWRIGHT.xlsx', header=0)
df_IQ = pd.read_excel('/home/cyrf0006/research/PeopleStuff/ColbourneStuff/AZMP_AIR_TEMP_COMPOSITE_IQALUIT.xlsx', header=0)
df_NK = pd.read_excel('/home/cyrf0006/research/PeopleStuff/ColbourneStuff/AZMP_AIR_TEMP_COMPOSITE_NUUK.xlsx', header=0)
df_SJ = pd.read_excel('/home/cyrf0006/research/PeopleStuff/ColbourneStuff/AZMP_AIR_TEMP_COMPOSITE_STJOHNS.xlsx', header=0)

# Rename columns
col_names = ['Year', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
df_BB.columns = col_names
df_CA.columns = col_names
df_IQ.columns = col_names
df_NK.columns = col_names
df_SJ.columns = col_names
# Set index Year
df_BB = df_BB.set_index('Year', drop=True)
df_CA = df_CA.set_index('Year', drop=True)
df_IQ = df_IQ.set_index('Year', drop=True)
df_NK = df_NK.set_index('Year', drop=True)
df_SJ = df_SJ.set_index('Year', drop=True)


# Stack months under Years (pretty cool!)
df_BB = df_BB.stack() 
df_CA = df_CA.stack() 
df_IQ = df_IQ.stack() 
df_NK = df_NK.stack() 
df_SJ = df_SJ.stack() 

# Transform to a series with values based the 15th of each month (had to convert years to string)
df_BB.index = pd.to_datetime('15-' + df_BB.index.get_level_values(1) + '-' + df_BB.index.get_level_values(0).values.astype(np.str))
df_CA.index = pd.to_datetime('15-' + df_CA.index.get_level_values(1) + '-' + df_CA.index.get_level_values(0).values.astype(np.str))
df_IQ.index = pd.to_datetime('15-' + df_IQ.index.get_level_values(1) + '-' + df_IQ.index.get_level_values(0).values.astype(np.str))
df_NK.index = pd.to_datetime('15-' + df_NK.index.get_level_values(1) + '-' + df_NK.index.get_level_values(0).values.astype(np.str))
df_SJ.index = pd.to_datetime('15-' + df_SJ.index.get_level_values(1) + '-' + df_SJ.index.get_level_values(0).values.astype(np.str))

# NEW FROM 2019 (replace Excel from DMI data):
df_NUUK = pd.read_pickle('Nuuk_air_temp.pkl')  

# Concatenate all timeseries
df = pd.concat([df_SJ, df_BB, df_CA, df_IQ, df_NUUK], axis=1)
df.columns = ['StJohns', 'Bonavista', 'Cartwright','Iqaluit', 'Nuuk']


## ---- Monthly anomalies for current year ---- ##
df_clim_period = df[(df.index.year>=clim_year[0]) & (df.index.year<=clim_year[1])]
df_monthly_stack = df_clim_period.groupby([(df_clim_period.index.year),(df_clim_period.index.month)]).mean()
df_monthly_clim = df_monthly_stack.mean(level=1)
df_monthly_std = df_monthly_stack.std(level=1)

df_current_year = df[df.index.year==current_year]
year_index = df_current_year.index # backup index
df_current_year.index=df_monthly_std.index # reset index
anom = df_current_year - df_monthly_clim
std_anom = (df_current_year - df_monthly_clim)/df_monthly_std
#std_anom.index = year_index.month # replace index
std_anom.index = year_index.strftime('%b') # replace index (by text)
anom.index = year_index.strftime('%b') # replace index (by text)


## ---- Annual anomalies ---- ##
df_annual = df.resample('As').mean()
df_annual = df_annual[df_annual.index.year>=1950]
clim = df_annual[(df_annual.index.year>=clim_year[0]) & (df_annual.index.year<=clim_year[1])].mean()
std = df_annual[(df_annual.index.year>=clim_year[0]) & (df_annual.index.year<=clim_year[1])].std()
std_anom_annual = (df_annual - clim)/std
std_anom_annual.index = std_anom_annual.index.year

# Save for scorecards
std_anom_annual.to_pickle('airT_std_anom.pkl')
df.to_pickle('airT_monthly.pkl')


## ---- plot monthly ---- ##
#std_anom.plot(kind='bar', stacked=True, cmap='YlGn')
#plt.grid('on')
#plt.show()

ax = anom.plot(kind='bar', stacked=True, cmap='YlGn')
plt.grid('on')
ax.set_ylabel(r'[$^{\circ}$C]')
ax.set_title(np.str(current_year) + ' Air temperature anomalies')
#ax.legend(loc='upper center')
plt.ylim([-6, 14])


fig = ax.get_figure()
fig.set_size_inches(w=9,h=6)
fig_name = 'air_temp_2019.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)

# Save French Figure
french_months = ['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D']
ax.set_title(' Anomalies des températures de l\'air - ' + np.str(current_year))
ax.set_xticklabels(french_months, rotation='horizontal')
fig_name = 'air_temp_2019_FR.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)


## ---- plot annual ---- ##
n = 5 # xtick every n years
ax = std_anom_annual.plot(kind='bar', stacked=True, cmap='YlGn')
ticks = ax.xaxis.get_ticklocs()
ticklabels = [l.get_text() for l in ax.xaxis.get_ticklabels()]
ax.xaxis.set_ticks(ticks[::n])
ax.xaxis.set_ticklabels(ticklabels[::n])
plt.grid('on')
ax.set_ylabel(r'Standardized anomaly')
ax.set_title('Annual air temperature anomalies')
fig = ax.get_figure()
fig.set_size_inches(w=12,h=8)
fig_name = 'air_temp_anom.png'
#plt.annotate('data source: www.ncdc.noaa.gov/teleconnections/', xy=(.58, .01), xycoords='figure fraction', annotation_clip=False, FontSize=12)
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)

# Save in French
ax.set_ylabel(r'Anomalie normalisée')
ax.set_title('Anomalies des températures de l\'air')
fig_name = 'air_temp_anom_FR.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)


keyboard

## ---- Tried to overlay NAO, but did not work. I give up for now. ---- ##
df_winter = pd.read_pickle('winterNAO_1951-2019.pkl')
df_roll = df_winter.rolling(window=5, center=True).mean()
df_roll = df_roll.rename(columns={"Value": "Winter NAO"})
df = pd.concat([std_anom_annual, df_roll])
df.index.name = 'year'
df.reset_index(inplace=True)

df['year'] = df['year'].astype("string") # Let them be strings!
fig, ax = plt.subplots(figsize = (15,8))
df.plot(x = 'year', y = ['Bonavista','Cartwright'], kind = 'bar', ax = ax)
df.plot(x = 'year', y = ['Winter NAO'], kind = 'line', ax = ax)

ax = df.plot(x = 'year', y = ['Bonavista','Cartwright'], kind = 'bar')#x='month', linestyle='-', marker='o')
df.plot(x='year', y='Winter NAO', kind='line', ax=ax)


fig, ax = plt.subplots(figsize = (15,8))
#std_anom_annual.plot(kind='bar', stacked=True, cmap='YlGn', ax=ax)
#df_roll.plot(kind='line', color='k', linewidth=4, ax=ax)

df.plot(y = ['Winter NAO'], kind = 'line', ax = ax)
df.plot(y= ['Bonavista','Cartwright'], kind = 'bar', ax = ax)

ticks = ax.xaxis.get_ticklocs()
ticklabels = [l.get_text() for l in ax.xaxis.get_ticklabels()]
ax.xaxis.set_ticks(ticks[::n])
ax.xaxis.set_ticklabels(ticklabels[::n])
plt.grid('on')
ax.set_ylabel(r'Standardized anomaly')
ax.set_title('Annual air temperature anomalies')
fig = ax.get_figure()
fig.set_size_inches(w=12,h=8)
fig_name = 'air_temp_anom_NAO.png'
#plt.annotate('data source: www.ncdc.noaa.gov/teleconnections/', xy=(.58, .01), xycoords='figure fraction', annotation_clip=False, FontSize=12)
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)

## Test only with summer data
df_summer = df[(df.index.month >= 6) & (df.index.month<=9)]
df_annual = df_summer.resample('As').mean()
df_annual = df_annual[df_annual.index.year>=1950]
clim = df_annual[(df_annual.index.year>=clim_year[0]) & (df_annual.index.year<=clim_year[1])].mean()
std = df_annual[(df_annual.index.year>=clim_year[0]) & (df_annual.index.year<=clim_year[1])].std()
std_anom_annual = (df_annual - clim)/std
std_anom_annual.index = std_anom_annual.index.year

ax = std_anom_annual['Nuuk'].plot()
std_anom_annual['Iqaluit'].plot(ax=ax)
plt.show()
