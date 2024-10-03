'''
# AZMP reporting - Air temperature from Colbourne's Excel sheets
#  check in /home/cyrf0006/research/AZMP_stateReports/2017

## ---- Link info ---- ##
# There are two types of product.
# Homogenized Temperature is the best since corrections are applied:
# ftp://ccrp.tor.ec.gc.ca/pub/AHCCD/Homog_monthly_mean_temp.zip

# Update 2022 with Gen.3, may want to use this address:
# https://www.canada.ca/en/environment-climate-change/services/climate-change/science-research-data/climate-trends-variability/adjusted-homogenized-canadian-data/surface-air-temperature-access.html
#
# Since these are available with delay, it is sometimes necessary to use standard monthly temperature:
# http://climate.weather.gc.ca/prods_servs/cdn_climate_summary_e.html
# http://climate.weather.gc.ca/prods_servs/cdn_climate_summary_report_e.html?intYear=2018&intMonth=2&prov=NL&dataFormat=csv&btnSubmit=Download+data
# !!! Need to write a shell script / crontab to update these files automatically !!!

8400601	BONAVISTA
8501106	CARTWRIGHT
8403505	ST_JOHN_S
8403603	ST_JOHN_WEST
2402592	IQALUIT

I took NUUK temperature here:
https://www.dmi.dk/publikationer/
https://www.dmi.dk/vejrarkiv/
using file 4250_2014_2018.csv
( this one ends in 2013: https://crudata.uea.ac.uk/cru/data/greenland/nuuk.dat)

I generated historical data from here (see azmp_dmi_nuukAirT.py):
JOURNAL OF GEOPHYSICAL RESEARCH, VOL. 111, D11105, doi:10.1029/2005JD006810, 2006
https://www.dmi.dk/fileadmin/user_upload/Rapporter/TR/2018/DMIRep18-05.zip
https://www.dmi.dk/publikationer/

** Note that NUUK Air temperature is also provided in ices/iroc by Boris**

Frederic.Cyr@dfo-mpo.gc.ca
Nov. 2020 (this is now the good version using AHCCD)


'''
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import datetime
import matplotlib.dates as mdates
from scipy.interpolate import griddata
import urllib.request
import zipfile

# Adjust fontsize/weight
font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 14}
plt.rc('font', **font)

   
#clim_year = [1981, 2010]
clim_year = [1991, 2020]
print('Enter the year of interest: ')
current_year = input()
current_year = int(current_year)
use_climate_summaries = True
update_files = True

#Download the data for the year
if update_files:
    for month in np.arange(1,12+1):
        urllib.request.urlretrieve(
            'http://climate.weather.gc.ca/prods_servs/cdn_climate_summary_report_e.html?intYear='+str(current_year)+'&intMonth='+str(month)+'&prov=NL&dataFormat=csv&btnSubmit=Download+data',
            os.path.expanduser('~/github/AZMP-NL/external_data/ECCC/climate_summaries/en_climate_summaries_NL_'+"%.2d" % month+'-'+str(current_year)+'.csv'))
        urllib.request.urlretrieve(
            'http://climate.weather.gc.ca/prods_servs/cdn_climate_summary_report_e.html?intYear='+str(current_year)+'&intMonth='+str(month)+'&prov=NU&dataFormat=csv&btnSubmit=Download+data',
            os.path.expanduser('~/github/AZMP-NL/external_data/ECCC/climate_summaries/en_climate_summaries_NU_'+"%.2d" % month+'-'+str(current_year)+'.csv'))

#Update the Homogenized Temperature Zip file
if update_files:
    urllib.request.urlretrieve(
        'http://crd-data-donnees-rdc.ec.gc.ca/CDAS/products/AHCCD/Homog_monthly_mean_temp_Gen3.zip',
        os.path.expanduser('~/github/AZMP-NL/external_data/ECCC/homog_monthly_mean_temp/Homog_monthly_mean_temp_Gen3.zip'))
    with zipfile.ZipFile(os.path.expanduser('~/github/AZMP-NL/external_data/ECCC/homog_monthly_mean_temp/Homog_monthly_mean_temp_Gen3.zip'), 'r') as zip_ref:
        zip_ref.extractall(os.path.expanduser('~/github/AZMP-NL/external_data/ECCC/homog_monthly_mean_temp/'))


## ---- If climate summaries are needed ---- ##
empty_frame_NL = pd.Series(np.array([np.nan,np.nan,np.nan]),index=['8400601','8403505','8501106'])
empty_frame_NL.index.name = 'Clim_ID'
empty_frame_NL.name = 'Tm'
empty_frame_NU = pd.Series(np.array([np.nan]),index=['2402592'])
empty_frame_NU.index.name = 'Clim_ID'
empty_frame_NU.name = 'Tm'


if use_climate_summaries:
    # NL
    NL_01 = pd.read_csv('~/github/AZMP-NL/external_data/ECCC/climate_summaries/en_climate_summaries_NL_01-' + str(current_year) + '.csv', index_col='Clim_ID')
    if NL_01.empty:
        NL_01 = empty_frame_NL
    else:
        NL_01 = NL_01.loc[['8400601','8403505','8501106']].Tm
    NL_02 = pd.read_csv('~/github/AZMP-NL/external_data/ECCC/climate_summaries/en_climate_summaries_NL_02-' + str(current_year) + '.csv', index_col='Clim_ID')
    if NL_02.empty:
        NL_02 = empty_frame_NL
    else:
        NL_02 = NL_02.loc[['8400601','8403505','8501106']].Tm
    NL_03 = pd.read_csv('~/github/AZMP-NL/external_data/ECCC/climate_summaries/en_climate_summaries_NL_03-' + str(current_year) + '.csv', index_col='Clim_ID')
    if NL_03.empty:
        NL_03 = empty_frame_NL
    else:
        NL_03 = NL_03.loc[['8400601','8403505','8501106']].Tm
    NL_04 = pd.read_csv('~/github/AZMP-NL/external_data/ECCC/climate_summaries/en_climate_summaries_NL_04-' + str(current_year) + '.csv', index_col='Clim_ID')
    if NL_04.empty:
        NL_04 = empty_frame_NL
    else:
        NL_04 = NL_04.loc[['8400601','8403505','8501106']].Tm
    NL_05 = pd.read_csv('~/github/AZMP-NL/external_data/ECCC/climate_summaries/en_climate_summaries_NL_05-' + str(current_year) + '.csv', index_col='Clim_ID')
    if NL_05.empty:
        NL_05 = empty_frame_NL
    else:
        NL_05 = NL_05.loc[['8400601','8403505','8501106']].Tm
    NL_06 = pd.read_csv('~/github/AZMP-NL/external_data/ECCC/climate_summaries/en_climate_summaries_NL_06-' + str(current_year) + '.csv', index_col='Clim_ID')
    if NL_06.empty:
        NL_06 = empty_frame_NL
    else:
        NL_06 = NL_06.loc[['8400601','8403505','8501106']].Tm
    NL_07 = pd.read_csv('~/github/AZMP-NL/external_data/ECCC/climate_summaries/en_climate_summaries_NL_07-' + str(current_year) + '.csv', index_col='Clim_ID')
    if NL_07.empty:
        NL_07 = empty_frame_NL
    else:
        NL_07 = NL_07.loc[['8400601','8403505','8501106']].Tm
    NL_08 = pd.read_csv('~/github/AZMP-NL/external_data/ECCC/climate_summaries/en_climate_summaries_NL_08-' + str(current_year) + '.csv', index_col='Clim_ID')
    if NL_08.empty:
        NL_08 = empty_frame_NL
    else:
        NL_08 = NL_08.loc[['8400601','8403505','8501106']].Tm
    NL_09 = pd.read_csv('~/github/AZMP-NL/external_data/ECCC/climate_summaries/en_climate_summaries_NL_09-' + str(current_year) + '.csv', index_col='Clim_ID')
    if NL_09.empty:
        NL_09 = empty_frame_NL
    else:
        NL_09 = NL_09.loc[['8400601','8403505','8501106']].Tm
    NL_10 = pd.read_csv('~/github/AZMP-NL/external_data/ECCC/climate_summaries/en_climate_summaries_NL_10-' + str(current_year) + '.csv', index_col='Clim_ID')
    if NL_10.empty:
        NL_10 = empty_frame_NL
    else:
        NL_10 = NL_10.loc[['8400601','8403505','8501106']].Tm
    NL_11 = pd.read_csv('~/github/AZMP-NL/external_data/ECCC/climate_summaries/en_climate_summaries_NL_11-' + str(current_year) + '.csv', index_col='Clim_ID')
    if NL_11.empty:
        NL_11 = empty_frame_NL
    else:
        NL_11 = NL_11.loc[['8400601','8403505','8501106']].Tm
    NL_12 = pd.read_csv('~/github/AZMP-NL/external_data/ECCC/climate_summaries/en_climate_summaries_NL_12-' + str(current_year) + '.csv', index_col='Clim_ID')      
    if NL_12.empty:
        NL_12 = empty_frame_NL
    else:
        NL_12 = NL_12.loc[['8400601','8403505','8501106']].Tm
    df_NL = pd.concat([NL_01,NL_02,NL_03,NL_04,NL_05,NL_06,NL_07,NL_08,NL_09,NL_10,NL_11,NL_12], axis=1).T
    months = pd.Series(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']) 
    df_NL.index = pd.to_datetime('15-' + months + '-' + str(current_year)) 
    # NU
    NU_01 = pd.read_csv('~/github/AZMP-NL/external_data/ECCC/climate_summaries/en_climate_summaries_NU_01-' + str(current_year) + '.csv', index_col='Clim_ID')
    if NU_01.empty:
        NU_01 = empty_frame_NU
    else:
        NU_01 = NU_01.loc[['2402592']].Tm
    NU_02 = pd.read_csv('~/github/AZMP-NL/external_data/ECCC/climate_summaries/en_climate_summaries_NU_02-' + str(current_year) + '.csv', index_col='Clim_ID')
    if NU_02.empty:
        NU_02 = empty_frame_NU
    else:
        NU_02 = NU_02.loc[['2402592']].Tm
    NU_03 = pd.read_csv('~/github/AZMP-NL/external_data/ECCC/climate_summaries/en_climate_summaries_NU_03-' + str(current_year) + '.csv', index_col='Clim_ID')
    if NU_03.empty:
        NU_03 = empty_frame_NU
    else:
        NU_03 = NU_03.loc[['2402592']].Tm
    NU_04 = pd.read_csv('~/github/AZMP-NL/external_data/ECCC/climate_summaries/en_climate_summaries_NU_04-' + str(current_year) + '.csv', index_col='Clim_ID')
    if NU_04.empty:
        NU_04 = empty_frame_NU
    else:
        NU_04 = NU_04.loc[['2402592']].Tm
    NU_05 = pd.read_csv('~/github/AZMP-NL/external_data/ECCC/climate_summaries/en_climate_summaries_NU_05-' + str(current_year) + '.csv', index_col='Clim_ID')
    if NU_05.empty:
        NU_05 = empty_frame_NU
    else:
        NU_05 = NU_05.loc[['2402592']].Tm
    NU_06 = pd.read_csv('~/github/AZMP-NL/external_data/ECCC/climate_summaries/en_climate_summaries_NU_06-' + str(current_year) + '.csv', index_col='Clim_ID')
    if NU_06.empty:
        NU_06 = empty_frame_NU
    else:
        NU_06 = NU_06.loc[['2402592']].Tm
    NU_07 = pd.read_csv('~/github/AZMP-NL/external_data/ECCC/climate_summaries/en_climate_summaries_NU_07-' + str(current_year) + '.csv', index_col='Clim_ID')
    if NU_07.empty:
        NU_07 = empty_frame_NU
    else:
        NU_07 = NU_07.loc[['2402592']].Tm
    NU_08 = pd.read_csv('~/github/AZMP-NL/external_data/ECCC/climate_summaries/en_climate_summaries_NU_08-' + str(current_year) + '.csv', index_col='Clim_ID')
    if NU_08.empty:
        NU_08 = empty_frame_NU
    else:
        NU_08 = NU_08.loc[['2402592']].Tm
    NU_09 = pd.read_csv('~/github/AZMP-NL/external_data/ECCC/climate_summaries/en_climate_summaries_NU_09-' + str(current_year) + '.csv', index_col='Clim_ID')
    if NU_09.empty:
        NU_09 = empty_frame_NU
    else:
        NU_09 = NU_09.loc[['2402592']].Tm
    NU_10 = pd.read_csv('~/github/AZMP-NL/external_data/ECCC/climate_summaries/en_climate_summaries_NU_10-' + str(current_year) + '.csv', index_col='Clim_ID')
    if NU_10.empty:
        NU_10 = empty_frame_NU
    else:
        NU_10 = NU_10.loc[['2402592']].Tm
    NU_11 = pd.read_csv('~/github/AZMP-NL/external_data/ECCC/climate_summaries/en_climate_summaries_NU_11-' + str(current_year) + '.csv', index_col='Clim_ID')
    if NU_11.empty:
        NU_11 = empty_frame_NU
    else:
        NU_11 = NU_11.loc[['2402592']].Tm
    NU_12 = pd.read_csv('~/github/AZMP-NL/external_data/ECCC/climate_summaries/en_climate_summaries_NU_12-' + str(current_year) + '.csv', index_col='Clim_ID')      
    if NU_12.empty:
        NU_12 = empty_frame_NU
    else:
        NU_12 = NU_12.loc[['2402592']].Tm
    df_NU = pd.concat([NU_01,NU_02,NU_03,NU_04,NU_05,NU_06,NU_07,NU_08,NU_09,NU_10,NU_11,NU_12], axis=1).T
    df_NU.index = pd.to_datetime('15-' + months + '-' + str(current_year)) 
    
## ---- Read 4 stations of interest (AHCCD) ---- ##
## 1. Bonavista - 8400601
# tmp file without blank space
with open(os.path.expanduser('~/github/AZMP-NL/external_data/ECCC/homog_monthly_mean_temp/mm8400601.txt'), 'r') as f:
    lines = f.readlines()
lines = [line.replace(' ', '') for line in lines]
with open('/tmp/tmp.txt', 'w') as f:
    f.writelines(lines)
df = pd.read_csv('/tmp/tmp.txt', header=2, usecols=['Year', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
df = df.drop(df.index[0]) # Drop French columns
# set Year as index
df = df.set_index('Year')  
# Remove white space in Year values
df = df.set_index(df.index.str.strip())
# Stack months under Years (pretty cool!)
df = df.stack() 
# Transform to a series with values based the 15th of each month
df.index = pd.to_datetime('15-' + df.index.get_level_values(1) + '-' + df.index.get_level_values(0))
# Transform series to numeric
df = pd.to_numeric(df)
# Replace missing value by NaNs
df = df.replace(-9999.9, np.nan)
df_BB = df.copy()
del df
# Append climate summaries if needed
if use_climate_summaries:
    df_BB = pd.concat([df_BB,df_NL['8400601']])


## 2. St. John's - 8403505
# tmp file without blank space
with open(os.path.expanduser('~/github/AZMP-NL/external_data/ECCC/homog_monthly_mean_temp/mm8403505.txt'), 'r') as f:
    lines = f.readlines()
lines = [line.replace(' ', '') for line in lines]
with open('/tmp/tmp.txt', 'w') as f:
    f.writelines(lines)
df = pd.read_csv('/tmp/tmp.txt', header=2, usecols=['Year', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
df = df.drop(df.index[0]) # Drop French columns
# set Year as index
df = df.set_index('Year')  
# Remove white space in Year values
df = df.set_index(df.index.str.strip())
# Stack months under Years (pretty cool!)
df = df.stack() 
# Transform to a series with values based the 15th of each month
df.index = pd.to_datetime('15-' + df.index.get_level_values(1) + '-' + df.index.get_level_values(0))
# Transform series to numeric
df = pd.to_numeric(df)
# Replace missing value by NaNs
df = df.replace(-9999.9, np.nan)
df_SJ = df.copy()
del df
# Append climate summaries if needed
if use_climate_summaries:
    df_SJ = pd.concat([df_SJ,df_NL['8403505']])
    
## 3. Cartwright - 8501106
# tmp file without blank space
with open(os.path.expanduser('~/github/AZMP-NL/external_data/ECCC/homog_monthly_mean_temp/mm8501106.txt'), 'r') as f:
    lines = f.readlines()
lines = [line.replace(' ', '') for line in lines]
with open('/tmp/tmp.txt', 'w') as f:
    f.writelines(lines)
df = pd.read_csv('/tmp/tmp.txt', header=2, usecols=['Year', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
df = df.drop(df.index[0]) # Drop French columns
# set Year as index
df = df.set_index('Year')  
# Remove white space in Year values
df = df.set_index(df.index.str.strip())
# Stack months under Years (pretty cool!)
df = df.stack() 
# Transform to a series with values based the 15th of each month
df.index = pd.to_datetime('15-' + df.index.get_level_values(1) + '-' + df.index.get_level_values(0))
# Transform series to numeric
df = pd.to_numeric(df)
# Replace missing value by NaNs
df = df.replace(-9999.9, np.nan)
df_CA = df.copy()
del df
# Append climate summaries if needed
if use_climate_summaries:
    df_CA = pd.concat([df_CA,df_NL['8501106']])
    
## 4. Iqaluit - 2402592
# tmp file without blank space
with open(os.path.expanduser('~/github/AZMP-NL/external_data/ECCC/homog_monthly_mean_temp/mm2402592.txt'), 'r') as f:
    lines = f.readlines()
lines = [line.replace(' ', '') for line in lines]
with open('/tmp/tmp.txt', 'w') as f:
    f.writelines(lines)
df = pd.read_csv('/tmp/tmp.txt', header=2, usecols=['Year', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
df = df.drop(df.index[0]) # Drop French columns
# set Year as index
df = df.set_index('Year')  
# Remove white space in Year values
df = df.set_index(df.index.str.strip())
# Stack months under Years (pretty cool!)
df = df.stack() 
# Transform to a series with values based the 15th of each month
df.index = pd.to_datetime('15-' + df.index.get_level_values(1) + '-' + df.index.get_level_values(0))
# Transform series to numeric
df = pd.to_numeric(df)
# Replace missing value by NaNs
df = df.replace(-9999.9, np.nan)
df_IQ = df.copy()
del df
# Append climate summaries if needed
if use_climate_summaries:
    df_IQ = pd.concat([df_IQ,df_NU['2402592']])
    
## 5. NUUK - see azmp_dmi_nuukAirT.py
df_NUUK = pd.read_pickle('Nuuk_air_temp.pkl')  

## ---- Concatenate all timeseries ---- ##
df = pd.concat([df_NUUK, df_IQ, df_CA, df_BB, df_SJ], axis=1)
df.columns = ['Nuuk', 'Iqaluit', 'Cartwright', 'Bonavista', 'StJohns']

## ---- Monthly anomalies for current year ---- ##
df_clim_period = df[(df.index.year>=clim_year[0]) & (df.index.year<=clim_year[1])]
df_monthly_stack = df_clim_period.groupby([(df_clim_period.index.year),(df_clim_period.index.month)]).mean()
df_monthly_clim = df_monthly_stack.groupby(level=1).mean()
df_monthly_std = df_monthly_stack.groupby(level=1).std()

df_current_year = df[df.index.year==current_year]
year_index = df_current_year.index # backup index
df_current_year.index=df_monthly_std.index # reset index
anom = df_current_year - df_monthly_clim
std_anom = (df_current_year - df_monthly_clim)/df_monthly_std
#std_anom.index = year_index.month # replace index
std_anom.index = year_index.strftime('%b') # replace index (by text)
anom.index = year_index.strftime('%b') # replace index (by text)

## ---- Annual anomalies (2021 version) ---- ##
# how to select one year: df_monthly_stack.xs(2015,level=0)
df_stack = df.groupby([(df.index.year),(df.index.month)]).mean() 
df_stack_anom = df_stack.sub(df_monthly_clim, level=1)
df_annual_anom = df_stack_anom.groupby(level=0).mean()
# Calculate df_annual* , a "fake" annual mean T
# (used as an annual time series from which anomalies are calculated)
df_annual_star =  df_annual_anom + df_monthly_clim.mean() 
clim = df_annual_star[(df_annual_star.index>=clim_year[0]) & (df_annual_star.index<=clim_year[1])].mean()
std = df_annual_star[(df_annual_star.index>=clim_year[0]) & (df_annual_star.index<=clim_year[1])].std()
anom_annual = (df_annual_star - clim)
std_anom_annual = anom_annual/std

# Save for scorecards
anom_annual.to_pickle('airT_anom.pkl') # for IROC
std_anom_annual.to_pickle('airT_std_anom.pkl')
df.to_pickle('airT_monthly.pkl')

#Load back in if starting here!
std_anom_annual = pd.read_pickle('airT_std_anom.pkl')

# restrict time for following
std_anom_annual = std_anom_annual[std_anom_annual.index>=1950]

## ---- Annual anomalies (old version) ---- ##
## df_annual = df.resample('As').mean()
## clim2 = df_annual[(df_annual.index.year>=clim_year[0]) & (df_annual.index.year<=clim_year[1])].mean()
## std2 = df_annual[(df_annual.index.year>=clim_year[0]) & (df_annual.index.year<=clim_year[1])].std()
## anom_annual2 = (df_annual - clim2)
## anom_annual2.index = anom_annual2.index.year
## std_anom_annual2 = anom_annual2/std2


## ---- plot monthly ---- ##
ax = anom.plot(kind='bar', stacked=True, cmap='YlGn')
plt.grid('on')
ax.set_ylabel(r'[$^{\circ}$C]')
ax.set_title(str(current_year) + ' Air temperature anomalies')
#ax.legend(loc='upper center')
plt.ylim([-18, 15])

fig = ax.get_figure()
fig.set_size_inches(w=9,h=6)
fig_name = 'air_temp_' + str(current_year) + '.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)

# Save French Figure
french_months = ['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D']
ax.set_title(' Anomalies des températures de l\'air - ' + str(current_year))
ax.set_xticklabels(french_months, rotation='horizontal')
fig_name = 'air_temp_' + str(current_year) + '_FR.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)


## ---- plot annual normalized---- ##
#std_anom_annual_norm = std_anom_annual/std_anom_annual.shape[1]
std_anom_annual_norm = std_anom_annual.divide((5 - std_anom_annual.isna().sum(axis=1)).values, axis=0)
n = 5 # xtick every n years
ax = std_anom_annual_norm.plot(kind='bar', stacked=True, cmap='YlGn')
ticks = ax.xaxis.get_ticklocs()
ticklabels = [l.get_text() for l in ax.xaxis.get_ticklabels()]
ax.xaxis.set_ticks(ticks[::n])
ax.xaxis.set_ticklabels(ticklabels[::n])
plt.fill_between([ticks[0]-1, ticks[-1]+1], [-.5, -.5], [.5, .5], facecolor='gray', alpha=.2)
plt.grid('on')
ax.set_ylabel(r'Standardized anomaly')
ax.set_title('Annual air temperature anomalies')
fig = ax.get_figure()
fig.set_size_inches(w=12,h=8)
fig_name = 'air_temp_anom.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)

# Save in French
ax.set_ylabel(r'Anomalie normalisée')
ax.set_title('Anomalies des températures de l\'air')
fig_name = 'air_temp_anom_FR.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)


#### ------------- Stacked airTemp - Red/Blue ---------------- ####
reds = plt.cm.Reds(np.linspace(0,1, num=6))
blues = plt.cm.Blues(np.linspace(0,1, num=6))
colors = np.vstack((blues[1:,:], reds[1:,:]))
colors_ordered = colors[[0,5,1,6,2,7,3,8,4,9],:]

air_stack = pd.concat([std_anom_annual.Nuuk[std_anom_annual.Nuuk<0], std_anom_annual.Nuuk[std_anom_annual.Nuuk>0],
                       std_anom_annual.Iqaluit[std_anom_annual.Iqaluit<0], std_anom_annual.Iqaluit[std_anom_annual.Iqaluit>0],
                       std_anom_annual.Cartwright[std_anom_annual.Cartwright<0], std_anom_annual.Cartwright[std_anom_annual.Cartwright>0],
                       std_anom_annual.Bonavista[std_anom_annual.Bonavista<0], std_anom_annual.Bonavista[std_anom_annual.Bonavista>0],
                       std_anom_annual.StJohns[std_anom_annual.StJohns<0], std_anom_annual.StJohns[std_anom_annual.StJohns>0]],
                       axis=1)/5

fig, ax = plt.subplots(nrows=1, ncols=1)
n = 5 # xtick every n years
sign=air_stack>0
air_stack.plot(kind='bar', stacked=True, color=colors_ordered, zorder=10, ax=ax, legend=False)
ticks = ax.xaxis.get_ticklocs()
ticklabels = [l.get_text() for l in ax.xaxis.get_ticklabels()]
ax.xaxis.set_ticks(ticks[::n])
ax.xaxis.set_ticklabels(ticklabels[::n])
ax.set_ylabel(r'Standardized anomaly')
ax.set_title('Annual air temperature anomalies')
plt.fill_between([ticks[0]-1, ticks[-1]+1], [-.5, -.5], [.5, .5], facecolor='gray', alpha=.2)
plt.grid()
plt.ylim([-2.3,2.3])
#ax.yaxis.set_ticks(np.arange(-2.5, 3, .5))
# Custom legend
import matplotlib.lines as mlines
legend_elements1 = [mlines.Line2D([],[], marker='s',linestyle='None', color=colors_ordered[0], markersize=10, label='\nNuuk'),
                    mlines.Line2D([],[], marker='s',linestyle='None', color=colors_ordered[1], markersize=10, label=''),
                    mlines.Line2D([],[], marker='s',linestyle='None', color=colors_ordered[2], markersize=10, label='\nIqaluit'),
                    mlines.Line2D([],[], marker='s',linestyle='None', color=colors_ordered[3], markersize=10, label=''),
                    mlines.Line2D([],[], marker='s',linestyle='None', color=colors_ordered[4], markersize=10, label='\nCartwright'),
                    mlines.Line2D([],[], marker='s',linestyle='None', color=colors_ordered[5], markersize=10, label='')
                    ]
legend_elements2 = [mlines.Line2D([],[], marker='s',linestyle='None', color=colors_ordered[6], markersize=10, label='\nBonavista'),
                    mlines.Line2D([],[], marker='s',linestyle='None', color=colors_ordered[7], markersize=10, label=''),
                    mlines.Line2D([],[], marker='s',linestyle='None', color=colors_ordered[8], markersize=10, label='\nSt. John''s'),
                    mlines.Line2D([],[], marker='s',linestyle='None', color=colors_ordered[9], markersize=10, label='')
                    ]

    
L2 = ax.legend(handles=legend_elements2, loc=9, bbox_to_anchor=(0.3, 1))
ax.legend(handles=legend_elements1, loc=2)
ax.add_artist(L2)


## ---- plot annual normalized (with scorecards) ---- ##
# preamble
from matplotlib.colors import from_levels_and_colors
# Build the colormap
vmin = -3.49
vmax = 3.49
midpoint = 0
levels = np.linspace(vmin, vmax, 15)
midp = np.mean(np.c_[levels[:-1], levels[1:]], axis=1)
colvals = np.interp(midp, [vmin, midpoint, vmax], [-1, 0., 1])
normal = plt.Normalize(-3.49, 3.49)
reds = plt.cm.Reds(np.linspace(0,1, num=7))
blues = plt.cm.Blues_r(np.linspace(0,1, num=7))
whites = [(1,1,1,1)]*2
colors = np.vstack((blues[0:-1,:], whites, reds[1:,:]))
colors = np.concatenate([[colors[0,:]], colors, [colors[-1,:]]], 0)
cmap, norm = from_levels_and_colors(levels, colors, extend='both')
cmap_r, norm_r = from_levels_and_colors(levels, np.flipud(colors), extend='both')

     
fig, ax = plt.subplots() 
n = 5 # xtick every n years
std_anom_annual_norm.plot(kind='bar', stacked=True, cmap='YlGn', ax=ax)
ticks = ax.xaxis.get_ticklocs()
ticklabels = [l.get_text() for l in ax.xaxis.get_ticklabels()]
ax.xaxis.set_ticks(ticks[::n])
ax.xaxis.set_ticklabels(ticklabels[::n])
plt.fill_between([ticks[0]-1, ticks[-1]+1], [-.5, -.5], [.5, .5], facecolor='gray', alpha=.2)
plt.grid('on')
ax.set_ylabel(r'Standardized anomaly')
ax.set_title('Annual air temperature anomalies')
colors = cmap(normal(np.nansum(std_anom_annual_norm.values, axis=1)))
cell_text = np.nansum(std_anom_annual_norm.values, axis=1).round(1).astype('str')
the_table = ax.table(cellText=[np.nansum(std_anom_annual_norm.values, axis=1).round(1)],
        rowLabels=['Air Temp. subindex'],
        colLabels=None,
        cellColours = [colors],
        cellLoc = 'center', rowLoc = 'center',
        loc='bottom', bbox=[0, -0.14, 1, 0.05])
the_table.auto_set_font_size (False)
the_table.set_fontsize(9)

for key, cell in the_table.get_celld().items():
    if key[1] == -1:
        cell.set_linewidth(0)
        cell.set_fontsize(10)
    else:
        cell._text.set_rotation(90)

fig.set_size_inches(w=13,h=9.5)
fig_name = 'air_temp_climate_index.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim -bordercolor White -border 10x10 ' + fig_name + ' ' + fig_name)


# Save in French
ax.set_ylabel(r'Anomalie normalisée')
ax.set_title('Anomalies des températures de l\'air')
fig_name = 'air_temp_climate_index_FR.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)
