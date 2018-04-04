# AZMP reporting - Air temperature from Colbourne's Excel sheets
#  check in /home/cyrf0006/research/AZMP_stateReports/2017

## ---- Link info ---- ##
# There are two types of product.
# Homogenized Temperature is the best since corrections are applied:
# ftp://ccrp.tor.ec.gc.ca/pub/AHCCD/Homog_monthly_mean_temp.zip
#
# Since these are available with delay, it is sometimes necessary to use standard monthly temperature:
# http://climate.weather.gc.ca/prods_servs/cdn_climate_summary_e.html
# http://climate.weather.gc.ca/prods_servs/cdn_climate_summary_report_e.html?intYear=2018&intMonth=2&prov=NL&dataFormat=csv&btnSubmit=Download+data
# !!! Need to write a shell script / crontab to update these files automatically !!!


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import datetime
import gsw
from seawater import extras as swx
import matplotlib.dates as mdates
from scipy.interpolate import griddata

# Adjust fontsize/weight
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 18}
plt.rc('font', **font)

## ----  Load data from Excel sheet ---- ##
## NEED TO EDIT BY HAND TO ADD A COMMA AFTER AUTUMN and AUTOMNE (!!! TRY TO AUTOMIZED THIS)
df = pd.read_csv('/home/cyrf0006/Data/AZMP/AirT/mm8403505.txt', header=2)

# Drop French columns (Canadian bilinguism even in datafile!)
df = df.drop(df.index[0]) # Drop French columns

# Drop dump white space in column names
df = df.rename(columns=lambda x: x.strip()) 

# set Year as index
df = df.set_index('Year')  

# Keep only months
df = df[['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']]

# Remove white space in Year values
df = df.set_index(df.index.str.strip())

# Stack months under Years (pretty cool!)
df = df.stack() 

# Transform to a series with values based the 15th of each month
df.index = pd.to_datetime('15-' + df.index.get_level_values(1) + '-' + df.index.get_level_values(0))

# Transform series to numeric
df = pd.to_numeric(df)

# Replace missing value by NaNs
df = df.replace(df[df<-9999], np.nan)
## ------------------------------------- ##


# raw plot
df.plot()
plt.show()

# Annual anomalies
df_annual = df.resample('A').mean()
df_annual = df_annual['1950':'2018']
df_anom = df_annual - df_annual['1981':'2010'].mean()
df_std_anom = (df_annual - df_annual['1981':'2010'].mean()) / df_annual['1981':'2010'].std()

df_std_anom.plot()
plt.grid()
plt.ylim([-2,2])
plt.show()


# Monthly annomalies
