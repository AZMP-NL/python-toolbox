'''

This is to extract historical air temperature data from NUUK based on :
JOURNAL OF GEOPHYSICAL RESEARCH, VOL. 111, D11105, doi:10.1029/2005JD006810, 2006
https://www.dmi.dk/fileadmin/user_upload/Rapporter/TR/2018/DMIRep18-05.zip
https://www.dmi.dk/publikationer/

Frederic.Cyr@dfo-mpo.gc.ca - February 2020

'''


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import datetime
import matplotlib.dates as mdates
from scipy.interpolate import griddata
import os


# Read historical data from Vinther et al. 2006 (JGR 111)
df = pd.read_csv('~/github/AZMP-NL/external_data/DMI/gr_monthly_all_1784_2017.csv', delimiter=";", decimal=",")
# keep only air Temp at Nuuk
df = df[(df.stat_no==4250) & (df.elem_no==101)]
df.index = df.year
df = df[['jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec']]

# Convert to timeseries index
df = df.stack() 
df_dates = '15-' + df.index.get_level_values(1) + '-' + df.index.get_level_values(0).values.astype(str)
df.index = [pd.to_datetime(i) for i in df_dates]

# For year 2021, I downloaded: https://www.dmi.dk/fileadmin/Rapporter/2021/DMIRep21-08_new_dataformat_2014_2020.zip
# But haven't analysed it.

# Read more recent measures from https://www.dmi.dk/vejrarkiv/ (select Gronland - Nuuk - Maneder - 2020 + download)
df_14 = pd.read_csv('~/github/AZMP-NL/external_data/DMI/nuuk-2014.csv', delimiter=";")
df_15 = pd.read_csv('~/github/AZMP-NL/external_data/DMI/nuuk-2015.csv', delimiter=";")
df_16 = pd.read_csv('~/github/AZMP-NL/external_data/DMI/nuuk-2016.csv', delimiter=";")
df_17 = pd.read_csv('~/github/AZMP-NL/external_data/DMI/nuuk-2017.csv', delimiter=";")
df_18 = pd.read_csv('~/github/AZMP-NL/external_data/DMI/nuuk-2018.csv', delimiter=";")
df_19 = pd.read_csv('~/github/AZMP-NL/external_data/DMI/nuuk-2019.csv', delimiter=";")
df_20 = pd.read_csv('~/github/AZMP-NL/external_data/DMI/nuuk-2020.csv', delimiter=";")
df_21 = pd.read_csv('~/github/AZMP-NL/external_data/DMI/nuuk-2021.csv', delimiter=";")
df_22 = pd.read_csv('~/github/AZMP-NL/external_data/DMI/nuuk-2022.csv', delimiter=";")
df_23 = pd.read_csv('~/github/AZMP-NL/external_data/DMI/nuuk-2023.csv', delimiter=",")

df_all = pd.concat([df_14.Middel, df_15.Middel, df_16.Middel, df_17.Middel, df_18.Middel, df_19.Middel, df_20.Middel, df_21.Middel, df_22.Middel, df_23.Middel], axis=1, keys=np.arange(2014, 2024))

df_recent = df_all.T
df_recent.columns = ['jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec']
df_recent = df_recent.stack()
df_recent_dates = '15-' + df_recent.index.get_level_values(1) + '-' + df_recent.index.get_level_values(0).values.astype(str)
df_recent.index = [pd.to_datetime(i) for i in df_recent_dates]


# Plot 2 timeseries
df.plot()
df_recent.plot()

# Merge and save
df_merged = pd.concat([df, df_recent], axis=1)
df_merged = df_merged.mean(axis=1)

df_merged.to_pickle('Nuuk_air_temp.pkl')
