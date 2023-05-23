'''
AZMP reporting - Headlands

Used as a test script.
Working parts were transfered into hl_module, e.g.:
-> hl.annual_plot_anom('Stock_Cove.nc', "Stock Cove", 2021) 


WORK IN PROGRESS

(script ran in /home/cyrf0006/research/Headlands)
* Should go in AZMP


Frederic.Cyr@dfo-mpo.gc.ca - February 2022


'''

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import hl_modules as hl  
import xarray as xr

# Adjust fontsize/weight
font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 10}
plt.rc('font', **font)

clim_years = [1991, 2020]
current_year = 2021
data_file = 'Arnolds_Cove.nc'
anom_months = [5,6,7,8,9,10]



## ---- Load and prepare the data ---- ##
ds = xr.open_dataset(data_file)   
df = ds.temperature.to_pandas()

## df_clim = df[(df.index.year>=clim_years[0]) & (df.index.year<=clim_years[1])]
## df_monthly_clim = df_clim.groupby(df_clim.index.month).mean() 
## df_monthly_std = df_clim.groupby(df_clim.index.month).std() 
## df_year = df[df.index.year == current_year] # >
## df_monthly = df_year.groupby(df_year.index.month).mean()

# monthly anom for scorecards
df_stack = df.groupby([(df.index.year),(df.index.month)]).mean().squeeze()
df_unstack = df_stack.unstack()
df_clim_period = df[(df.index.year>=clim_years[0]) & (df.index.year<=clim_years[1])]
df_clim_stack = df_clim_period.groupby([(df_clim_period.index.year),(df_clim_period.index.month)]).mean().squeeze()
df_monthly_clim = df_clim_stack.mean(level=1)
df_monthly_std = df_clim_stack.std(level=1)
monthly_anom = df_unstack - df_monthly_clim 
monthly_stdanom = (df_unstack - df_monthly_clim) /  df_monthly_std
# Restricts months
monthly_stdanom = monthly_stdanom[anom_months]
# Threshold on no. of months to consider anomaly valid (50%)
monthly_stdanom.loc[monthly_stdanom.count(axis=1)<np.ceil(len(anom_months)/2)]=np.nan
# annual mean std anomaly
anom_std = monthly_stdanom[anom_months].mean(axis=1)

# Add mean and std
anom_std.at['MEAN'] = df_monthly_clim.loc[anom_months].mean()
anom_std.at['SD'] = df_monthly_std.loc[anom_months].mean()
