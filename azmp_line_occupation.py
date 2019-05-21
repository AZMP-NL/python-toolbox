'''
This is a script used to extract long-term data from Station 27

Tried this:
'stn27' # 494
'S27' # 2
's27' # 0
'S27-01' # 20
's27-01' # 0
'STN27' # 528
'STN-27' # 102
'stn-27' # 1
'station27' # 140
'Station27' # 71
'STATION27' # 207
'Stat27' # 1
'stat27' # 1
'STAT27' # 4
'STA27' # 2
'Sta27' # 0
'sta27' # 3

'''

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr
import datetime
import os

#import water_masses as wm

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 14}
plt.rc('font', **font)


ds = xr.open_mfdataset('/home/cyrf0006/data/dev_database/*.nc')

# Select a depth range
ds = ds.sel(level=ds['level']<200)
ds = ds.sel(level=ds['level']>0)


da = ds['comments']
df = da.to_pandas()

df = df[df.index.year>=1950]

df_MB = df[df.values == 'MB-01']
df_BI = df[df.values == 'BI-02']

