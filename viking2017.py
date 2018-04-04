'''
This is a test with xarray, see if I can manipulate a lot of netCDF files.
Try this script in /home/cyrf0006/research/AZMP_database/2017_data
'''

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr
import datetime
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from scipy.interpolate import griddata # here should remove nans or empty profiles
import netCDF4
#import matplotlib
#matplotlib.interactive(True)

# For plots
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 14}
plt.rc('font', **font)

# This is a dataset
# in /home/cyrf0006/research/AZMP_database/2017_data
ds = xr.open_mfdataset('2017_viking.nc')

# Some utils:
# np.unique(ds['instrument_ID'].values)
# np.unique(ds['survey_ID'].values)

# Select Temperature
df = ds['temperature'].to_dataframe()
df = df.unstack()
plt.contourf(df.index, np.array(df.columns.levels[1], dtype=float), df.values.T)
plt.gca().invert_yaxis()
plt.show()
