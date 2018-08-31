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

#da_sal = ds_sub['salinity']
#df_sal = da_sal.to_pandas()

# Select Temperature
df = ds['temperature'].to_dataframe()
df = df.unstack()
df = df.resample('1D').mean()
df = df.dropna(how='all')

fig, ax = plt.subplots(nrows=1, ncols=1)
c = plt.contourf(df.index, np.array(df.columns.levels[1], dtype=float), df.values.T, extend='both')
plt.plot(df.index, np.repeat(0, df.index.size), '|k', markersize=20)
plt.ylim([0, 175])

plt.grid('on')
plt.xlabel('Time', fontsize=15, fontweight='bold')
plt.ylabel('Depth (m)', fontsize=15, fontweight='bold')
plt.ylim([0, 175])
plt.gca().invert_yaxis()

cax = fig.add_axes([0.91, .15, 0.01, 0.7])
cb = plt.colorbar(c, cax=cax, orientation='vertical')
cb.set_label(r'$\rm T(^{\circ}C)$', fontsize=12, fontweight='normal')

# Save Figure
fig.set_size_inches(w=12, h=6)
fig.set_dpi(200)
outfile = 'Viking2017_temp.png'
fig.savefig(outfile)
#os.system('convert -trim ' + outfile + ' ' + outfile)

