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
ds = xr.open_dataset('2017_viking.nc')

# Some utils:
# np.unique(ds['instrument_ID'].values)
# np.unique(ds['survey_ID'].values)

#da_sal = ds_sub['salinity']
#df_sal = da_sal.to_pandas()

## ---- Select variables ---- ##
df_temp = ds['temperature'].to_dataframe()
df_temp = df_temp.unstack()
df_temp = df_temp.resample('1D').mean()
df_temp = df_temp.dropna(how='all')
df_sal = ds['salinity'].to_dataframe()
df_sal = df_sal.unstack()
df_sal = df_sal.resample('1D').mean()
df_sal = df_sal.dropna(how='all')
df_sig = ds['sigma-t'].to_dataframe()
df_sig = df_sig.unstack()
df_sig = df_sig.resample('1D').mean()
df_sig = df_sig.dropna(how='all')

Vsig = np.arange(21,27)
Vtemp = np.arange(-2, 20, 2)
Vsal = np.arange(29.5, 34, .5)


## ---- plot temperature ---- ##
fig, ax = plt.subplots(nrows=1, ncols=1)
c = plt.contourf(df_temp.index, np.array(df_temp.columns.levels[1], dtype=float), df_temp.values.T, Vtemp, extend='both', cmap=plt.cm.RdBu_r)
cc = plt.contour(df_sig.index, np.array(df_sig.columns.levels[1], dtype=float), df_sig.values.T, Vsig, colors='k')
plt.plot(df_temp.index, np.repeat(0, df_temp.index.size), '|k', markersize=20)
plt.ylim([0, 175])
plt.clabel(cc, inline=1, fontsize=10, colors='k', fmt='%1.1f')

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


## ---- plot Salinity ---- ##
fig, ax = plt.subplots(nrows=1, ncols=1)
c = plt.contourf(df_sal.index, np.array(df_sal.columns.levels[1], dtype=float), df_sal.values.T, Vsal, extend='both', cmap=plt.cm.RdBu_r)
cc = plt.contour(df_sig.index, np.array(df_sig.columns.levels[1], dtype=float), df_sig.values.T, Vsig, colors='k')
plt.plot(df_sal.index, np.repeat(0, df_sal.index.size), '|k', markersize=20)
plt.ylim([0, 175])
#plt.clabel(cc, inline=1, fontsize=10, colors='k', fmt='%1.1f')

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
outfile = 'Viking2017_sal.png'
fig.savefig(outfile)




#os.system('convert -trim ' + outfile + ' ' + outfile)
