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
from scipy.interpolate import griddata # here should remove nans or empty profiles
import netCDF4
import os
from matplotlib import dates

#import matplotlib
#matplotlib.interactive(True)
 
# For plots
font = {'family' : 'times',
        'weight' : 'normal',
        'size'   : 14}
plt.rc('font', **font)

# This is a dataset
# in /home/cyrf0006/research/AZMP_database/2017_data
ds = xr.open_dataset('/home/cyrf0006/data/dev_database/viking_nc/2019_viking.nc')

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

df_fluo = ds['fluorescence'].to_dataframe()
df_fluo = df_fluo.unstack()
df_fluo = df_fluo.resample('1D').mean()
df_fluo = df_fluo.dropna(how='all')

Vsig = np.arange(21,27)
Vtemp = np.arange(-2, 20, .5)
Vsal = np.arange(29.5, 34, .25)
XLIM = [datetime.date(2019, 6, 1), datetime.date(2019, 11, 1)]


## ---- plot temperature ---- ##
fig, ax = plt.subplots(nrows=1, ncols=1)
c = plt.contourf(df_temp.index, np.array(df_temp.columns.levels[1], dtype=float), df_temp.values.T, Vtemp, extend='both', cmap=plt.cm.RdBu_r)
cc = plt.contour(df_sig.index, np.array(df_sig.columns.levels[1], dtype=float), df_sig.values.T, Vsig, colors='k')
plt.plot(df_temp.index, np.repeat(0, df_temp.index.size), '|k', markersize=20)
plt.ylim([0, 175])
plt.xlim([XLIM[0], XLIM[1]])
plt.clabel(cc, inline=1, fontsize=10, colors='k', fmt='%1.1f')

plt.grid('on')
#plt.xlabel('Time', fontsize=15, fontweight='bold')
plt.ylabel('Depth (m)', fontsize=15, fontweight='bold')
plt.ylim([0, 175])
plt.gca().invert_yaxis()

cax = fig.add_axes([0.91, .15, 0.01, 0.7])
cb = plt.colorbar(c, cax=cax, orientation='vertical')
cb.set_label(r'$\rm T(^{\circ}C)$', fontsize=12, fontweight='normal')

# Save Figure
outfile = 'Viking2019_temp.png'
fig.savefig(outfile, dpi=200)
os.system('convert -trim ' + outfile + ' ' + outfile)

# Save Figure in French
outfile = 'Viking2019_temp_FR.png'
plt.ylabel('Profondeur (m)', fontsize=15, fontweight='bold')
fig.savefig(outfile, dpi=200)
os.system('convert -trim ' + outfile + ' ' + outfile)


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
plt.xlim([XLIM[0], XLIM[1]])
plt.gca().invert_yaxis()

cax = fig.add_axes([0.91, .15, 0.01, 0.7])
cb = plt.colorbar(c, cax=cax, orientation='vertical')
cb.set_label(r'$\rm S$', fontsize=12, fontweight='normal')

# Save Figure
fig.set_size_inches(w=12, h=6)
fig.set_dpi(200)
outfile = 'Viking2019_sal.png'
fig.savefig(outfile)
os.system('convert -trim ' + outfile + ' ' + outfile)

# Save Figure in French
outfile = 'Viking2019_sal_FR.png'
plt.ylabel('Profondeur (m)', fontsize=15, fontweight='bold')
fig.savefig(outfile, dpi=200)
os.system('convert -trim ' + outfile + ' ' + outfile)

## ----  plot both as suubplots ---- ##

fig = plt.figure()
# ax1
ax = plt.subplot2grid((2, 1), (0, 0))
c = plt.contourf(df_temp.index, np.array(df_temp.columns.levels[1], dtype=float), df_temp.values.T, Vtemp, extend='both', cmap=plt.cm.RdBu_r)
cc = plt.contour(df_sig.index, np.array(df_sig.columns.levels[1], dtype=float), df_sig.values.T, Vsig, colors='k')
plt.plot(df_temp.index, np.repeat(0, df_temp.index.size), '|k', markersize=20)
plt.ylim([0, 175])
plt.xlim([XLIM[0], XLIM[1]])
plt.clabel(cc, inline=1, fontsize=10, colors='k', fmt='%1.1f')
plt.grid('on')
plt.ylabel('Depth (m)', fontsize=15, fontweight='bold')
plt.ylim([0, 175])
plt.gca().invert_yaxis()
cax = fig.add_axes([0.91, .55, 0.01, 0.32])
cb = plt.colorbar(c, cax=cax, orientation='vertical')
ax.xaxis.label.set_visible(False)
ax.tick_params(labelbottom='off')
cb.set_label(r'$\rm T(^{\circ}C)$', fontsize=12, fontweight='normal')

# ax2
ax2 = plt.subplot2grid((2, 1), (1, 0))
c = plt.contourf(df_sal.index, np.array(df_sal.columns.levels[1], dtype=float), df_sal.values.T, Vsal, extend='both', cmap=plt.cm.RdBu_r)
cc = plt.contour(df_sig.index, np.array(df_sig.columns.levels[1], dtype=float), df_sig.values.T, Vsig, colors='k')
plt.plot(df_sal.index, np.repeat(0, df_sal.index.size), '|k', markersize=20)
plt.ylim([0, 175])
plt.grid('on')
#plt.xlabel('Time', fontsize=15, fontweight='bold')
plt.ylabel('Depth (m)', fontsize=15, fontweight='bold')
plt.ylim([0, 175])
plt.xlim([XLIM[0], XLIM[1]])
plt.gca().invert_yaxis()
cax2 = fig.add_axes([0.91, .12, 0.01, 0.32])
cb = plt.colorbar(c, cax=cax2, orientation='vertical')
cb.set_label(r'$\rm S$', fontsize=12, fontweight='normal')

# Save Figure
fig.set_size_inches(w=12, h=12)
outfile = 'Viking2019.png'
fig.savefig(outfile, dpi=200)
os.system('convert -trim ' + outfile + ' ' + outfile)

# French Figure
fig.set_size_inches(w=12, h=12)
ax.set_ylabel('Profondeur (m)', fontsize=15, fontweight='bold')
ax2.set_ylabel('Profondeur (m)', fontsize=15, fontweight='bold')
outfile = 'Viking2019_FR.png'
fig.savefig(outfile, dpi=200)
os.system('convert -trim ' + outfile + ' ' + outfile)

