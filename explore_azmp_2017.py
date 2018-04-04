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
# in /home/cyrf0006/research/AZMP_database/netCDF_first_set
ds = xr.open_mfdataset('2017*.nc')
zVec = ds.level.to_masked_array()

# Some utils:
# np.unique(ds['instrument_ID'].values)
# np.unique(ds['survey_ID'].values)

# Select Survey and switch to DataFrame
ds_tel173 = ds.sel(time=ds['survey_ID'].values=='39173')
df_tel173 = ds_tel173.to_dataframe()

ds_tel176 = ds.sel(time=ds['survey_ID'].values=='39176')
df_tel176 = ds_tel176.to_dataframe()

## ----  plot FC-line ---- ##
df_tel173_FC = df_tel173[df_tel173['comments'].str.contains("FC-")]
sr_temp = df_tel173_FC['temperature']
sr_lon = df_tel173_FC['longitude']
sr_lat = df_tel173_FC['latitude']
sr_stn = df_tel173_FC['comments']

df_stn = sr_stn.unstack()
df_temp = sr_temp.unstack()
df_lon = sr_lon.unstack()

df_temp.columns = df_lon.values[0,]

plt.contourf(df_temp.columns, df_temp.index, df_temp, 20, cmap=plt.cm.RdBu_r)
plt.ylim([0, 500])
plt.gca().invert_yaxis()
plt.show()

# add bathymetry
bb_bathy = '/home/cyrf0006/github/AZMP-NL/bathymetry/bottom_profiles/bbline.txt'
from numpy import loadtxt
bathy = loadtxt(bb_bathy, delimiter=",", unpack=False)
bathy_x = bathy[:,0]/1000.0
bathy_y = np.abs(bathy[:,1])

bathy_x_close = np.append(bathy_x, np.array([bathy_x[-1], bathy_x[0], bathy_x[0]]))
bathy_y_close = np.append(bathy_y ,np.array([bathy_y.max(), bathy_y.max(), bathy_y[0]]))
bathymetry = zip(bathy_x_close, bathy_y_close)

# Check maximum depth
cast_depth = []
for i in distance:
    min_idx = np.argmin(np.abs(i-bathy_x))
    cast_depth = np.append(cast_depth, bathy_y[min_idx])


Bgon = plt.Polygon(bathymetry,color=np.multiply([1,.9333,.6667],.4), alpha=1)
ax.add_patch(Bgon)
for i in range(0,len(distance)):
    plt.plot(np.array([distance[i], distance[i]]), np.array([Pbin[0], cast_depth[i]]), '--k', linewidth=0.1)


    

fig.set_size_inches(w=9,h=9)
fig_name = 'AZMP_surveys_2017.png'
fig.set_dpi(300)
fig.savefig(fig_name)



