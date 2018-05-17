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
ds = xr.open_mfdataset('/home/cyrf0006/Data/AZMP/dev_database/2017.nc')
zVec = ds.level.to_masked_array()

# Selection of a subset region
ds = ds.where((ds.longitude>-70) & (ds.longitude<-40), drop=True)
ds = ds.where((ds.latitude>40) & (ds.latitude<70), drop=True)

# Sort time dimension (this takes time to display!!)
#ds = ds.isel(time=np.argsort(ds.time))

# Selection according to season (xarray seasons)
ds_spring = ds.sel(time=ds['time.season']=='MAM')
df_spring = ds_spring.to_dataframe()
ds_summer = ds.sel(time=ds['time.season']=='JJA')
df_summer = ds_summer.to_dataframe()
ds_fall = ds.sel(time=ds['time.season']=='SON')
df_fall = ds_fall.to_dataframe()

df = ds.to_dataframe()

## ---- Bathymetry ---- ####
dataFile = '/home/cyrf0006/Data/GEBCO/GRIDONE_1D.nc'
decim_scale = 4
v = np.linspace(0, 4000, 9)

# Load data
dataset = netCDF4.Dataset(dataFile)

# Extract variables
x = dataset.variables['x_range']
y = dataset.variables['y_range']
spacing = dataset.variables['spacing']

# Compute Lat/Lon
nx = int((x[-1]-x[0])/spacing[0]) + 1  # num pts in x-dir
ny = int((y[-1]-y[0])/spacing[1]) + 1  # num pts in y-dir

lonz = np.linspace(x[0],x[-1],nx)
latz = np.linspace(y[0],y[-1],ny)

# Reshape data
zz = dataset.variables['z']
Z = zz[:].reshape(ny, nx)

# Reduce data according to Region params
lonz = lonz[::decim_scale]
latz = latz[::decim_scale]
Z = Z[::decim_scale, ::decim_scale]


## ---- Plot stations position on map ---- ##
lat_spring = df_spring['latitude'].values
lon_spring = df_spring['longitude'].values
lat_summer = df_summer['latitude'].values
lon_summer = df_summer['longitude'].values
lat_fall = df_fall['latitude'].values
lon_fall = df_fall['longitude'].values

lat_all = df['latitude'].values
lon_all = df['longitude'].values

lon_0 = -50
lat_0 = 50

lonLims = [-70, -40]
latLims = [40, 65]
proj = 'merc'

fig = plt.figure()

m = Basemap(projection='merc',lon_0=lon_0,lat_0=lat_0, llcrnrlon=lonLims[0],llcrnrlat=latLims[0],urcrnrlon=lonLims[1],urcrnrlat=latLims[1], resolution='l')

# Bathymetry
x,y = m(*np.meshgrid(lonz,latz))
v = np.linspace(-4000, 0, 9)
c = m.contourf(x, y, np.flipud(Z), v, cmap=plt.cm.PuBu_r, extend="min");
m.fillcontinents(color='grey');

# Add Colorbar
cb = plt.colorbar(c, orientation='vertical')
cb.set_label('Depth (m)', fontSize=15, fontWeight='bold')

# Add casts identification
x, y = m(lon_spring, lat_spring)
pts1 = m.scatter(x,y, s=50, marker='o',color='seagreen', label='Spring')
#pts.set_facecolor([(.18039, .5451, .3412, .3)])
## pts.set_facecolor([(0, 145/255., 106/255., .5)])
x, y = m(lon_summer, lat_summer)
pts2 = m.scatter(x,y,s=25, marker='o',color='orange', label='Summer')
#pts.set_facecolor([(1.0, 0.6470588235294118, 0.0, .3)])
x, y = m(lon_fall, lat_fall)
pts3 = m.scatter(x,y,s=10, marker='o',color='red', label='Fall')

plt.legend(handles=[pts1, pts2, pts3])

# Add Grid Lines
m.drawparallels(np.arange(latLims[0], latLims[1], 10.), labels=[1,0,0,0], fontsize=10)
m.drawmeridians(np.arange(lonLims[0], lonLims[1], 10.), labels=[0,0,0,1], fontsize=10)

# Add Coastlines, States, and Country Boundaries
m.drawcoastlines()
m.drawstates()
m.drawcountries()

# Add Title
plt.title('Available data 2017 (multi-sources)', fontweight='bold')


fig.set_size_inches(w=9,h=9)
fig_name = 'AZMP_2017_all.png'
fig.set_dpi(300)
fig.savefig(fig_name)



