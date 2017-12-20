'''
This is a test with xarray, see if I can manipulate a lot of netCDF files.
Try this script in /home/cyrf0006/research/AZMP_database
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

# This is a dataset
ds = xr.open_mfdataset('*.nc')


# Select a depth range
## ds = ds.sel(level=ds['level']<500)
## ds = ds.sel(level=ds['level']>10)

## # Selection of a subset region
ds = ds.where((ds.longitude>-58) & (ds.longitude<-46), drop=True)
## #ds = ds.where((ds.latitude>50) & (ds.latitude<55), drop=True)
ds = ds.where((ds.latitude>42) & (ds.latitude<56), drop=True)

# Select only one year
ds = ds.sel(time=ds['time.year']==2008)
# Select only summer
ds = ds.sel(time=ds['time.season']=='SON')


## --- Map temperature --- ##

# Temperature to pandas
da_temp = ds.temperature
df_temp = da_temp.to_pandas()

# Find depth range to averae
idx = np.squeeze(np.where((df_temp.columns>=50) & (df_temp.columns<=150)))
#idx = np.squeeze(np.where((df_temp.columns>=10) & (df_temp.columns<=50)))

# All 3 arrays
temp = np.array(df_temp[df_temp.columns[idx]].mean(axis=1))
lons = np.array(ds.longitude)
lats = np.array(ds.latitude)

# Meshgrid 1D data (after removing NaNs)
idx = np.argwhere(~np.isnan(temp))
x = np.arange(np.min(lons[idx]), np.max(lons[idx]), .1)
y = np.arange(np.min(lats[idx]), np.max(lats[idx]), .1)  
lon, lat = np.meshgrid(x,y)

# griddata
LN = np.squeeze(lons[idx])
LT = np.squeeze(lats[idx])
TT = np.squeeze(temp[idx])
temp_itp = griddata((LN, LT), TT, (lon, lat), method='linear')

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


## ---- Plot ---- ## 
lon_0 = lons.mean()
lat_0 = lats.mean()
lon_0 = -50
lat_0 = 50
lonLims = [-58, -46]
latLims = [42, 56]
proj = 'merc'
m = Basemap(projection='merc',lon_0=lon_0,lat_0=lat_0, llcrnrlon=lonLims[0],llcrnrlat=latLims[0],urcrnrlon=lonLims[1],urcrnrlat=latLims[1], resolution='l')

x,y = m(*np.meshgrid(lonz,latz))
c = m.contour(x, y, np.flipud(-Z), v, colors='darkgrey');
#c = m.contourf(x, y, np.flipud(Z), v, cmap=plt.cm.PuRd_r, extend="min");
m.fillcontinents(color='grey');

#v = np.arange(np.floor(np.min(tmax)), np.ceil(np.max(tmax))+1)
v = 20
xi, yi = m(lon, lat)
lon_casts, lat_casts = m(lons[idx], lats[idx])
cs = m.contourf(xi,yi,temp_itp, v)

m.fillcontinents(color='grey');

# Add Colorbar
cbar = m.colorbar(cs, location='bottom')
#cbar.set_label(tmax_units)

# Add casts identification
m.plot(lon_casts, lat_casts, '.k')

# Add Grid Lines
m.drawparallels(np.arange(latLims[0], latLims[1], 1.), labels=[1,0,0,0], fontsize=10)
m.drawmeridians(np.arange(lonLims[0], lonLims[1], 1.), labels=[0,0,0,1], fontsize=10)

# Add Coastlines, States, and Country Boundaries
m.drawcoastlines()
m.drawstates()
m.drawcountries()

# Add Title
plt.title('Temperature')

plt.show()





plt.contourf(df_temp.index, df_temp.columns, df_temp.T, 20)


# HERE!!


# this is now a dataArray
da_temp = ds['temperature']
da_lat = ds['latitude']
da_lon = ds['longitude']


# Averrage profile of the whole timeseries
A = da_temp.mean(dim='time')

# Selection of only summer data
summer_temp = temp.sel(time=temp['time.season']=='JJA')
summer_lat = lat.sel(time=lat['time.season']=='JJA')
summer_lon = lon.sel(time=lon['time.season']=='JJA')




ds = xr.Dataset()
ds['data1'] = xr.DataArray(np.arange(100), coords={'t1': np.linspace(0, 1, 100)})
ds['data1b'] = xr.DataArray(np.arange(100, 200), coords={'t1': np.linspace(0, 1, 100)})
ds['data2'] = xr.DataArray(np.arange(200, 500), coords={'t2': np.linspace(0, 1, 300)})
ds['data2b'] = xr.DataArray(np.arange(600, 900), coords={'t2': np.linspace(0, 1, 300)})
ds.where(ds.data1 < 50, drop=True)



