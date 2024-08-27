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

# For plots
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 14}
plt.rc('font', **font)

# This is a dataset
# in /home/cyrf0006/research/AZMP_database/netCDF_first_set
ds = xr.open_mfdataset('2017.nc')


# Select only summer
ds = ds.sel(time=ds['time.season']=='JJA')


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

lonLims = [-70, -40]
latLims = [40, 65]
proj = 'merc'



fig = plt.figure()

m = Basemap(projection='merc',lon_0=lon_0,lat_0=lat_0, llcrnrlon=lonLims[0],llcrnrlat=latLims[0],urcrnrlon=lonLims[1],urcrnrlat=latLims[1], resolution='l')

x,y = m(*np.meshgrid(lonz,latz))
c = m.contour(x, y, np.flipud(-Z), v, colors='darkgrey');
#c = m.contourf(x, y, np.flipud(Z), v, cmap=plt.cm.PuRd_r, extend="min");
m.fillcontinents(color='grey');

#v = np.arange(np.floor(np.min(tmax)), np.ceil(np.max(tmax))+1)
v = np.linspace(-2, 24, 56)
xi, yi = m(lon, lat)
lon_casts, lat_casts = m(lons[idx], lats[idx])
cs = m.contourf(xi,yi,temp_itp, v, cmap=plt.cm.RdBu_r)

m.fillcontinents(color='grey');

# Add Colorbar
cbar = m.colorbar(cs, ticks=np.linspace(0, 22, 12), location='right')
#cbar.set_label(tmax_units)

# Add casts identification
x, y = m(lons, lats)
m.scatter(x,y,10,marker='o',color='k', alpha=.5)
#m.plot(x, y, '.k')

# Add Grid Lines
m.drawparallels(np.arange(latLims[0], latLims[1], 10.), labels=[1,0,0,0], fontsize=10)
m.drawmeridians(np.arange(lonLims[0], lonLims[1], 10.), labels=[0,0,0,1], fontsize=10)

# Add Coastlines, States, and Country Boundaries
m.drawcoastlines()
m.drawstates()
m.drawcountries()

# Add Title
plt.title('Temperature (50-150m)', fontweight='bold')



fig.set_size_inches(w=9,h=9)
fig_name = 'TCIL_2017.png'
fig.savefig(fig_name)



