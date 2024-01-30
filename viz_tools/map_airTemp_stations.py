'''
A map for AZMP-NL air temperature stations

'''


import netCDF4
from mpl_toolkits.basemap import Basemap
import numpy as  np
import matplotlib.pyplot as plt
import openpyxl, pprint
import pandas as pd
import xarray as xr
import cmocean
import os

## ---- Region parameters ---- ##
dataFile = '/home/cyrf0006/data/GEBCO/GRIDONE_1D.nc'
AZMP_airTemp_file = '/home/cyrf0006/github/AZMP-NL/utils/airTemp_sites.xlsx'
lon_0 = -50
lat_0 = 50
lonLims = [-72, -35]
latLims = [45, 66]
#lonLims = [-80, -30]
#latLims = [30, 61]
proj= 'merc'
decim_scale = 10
fig_name = 'map_air_temp.png'

## ---- Bathymetry ---- ####
#v = np.linspace(-3500, 0, 36)
v = np.linspace(0, 3500, 36)


# Load data
dataset = netCDF4.Dataset(dataFile)

# Extract variables
x = dataset.variables['x_range']
y = dataset.variables['y_range']
spacing = dataset.variables['spacing']

# Compute Lat/Lon
nx = int((x[-1]-x[0])/spacing[0]) + 1  # num pts in x-dir
ny = int((y[-1]-y[0])/spacing[1]) + 1  # num pts in y-dir

lon = np.linspace(x[0],x[-1],nx)
lat = np.linspace(y[0],y[-1],ny)

# Reshape data
zz = dataset.variables['z']
Z = zz[:].reshape(ny, nx)

# Reduce data according to Region params and decim scale
idx_lon = np.where((lon>=lonLims[0]) & (lon<=lonLims[1]))
idx_lat = np.where((lat>=latLims[0]) & (lat<=latLims[1]))
Z = Z[idx_lat[0][0]:idx_lat[0][-1]+1, idx_lon[0][0]:idx_lon[0][-1]+1]
lon = lon[idx_lon[0]]
lat = lat[idx_lat[0]]

lon = lon[::decim_scale]
lat = lat[::decim_scale]
Z = Z[::decim_scale, ::decim_scale]

## ---- AZMP air temp Stations  ---- ##
df = pd.read_excel(AZMP_airTemp_file)
air_stationLat = df['latitude'].values
air_stationLon = df['longitude'].values

## ---- plot ---- ##
fig, ax = plt.subplots(nrows=1, ncols=1)
m = Basemap(ax=ax, projection='merc',lon_0=lon_0,lat_0=lat_0, llcrnrlon=lonLims[0],llcrnrlat=latLims[0],urcrnrlon=lonLims[1],urcrnrlat=latLims[1], resolution='l')

x,y = m(*np.meshgrid(lon,lat))
c = m.contourf(x, y, np.flipud(-Z), v, cmap=cmocean.cm.deep, extend="max");
cc = m.contour(x, y, np.flipud(-Z), [100, 500, 1000, 3000, 4000], colors='lightgrey', linewidths=.5);
plt.clabel(cc, inline=1, fontsize=10, colors='gray', fmt='%d')
ccc = m.contour(x, y, np.flipud(-Z), [0], colors='black');
m.fillcontinents(color='peru');
m.drawparallels(np.arange(10,70,10), labels=[1,0,0,0], fontsize=12, fontweight='bold');
m.drawmeridians(np.arange(-80, 5, 10), labels=[0,0,0,1], fontsize=12, fontweight='bold');
#plt.title('Standard air temperature sites')


# add stations
x, y = m(air_stationLon,air_stationLat)
plt.text(x[0], y[0], '  St. John''s', horizontalalignment='left', verticalalignment='center', fontsize=14, color='r', fontweight='bold')
plt.text(x[1], y[1], '  Bonavista', horizontalalignment='left', verticalalignment='center', fontsize=14, color='r', fontweight='bold')
plt.text(x[2], y[2], '   Cartwright', horizontalalignment='left', verticalalignment='center', fontsize=14, color='r', fontweight='bold')
plt.text(x[3], y[3], '  Iqaluit', horizontalalignment='left', verticalalignment='center', fontsize=14, color='r', fontweight='bold')
plt.text(x[4], y[4], 'Nuuk  ', horizontalalignment='right', verticalalignment='center', fontsize=14, color='r', fontweight='bold')

#plt.annotate('St. John''s', xy=(x[0], y[0]),  xycoords='data', xytext=(x[0], y[0]), color='r')

m.scatter(x,y,100,marker='*',color='r', zorder=10)


#### ---- Save Figure ---- ####
fig.set_size_inches(w=9, h=9)
#fig.tight_layout() 
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)

