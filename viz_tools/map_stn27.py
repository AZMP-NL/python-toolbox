'''
   map of the North Atlantic. Originally realized for MOPGA project.
'''

import netCDF4
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import numpy as  np
import matplotlib.pyplot as plt
import openpyxl, pprint


## ---- Region parameters ---- ##
dataFile = '/home/cyrf0006/Data/GEBCO/GRIDONE_1D.nc'
lonLims = [-55, -50]
latLims = [45, 50]
lon_0 = np.mean(lonLims)
lat_0 = np.mean(latLims)
proj = 'merc'
decim_scale = 1
fig_name = 'map_stn27.png'

## ---- Bathymetry ---- ####
v = np.linspace(0, 500, 11)
v2 = np.array([200,300])

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
Z = np.flipud(Z) # <------------ important!!!


# Reduce data according to Region params and decim scale
idx_lon = np.where((lon>=lonLims[0]) & (lon<=lonLims[1]))
idx_lat = np.where((lat>=latLims[0]) & (lat<=latLims[1]))

Z = Z[idx_lat[0][0]:idx_lat[0][-1]+1, idx_lon[0][0]:idx_lon[0][-1]+1]
lon = lon[idx_lon[0]]
lat = lat[idx_lat[0]]

lon = lon[::decim_scale]
lat = lat[::decim_scale]
Z = Z[::decim_scale, ::decim_scale]

## ---- plot ---- ##
## ----- plot with inset ---- ##
fig = plt.figure()
ax = fig.add_subplot(111)
m = Basemap(projection='merc',lon_0=lon_0,lat_0=lat_0, llcrnrlon=lonLims[0],llcrnrlat=latLims[0],urcrnrlon=lonLims[1],urcrnrlat=latLims[1], resolution='h')
x,y = m(*np.meshgrid(lon,lat))
c = m.contourf(x, y, -Z, v, cmap=plt.cm.PuBu_r, extend="min");
#c = m.contour(x, y, -Z, v2, colors='k');

m.fillcontinents(color='grey');
m.drawparallels(np.arange(latLims[0],latLims[1],1), labels=[1,0,0,0], fontsize=12, fontweight='bold');
m.drawmeridians(np.arange(lonLims[0], lonLims[1], 1), labels=[0,0,0,1], fontsize=12, fontweight='bold');

## Text on map
lat_27 = 47.5373
lon_27 = -52.5838

x, y = m(lon_27, lat_27)
m.scatter(x,y,5,marker='p',color='red')
plt.text(x, y, ' Stn 27', horizontalalignment='left', verticalalignment='center', fontsize=10, color='red', fontweight='bold')


#### ---- Save Figure ---- ####
fig.set_size_inches(w=5, h=6)
fig.set_dpi(300)
fig.savefig(fig_name)
#plt.show()

