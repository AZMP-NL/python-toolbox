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
lonLims = [-65, -40]
latLims = [40, 60]
lon_0 = np.mean(lonLims)
lat_0 = np.mean(latLims)
proj = 'merc'
decim_scale = 1
fig_name = 'map_stn27_inset.png'

## ---- Bathymetry ---- ####
v = np.linspace(0, 5000, 21)
v1 = np.linspace(0, 4000, 5)
v2 = np.linspace(0, 500, 21)

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

map = Basemap(projection='cyl',lon_0=lon_0,lat_0=lat_0, llcrnrlon=lonLims[0],llcrnrlat=latLims[0],urcrnrlon=lonLims[1],urcrnrlat=latLims[1], resolution='h')
x,y = map(*np.meshgrid(lon,lat))
c = map.contourf(x, y, -Z, v, cmap=plt.cm.PuBu_r);
c2 = map.contour(x, y, -Z, v1, colors='k');
map.fillcontinents(color='grey');
map.drawparallels(np.arange(10,70,10), labels=[1,0,0,0], fontsize=10, fontweight='normal');
map.drawmeridians(np.arange(-80, 5, 10), labels=[0,0,0,1], fontsize=10, fontweight='normal');
cb = plt.colorbar(c, orientation='horizontal')
cb.set_label('Depth (m)', fontsize=11, fontweight='normal')

map.drawmapboundary()
#map.drawcoastlines()


## ---- Inset ---- ##
axins = zoomed_inset_axes(ax, 5, loc=1)
axins.set_xlim(-53.75, -51.75)
axins.set_ylim(46, 48)
#axins.set_edgecolor('r')
plt.xticks(visible=False)
plt.yticks(visible=False)
map2 = Basemap(llcrnrlon=-53.75, llcrnrlat=46, urcrnrlon=-51.75,urcrnrlat=48, ax=axins, resolution='f')
x,y = map2(*np.meshgrid(lon,lat))
#c = map2.contour(x, y, -Z, v2, colors='k');
c = map2.contourf(x, y, -Z, v2, cmap=plt.cm.PuBu_r);
map2.fillcontinents(color='grey');
map2.drawparallels(np.arange(-46,48,1))
map2.drawmeridians(np.arange(-54, -52, 1))

## Stn27 on map
lat_27 = 47.5373
lon_27 = -52.5838
x, y = map2(lon_27, lat_27)
map2.scatter(x,y,5,marker='p',color='red')
plt.text(x, y, ' Stn 27', horizontalalignment='left', verticalalignment='center', fontsize=10, color='red', fontweight='bold')

# Plot inset lines
mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="r")


#### ---- Save Figure ---- ####
fig.set_size_inches(w=5, h=5)
fig.set_dpi(300)
fig.savefig(fig_name)
#plt.show()

