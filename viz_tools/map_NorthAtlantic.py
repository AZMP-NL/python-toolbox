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
lon_0 = -50
lat_0 = 50
lonLims = [-65, 10]
latLims = [40, 65]
proj = 'merc'
decim_scale = 4
stationFile = '/home/cyrf0006/research/AZMP_surveys/STANDARD_SECTIONS.xlsx'
fig_name = 'map_NorthAtlantic.png'

## ---- Bathymetry ---- ####
v = np.linspace(-4000, 0, 9)

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

# Reduce data according to Region params
lon = lon[::decim_scale]
lat = lat[::decim_scale]
Z = Z[::decim_scale, ::decim_scale]

## ---- plot ---- ##
fig = plt.figure(1)
#m = Basemap(projection='ortho',lon_0=lon_0,lat_0=lat_0,resolution=None)
m = Basemap(projection='merc',lon_0=lon_0,lat_0=lat_0, llcrnrlon=lonLims[0],llcrnrlat=latLims[0],urcrnrlon=lonLims[1],urcrnrlat=latLims[1], resolution='l')
x,y = m(*np.meshgrid(lon,lat))
c = m.contourf(x, y, np.flipud(Z), v, cmap=plt.cm.PuBu_r, extend="min");
#c = m.contourf(x, y, np.flipud(Z), v, cmap=plt.cm.PuRd_r, extend="min");
m.fillcontinents(color='grey');
m.drawparallels(np.arange(10,70,10), labels=[1,0,0,0], fontsize=12, fontweight='bold');
m.drawmeridians(np.arange(-80, 5, 10), labels=[0,0,0,1], fontsize=12, fontweight='bold');
cb = plt.colorbar(c, orientation='horizontal')
for l in cb.ax.yaxis.get_ticklabels():
    l.set_weight("bold")
    l.set_fontsize(10)
cb.set_label('Depth (m)', fontsize=13, fontweight='bold')



## Text on map
lat_Brest = 48.4 
lon_Brest = -4.5
lat_SJ = 47.6 
lon_SJ = -52.7
x, y = m(lon_SJ, lat_SJ)
m.scatter(x,y,5,marker='p',color='red')
plt.text(x, y, ' St. John\'s', horizontalalignment='left', verticalalignment='center', fontsize=10, color='red', fontweight='bold')
x, y = m(lon_Brest, lat_Brest)
m.scatter(x,y,5,marker='p',color='red')
plt.text(x, y, 'Brest ', horizontalalignment='right', verticalalignment='center', fontsize=10, color='red', fontweight='bold')




#### ---- Save Figure ---- ####
fig.set_size_inches(w=9, h=6)
fig.set_dpi(200)
fig.savefig(fig_name)
#plt.show()

