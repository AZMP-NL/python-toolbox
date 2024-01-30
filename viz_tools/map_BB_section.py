'''
Map of Bonavista section used to plan glider deployment.

'''


import netCDF4
from mpl_toolkits.basemap import Basemap
import numpy as  np
import matplotlib.pyplot as plt
import openpyxl, pprint

## ---- Region parameters ---- ##
dataFile = '/home/cyrf0006/data/GEBCO/GRIDONE_1D.nc'
lon_0 = -50
lat_0 = 50
lonLims = [-55, -45]
latLims = [46, 52]
proj = 'merc'
decim_scale = 4
stationFile = '/home/cyrf0006/github/AZMP-NL/data/STANDARD_SECTIONS.xlsx'
fig_name = 'map_BB_section.png'

## ---- Bathymetry ---- ####
v = np.linspace(-2500, 0, 21)

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

## ---- Station info ---- ##
import pandas as pd
df = pd.read_excel(stationFile)
#print the column names
print df.columns
#get the values for a given column
sections = df['SECTION'].values
stations = df['STATION'].values
stationLat = df['LAT'].values
stationLon = df['LONG.1'].values

index_BB = df.SECTION[df.SECTION=="BONAVISTA"].index.tolist()

## ---- plot ---- ##
fig = plt.figure(1)
plt.clf()
#m = Basemap(projection='ortho',lon_0=lon_0,lat_0=lat_0,resolution=None)
m = Basemap(projection='merc',lon_0=lon_0,lat_0=lat_0, llcrnrlon=lonLims[0],llcrnrlat=latLims[0],urcrnrlon=lonLims[1],urcrnrlat=latLims[1], resolution='f')
x,y = m(*np.meshgrid(lon,lat))
c = m.contourf(x, y, np.flipud(Z), v, cmap=plt.cm.PuBu_r, extend="min");
#c = m.contourf(x, y, np.flipud(Z), v, cmap=plt.cm.PuRd_r, extend="min");
m.fillcontinents(color='grey');
m.drawparallels(np.arange(10,70,10), labels=[1,0,0,0], fontsize=12, fontweight='bold');
m.drawmeridians(np.arange(-80, 5, 10), labels=[0,0,0,1], fontsize=12, fontweight='bold');
cb = plt.colorbar(c)
for l in cb.ax.yaxis.get_ticklabels():
    l.set_weight("bold")
    l.set_fontsize(10)
cb.set_label('Depth (m)', fontsize=13, fontweight='bold')
plt.title("AZMP-NL BB Section", fontsize=13, fontweight='bold')

# plot stations
x, y = m(stationLon[index_BB],stationLat[index_BB])
m.scatter(x,y,3,marker='o',color='r')
plt.text(x[-1], y[-1], ' BB', horizontalalignment='left', verticalalignment='center', fontsize=10, color='r', fontweight='bold')


#### ---- Save Figure ---- ####
fig.set_size_inches(w=8, h=6)
fig.set_dpi(200)
fig.savefig(fig_name)
#plt.show()

