

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
lonLims = [-70, 0]
latLims = [30, 80]
proj = 'merc'
decim_scale = 4
stationFile = '/home/cyrf0006/research/AZMP_surveys/STANDARD_SECTIONS.xlsx'
fig_name = 'map_NorthAtlantic.png'

## ---- Bathymetry ---- ####
v = np.linspace(-6000, 0, 11)

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

index_SEGB = df.SECTION[df.SECTION=="SOUTHEAST GRAND BANK"].index.tolist()
index_FC = df.SECTION[df.SECTION=="FLEMISH CAP"].index.tolist()
index_BB = df.SECTION[df.SECTION=="BONAVISTA"].index.tolist()
index_WB = df.SECTION[df.SECTION=="WHITE BAY"].index.tolist()
index_SI = df.SECTION[df.SECTION=="SEAL ISLAND"].index.tolist()
index_MB = df.SECTION[df.SECTION=="MAKKOVIK BANK"].index.tolist()
index_BI = df.SECTION[df.SECTION=="BEACH ISLAND"].index.tolist()
index_FI = df.SECTION[df.SECTION=="FUNK ISLAND"].index.tolist()
index_S27 = df.SECTION[df.SECTION=="STATION 27"].index.tolist()
index_SESPB = df.SECTION[df.SECTION=="SOUTHEAST ST PIERRE BANK"].index.tolist()
index_SWSPB = df.SECTION[df.SECTION=="SOUTHWEST ST PIERRE BANK"].index.tolist()
index_SS = df.SECTION[df.SECTION=="SMITH SOUND"].index.tolist()



## ----- plot with inset ---- ##
fig = plt.figure()
ax = fig.add_subplot(111)

map = Basemap(projection='cyl',lon_0=lon_0,lat_0=lat_0, llcrnrlon=lonLims[0],llcrnrlat=latLims[0],urcrnrlon=lonLims[1],urcrnrlat=latLims[1], resolution='l')
x,y = map(*np.meshgrid(lon,lat))
c = map.contourf(x, y, np.flipud(Z), v, cmap=plt.cm.PuBu_r);
map.fillcontinents(color='grey');
map.drawparallels(np.arange(10,70,10), labels=[1,0,0,0], fontsize=10, fontweight='normal');
map.drawmeridians(np.arange(-80, 5, 10), labels=[0,0,0,1], fontsize=10, fontweight='normal');
#cax = plt.axes([0.9,0.1,0.1,0.8])
#plt.colorbar(c, cax=cax)
cb = plt.colorbar(c, orientation='horizontal')
## for l in cb.ax.yaxis.get_ticklabels():
##     l.set_weight("bold")
##     l.set_fontsize(10)
cb.set_label('Depth (m)', fontsize=11, fontweight='normal')

map.drawmapboundary()
#map.drawcoastlines()

## Text on map
lat_Brest = 48.4 
lon_Brest = -4.5
lat_SJ = 47.6 
lon_SJ = -52.7
x, y = map(lon_SJ, lat_SJ)
map.scatter(x,y,10,marker='o',color='red')
plt.text(x, y, ' St. John\'s', horizontalalignment='left', verticalalignment='center', fontsize=10, color='red', fontweight='bold')
x, y = map(lon_Brest, lat_Brest)
map.scatter(x,y,10,marker='o',color='red')
plt.text(x, y, 'Brest ', horizontalalignment='right', verticalalignment='center', fontsize=10, color='red', fontweight='bold')


## ## Fancy arrows (To highlight GS and LC)
## x_text, y_text = map(-90,-9)
## x_gs, y_gs = map(-50, 41)
## plt.annotate(' ',
##             xy=(x_gs, y_gs), xycoords='data',
##             xytext=(x_text, y_text), textcoords='offset points',
##             size=20,
##             # bbox=dict(boxstyle="round", fc="0.8"),
##             arrowprops=dict(arrowstyle="fancy",
##                             #fc="0.6", ec="none",
##                             facecolor='orange',
##                             connectionstyle="angle3,angleA=25,angleB=-55"))

## x_text, y_text = map(-18,60)
## x_lc, y_lc = map(-48, 50)
## plt.annotate(' ',
##             xy=(x_lc, y_lc), xycoords='data',
##             xytext=(x_text, y_text), textcoords='offset points',
##             size=20,
##             # bbox=dict(boxstyle="round", fc="0.8"),
##             arrowprops=dict(arrowstyle="fancy",
##                             #fc="0.6", ec="none",
##                             facecolor='cornflowerblue',
##                             connectionstyle="angle3,angleA=180,angleB=-750"))

## x, y = map(-60, 38)
## plt.text(x, y, 'GS', horizontalalignment='right', verticalalignment='center', fontsize=10, color='orange', fontweight='bold')
## x, y = map(-53, 57)
## plt.text(x, y, 'LC', horizontalalignment='right', verticalalignment='center', fontsize=10, color='black', fontweight='bold')

## Inset
axins = zoomed_inset_axes(ax, 1, loc=9)
axins.set_xlim(-63, -40)
axins.set_ylim(40, 60)
#axins.set_edgecolor('r')
plt.xticks(visible=False)
plt.yticks(visible=False)

#map2 = Basemap(projection='cyl', llcrnrlon=-70, llcrnrlat=40, urcrnrlon=-40,urcrnrlat=65, ax=axins)
map2 = Basemap(llcrnrlon=-63, llcrnrlat=40, urcrnrlon=-40,urcrnrlat=60, ax=axins, resolution='h')
x,y = map2(*np.meshgrid(lon,lat))
c = map2.contourf(x, y, np.flipud(Z), v, cmap=plt.cm.PuBu_r, extend="min");
#c = map.contourf(x, y, np.flipud(Z), v, cmap=plt.cm.PuRd_r, extend="min");
map2.fillcontinents(color='grey');
map2.drawparallels(np.arange(10,70,10))
map2.drawmeridians(np.arange(-80, 5, 10))


## Stations
x, y = map2(stationLon[index_SEGB],stationLat[index_SEGB])
map2.scatter(x,y,3,marker='o',color='r')
x, y = map2(stationLon[index_FC],stationLat[index_FC])
map2.scatter(x,y,3,marker='o',color='r')
x, y = map2(stationLon[index_BB],stationLat[index_BB])
map2.scatter(x,y,3,marker='o',color='r')
x, y = map2(stationLon[index_WB],stationLat[index_WB])
map2.scatter(x,y,3,marker='o',color='r')
x, y = map2(stationLon[index_SI],stationLat[index_SI])
map2.scatter(x,y,3,marker='o',color='r')
x, y = map2(stationLon[index_MB],stationLat[index_MB])
map2.scatter(x,y,3,marker='o',color='r')
x, y = map2(stationLon[index_BI],stationLat[index_BI])
map2.scatter(x,y,3,marker='o',color='r')
x, y = map2(stationLon[index_FI],stationLat[index_FI])
map2.scatter(x,y,3,marker='o',color='r')
x, y = map2(stationLon[index_SWSPB],stationLat[index_SWSPB])
map2.scatter(x,y,3,marker='o',color='r')


mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="r")


#### ---- Save Figure ---- ####
fig.set_size_inches(w=6, h=6)
fig.set_dpi(200)
fig.savefig(fig_name)
#plt.show()

