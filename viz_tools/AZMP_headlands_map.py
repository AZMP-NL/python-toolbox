'''
Script created to plot a map with SST boxes on it.

To extract the box lat/lon, I did this:
df = pd.read_csv('/home/cyrf0006/github/AZMP-NL/utils/SST_boxes_coordinates.csv', header=None)
# from ftp://ftp.dfo-mpo.gc.ca/bometrics/noaa/noaa-dat_stats.in.all (removed 1st line)
box_name = df.iloc[::2]
box_coord = df.iloc[1::2]
box_name = box_name.reset_index(drop=True)
box_coord = box_coord.reset_index(drop=True)
df = pd.concat([box_name, box_coord], axis=1)
df.to_csv('SST_boxes_coords.csv')
box_name.to_csv('SST_boxes_coords_name.csv')  # <---- Easier to play with
box_coord.to_csv('SST_boxes_coords_coords.csv')

Then I opened it and copy-pasted it in /home/cyrf0006/github/AZMP-NL/utils/SST_boxes.xslx

'''


import netCDF4
from mpl_toolkits.basemap import Basemap
import numpy as  np
import matplotlib.pyplot as plt
import openpyxl, pprint
import pandas as pd
import os
import cmocean

## ---- Region parameters ---- ##
dataFile = '/home/cyrf0006/data/GEBCO/GRIDONE_1D.nc'
lon_0 = -50
lat_0 = 50
#lonLims = [-60, -52]
#latLims = [46, 52]
lonLims = [-61.5, -50.25]
latLims = [45, 52.25]
proj = 'merc'
decim_scale = 4
stationFile = '/home/cyrf0006/github/AZMP-NL/utils/Headlands_sites.xlsx'
fig_name = 'map_headlands_sites.png'

## ---- Load SST boxes ---- ##
df_box = pd.read_excel('/home/cyrf0006/github/AZMP-NL/utils/SST_boxes.xslx')
#df_box = df_box[df_box.region=='NL']
df_box = df_box[(df_box.region=='NL') | (df_box.region=='GSL')]


## ---- Bathymetry ---- ####
v = np.linspace(0, 500, 11)

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
## idx_lon = np.where((lon>=lonLims[0]) & (lon<=lonLims[1]))
## idx_lat = np.where((lat>=latLims[0]) & (lat<=latLims[1]))
## Z = Z[idx_lat[0][0]:idx_lat[0][-1]+1, idx_lon[0][0]:idx_lon[0][-1]+1]
## lon = lon[idx_lon[0]]
## lat = lat[idx_lat[0]]

lon = lon[::decim_scale]
lat = lat[::decim_scale]
Z = Z[::decim_scale, ::decim_scale]

## ---- Station info ---- ##
import pandas as pd
df = pd.read_excel(stationFile)
#print the column names
print df.columns
#get the values for a given column
station = df['Site'].values
stationLat = df['Lat'].values
stationLon = df['Lon'].values


## ---- plot ---- ##
fig = plt.figure(1)
#m = Basemap(projection='ortho',lon_0=lon_0,lat_0=lat_0,resolution=None)
m = Basemap(projection='merc',lon_0=lon_0,lat_0=lat_0, llcrnrlon=lonLims[0],llcrnrlat=latLims[0],urcrnrlon=lonLims[1],urcrnrlat=latLims[1], resolution='h')
x,y = m(*np.meshgrid(lon,lat))
c = m.contourf(x, y, np.flipud(-Z), v, cmap=cmocean.cm.deep, extend="both");
cc = m.contour(x, y, np.flipud(-Z), [25, 50, 100, 200, 400], colors='lightgrey', linewidths=.5);
plt.clabel(cc, inline=1, fontsize=10, colors='gray', fmt='%d')
#ccc = m.contour(x, y, np.flipud(-Z), [0], colors='black');
m.fillcontinents(color='peru');
m.drawparallels(np.arange(10,70,2), labels=[1,0,0,0], fontsize=12, fontweight='bold');
m.drawmeridians(np.arange(-80, 5, 2), labels=[0,0,0,1], fontsize=12, fontweight='bold');

# plot Headlands
sites = df.Site.values
lats = df.Lat.values
lons = df.Lon.values

for idx, name in enumerate(sites):
    x, y = m(lons[idx], lats[idx])
    m.plot(x,y,'p',color='k')
    #plt.text(x, y, np.str(' '+name), horizontalalignment='left', verticalalignment='center', fontsize=15, color='k', fontweight='bold')
    plt.text(x, y, np.str(name), horizontalalignment='center', verticalalignment='top', fontsize=15, color='k', fontweight='bold')

# plot SST_boxes
abbr = df_box.abbr.values
lat_min = df_box.lat_min.values
lat_max = df_box.lat_max.values
lon_min = df_box.lon_min.values
lon_max = df_box.lon_max.values
for idx, name in enumerate(abbr):
    xbox = np.array([lon_min[idx], lon_max[idx], lon_max[idx], lon_min[idx], lon_min[idx]])
    ybox = np.array([lat_min[idx], lat_min[idx], lat_max[idx], lat_max[idx], lat_min[idx]])
    x, y = m(xbox, ybox)
    if (name=='NEGSL') | (name=='SAB') | (name=='NENL') | (name=='AC') | (name=='SPB'):
        m.plot(x,y,color='k')
        plt.text(x.mean(), y.mean(), name, horizontalalignment='center', verticalalignment='center', fontsize=10, color='k', fontweight='bold')
    elif (name=='CS'):
        m.plot(x,y,color='w')
        plt.text(x.mean(), y.mean(), name, horizontalalignment='center', verticalalignment='center', fontsize=10, color='w', fontweight='bold')


#### ---- Save Figure ---- ####
fig.set_size_inches(w=10, h=12)
fig.savefig(fig_name, dpi=200)
os.system('convert -trim ' + fig_name + ' ' + fig_name)
#plt.show()

