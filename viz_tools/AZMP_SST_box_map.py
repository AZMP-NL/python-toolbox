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
lonLims = [-60, -50]
latLims = [45, 55]
lonLims = [-70, -40]
latLims = [40, 65]
proj = 'merc'
decim_scale = 4
stationFile = '/home/cyrf0006/github/AZMP-NL/data/STANDARD_SECTIONS.xlsx'
fig_name = 'AZMP_SST_boxes.png'

## ---- Load SST boxes ---- ##
df_box = pd.read_excel('/home/cyrf0006/github/AZMP-NL/utils/SST_boxes.xslx')
#df_box = df_box[df_box.region=='NL']
#df_box = df_box[(df_box.region=='NL') | (df_box.region=='GSL')]
df_box = df_box[(df_box.region=='NL') | (df_box.region=='GSL') | (df_box.region=='Nlx')]

## ---- Bathymetry ---- ####
v = np.linspace(0, 3500, 36)
#v = np.linspace(-4000, 0, 9)

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

## ---- plot ---- ##
fig = plt.figure(1)
#m = Basemap(projection='ortho',lon_0=lon_0,lat_0=lat_0,resolution=None)
m = Basemap(projection='merc',lon_0=lon_0,lat_0=lat_0, llcrnrlon=lonLims[0],llcrnrlat=latLims[0],urcrnrlon=lonLims[1],urcrnrlat=latLims[1], resolution='l')
x,y = m(*np.meshgrid(lon,lat))
c = m.contourf(x, y, np.flipud(-Z), v, cmap=cmocean.cm.deep, extend="max");
cc = m.contour(x, y, np.flipud(-Z), [100, 500, 1000, 3000, 4000], colors='lightgrey', linewidths=.5);
plt.clabel(cc, inline=1, fontsize=10, colors='gray', fmt='%d')
ccc = m.contour(x, y, np.flipud(-Z), [0], colors='black');
m.fillcontinents(color='peru');
m.drawparallels(np.arange(10,70,10), labels=[1,0,0,0], fontsize=12, fontweight='bold');
m.drawmeridians(np.arange(-80, 5, 10), labels=[0,0,0,1], fontsize=12, fontweight='bold');
#plt.title('Standard air temperature sites')

## x,y = m(*np.meshgrid(lon,lat))
## c = m.contourf(x, y, np.flipud(Z), v, cmap=plt.cm.PuBu_r, extend="min");
## #c = m.contourf(x, y, np.flipud(Z), v, cmap=plt.cm.PuRd_r, extend="min");
## m.fillcontinents(color='grey');
## m.drawparallels(np.arange(10,70,10), labels=[1,0,0,0], fontsize=12, fontweight='bold');
## m.drawmeridians(np.arange(-80, 5, 10), labels=[0,0,0,1], fontsize=12, fontweight='bold');
## cb = plt.colorbar(c)
## for l in cb.ax.yaxis.get_ticklabels():
##     l.set_weight("bold")
##     l.set_fontsize(10)
## cb.set_label('Depth (m)', fontsize=13, fontweight='bold')
## #plt.title("AZMP-NL standard lines", fontsize=13, fontweight='bold')

# plot AZMP stations
x, y = m(stationLon[index_SEGB],stationLat[index_SEGB])
m.scatter(x,y,3,marker='o',color='r')
plt.text(x[-1], y[-1], ' SEGB', horizontalalignment='left', verticalalignment='center', fontsize=10, color='r', fontweight='bold')
x, y = m(stationLon[index_FC],stationLat[index_FC])
m.scatter(x,y,3,marker='o',color='r')
plt.text(x[-1], y[-1], ' FC', horizontalalignment='left', verticalalignment='center', fontsize=10, color='r', fontweight='bold')
x, y = m(stationLon[index_BB],stationLat[index_BB])
m.scatter(x,y,3,marker='o',color='r')
plt.text(x[-1], y[-1], ' BB', horizontalalignment='left', verticalalignment='center', fontsize=10, color='r', fontweight='bold')
x, y = m(stationLon[index_WB],stationLat[index_WB])
m.scatter(x,y,3,marker='o',color='r')
plt.text(x[-1], y[-1], ' WB', horizontalalignment='left', verticalalignment='center', fontsize=10, color='r', fontweight='bold')
x, y = m(stationLon[index_SI],stationLat[index_SI])
m.scatter(x,y,3,marker='o',color='r')
plt.text(x[-1], y[-1], ' SI', horizontalalignment='left', verticalalignment='center', fontsize=10, color='r', fontweight='bold')
x, y = m(stationLon[index_MB],stationLat[index_MB])
m.scatter(x,y,3,marker='o',color='r')
plt.text(x[-1], y[-1], ' MB', horizontalalignment='left', verticalalignment='center', fontsize=10, color='r', fontweight='bold')
x, y = m(stationLon[index_BI],stationLat[index_BI])
m.scatter(x,y,3,marker='o',color='r')
plt.text(x[-1], y[-1], ' BI', horizontalalignment='left', verticalalignment='center', fontsize=10, color='r', fontweight='bold')
x, y = m(stationLon[index_FI],stationLat[index_FI])
m.scatter(x,y,3,marker='o',color='lightcoral')
plt.text(x[-1], y[-1], ' FI', horizontalalignment='left', verticalalignment='center', fontsize=10, color='lightcoral', fontweight='bold')
x, y = m(stationLon[index_S27],stationLat[index_S27])
m.scatter(x,y,3,marker='o',color='lightcoral')
#m.scatter(x[0],y[0],25,marker='p',color='r')
plt.text(x[-1], y[-1], ' S27', horizontalalignment='left', verticalalignment='center', fontsize=10, color='lightcoral', fontweight='bold')
x, y = m(stationLon[index_SESPB],stationLat[index_SESPB])
m.scatter(x,y,3,marker='o',color='r')
plt.text(x[-1], y[-1], 'SESPB ', horizontalalignment='right', verticalalignment='center', fontsize=10, color='r', fontweight='bold')
x, y = m(stationLon[index_SWSPB],stationLat[index_SWSPB])
m.scatter(x,y,3,marker='o',color='r')
plt.text(x[-1], y[-1], 'SWSPB ', horizontalalignment='right', verticalalignment='center', fontsize=10, color='r', fontweight='bold')
x, y = m(stationLon[index_SS],stationLat[index_SS])
m.scatter(x,y,3,marker='o',color='lightcoral')
plt.text(x[0], y[0], 'SS ', horizontalalignment='right', verticalalignment='center', fontsize=10, color='lightcoral', fontweight='bold')

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
    if name=='CLS':
        m.plot(x,y,color='w')
        x, y = m(xbox, ybox+1.5)
        plt.text(x.mean(), y.mean(), name, horizontalalignment='center', verticalalignment='center', fontsize=10, color='w', fontweight='bold')
    elif (name=='BRA') | (name=='NCLS') | (name=='OK'):
        m.plot(x,y,color='w')
        plt.text(x.mean(), y.mean(), name, horizontalalignment='center', verticalalignment='center', fontsize=10, color='w', fontweight='bold')
    else:
        m.plot(x,y,color='k')
        plt.text(x.mean(), y.mean(), name, horizontalalignment='center', verticalalignment='center', fontsize=10, color='k', fontweight='bold')




#### ---- Save Figure ---- ####
fig.set_size_inches(w=10, h=12)
fig.savefig(fig_name, dpi=200)
os.system('convert -trim ' + fig_name + ' ' + fig_name)
#plt.show()

