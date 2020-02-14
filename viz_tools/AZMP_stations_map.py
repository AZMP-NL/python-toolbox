## from netCDF4 import Dataset, num2date
## import time, calendar, datetime, numpy
## from mpl_toolkits.basemap import Basemap
## import matplotlib.pyplot as plt
## import urllib, os
## # data downloaded from the form at
## # http://coastwatch.pfeg.noaa.gov/erddap/tabledap/apdrcArgoAll.html
## filename, headers = urllib.urlretrieve('http://coastwatch.pfeg.noaa.gov/erddap/tabledap/apdrcArgoAll.nc?longitude,latitude,time&longitude>=0&longitude<=360&latitude>=-90&latitude<=90&time>=2010-01-01&time<=2010-01-08&distinct()')
## dset = Dataset(filename)
## lats = dset.variables['latitude'][:]
## lons = dset.variables['longitude'][:]
## time = dset.variables['time']
## times = time[:]
## t1 = times.min(); t2 = times.max()
## date1 = num2date(t1, units=time.units)
## date2 = num2date(t2, units=time.units)
## dset.close()
## os.remove(filename)
## # draw map with markers for float locations
## m = Basemap(projection='hammer',lon_0=180)
## x, y = m(lons,lats)
## m.drawmapboundary(fill_color='#99ffff')
## m.fillcontinents(color='#cc9966',lake_color='#99ffff')
## m.scatter(x,y,3,marker='o',color='k')
## plt.title('Locations of %s ARGO floats active between %s and %s' %\
##         (len(lats),date1,date2),fontsize=12)
## plt.show()


import netCDF4
from mpl_toolkits.basemap import Basemap
import numpy as  np
import matplotlib.pyplot as plt
import openpyxl, pprint

## ---- Region parameters ---- ##
dataFile = '/home/cyrf0006/data_orwell/GEBCO_orwell/GRIDONE_1D.nc'
lon_0 = -50
lat_0 = 50
lonLims = [-70, -40]
latLims = [40, 65]
proj = 'merc'
decim_scale = 4
stationFile = '/home/cyrf0006/github/AZMP-NL/data/STANDARD_SECTIONS.xlsx'
fig_name = 'AZMP_sections.png'
ephem = 'ephem_calval.txt'
swot_kml = 'SWOT_Science_sept2015_Swath_10_60.kml'

# -> for AZMP database anaylis map:
lon_box = [-55, -50, -50, -55, -55]
lat_box = [45, 45, 50, 50, 45]


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

## ---- Station info ---- ##
import pandas as pd
df = pd.read_excel(stationFile)
#print the column names
print(df.columns)
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

## ---- Ephemerides ---- ##
eph = np.genfromtxt(ephem, dtype=str)
swot_lon = eph[:,1].astype(np.float)-180.0
swot_lat = eph[:,2].astype(np.float)

asign = np.sign(swot_lon)
signchange = ((np.roll(asign, 1) - asign) != 0).astype(int)
signchange = np.squeeze(np.array(np.where(signchange==1), int))

swot_segment_lat = []
swot_segment_lon = []
idx_end = 0
for g in signchange:
    idx_beg = idx_end+1
    idx_end = g
    swot_segment_lat.append(swot_lat[idx_beg:idx_end])
    swot_segment_lon.append(swot_lon[idx_beg:idx_end])
swot_segment_lat.append(swot_lat[idx_end:-1])
swot_segment_lon.append(swot_lon[idx_end:-1])
    
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
cb = plt.colorbar(c)
for l in cb.ax.yaxis.get_ticklabels():
    l.set_weight("bold")
    l.set_fontsize(10)
cb.set_label('Depth (m)', fontsize=13, fontweight='bold')
plt.title("AZMP-NL Standard Sections", fontsize=13, fontweight='bold')

# plot stations
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
m.scatter(x[0],y[0],25,marker='p',color='r')
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



#### ---- Save Figure ---- ####
fig.set_size_inches(w=8, h=9)
fig.set_dpi(200)
fig.savefig(fig_name)
#plt.show()

