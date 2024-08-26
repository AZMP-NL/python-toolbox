'''
This script generates the map used in CSAS ResDoc.
Inspired after ~/github/python/nafc/map_NL_index.py
used for the NL climate index.

Run in /home/cyrf0006/AZMP/utils

Frederic.Cyr@dfo-mpo.gc.ca
March 2021
'''


import os
import h5py                                                                
import netCDF4
os.environ['PROJ_LIB'] = '/home/cyrf0006/anaconda3/share/proj'
from mpl_toolkits.basemap import Basemap
import numpy as  np
import matplotlib.pyplot as plt
import openpyxl, pprint
import pandas as pd
import cmocean
import azmp_utils as azu
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

## ---- Region parameters ---- ##
dataFile = '/home/cyrf0006/data/GEBCO/GEBCO_2014_1D.nc'
lon_0 = -50
lat_0 = 50
#lonLims = [-70, -35]
#latLims = [38, 68]
lonLims = [-77, -35]
latLims = [35, 68]
proj = 'merc'
decim_scale = 4
stationFile = '/home/cyrf0006/github/AZMP-NL/utils/STANDARD_SECTIONS.xlsx'
fig_name = 'map_WGOH.png'
#fig_name = 'map_WGOH_allairT.png'
AZMP_airTemp_file = '/home/cyrf0006/github/AZMP-NL/utils/airTemp_sites.xlsx'
lightdeep = cmocean.tools.lighten(cmocean.cm.deep, .9)
v = np.linspace(0, 3500, 36)
#v = np.linspace(-4000, 0, 9)
v = np.linspace(0, 5400, 55) # Olivia's (good one, but loose channels on the shelf)
#v = np.linspace(0, 4000, 41) # t

## ---- Ice region ---- ##
ice_region = azu.get_ice_regions()

# Some polygons for regions
# NEC:
xlon = [-66.80, -66.00, -66.00, -66.80, -66.80]
xlat = [42.20, 42.20, 42.60, 42.60, 42.20]
NEC = {'lat' : xlat, 'lon' : xlon}
# EGOM:
xlon = [-68.00, -67.50, -67.00, -66.00, -66.00, -66.50, -67.00, -68.00, -68.00]
xlat = [41.80, 42.10, 42.10, 42.10, 43.00, 44.20, 44.20, 44.20, 41.80]
EGOM = {'lat' : xlat, 'lon' : xlon}
US_poly = {}
US_poly['NEC'] = NEC
US_poly['EGOM'] = EGOM
# CLS:
xlon = [-53.1770, -53.8681, -54.6382, -54.7567, -53.2856, -52.4069, -51.1432, -49.6523, -48.4676, -48.3590, -48.4500, -49.2000, -50.4126, -52.0614, -52.6044, -53.1770]
xlat = [55.9918, 56.9547, 57.7997, 58.9001, 59.9219, 60.1184, 59.8237, 59.2538, 58.5071, 57.6032, 57.1217, 56.0999, 55.4318, 55.5791,  55.9034, 55.9918]
CLS = {'lat' : xlat, 'lon' : xlon}
CLS_poly = {}
CLS_poly['CLS'] = CLS

sNEC = [-66, 42.2]
sGB = [-68, 41.5]
sEGOM = [-67.5, 43.5]
sWGOM = [-69.5, 42.5]
sNMAB = [-71.5, 40.5]
sSMAB = [-74.5, 38.5]


## ---- Bathymetry ---- ##
print('Load and grid bathymetry')
# h5 file
h5_outputfile = 'wgoh_east_northam_shelf_bathymetry.h5'
if os.path.isfile(h5_outputfile):
     print([h5_outputfile + ' exists! Reading directly'])
     h5f = h5py.File(h5_outputfile,'r')
     lon = h5f['lon'][:]
     lat = h5f['lat'][:]
     Z = h5f['Z'][:]
     h5f.close()

else:
    # Extract variables
    dataset = netCDF4.Dataset(dataFile)
    x = [-179-59.75/60, 179+59.75/60] # to correct bug in 30'' dataset?
    y = [-89-59.75/60, 89+59.75/60]
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
    # Reduce data according to Region params
    idx_lon = np.where((lon>=lonLims[0]) & (lon<=lonLims[1]))
    idx_lat = np.where((lat>=latLims[0]) & (lat<=latLims[1]))
    Z = Z[idx_lat[0][0]:idx_lat[0][-1]+1, idx_lon[0][0]:idx_lon[0][-1]+1]
    lon = lon[idx_lon[0]]
    lat = lat[idx_lat[0]]

    # Save data for later use
    h5f = h5py.File(h5_outputfile, 'w')
    h5f.create_dataset('lon', data=lon)
    h5f.create_dataset('lat', data=lat)
    h5f.create_dataset('Z', data=Z)
    h5f.close()
    print(' -> Done!')


## ---- AZMP air temp Stations  ---- ##
df = pd.read_excel(AZMP_airTemp_file)
air_stationLat = df['latitude'].values
air_stationLon = df['longitude'].values

## ---- NAFO divisions ---- ##
#nafo_div = azu.get_nafo_divisions()

## ---- Ice region ---- ##
#ice_region = azu.get_ice_regions()

## ---- Load SST boxes ---- ##
#df_box = pd.read_excel('/home/cyrf0006/github/AZMP-NL/utils/SST_boxes.xslx')
#df_box = df_box[(df_box.region=='NL') | (df_box.region=='GSL') | (df_box.region=='NS')]

## ---- Station info ---- ##
import pandas as pd
df = pd.read_excel(stationFile)
sections = df['SECTION'].values
stations = df['STATION'].values
stationLat = df['LAT'].values
stationLon = df['LONG.1'].values

index_BI = df.SECTION[df.SECTION=="BEACH ISLAND"].index.tolist()
index_MB = df.SECTION[df.SECTION=="MAKKOVIK BANK"].index.tolist()
index_SI = df.SECTION[df.SECTION=="SEAL ISLAND"].index.tolist()
index_WB = df.SECTION[df.SECTION=="WHITE BAY"].index.tolist()
index_BB = df.SECTION[df.SECTION=="BONAVISTA"].index.tolist()
index_FC = df.SECTION[df.SECTION=="FLEMISH CAP"].index.tolist()
index_SEGB = df.SECTION[df.SECTION=="SOUTHEAST GRAND BANK"].index.tolist()
index_SESPB = df.SECTION[df.SECTION=="SOUTHEAST ST PIERRE BANK"].index.tolist()
index_SWSPB = df.SECTION[df.SECTION=="SOUTHWEST ST PIERRE BANK"].index.tolist()

index_S27 = df.SECTION[df.SECTION=="STATION 27"].index.tolist()

## ---- plot ---- ##
fig, ax = plt.subplots(nrows=1, ncols=1)
m = Basemap(projection='merc',lon_0=lon_0,lat_0=lat_0, llcrnrlon=lonLims[0],llcrnrlat=latLims[0],urcrnrlon=lonLims[1],urcrnrlat=latLims[1], resolution='l')
x,y = m(*np.meshgrid(lon,lat))
c = plt.contourf(x, y, -Z, levels=[10,100,150,200,300,500,1000,2000,3000,4000], cmap=lightdeep, extend='max', zorder=1)
cc = plt.contour(x, y, -Z, [100, 500, 1000, 2000, 3000, 4000], colors='lightgrey', linewidths=.5, zorder=1)
ccc = plt.contour(x, y, -Z, [1000], colors='darkgrey', linewidths=2, zorder=1)
plt.clabel(cc, inline=1, fontsize=10, colors='lightgrey', fmt='%d')

m.fillcontinents(color='tan', zorder=50);
m.drawparallels(np.arange(10,70,5), labels=[1,0,0,0], fontsize=12, fontweight='bold');
m.drawmeridians(np.arange(-80, 5, 5), labels=[0,0,0,1], fontsize=12, fontweight='bold');


# plot AZMP stations
x, y = m(stationLon[index_SI],stationLat[index_SI])
#m.scatter(x,y,10,marker='o',color='r', zorder=10)
m.plot(x,y,color='r', zorder=10)
plt.text(x[-1], y[-1], ' Seal Island', horizontalalignment='left', verticalalignment='center', fontsize=13, color='r', fontweight='bold', zorder=100)

x, y = m(stationLon[index_BB],stationLat[index_BB])
#m.scatter(x,y,10,marker='o',color='r', zorder=10)
m.plot(x,y,color='r', zorder=10)
plt.text(x[-1], y[-1], ' Bonavista', horizontalalignment='left', verticalalignment='center', fontsize=13, color='r', fontweight='bold', zorder=100)

x, y = m(stationLon[index_FC],stationLat[index_FC])
#m.scatter(x,y,10,marker='o',color='r', zorder=100)
m.plot(x,y,color='r', zorder=100)
plt.text(x[-6], y[-1]-50000, ' Flemish Cap', horizontalalignment='center', verticalalignment='top', fontsize=13, color='r', fontweight='bold', zorder=100)


# Add high-res stations
x, y = m(stationLon[index_S27],stationLat[index_S27])
m.scatter(x[0],y[0],100,marker='*',color='r', zorder=100)
plt.text(x[0], y[0]+50000, 'S27', horizontalalignment='left', verticalalignment='bottom', fontsize=13, color='r', fontweight='bold', zorder=100)
x, y = m(-63.317, 44.267)
m.scatter(x,y,100,marker='*',color='m', zorder=100)
plt.text(x, y, '  HLX-2', horizontalalignment='right', verticalalignment='bottom', fontsize=10, color='m', fontweight='bold', zorder=100)
x, y = m(-66.85, 44.93)
m.scatter(x,y,100,marker='*',color='m', zorder=100)
plt.text(x, y, ' Prince 5  ', horizontalalignment='right', verticalalignment='bottom', fontsize=10, color='m', fontweight='bold', zorder=100)
# Add other stations
## x, y = m(-53.37, 63.88)
## m.scatter(x,y,100,marker='*',color='r', zorder=100)
## plt.text(x, y, '  FB4', horizontalalignment='left', verticalalignment='center', fontsize=10, color='r', fontweight='bold', zorder=100)
## x, y = m(-50, 60.47)
## m.scatter(x,y,100,marker='*',color='r', zorder=100)
## plt.text(x, y, 'CD3  ', horizontalalignment='right', verticalalignment='center', fontsize=10, color='r', fontweight='bold', zorder=100)
# Regions
x, y = m(-59, 45)
m.scatter(x,y,100,marker='*',color='m', zorder=10)
plt.text(x, y, ' Misaine ', horizontalalignment='right', verticalalignment='bottom', fontsize=10, color='m', fontweight='bold', zorder=100)
x, y = m(-63, 44)
m.scatter(x,y,100,marker='*',color='m', zorder=100)
plt.text(x, y, '  Emerald', horizontalalignment='left', verticalalignment='top', fontsize=10, color='m', fontweight='bold', zorder=100)
## x, y = m(np.mean(US_poly['NEC']['lon']), np.mean(US_poly['NEC']['lat'])) 
## m.scatter(x,y,100,marker='*',color='m', zorder=100)
## plt.text(x, y, 'NEC ', horizontalalignment='right', verticalalignment='top', fontsize=10, color='m', fontweight='bold', zorder=100)
## x, y = m(np.mean(US_poly['EGOM']['lon']), np.mean(US_poly['EGOM']['lat'])) 
## m.scatter(x,y,100,marker='*',color='m', zorder=100)
## plt.text(x, y, 'EGOM ', horizontalalignment='right', verticalalignment='bottom', fontsize=10, color='m', fontweight='bold', zorder=100)

# add air temperature stations
x, y = m(air_stationLon,air_stationLat)
plt.scatter(x[2],y[2],100, marker='.', color='red', zorder=100)
## plt.scatter(x,y,100, marker='.', color='red', zorder=100)
## plt.text(x[0], y[0], '  St. John''s ', horizontalalignment='right', verticalalignment='center', fontsize=13, color='red', fontweight='bold', zorder=100)
## plt.text(x[1], y[1], '  Bonavista ', horizontalalignment='right', verticalalignment='center', fontsize=13, color='red', fontweight='bold', zorder=100)
plt.text(x[2], y[2], '   Cartwright  ', horizontalalignment='right', verticalalignment='center', fontsize=13, color='red', fontweight='bold', zorder=100)
## plt.text(x[3], y[3], '  Iqaluit  ', horizontalalignment='left', verticalalignment='center', fontsize=13, color='red', fontweight='bold', zorder=100)
## plt.text(x[4], y[4], ' Nuuk  ', horizontalalignment='left', verticalalignment='center', fontsize=13, color='red', fontweight='bold', zorder=100)

# Add USA stations:
x, y = m(sNEC[0], sNEC[1])
plt.scatter(x,y,100, marker='*', color='tab:orange', zorder=100)
plt.text(x, y, ' NEC ', horizontalalignment='left', verticalalignment='top', fontsize=10, color='tab:orange', fontweight='bold', zorder=100)
x, y = m(sGB[0], sGB[1])
plt.scatter(x,y,100, marker='*', color='tab:orange', zorder=100)
plt.text(x, y, 'GB ', horizontalalignment='right', verticalalignment='top', fontsize=10, color='tab:orange', fontweight='bold', zorder=100)
x, y = m(sEGOM[0], sEGOM[1])
plt.scatter(x,y,100, marker='*', color='tab:orange', zorder=100)
plt.text(x, y, 'EGOMs ', horizontalalignment='right', verticalalignment='bottom', fontsize=10, color='tab:orange', fontweight='bold', zorder=100)
x, y = m(sWGOM[0], sWGOM[1])
plt.scatter(x,y,100, marker='*', color='tab:orange', zorder=100)
plt.text(x, y, 'WGOM ', horizontalalignment='right', verticalalignment='bottom', fontsize=10, color='tab:orange', fontweight='bold', zorder=100)
x, y = m(sNMAB[0], sNMAB[1])
plt.scatter(x,y,100, marker='*', color='tab:orange', zorder=100)
plt.text(x, y, 'NMAB ', horizontalalignment='right', verticalalignment='top', fontsize=10, color='tab:orange', fontweight='bold', zorder=100)
x, y = m(sSMAB[0], sSMAB[1])
plt.scatter(x,y,100, marker='*', color='tab:orange', zorder=100)
plt.text(x, y, 'SMAB ', horizontalalignment='right', verticalalignment='top', fontsize=10, color='tab:orange', fontweight='bold', zorder=100)


# Add NAFO divisions
## div_toplot = ['0B', '1C', '1D', '1E', '1F', '2G', '2H', '2J', '3K', '3L', '3N', '3M', '3O', '3Ps', '4Vn', '4Vs', '4W', '4X']
## for div in div_toplot:
##     div_lon, div_lat = m(nafo_div[div]['lon'], nafo_div[div]['lat'])
##     m.plot(div_lon, div_lat, 'dimgray', linewidth=2, zorder=10)
##     if (div == '3M') | (div == '4W'):
##         ax.text(np.mean(div_lon), np.mean(div_lat), div, fontsize=12, color='dimgray', fontweight='bold', zorder=100, verticalalignment='top', horizontalalignment='right')    
##     elif (div == '4X'):
##         ax.text(np.mean(div_lon), np.mean(div_lat), div, fontsize=12, color='dimgray', fontweight='bold', zorder=100, verticalalignment='top', horizontalalignment='left')
##     elif (div == '4Vn'):
##         ax.text(np.mean(div_lon), np.mean(div_lat), div, fontsize=12, color='dimgray', fontweight='bold', zorder=100, verticalalignment='center', horizontalalignment='right')
##     elif (div == '3M'):
##         ax.text(np.mean(div_lon), np.mean(div_lat), div, fontsize=12, color='dimgray', fontweight='bold', zorder=100, verticalalignment='top', horizontalalignment='center')
##     elif (div == '3Ps') | (div == '0B'):
##         ax.text(np.mean(div_lon), np.mean(div_lat), div, fontsize=12, color='dimgray', fontweight='bold', zorder=100, verticalalignment='bottom', horizontalalignment='right')
##     else:
##         ax.text(np.mean(div_lon), np.mean(div_lat), div, fontsize=12, color='dimgray', fontweight='bold', zorder=100)    

# Add USA polygons:
## poly_lon, poly_lat = m(US_poly['NEC']['lon'], US_poly['NEC']['lat'])
## m.plot(poly_lon, poly_lat, 'm', linewidth=2, zorder=10)
# Add CLS polygons:
## poly_lon, poly_lat = m(CLS_poly['CLS']['lon'], CLS_poly['CLS']['lat'])
## m.plot(poly_lon, poly_lat, 'm', linewidth=2, zorder=10)

## patches = []
## NEC = np.array([[-66.80, 42.20], [-66.00, 42.20], [-66.00, 42.60], [-66.80, 42.60], [-66.80, 42.20]])
## EGOM = np.array([[-68.00, 41.80], [-67.50, 42.10], [-67.00, 42.10], [-66.00, 42.10], [-66.00, 43.00], [-66.50, 44.20], [-67.00, 44.20], [-68.00, 44.20], [-68.00, 41.80]])
## patches.append(Polygon(NEC))
## patches.append(Polygon(EGOM))
## ax.add_collection(PatchCollection(patches, facecolor='lightgreen', edgecolor='k', linewidths=0.5))

        
# Add Ice regions
reg_toplot = ['NLab', 'SLab', 'Nfld']
for reg in reg_toplot:
     reg_lon, reg_lat = m(ice_region[reg]['lon'], ice_region[reg]['lat'])
     m.plot(reg_lon, reg_lat, 'red', linestyle='--', linewidth=2, zorder=20)
     plt.text(np.max(reg_lon), np.max(reg_lat), ' '+reg, horizontalalignment='left', verticalalignment='bottom', fontsize=10, color='r', fontweight='bold', zorder=100)

    
     
## # Add Icebergs regions
## berg_lon, berg_lat = m([-60, -36, -36, -60, -60], [48, 48, 27, 27, 48])
## m.plot(berg_lon, berg_lat, 'dimgray', linestyle='--', linewidth=2, zorder=20)

    
#### ---- Save Figure ---- ####
fig.set_size_inches(w=10, h=12)
fig.savefig(fig_name, dpi=200)
os.system('convert -trim ' + fig_name + ' ' + fig_name)
#plt.show()

