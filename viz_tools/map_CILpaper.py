'''
This script generates the map used in Heather's CIL paper

Run in /home/cyrf0006/AZMP/utils

Frederic.Cyr@dfo-mpo.gc.ca
January 2024
'''


import os
import h5py                                                                
import netCDF4
import numpy as  np
import matplotlib.pyplot as plt
import openpyxl, pprint
import pandas as pd
import cmocean as cmo
#import azmp_utils as azu
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
# Cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature

## ---- Region parameters ---- ##
dataFile = '/home/cyrf0006/data/GEBCO/GEBCO_2014_1D.nc'
lon_0 = -50
lat_0 = 50
lonLims = [-70, -35]
latLims = [38, 68]
proj = 'merc'
decim_scale = 4
stationFile = '/home/cyrf0006/github/AZMP-NL/data/STANDARD_SECTIONS.xlsx'
fig_name = 'map_CILpaper.png'
lightdeep = cmo.tools.lighten(cmo.cm.deep, .9)
v = np.linspace(0, 3500, 36)
#v = np.linspace(-4000, 0, 9)
v = np.linspace(0, 5400, 55) # Olivia's (good one, but loose channels on the shelf)
#v = np.linspace(0, 4000, 41) # t

## ---- Bathymetry ---- ##
print('Load and grid bathymetry')
# h5 file
h5_outputfile = 'nl_stacfis_bathymetry.h5'
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
fig = plt.figure(figsize=(10,12))

#Set up the map
land_10m = cfeature.NaturalEarthFeature('physical','land','50m',facecolor='tan')
ax = plt.subplot(projection=ccrs.Mercator())

#Plot the coastline
ax.set_facecolor('white')
ax.add_feature(land_10m,zorder=4)
ax.set_extent([lonLims[0],lonLims[1],latLims[0],65])

#Plot the bathymetry

plt.contour(
     lon,
     lat,
     -Z,
     levels=[100,200,500,1000,2000,3000,4000],
     colors='grey',
     linewidths=0.75,
     zorder=3,
     transform=ccrs.PlateCarree())

depth = plt.contourf(
     lon,
     lat,
     -Z,
     levels=[10,100,150,200,300,500,1000,2000,3000,4000],
     cmap=lightdeep,
     extend='max',
     zorder=2,
     transform=ccrs.PlateCarree())

# Add Stations
x, y = stationLon[index_BI],stationLat[index_BI]
plt.plot(x,y, markersize=1, marker='o',color='r', zorder=10, transform=ccrs.PlateCarree())
plt.text(x[-1], y[-1], ' Beachy Island', horizontalalignment='left', verticalalignment='center', fontsize=13, color='r', fontweight='bold',zorder=6,transform=ccrs.PlateCarree())

x, y = stationLon[index_MB],stationLat[index_MB]
plt.plot(x,y, markersize=1, marker='o',color='r', zorder=10, transform=ccrs.PlateCarree())
plt.text(x[-1], y[-1], ' Makkovik Bank', horizontalalignment='left', verticalalignment='center', fontsize=13, color='r', fontweight='bold',zorder=6,transform=ccrs.PlateCarree())

x, y = stationLon[index_SI],stationLat[index_SI]
plt.plot(x,y, markersize=1, marker='o',color='r', zorder=10, transform=ccrs.PlateCarree())
plt.text(x[-1], y[-1], ' Seal Island', horizontalalignment='left', verticalalignment='center', fontsize=13, color='r', fontweight='bold',zorder=6,transform=ccrs.PlateCarree())

x, y = stationLon[index_WB],stationLat[index_WB]
plt.plot(x,y, markersize=1, marker='o',color='r', zorder=10, transform=ccrs.PlateCarree())
plt.text(x[-1], y[-1], ' White Bay', horizontalalignment='left', verticalalignment='center', fontsize=13, color='r', fontweight='bold', zorder=6,transform=ccrs.PlateCarree())

x, y = stationLon[index_BB],stationLat[index_BB]
plt.plot(x,y, markersize=1, marker='o',color='r', zorder=10, transform=ccrs.PlateCarree())
plt.text(x[-1], y[-1], ' Bonavista', horizontalalignment='left', verticalalignment='center', fontsize=13, color='r', fontweight='bold', zorder=6,transform=ccrs.PlateCarree())

x, y = stationLon[index_FC],stationLat[index_FC]
plt.plot(x,y, markersize=1, marker='o',color='r', zorder=10, transform=ccrs.PlateCarree())
plt.text(x[-6], y[-1], ' Flemish Cap', horizontalalignment='center', verticalalignment='bottom', fontsize=13, color='r', fontweight='bold', zorder=6,transform=ccrs.PlateCarree())

## x, y = m(stationLon[index_SEGB],stationLat[index_SEGB])
## plt.plot(x,y, markersize=1, marker='o',color='r', zorder=10, transform=ccrs.PlateCarree())
## plt.text(x[-1], y[-1], ' SE Grand Bank', horizontalalignment='left', verticalalignment='center', fontsize=13, color='r', fontweight='bold', zorder=6,transform=ccrs.PlateCarree())

## x, y = m(stationLon[index_SESPB],stationLat[index_SESPB])
## plt.plot(x,y, markersize=1, marker='o',color='r', zorder=10, transform=ccrs.PlateCarree())
## plt.text(x[-1], y[-1], ' SE St. Pierre Bank', horizontalalignment='right', verticalalignment='top', fontsize=13, color='r', fontweight='bold', zorder=6,transform=ccrs.PlateCarree())

## x, y = m(stationLon[index_SWSPB],stationLat[index_SWSPB])
## plt.plot(x,y, markersize=1, marker='o',color='r', zorder=10, transform=ccrs.PlateCarree())
## plt.text(x[-1], y[-1], ' SW St. Pierre Bank', horizontalalignment='right', verticalalignment='top', fontsize=13, color='r', fontweight='bold', zorder=6,transform=ccrs.PlateCarree())

# Add high-res stations
x, y = stationLon[index_S27],stationLat[index_S27]
plt.plot(x[0],y[0], markersize=10, marker='*',color='r', zorder=100, transform=ccrs.PlateCarree())
plt.text(x[0], y[0], ' S27', horizontalalignment='left', verticalalignment='center', fontsize=13, color='r', fontweight='bold', zorder=100, transform=ccrs.PlateCarree())

#Latitude and Longitude
gl = ax.gridlines(draw_labels=['left','bottom'],xlocs=[-40,-45,-50,-55,-60,-65,-70],ylocs=[40,45,50,55,60,65],dms=True,x_inline=False,y_inline=False,linestyle='--')
gl.xlabel_style = {'size': 8}
gl.ylabel_style = {'size': 8}

#Colourbar
cax = fig.add_axes([ax.get_position().xmin,0.05,ax.get_position().xmax-ax.get_position().xmin,0.02])
cb = plt.colorbar(depth, cax=cax, orientation='horizontal')
cb.set_label('Depth (m)', fontsize=8, fontweight='normal')
cb.ax.tick_params(labelsize=8)

    
#### ---- Save Figure ---- ####
fig.savefig(fig_name, dpi=200)
os.system('convert -trim ' + fig_name + ' ' + fig_name)
#plt.show()

