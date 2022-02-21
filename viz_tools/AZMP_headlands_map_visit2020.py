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

import os
import netCDF4
import h5py                                                                
os.environ['PROJ_LIB'] = '/home/cyrf0006/anaconda3/share/proj'
from mpl_toolkits.basemap import Basemap
import numpy as  np
import matplotlib.pyplot as plt
import openpyxl, pprint
import pandas as pd
import cmocean
import cmocean.cm as cmo

import cartopy.crs as ccrs
import cartopy.feature as cpf
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

## ---- Region parameters ---- ##
dataFile = '/home/cyrf0006/data/GEBCO/GEBCO_2014_1D.nc'
lon_0 = -50
lat_0 = 50
#lonLims = [-60, -52]
#latLims = [46, 52]
lonLims = [-60, -52]
latLims = [46, 52]
proj = 'merc'
decim_scale = 4
stationFile = '/home/cyrf0006/github/AZMP-NL/utils/Headlands_sites.xlsx'
fig_name = 'map_LTTMP_visit2020_admiral.png'
v = np.linspace(0, 500, 11)
newStationFile = '/home/cyrf0006/AZMP/Headlands/LTTMPup_stations_20200928.csv'

sites_2020 = ['Hampden', 'South Brook', 'Lewisporte', 'Comfort Cove', 'Cape Freels', 'Stock Cove', 'Bonavista', 'Lower Lance Cove', 'Dildo', 'Holyrood', 'St Johns', 'Admiral Beach']
#sites_2020 = ['Springdale', 'Lewisporte', 'Comfort Cove', 'Cape Freels', 'Bonavista', 'Lower Lance Cove', 'Dildo', 'Holyrood', 'St Johns', 'Admiral Beach']



## ---- Bathymetry ---- ##
print('Load and grid bathymetry')
# h5 file
h5_outputfile = 'hl_bathymetry.h5'
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
df = pd.read_excel(stationFile)
dfn = pd.read_csv(newStationFile)
# Keep only 2020 visit
dfn = dfn[dfn['Station'].isin(sites_2020)]
df = df[df['Site'].isin(sites_2020)]

## ---- plot ---- ##
#plt.title('Standard air temperature sites')
print('--- Now plot ---')
fig = plt.figure(figsize=(7,5))
ax = fig.add_subplot(111, projection=ccrs.Mercator())
ax.set_extent([lonLims[0], lonLims[1], latLims[0], latLims[1]], crs=ccrs.PlateCarree())
ax.add_feature(cpf.NaturalEarthFeature('physical', 'coastline', '50m', edgecolor='k', alpha=0.7, linewidth=0.6, facecolor='black'), zorder=1)
m=ax.gridlines(linewidth=0.5, color='black', draw_labels=True, alpha=0.5)
m.xlabels_top=False
m.ylabels_right=False
m.xlocator = mticker.FixedLocator([-60, -58, -56, -54, -52])
m.ylocator = mticker.FixedLocator([46, 48, 50, 52])
m.xformatter = LONGITUDE_FORMATTER
m.yformatter = LATITUDE_FORMATTER
m.ylabel_style = {'size': 7, 'color': 'black', 'weight':'bold'}
m.xlabel_style = {'size': 7, 'color': 'black', 'weight':'bold'}
lightdeep = cmocean.tools.lighten(cmo.deep, 0.5)
c = plt.contourf(lon, lat, -Z, v, transform=ccrs.PlateCarree(), cmap=lightdeep, extend='max', zorder=5)
cc = plt.contour(lon, lat, -Z, [100, 500], colors='silver', linewidths=.5, transform=ccrs.PlateCarree(), zorder=10)
plt.clabel(cc, inline=True, fontsize=7, fmt='%i')

# plot Headlands
sites = df.Site.values
lats = df.Lat.values
lons = df.Lon.values
ax.plot(lons, lats, '.', color='palevioletred', transform=ccrs.PlateCarree(), markersize=15, zorder=10)

# Plot other sites
df_smart = dfn[dfn.comments.str.contains('smartAtlantic')]
ax.plot(df_smart['Longitude(DecDeg)'], df_smart['Latitude(DecDeg)'], '.', color='green', transform=ccrs.PlateCarree(), markersize=15, zorder=10)

df_new = dfn[dfn.comments.str.contains('new site')]
ax.plot(df_new['Longitude(DecDeg)'], df_new['Latitude(DecDeg)'], '.', color='cyan', transform=ccrs.PlateCarree(), markersize=15, zorder=10)

## df_aq = dfn[dfn.comments.str.contains('Aquaculture')]
## ax.plot(df_aq['Longitude(DecDeg)'], df_aq['Latitude(DecDeg)'], '.', color='salmon', transform=ccrs.PlateCarree(), markersize=15, zorder=10)

df_chs = dfn[dfn.comments.str.contains('CHS')]
ax.plot(df_chs['Longitude(DecDeg)'], df_chs['Latitude(DecDeg)'], '.', color='red', transform=ccrs.PlateCarree(), markersize=15, zorder=10)

# add text Headlands
for idx, name in enumerate(sites):
    if (name == "Bristol's Hope") | (name == "Upper Gullies") | (name == "Comfort Cove") | (name == "Cape Freels") | (name == "Stock Cove"):
        ax.text(lons[idx], lats[idx], name+' ', horizontalalignment='right', verticalalignment='top', fontsize=10, color='palevioletred', fontweight='bold', transform=ccrs.PlateCarree(), zorder=10)
    else:
        ax.text(lons[idx], lats[idx], name+' ', horizontalalignment='right', verticalalignment='bottom', fontsize=10, color='palevioletred', fontweight='bold', transform=ccrs.PlateCarree(), zorder=10)

# Station name SmartAtl.
for idx, name in enumerate(df_smart.Station):
    if (name == "St Johns"):
        ax.text(df_smart['Longitude(DecDeg)'].iloc[idx], df_smart['Latitude(DecDeg)'].iloc[idx], ' '+name, horizontalalignment='left', verticalalignment='bottom', fontsize=10, color='green', fontweight='bold', transform=ccrs.PlateCarree(), zorder=10)
    else:
        ax.text(df_smart['Longitude(DecDeg)'].iloc[idx], df_smart['Latitude(DecDeg)'].iloc[idx], name+' ', horizontalalignment='right', verticalalignment='top', fontsize=10, color='green', fontweight='bold', transform=ccrs.PlateCarree(), zorder=10)
  
## # Station name Aquaculture.
## for idx, name in enumerate(df_aq.Station):
##     ax.text(df_aq['Longitude(DecDeg)'].iloc[idx], df_aq['Latitude(DecDeg)'].iloc[idx], name+' ', horizontalalignment='right', verticalalignment='bottom', fontsize=10, color='salmon', fontweight='bold', transform=ccrs.PlateCarree(), zorder=10)

# Station name New sites.
for idx, name in enumerate(df_new.Station):
    if (name == "Port au Choix"):
        ax.text(df_new['Longitude(DecDeg)'].iloc[idx], df_new['Latitude(DecDeg)'].iloc[idx], ' '+name, horizontalalignment='left', verticalalignment='center', fontsize=10, color='cyan', fontweight='bold', transform=ccrs.PlateCarree(), zorder=10)
    else:
        ax.text(df_new['Longitude(DecDeg)'].iloc[idx], df_new['Latitude(DecDeg)'].iloc[idx], name+'  ', horizontalalignment='right', verticalalignment='center', fontsize=10, color='cyan', fontweight='bold', transform=ccrs.PlateCarree(), zorder=10)
    
# Station name CHS
for idx, name in enumerate(df_chs.Station):
    ax.text(df_chs['Longitude(DecDeg)'].iloc[idx], df_chs['Latitude(DecDeg)'].iloc[idx], name+'  ', horizontalalignment='right', verticalalignment='top', fontsize=10, color='red', fontweight='bold', transform=ccrs.PlateCarree(), zorder=10)

    
# Custom legend
import matplotlib.lines as mlines
legend_elements = [mlines.Line2D([],[], marker='.',linestyle='None', color='palevioletred', markersize=15, label='LTTMP'),
                   mlines.Line2D([],[], marker='.',linestyle='None', color='green', markersize=15, label='Smart Atl.'),
                   ## mlines.Line2D([],[], marker='.',linestyle='None', color='salmon', markersize=15, label='Aquaculture'),
                   mlines.Line2D([],[], marker='.',linestyle='None', color='red', markersize=15, label='CHS'),
                   mlines.Line2D([],[], marker='.',linestyle='None', color='cyan', markersize=15, label='New sites')]
ax.legend(handles=legend_elements)
      
# Save    
fig.set_size_inches(w=10, h=12)
fig.savefig(fig_name, dpi=200)
os.system('convert -trim ' + fig_name + ' ' + fig_name)


