'''
This script generates the map used in CSAS ResDoc.
It is similar to another script used for NAFO SCR doc: nafo_map.py

Check this for nafo boxes:
/home/cyrf0006/AZMP/state_reports/2018/nafo.html

Frederic.Cyr@dfo-mpo.gc.ca
June 2019
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

## ---- Region parameters ---- ##
dataFile = '/home/cyrf0006/data/GEBCO/GEBCO_2014_1D.nc'
lon_0 = -50
lat_0 = 50
#lonLims = [-60, -40]
#latLims = [38, 56]
lonLims = [-70, -39.9]
latLims = [38, 65]
proj = 'merc'
decim_scale = 4
stationFile = '/home/cyrf0006/github/AZMP-NL/data/STANDARD_SECTIONS.xlsx'
AZMP_airTemp_file = '/home/cyrf0006/github/AZMP-NL/utils/airTemp_sites.xlsx'
lightdeep = cmocean.tools.lighten(cmocean.cm.deep, .9)
v = np.linspace(0, 3500, 36)
#v = np.linspace(-4000, 0, 9)
v = np.linspace(0, 5400, 55) # Olivia's
LANGUAGE = 'french'
if LANGUAGE == 'english':
    print('In English!')
    fig_name = 'map_csas_nlci.png'
elif LANGUAGE == 'french':
    print('In French!')
    fig_name = 'map_csas_nlci_FR.png'
    
## ---- Bathymetry ---- ##
print('Load and grid bathymetry')
# h5 file
h5_outputfile = 'nl_climate_bathymetry.h5'
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
nafo_div = azu.get_nafo_divisions()

## ---- Ice region ---- ##
ice_region = azu.get_ice_regions()

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
m = Basemap(projection='merc',lon_0=lon_0,lat_0=lat_0, llcrnrlon=lonLims[0],llcrnrlat=latLims[0],urcrnrlon=lonLims[1],urcrnrlat=latLims[1], resolution='h')
x,y = m(*np.meshgrid(lon,lat))
c = plt.contourf(x, y, -Z, v, cmap=lightdeep, extend='max', zorder=1)
cc = plt.contour(x, y, -Z, [100, 500, 1000, 2000, 3000, 4000], colors='lightgrey', linewidths=.5, zorder=1)
ccc = plt.contour(x, y, -Z, [1000], colors='darkgrey', linewidths=2, zorder=1)
plt.clabel(cc, inline=1, fontsize=10, colors='lightgrey', fmt='%d')

m.fillcontinents(color='tan', zorder=50);
m.drawparallels(np.arange(10,70,5), labels=[1,0,0,0], fontsize=12, fontweight='bold');
m.drawmeridians(np.arange(-80, 5, 5), labels=[0,0,0,1], fontsize=12, fontweight='bold');

## plot AZMP stations
x, y = m(stationLon[index_BI],stationLat[index_BI])
m.scatter(x,y,10,marker='o',color='r', zorder=10)
plt.text(x[-1], y[-1], ' Beachy Island', horizontalalignment='left', verticalalignment='center', fontsize=13, color='r', fontweight='bold', zorder=100)

x, y = m(stationLon[index_MB],stationLat[index_MB])
m.scatter(x,y,10,marker='o',color='r', zorder=10)
plt.text(x[-1], y[-1], ' Makkovik', horizontalalignment='left', verticalalignment='center', fontsize=13, color='r', fontweight='bold', zorder=100)

x, y = m(stationLon[index_SI],stationLat[index_SI])
m.scatter(x,y,10,marker='o',color='r', zorder=10)
plt.text(x[-1], y[-1], ' Seal Island', horizontalalignment='left', verticalalignment='center', fontsize=13, color='r', fontweight='bold', zorder=100)

x, y = m(stationLon[index_WB],stationLat[index_WB])
m.scatter(x,y,10,marker='o',color='r', zorder=10)
plt.text(x[-1], y[-1], ' White Bay', horizontalalignment='left', verticalalignment='top', fontsize=13, color='r', fontweight='bold', zorder=100)

x, y = m(stationLon[index_BB],stationLat[index_BB])
m.scatter(x,y,10,marker='o',color='r', zorder=10)
plt.text(x[-1], y[-1], ' Bonavista', horizontalalignment='left', verticalalignment='center', fontsize=13, color='r', fontweight='bold', zorder=100)

x, y = m(stationLon[index_FC],stationLat[index_FC])
m.scatter(x,y,10,marker='o',color='r', zorder=100)
plt.text(x[-6], y[-1]-50000, ' Flemish Cap', horizontalalignment='center', verticalalignment='top', fontsize=13, color='r', fontweight='bold', zorder=100)

x, y = m(stationLon[index_SEGB],stationLat[index_SEGB])
m.scatter(x,y,10,marker='o',color='r', zorder=10)
plt.text(x[-1], y[-1], ' SE Grand Bank', horizontalalignment='left', verticalalignment='center', fontsize=13, color='r', fontweight='bold', zorder=100)

x, y = m(stationLon[index_SESPB],stationLat[index_SESPB])
m.scatter(x,y,10,marker='o',color='r', zorder=10)
plt.text(x[-1], y[-1], ' SE St. Pierre Bank', horizontalalignment='right', verticalalignment='top', fontsize=13, color='r', fontweight='bold', zorder=100)

x, y = m(stationLon[index_SWSPB],stationLat[index_SWSPB])
m.scatter(x,y,10,marker='o',color='r', zorder=10)
plt.text(x[-1], y[-1], ' SW St. Pierre Bank', horizontalalignment='right', verticalalignment='top', fontsize=13, color='r', fontweight='bold', zorder=100)


x, y = m(stationLon[index_S27],stationLat[index_S27])
m.scatter(x[0],y[0],100,marker='*',color='r', zorder=100)
plt.text(x[0], y[0]+50000, 'S27', horizontalalignment='left', verticalalignment='bottom', fontsize=13, color='r', fontweight='bold', zorder=100)

# add air temperature stations
x, y = m(air_stationLon,air_stationLat)
plt.scatter(x,y,100, marker='.', color='saddlebrown', zorder=100)
plt.text(x[0], y[0], '  St. John''s ', horizontalalignment='right', verticalalignment='center', fontsize=13, color='saddlebrown', fontweight='bold', zorder=100)
plt.text(x[1], y[1], '  Bonavista ', horizontalalignment='right', verticalalignment='center', fontsize=13, color='saddlebrown', fontweight='bold', zorder=100)
plt.text(x[2], y[2], '   Cartwright  ', horizontalalignment='right', verticalalignment='center', fontsize=13, color='saddlebrown', fontweight='bold', zorder=100)
plt.text(x[3], y[3], '  Iqaluit  ', horizontalalignment='left', verticalalignment='center', fontsize=13, color='saddlebrown', fontweight='bold', zorder=100)
plt.text(x[4], y[4], ' Nuuk  ', horizontalalignment='left', verticalalignment='center', fontsize=13, color='saddlebrown', fontweight='bold', zorder=100)


# Add NAFO divisions
div_toplot = ['2G', '2H', '2J', '3K', '3L', '3N', '3M', '3O', '3Ps']
for div in div_toplot:
    div_lon, div_lat = m(nafo_div[div]['lon'], nafo_div[div]['lat'])
    m.plot(div_lon, div_lat, 'black', linewidth=2, zorder=10)
    if (div == '3M') | (div == '4X'):
        ax.text(np.mean(div_lon), np.mean(div_lat), div, fontsize=12, color='black', fontweight='bold', zorder=100, verticalalignment='top', horizontalalignment='right')    
    elif (div == '3M'):
        ax.text(np.mean(div_lon), np.mean(div_lat), div, fontsize=12, color='black', fontweight='bold', zorder=100, verticalalignment='top', horizontalalignment='center')
    elif (div == '3Ps'):
        ax.text(np.mean(div_lon), np.mean(div_lat), div, fontsize=12, color='black', fontweight='bold', zorder=100, verticalalignment='bottom', horizontalalignment='right')
    else:
        ax.text(np.mean(div_lon), np.mean(div_lat), div, fontsize=12, color='black', fontweight='bold', zorder=100)    

# Add Ice regions
reg_toplot = ['NLab', 'SLab', 'Nfld']

for reg in reg_toplot:
    reg_lon, reg_lat = m(ice_region[reg]['lon'], ice_region[reg]['lat'])
    m.plot(reg_lon, reg_lat, 'darkmagenta', linestyle='--', linewidth=2, zorder=20)

    if LANGUAGE == 'english':
        if reg == 'NLab':
            ax.text(np.max(reg_lon), np.max(reg_lat), 'North. Labrador', horizontalalignment='right', verticalalignment='bottom', fontsize=13, color='darkmagenta', fontweight='bold', zorder=100)
        elif reg == 'SLab':
            ax.text(np.max(reg_lon), np.max(reg_lat), 'South. Labrador', horizontalalignment='right', verticalalignment='bottom', fontsize=13, color='darkmagenta', fontweight='bold', zorder=100)
        elif reg == 'Nfld':
            ax.text(np.max(reg_lon), np.max(reg_lat), 'Newfoundland', horizontalalignment='right', verticalalignment='bottom', fontsize=13, color='darkmagenta', fontweight='bold', zorder=100)
    elif LANGUAGE == 'french':
        if reg == 'NLab':
            ax.text(np.max(reg_lon), np.max(reg_lat), 'Labrador Nord', horizontalalignment='right', verticalalignment='bottom', fontsize=13, color='darkmagenta', fontweight='bold', zorder=100)
        elif reg == 'SLab':
            ax.text(np.max(reg_lon), np.max(reg_lat), 'Labrador Sud', horizontalalignment='right', verticalalignment='bottom', fontsize=13, color='darkmagenta', fontweight='bold', zorder=100)
        elif reg == 'Nfld':
            ax.text(np.max(reg_lon), np.max(reg_lat), 'Terre-Neuve', horizontalalignment='right', verticalalignment='bottom', fontsize=13, color='darkmagenta', fontweight='bold', zorder=100)      
# Add Icebergs regions
berg_lon, berg_lat = m([-60, -36, -36, -60, -60], [48, 48, 27, 27, 48])
m.plot(berg_lon, berg_lat, 'cyan', linestyle='--', linewidth=2, zorder=20)

    
#### ---- Save Figure ---- ####
fig.set_size_inches(w=10, h=12)
fig.savefig(fig_name, dpi=200)
os.system('convert -trim ' + fig_name + ' ' + fig_name)
#plt.show()

