'''
A map of G. Han's transport section
'''



import netCDF4
from mpl_toolkits.basemap import Basemap
import numpy as  np
import pandas as pd
import matplotlib.pyplot as plt
import openpyxl, pprint
import shapefile 
import cmocean
import os

## ---- Region parameters ---- ##
dataFile = '/home/cyrf0006/data/GEBCO/GRIDONE_1D.nc'
lon_0 = -50
lat_0 = 50
lonLims = [-70, -40]
latLims = [40, 65]
proj = 'merc'
decim_scale = 4
stationFile = '/home/cyrf0006/github/AZMP-NL/data/STANDARD_SECTIONS.xlsx'
fig_name = 'AZMP_altimetry_sections.png'


## ---- Altimetry Sections ---- ##
infile_sections = '/home/cyrf0006/AZMP/altimetry/LC_transport_sections.csv'
df = pd.read_csv(infile_sections, header=None)

## ---- Bathymetry ---- ####
v = np.linspace(0, 3500, 36)
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


## ## ---- NAFO divisions ---- ##
## myshp = open('/home/cyrf0006/AZMP/utils/NAFO_divisions/Divisions/Divisions.shp', 'rb')
## mydbf = open('/home/cyrf0006/AZMP/utils/NAFO_divisions/Divisions/Divisions.dbf', 'rb')
## r = shapefile.Reader(shp=myshp, dbf=mydbf)
## records = r.records()
## shapes = r.shapes()



## ---- plot ---- ##
fig = plt.figure(1)
m = Basemap(projection='merc',lon_0=lon_0,lat_0=lat_0, llcrnrlon=lonLims[0],llcrnrlat=latLims[0],urcrnrlon=lonLims[1],urcrnrlat=latLims[1], resolution='h')
x,y = m(*np.meshgrid(lon,lat))
c = m.contourf(x, y, np.flipud(-Z), v, cmap=cmocean.cm.deep, extend="max", alpha=.5);
cc = m.contour(x, y, np.flipud(-Z), [100, 500, 1000, 3000, 4000], colors='lightgrey', linewidths=.5);
plt.clabel(cc, inline=1, fontsize=10, colors='gray', fmt='%d')
ccc = m.contour(x, y, np.flipud(-Z), [0], colors='black');
#c = m.contourf(x, y, np.flipud(Z), v, cmap=plt.cm.PuBu_r, extend="min");
m.fillcontinents(color='peru');
m.drawparallels(np.arange(10,70,10), labels=[1,0,0,0], fontsize=12, fontweight='bold');
m.drawmeridians(np.arange(-80, 5, 10), labels=[0,0,0,1], fontsize=12, fontweight='bold');

# plot section
for idx in df.index:
    x, y = m(np.array([df.iloc[idx][1],df.iloc[idx][3]]),np.array([df.iloc[idx][2],df.iloc[idx][4]]))  
    m.plot(x,y,'-r')
    # track names
    plt.text(x[-1], y[-1], '  ' + df.iloc[idx][0], horizontalalignment='left', verticalalignment='center', fontsize=10, color='r', fontweight='bold')
    

## for idx, rec in enumerate(records):
##     if rec[-1] == '':
##         continue
##     else:
##         coords = np.array(shapes[idx].points)
##         x,y = m(coords[:,0], coords[:,1])
##         m.plot(x,y,color='black')
##         bbox = np.array(shapes[idx].bbox)
##         x,y = m(np.mean([bbox[0], bbox[2]]), np.mean([bbox[1], bbox[3]]))                                                    
##         plt.text(x, y, rec[-1], fontsize=10, color='black', fontweight='bold')




#### ---- Save Figure ---- ####
fig.set_size_inches(w=12, h=12)
fig.set_dpi(200)
fig.savefig(fig_name, dpi=200)
os.system('convert -trim ' + fig_name + ' ' + fig_name)

