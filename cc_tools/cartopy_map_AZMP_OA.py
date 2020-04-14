# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 09:34:09 2019
add some comments
@author: gibbo
"""
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.colors import from_levels_and_colors
import matplotlib as mpl
import numpy as np
import seaborn as sns
import datetime
from scipy import stats
import water_masses as wm
#import seawater as swx
import cmocean
import cmocean.cm as cmo
import cartopy. crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.feature as cpf
from cartopy.mpl.geoaxes import GeoAxes
import matplotlib.ticker as mticker
import os
import netCDF4
import io
import h5py
import sys
import cc_variable_list as vl

  
#### Load Carbonate data into a Pandas DataFrame ######
df = pd.read_excel('C:\Users\gibbo\Documents\data\AZMP_OA\AZMP_OA_CO2stats.xlsx')
pd.set_option('display.max_rows', 500)

df = df.set_index('timestamp', drop=False)

#####################to remove unwanted sections or stations##################
#SSdrop = ('L3')#, 'YL', 'PS', 'BANQ', 'PL', 'SG', 'SPB', 'BP', 'LCC', 'LCM', 'VB')
#for i in SSdrop:
#    df.drop((df.loc[df['StationID'].str.contains(i, na=False)].index), inplace=True)
    
############################set data parameters here#######################
my_region=''
my_year='2015'
my_season='fall'
my_depth='bottom'
my_variable='Omega_A'

print (my_year+' '+my_season+' '+my_depth+' '+my_variable)

###############to extract the variable settings for plotting##################
v = vl.variable_parameters(my_variable)
num_levels = v[0]
vmin = v[1]
vmax = v[2]
midpoint = v[3]
colors = v[4]
ticks = v[5]
axis_label = v[6]
extent = v[7]


df = df.loc[my_year] ###locates year
df=df[df.Region.str.contains(my_region)] ###locates region
    
#######locates season##########
if my_season == 'spring':
    df = df[(df.index.month>=3) & (df.index.month<=6)]
elif my_season == 'summer':
    df = df[~df.Region.str.contains('SS')]
    df = df.assign(x=df.index.strftime('%m-%d')).query("'07-01' <= x <= '10-14'").drop('x',1)
elif my_season == 'fall':
    dfng = df[~df.Region.str.contains('SS')]
    dfng = dfng.assign(x=dfng.index.strftime('%m-%d')).query("'10-15' <= x <= '12-31'").drop('x',1)
    dfss = df[df.Region.str.contains('SS')]
    dfss = dfss[(dfss.index.month>=9) & (dfss.index.month<=12)]
    df = pd.concat([dfng, dfss], axis=0)
else:        
    print 'All seasons selected'


if df[my_variable].isna().values.all():# or df.size == 0:
    print ('!!! no data for this season !!!')
    sys.exit()

df.dropna(subset=[my_variable], axis=0, inplace=True)
df = df.reset_index(drop=True)

#######locates depth -- either surface, bottom, or within a range of depths############
if my_depth == 'surface':
    df = df.loc[df.groupby('StationID')['depth'].idxmin()] #group by station then pull "min or max depth"
    df = df.loc[df.depth <20] #take all depths >10m (for bottom) to eliminate lone surface samples, all depths <20m (for surface) to eliminate lone deep sample
if my_depth == 'bottom':
    df = df.loc[df.groupby('StationID')['depth'].idxmax()] #group by station then pull "min or max depth"
    df = df.loc[df.depth >10] #take all depths >10m (for bottom) to eliminate lone surface samples, all depths <20m (for surface) to eliminate lone deep sample
if my_depth == 'range':
    df = df.loc[(df.depth >=135) & (df.depth <=165)]
    df = df.loc[df.groupby('StationID')['depth'].idxmax()] #group by station then pull "min or max depth"


############create bathymetry################
os.environ["CARTOPY_USER_BACKGROUNDS"] = "C:\ProgramData\Anaconda2\Lib\site-packages\cartopy\BG"
dataFile = 'C:\Users\gibbo\Documents\data\GRIDONE_1D.nc'
###map boundaries###
lonLims = [-72, -41.5]  
latLims = [40.5, 58.4]

##### this is the dataset that will be plotted##########
lat_data = np.array(df.latitude)
lon_data = np.array(df.longitude)
data = np.array(df[my_variable])
lon_data = lon_data[~np.isnan(data)]
lat_data = lat_data[~np.isnan(data)]
data = data[~np.isnan(data)]

print('Load and grid bathymetry')
##Load data
h5_outputfile = 'cc_bathymetry.h5'

if os.path.isfile(h5_outputfile):
     print [h5_outputfile + ' exists! Reading directly']
     h5f = h5py.File(h5_outputfile,'r')
     lon = h5f['lon'][:]
     lat = h5f['lat'][:]
     Z = h5f['Z'][:]
     h5f.close()

#else:
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
## interpolate data on regular grid (temperature grid)
## Reshape data
zz = dataset.variables['z']
Z = zz[:].reshape(ny, nx)
Z = np.flipud(Z) # <------------ important!!!
# Reduce data according to Region params
idx_lon = np.where((lon>=lonLims[0]) & (lon<=lonLims[1]))
idx_lat = np.where((lat>=latLims[0]) & (lat<=latLims[1]))
Z = Z[idx_lat[0][0]:idx_lat[0][-1]+1, idx_lon[0][0]:idx_lon[0][-1]+1]
lon = lon[idx_lon[0]]
lat = lat[idx_lat[0]]
print(' -> Done!')

 # Save data for later use
if np.size(h5_outputfile):
       h5f = h5py.File(h5_outputfile, 'w')
       h5f.create_dataset('lon', data=lon)
       h5f.create_dataset('lat', data=lat)
       h5f.create_dataset('Z', data=Z)
       h5f.close()

#############draw map with bathymetry########################################################
fig = plt.figure(figsize=(7,5))
ax = fig.add_subplot(111, projection=ccrs.Mercator())
ax.set_extent([-72, -41.5, 40.5, 58.4], crs=ccrs.PlateCarree())
ax.add_feature(cpf.NaturalEarthFeature('physical', 'coastline', '50m', edgecolor='k', alpha=0.7, linewidth=0.6, facecolor='black'), zorder=1)#cpf.COLORS['land']))
m=ax.gridlines(linewidth=0.5, color='black', draw_labels=True, alpha=0.5)
m.xlabels_top=False
m.ylabels_right=False
m.xlocator = mticker.FixedLocator([-75, -70, -60, -50, -40])
m.ylocator = mticker.FixedLocator([40, 45, 50, 55, 60, 65])
m.xformatter = LONGITUDE_FORMATTER
m.yformatter = LATITUDE_FORMATTER
m.ylabel_style = {'size': 7, 'color': 'black', 'weight':'bold'}
m.xlabel_style = {'size': 7, 'color': 'black', 'weight':'bold'}
lightdeep = cmocean.tools.lighten(cmo.deep, 0.5)
ls = np.linspace(0, 5500, 20)
c = plt.contourf(lon, lat, -Z, ls, transform=ccrs.PlateCarree(), cmap=lightdeep, extend='max', zorder=5)
cc = plt.contour(lon, lat, -Z, [100, 500, 1000, 2000, 3000, 4000, 5000], colors='silver', linewidths=.5, transform=ccrs.PlateCarree(), zorder=10)
plt.clabel(cc, inline=True, fontsize=7, fmt='%i')

################to adjust the colorbar###################################
levels = np.linspace(vmin, vmax, num_levels)
midp = np.mean(np.c_[levels[:-1], levels[1:]], axis=1)
vals = np.interp(midp, [vmin, midpoint, vmax], [0, 0.5, 1])
colors = colors(vals)
colors=np.concatenate([[colors[0,:]], colors, [colors[-1,:]]],0)
cmap, norm = from_levels_and_colors(levels, colors, extend=extent)

#######################plot data onto map#################################
s = ax.scatter(lon_data,lat_data, c=data, s=40, lw=0.3, edgecolor='black', vmin=vmin, vmax=vmax, cmap=cmap, transform=ccrs.Geodetic(), zorder=20)
cax = plt.axes([0.83,0.125,0.03,0.756])
cb = plt.colorbar(s, cax=cax, extend=extent, ticks=ticks)
cb.ax.tick_params(labelsize=8)
cb.set_label(axis_label, fontsize=12, fontweight='normal')


#fig.savefig('C:\Users\gibbo\Documents\AZMP_OA_'+my_year+'_'+my_season+'_'+my_variable+'_'+my_depth+'.png', 
#            format='png', dpi=500, bbox_inches='tight')
#fig.savefig('C:\Users\gibbo\Documents\carbonates\other\maps\AZMP_empty.svg', format='svg', dpi=1000, bbox_inches='tight')