# -*- coding: utf-8 -*-
"""

To produce AZMP CC SAR figure
    
"""
import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.colors import from_levels_and_colors
import numpy as np
import h5py
import cmocean
import cmocean.cm as cmo
import cartopy. crs as ccrs
import cartopy.feature as cpf
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cc_variable_list_NL as vl

YEAR = 2019
SEASON = 'summer'
VARIABLE = 'Omega_Aragonite_(unitless)'
#VARIABLE = 'pH_Total_(total_scale)'
#VARIABLE = 'Oxygen_Saturation_(%)'
#VARIABLE = 'Dissolved_Oxygen_(mL/L)'
DEPTH='bottom'
#DEPTH='bottom'
# For colorbar:
v = vl.variable_parameters(VARIABLE)
num_levels = v[0]
vmin = v[1]
vmax = v[2]
midpoint = v[3]
colors = v[4]
ticks = v[5]
axis_label = v[6]
extent = v[7]

# Figure name
if VARIABLE == 'Omega_Aragonite_(unitless)':
    FIG_VAR = 'OmegaA'
elif VARIABLE == 'pH_Total_(total_scale)':
    FIG_VAR = 'pH'
elif VARIABLE == 'Oxygen_Saturation_(%)':
    FIG_VAR = 'DO_perc'    
elif VARIABLE == 'Dissolved_Oxygen_(mL/L)':
    FIG_VAR = 'DO'

# Read the entire AZMP dataset
df = pd.read_csv('/home/cyrf0006/github/AZMP-NL/datasets/carbonates/AZMP_carbon_data.csv', delimiter=',')

# Set index
df.set_index('Timestamp', inplace=True)
df.index = pd.to_datetime(df.index)

# Only year in review
df = df[df.index.year==YEAR]

# Only NL region
df = df[df.Region=='NL']

# Season
if SEASON == 'spring':
    df = df[(df.index.month>=3) & (df.index.month<=6)]
elif SEASON == 'summer':
    #df = df[~df.Region.str.contains('MAR')]
    #df = df.assign(x=df.index.strftime('%m-%d')).query("'07-01' <= x <= '10-14'").drop('x',1)
    df = df[(df.index.month>=6) & (df.index.month<=10)]
elif SEASON == 'fall':
    dfng = df[~df.Region.str.contains('MAR')]
    dfng = dfng.assign(x=dfng.index.strftime('%m-%d')).query("'10-15' <= x <= '12-31'").drop('x',1)
    dfss = df[df.Region.str.contains('MAR')]
    dfss = dfss[(dfss.index.month>=9) & (dfss.index.month<=12)]
    df = pd.concat([dfng, dfss], axis=0)
else:        
    print('All seasons selected')


if df[VARIABLE].isna().values.all():# or df.size == 0:
    print ('!!! no data for this season !!!')
    sys.exit()

df.dropna(subset=[VARIABLE, 'pH_Total_(total_scale)', 'Omega_Aragonite_(unitless)', 'Station_Name', 'Depth_(dbar)'], axis=0, inplace=True)
df = df.reset_index(drop=True)
# Set depth as float (had some problem with some data set)
df = df.astype({'Depth_(dbar)':'float', VARIABLE:'float'})  
#######locates depth -- either surface, bottom, or within a range of depths############
if DEPTH == 'surface':
    df = df.loc[df.groupby('Station_Name')['Depth_(dbar)'].idxmin()] #group by station then pull "min or max depth"
    df = df.loc[df['Depth_(dbar)'] <20] #take all depths >10m (for bottom) to eliminate lone surface samples, all depths <20m (for surface) to eliminate lone deep sample
if DEPTH == 'bottom':
    df = df.loc[df.groupby('Station_Name')['Depth_(dbar)'].idxmax()] #group by station then pull "min or max depth"
    df = df.loc[df['Depth_(dbar)'] >10] #take all depths >10m (for bottom) to eliminate lone surface samples

if DEPTH == 'range':
    df = df.loc[(df['Depth_(dbar)'] >=135) & (df.depth <=165)]
    df = df.loc[df.groupby('Station_Name')['Depth_(dbar)'].idxmax()] #group by station then pull "min or max depth"


############create bathymetry################
dataFile = '/home/cyrf0006/data/GEBCO/GEBCO_2014_1D.nc'
###map boundaries###
if SEASON == 'summer':
    lonLims = [-65, -41.5]  
    latLims = [45, 58.4]
else:
    lonLims = [-58, -41.5]  
    latLims = [40.5, 53]   

##### this is the dataset that will be plotted ##########
lat_data = np.array(df['Latitude_(degNorth)'])
lon_data = np.array(df['Longitude_(degEast)'])
data = np.array(df[VARIABLE])
lon_data = lon_data[~np.isnan(data)]
lat_data = lat_data[~np.isnan(data)]
data = data[~np.isnan(data)]

print('Load and grid bathymetry')
# h5 file
h5_outputfile = '/home/cyrf0006/AZMP/oa/cc_bathymetry.h5'
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

    # Save data for later use
    h5f = h5py.File(h5_outputfile, 'w')
    h5f.create_dataset('lon', data=lon)
    h5f.create_dataset('lat', data=lat)
    h5f.create_dataset('Z', data=Z)
    h5f.close()
    print(' -> Done!')

    
#############draw map with bathymetry########################################################


print('--- Now plot ---')
fig = plt.figure(figsize=(7,5))
ax = fig.add_subplot(111, projection=ccrs.Mercator())
ax.set_extent([lonLims[0], lonLims[1], latLims[0], latLims[1]], crs=ccrs.PlateCarree())
#ax.set_extent([-72, -41.5, 40.5, 55.1], crs=ccrs.PlateCarree())
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

# adjust the colorbar
levels = np.linspace(vmin, vmax, num_levels)
midp = np.mean(np.c_[levels[:-1], levels[1:]], axis=1)
vals = np.interp(midp, [vmin, midpoint, vmax], [0, 0.5, 1])
colors = colors(vals)
colors=np.concatenate([[colors[0,:]], colors, [colors[-1,:]]],0)
cmap, norm = from_levels_and_colors(levels, colors, extend=extent)

# plot data onto map
s = ax.scatter(lon_data,lat_data, c=data, s=20, lw=0.3, edgecolor='black', vmin=vmin, vmax=vmax, cmap=cmap, transform=ccrs.Geodetic(), zorder=20)
cax = plt.axes([0.83,0.125,0.03,0.756])
cb = plt.colorbar(s, cax=cax, extend=extent, ticks=ticks)
cb.ax.tick_params(labelsize=8)
cb.set_label(axis_label, fontsize=12, fontweight='normal')

# add text
if SEASON == 'summer':
    ax.text(-42, 57.8+.05, str(YEAR), horizontalalignment='right', color='black', zorder=200,
            fontsize=12, fontweight='bold', transform=ccrs.PlateCarree())
    ax.text(-42, 57.8, '('+SEASON+')', horizontalalignment='right', color='black', verticalalignment='top', zorder=200,
            fontsize=12, fontweight='bold', transform=ccrs.PlateCarree())    
else:
     ax.text(-42, 52+.05, str(YEAR), horizontalalignment='right', color='black', zorder=200,
            fontsize=12, fontweight='bold', transform=ccrs.PlateCarree())   
     ax.text(-42, 52, '('+SEASON+')', horizontalalignment='right', color='black', verticalalignment='top', zorder=200,
            fontsize=12, fontweight='bold', transform=ccrs.PlateCarree())   
        
fig_name = 'NL_OA_'+str(YEAR)+'_'+SEASON+'_'+FIG_VAR+'_'+DEPTH+'.png'

# Save figure
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)
#lt.close('all') 




# for 2020:
#montage AZMP_OA_2020_summer_OmegaA_surface.png AZMP_OA_2020_summer_pH_surface.png AZMP_OA_2020_summer_DO_perc_surface.png -tile 1x3 -geometry +10+10  -background white AZMP_OA_2020_summer_surface.png

# for 2021:
#montage AZMP_OA_2021_summer_OmegaA_surface.png AZMP_OA_2021_summer_pH_surface.png AZMP_OA_2021_summer_DO_perc_surface.png -tile 1x3 -geometry +10+10  -background white AZMP_OA_2021_summer_surface.png

#montage AZMP_OA_2020_summer_surface.png  AZMP_OA_2021_summer_surface.png -tile 2x1 -geometry +10+10  -background white AZMP_OA_2020-2021_summer_surface.png

#montage NL_OA_2019_summer_pH_bottom.png NL_OA_2019_summer_OmegaA_bottom.png NL_OA_2020_summer_pH_bottom.png NL_OA_2020_summer_OmegaA_bottom.png -tile 2x2 -geometry +10+10  -background white NL_OA_summer_2019-2020.png

#montage NL_OA_2019_summer_pH_bottom.png NL_OA_2019_fall_pH_bottom.png NL_OA_2020_summer_pH_bottom.png NL_OA_2020_fall_pH_bottom.png -tile 2x2 -geometry +10+10  -background white NL_pH_summer-fall_2019-2020.png

# Montage for same year summer and fall:
# montage NL_OA_2019_summer_pH_bottom.png NL_OA_2019_fall_pH_bottom.png NL_OA_2019_summer_OmegaA_bottom.png NL_OA_2019_fall_OmegaA_bottom.png -tile 2x2 -geometry +10+10  -background white NL_OA_2019.png

# montage NL_OA_2020_summer_pH_bottom.png NL_OA_2020_fall_pH_bottom.png NL_OA_2020_summer_OmegaA_bottom.png NL_OA_2020_fall_OmegaA_bottom.png -tile 2x2 -geometry +10+10  -background white NL_OA_2020.png


# montage NL_OA_2019_summer_DO_perc_bottom.png NL_OA_2019_fall_DO_perc_bottom.png NL_OA_2020_summer_DO_perc_bottom.png NL_OA_2020_fall_DO_perc_bottom.png-tile 2x2 -geometry +10+10  -background white NL_DO_perc_2019-2020.png


