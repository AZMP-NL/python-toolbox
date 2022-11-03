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
import cc_variable_list2 as vl

YEAR = 2019
SEASON = 'fall'
#VARIABLE = 'Omega_Aragonite_(--)'
#VARIABLE = 'pH_tot'
VARIABLE = 'Oxygen_Saturation_(%)'
#VARIABLE = 'Dissolved_Oxygen_(mL/L)'
DEPTH='bottom'
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


# Read the entire AZMP dataset
df = pd.read_csv('/home/cyrf0006/github/AZMP-NL/datasets/carbonates/AZMP_carbon_data.csv', delimiter=',')

# Set index
df.set_index('Timestamp', inplace=True)
df.index = pd.to_datetime(df.index)

# Only year in revieww
df = df[df.index.year==YEAR]

# Season
if SEASON == 'spring':
    df = df[(df.index.month>=3) & (df.index.month<=6)]
elif SEASON == 'summer':
    df = df[~df.Region.str.contains('MAR')]
    df = df.assign(x=df.index.strftime('%m-%d')).query("'07-01' <= x <= '10-14'").drop('x',1)
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

df.dropna(subset=[VARIABLE, 'pH_tot', 'Omega_Aragonite_(--)', 'Station_Name', 'Depth_(dbar)'], axis=0, inplace=True)
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
    if YEAR == 2019: # Some stations to be flagged (deep bottle not bottom)
        df.drop(df[df.Station_Name=='GULD_04'].index, inplace=True)
        df.drop(df[df.Station_Name=='HL_07'].index, inplace=True)
        df.drop(df[df.Station_Name=='LL_08'].index, inplace=True)

if DEPTH == 'range':
    df = df.loc[(df['Depth_(dbar)'] >=135) & (df.depth <=165)]
    df = df.loc[df.groupby('Station_Name')['Depth_(dbar)'].idxmax()] #group by station then pull "min or max depth"


############create bathymetry################
dataFile = '/home/cyrf0006/data/GEBCO/GEBCO_2014_1D.nc'
###map boundaries###
lonLims = [-72, -41.5]  
latLims = [40.5, 58.4]
#latLims = [40.5, 55.1]

##### this is the dataset that will be plotted ##########
lat_data = np.array(df['Latitude_(degNorth)'])
lon_data = np.array(df['Longitude_(degEast)'])
data = np.array(df[VARIABLE])
lon_data = lon_data[~np.isnan(data)]
lat_data = lat_data[~np.isnan(data)]
data = data[~np.isnan(data)]

print('Load and grid bathymetry')
# h5 file
h5_outputfile = 'cc_bathymetry.h5'
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
ax.set_extent([-72, -41.5, 40.5, 58.4], crs=ccrs.PlateCarree())
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
if VARIABLE == 'Omega_Aragonite_(--)':
    ax.text(-69, 56, str(YEAR), horizontalalignment='left', color='white',
            fontsize=12, fontweight='bold', transform=ccrs.PlateCarree())

# Save figure
fig.savefig('AZMP_OA_'+str(YEAR)+'_'+SEASON+'_'+VARIABLE+'_'+DEPTH+'.png', format='png', dpi=300, bbox_inches='tight')
#lt.close('all') 

keyboard


# Montage in Linux
#gm convert -trim AZMP_OA_2019_spring_Omega_A_bottom.png  AZMP_OA_2019_spring_Omega_A_bottom.png
#gm convert -trim AZMP_OA_2019_spring_satO2_perc_bottom.png AZMP_OA_2019_spring_satO2_perc_bottom.png
#gm convert -trim AZMP_OA_2019_spring_pH_bottom.png AZMP_OA_2019_spring_pH_bottom.png

# for 2020:
#montage AZMP_OA_2019_fall_Omega_Aragonite_\(--\)_bottom.png AZMP_OA_2019_fall_pH_tot_bottom.png AZMP_OA_2019_fall_Oxygen_Saturation_\(%\)_bottom.png -tile 1x3 -geometry +10+10  -background white AZMP_OA_2019_fall.png

#montage AZMP_OA_2020_fall_Omega_Aragonite_\(--\)_bottom.png AZMP_OA_2020_fall_pH_tot_bottom.png AZMP_OA_2020_fall_Oxygen_Saturation_\(%\)_bottom.png -tile 1x3 -geometry +10+10  -background white AZMP_OA_2020_fall.png

#montage AZMP_OA_2019_fall.png  AZMP_OA_2020_fall.png -tile 2x1 -geometry +10+10  -background white AZMP_OA_2019-2020_fall.png


## ## ---- Test with Bokeh plot ---- ##
from bokeh.plotting import figure, output_file, show
from bokeh.models import ColumnDataSource
from bokeh.models.tools import HoverTool
from bokeh.models import DatetimeTickFormatter
from bokeh.transform import factor_cmap, factor_mark
# For colorbar
from bokeh.models import ColorBar
from bokeh.palettes import Spectral8, Viridis7, Spectral7
from bokeh.transform import linear_cmap

## output_file('OA_2019.html')

## ## def datetime(x):
## ##     return np.array(x, dtype=np.datetime64)

## source = ColumnDataSource(df)
## #p = figure(x_axis_type="datetime", plot_width=900, plot_height=600)
## p = figure(plot_width=900, plot_height=600)
## p.scatter('Omega_A', 'pH', source=source, legend='Region', size=8,
##           color=factor_cmap('Region', 'Category10_3', df.Region.unique()))
## p.title.text = 'Carbonate parameters in the Atlantic Zone - Spring 2019, bottom'
## p.xaxis.axis_label = 'Saturation State (Aragonite)'
## p.yaxis.axis_label = 'pH'
## p.add_tools(HoverTool(
##     tooltips=[
##         ('pH', '@pH'),
##         ('Saturation Aragonite', '@Omega_A'),
##         ('Saturation Calcite', '@Omega_C'),
##         ('Total Alkalinity', '@TA'),
##         ('Inorg. Carbon', '@TIC'),
##         ('O2 sat (%)', '@satO2_perc'),
##         ('Region', '@Region'),
##         ('Station', '@Station_Name'),
##         ('Depth (m)', '@depth'),
##         ('Temperature', '@temperature'),
##         ('Salinity', '@salinity'),
##         ('Latitude', '@latitude'),
##         ('Longitude', '@longitude')      
##     ],

##     ## formatters={
##     ##     'time'      : 'datetime', # use 'datetime' formatter for 'date' field
##     ## },

##     # display a tooltip when cursor is vertically ('vline') aligned, horizontally ('hline') or above ('mouse') a glyph
##     mode='mouse'
## ))

## show(p)


# Rename columns for Bokeh
df = df.rename(columns={'Latitude_(degNorth)' : 'Latitude'})
df = df.rename(columns={'Longitude_(degEast)' : 'Longitude'})
df = df.rename(columns={'Omega_Aragonite_(--)' : 'Omega_Aragonite'})
df = df.rename(columns={'Omega_Calcite_(--)' : 'Omega_Calcite'})
df = df.rename(columns={'Total_Alkalinity_(umol/kg)' : 'Total_Alkalinity'})
df = df.rename(columns={'Inorganic_Carbon_(umol/kg)' : 'Inorganic_Carbon'})
df = df.rename(columns={'Oxygen_Saturation_(%)' : 'Oxygen_Saturation'})
df = df.rename(columns={'Dissolved_Oxygen_(mL/L)' : 'Dissolved_Oxygen'})
df = df.rename(columns={'Depth_(dbar)' : 'Depth'})
df = df.rename(columns={'Temperature_(degC)' : 'Temperature'})
df = df.rename(columns={'Salinity_(psu)' : 'Salinity'})

## # Another one in map-like display
# pH
Viridis7.reverse()
Spectral7.reverse()
mapper = linear_cmap(field_name='pH_tot', palette=Spectral7 ,low=7.5 ,high=8.2)
output_file('OA_2019_fall_pH_map.html')
source = ColumnDataSource(df)
p = figure(plot_width=900, plot_height=600)
p.scatter('Longitude', 'Latitude', source=source, size=8,
          line_color=mapper,color=mapper)
p.title.text = 'Carbonate parameters in the Atlantic Zone - Fall 2019, bottom'
p.xaxis.axis_label = 'Longitude'
p.yaxis.axis_label = 'Latitude'
p.add_tools(HoverTool(
    tooltips=[
        ('pH', '@pH_tot'),
        ('Saturation Aragonite', '@Omega_Aragonite'),
        ('Saturation Calcite', '@Omega_Calcite'),
        ('Total Alkalinity', '@Total_Alkalinity'),
        ('Inorg. Carbon', '@Inorganic_Carbon'),
        ('O2 sat (%)', '@Oxygen_Saturation'),
        ('O2 dissolved (ml/l)', '@Dissolved_Oxygen'),
        ('Region', '@Region'),
        ('Station', '@Station_Name'),
        ('Depth (m)', '@Depth'),
        ('Temperature', '@Temperature'),
        ('Salinity', '@Salinity'),
        ('Latitude', '@Latitude'),
        ('Longitude', '@Longitude')      
    ],

    mode='mouse'
))
color_bar = ColorBar(color_mapper=mapper['transform'], width=8,  location=(0,0),title="pH")
p.add_layout(color_bar, 'right')
show(p)

# 2. Omega_A
Spectral8.reverse()
mapper = linear_cmap(field_name='Omega_Aragonite', palette=Spectral8 ,low=0 ,high=2)
output_file('OA_2019_fall_omegaA_map.html')
source = ColumnDataSource(df)
p = figure(plot_width=900, plot_height=600)
p.scatter('Longitude', 'Latitude', source=source, size=8,
          line_color=mapper,color=mapper)
p.title.text = 'Carbonate parameters in the Atlantic Zone - Fall 2019, bottom'
p.xaxis.axis_label = 'Longitude'
p.yaxis.axis_label = 'Latitude'
p.add_tools(HoverTool(
    tooltips=[
        ('pH', '@pH_tot'),
        ('Saturation Aragonite', '@Omega_Aragonite'),
        ('Saturation Calcite', '@Omega_Calcite'),
        ('Total Alkalinity', '@Total_Alkalinity'),
        ('Inorg. Carbon', '@Inorganic_Carbon'),
        ('O2 sat (%)', '@Oxygen_Saturation'),
        ('O2 dissolved (ml/l)', '@Dissolved_Oxygen'),
        ('Region', '@Region'),
        ('Station', '@Station_Name'),
        ('Depth (m)', '@Depth'),
        ('Temperature', '@Temperature'),
        ('Salinity', '@Salinity'),
        ('Latitude', '@Latitude'),
        ('Longitude', '@Longitude')
    ],

    mode='mouse'
))
color_bar = ColorBar(color_mapper=mapper['transform'], width=8,  location=(0,0),title="OmegaA")
p.add_layout(color_bar, 'right')
show(p)


# 2. Omega_C
mapper = linear_cmap(field_name='Omega_Calcite', palette=Spectral8 ,low=0 ,high=2)
output_file('OA_2019_fall_omegaC_map.html')
source = ColumnDataSource(df)
p = figure(plot_width=900, plot_height=600)
p.scatter('Longitude', 'Latitude', source=source, size=8,
          line_color=mapper,color=mapper)
p.title.text = 'Carbonate parameters in the Atlantic Zone - Fall 2019, bottom'
p.xaxis.axis_label = 'Longitude'
p.yaxis.axis_label = 'Latitude'
p.add_tools(HoverTool(
    tooltips=[
        ('pH', '@pH_tot'),
        ('Saturation Aragonite', '@Omega_Aragonite'),
        ('Saturation Calcite', '@Omega_Calcite'),
        ('Total Alkalinity', '@Total_Alkalinity'),
        ('Inorg. Carbon', '@Inorganic_Carbon'),
        ('O2 sat (%)', '@Oxygen_Saturation'),
        ('O2 dissolved (ml/l)', '@Dissolved_Oxygen'),
        ('Region', '@Region'),
        ('Station', '@Station_Name'),
        ('Depth (m)', '@Depth'),
        ('Temperature', '@Temperature'),
        ('Salinity', '@Salinity'),
        ('Latitude', '@Latitude'),
        ('Longitude', '@Longitude')
    ],

    mode='mouse'
))
color_bar = ColorBar(color_mapper=mapper['transform'], width=8,  location=(0,0),title="OmegaC")
p.add_layout(color_bar, 'right')
show(p)


## # 3. T-S
## mapper = linear_cmap(field_name='depth', palette=cmo.cm.deep ,low=0 ,high=300)
## output_file('OA_T-S.html')
## source = ColumnDataSource(df)
## p = figure(plot_width=900, plot_height=600)
## p.scatter('longitude', 'latitude', source=source, size=8,
##           line_color=mapper,color=mapper)
## p.title.text = 'Carbonate parameters in the Atlantic Zone - 2014-2019'
## p.xaxis.axis_label = 'Longitude'
## p.yaxis.axis_label = 'Latitude'
## p.add_tools(HoverTool(
##     tooltips=[
##         ('Temperature', '@temperature'),
##         ('Salinity', '@salinity'),
##         ('Depth (m)', '@depth'),
##         ('pH', '@pH'),
##         ('Saturation Aragonite', '@Omega_A'),
##         ('Saturation Calcite', '@Omega_C'),
##         ('Total Alkalinity', '@TA'),
##         ('Inorg. Carbon', '@TIC'),
##         ('O2 sat (%)', '@satO2_perc'),
##         ('O2 dissolved (ml/l)', '@O2_ml_l'),
##         ('Region', '@Region'),
##         ('Station', '@Station_Name'),
##         ('Timestamp', '@timestamp'),
##         ('Latitude', '@latitude'),
##         ('Longitude', '@longitude')      
##     ],

##     mode='mouse'
## ))
## color_bar = ColorBar(color_mapper=mapper['transform'], width=8,  location=(0,0),title="T-S diagram")
## p.add_layout(color_bar, 'right')
## show(p)


## #Seasonal analysis
## df = df[df.depth<=50]
## df_NL = df[df.Region == 'NL']
## df_IML = df[df.Region == 'GSL']
## df_MAR = df[df.Region == 'SS']

## NL = df_NL.groupby([(df_NL.index.year),(df_NL.index.month)]).mean()
## NL = NL.mean(level=1)
## IML = df_IML.groupby([(df_IML.index.year),(df_IML.index.month)]).mean()
## IML = IML.mean(level=1)
## MAR = df_MAR.groupby([(df_MAR.index.year),(df_MAR.index.month)]).mean()
## MAR = MAR.mean(level=1)

## plt.plot(NL.sort_index().Omega_A)
## plt.plot(IML.sort_index().Omega_A)
## plt.plot(MAR.sort_index().Omega_A)
## plt.legend(['NL', 'GSL', 'SS'])

