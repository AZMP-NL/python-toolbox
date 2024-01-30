'''
This is a modified copy of cs_section_plot.py
which generates figures for the manuscript (Figure 9&10)
 
Refers to the script above for other usage.

Frederic.Cyr@dfo-mpo.gc.ca - July 2023
'''

import os
import pandas as pd
#matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.colors import from_levels_and_colors
import numpy as np
from math import radians, cos, sin, asin, sqrt
import azmp_sections_tools as azst
import cmocean.cm as cmo
#import cc_variable_list_NL as vl
# For NL work:
import cc_variable_list2 as vl


## ----- Parameters to edit ---- ##
# Variable you want to plot
#VAR = 'Omega_Aragonite_(unitless)'
#VAR = 'pH_Total_(total_scale)'
VAR = 'Oxygen_Saturation_(%)'
#VAR = 'Dissolved_Oxygen_(mL/L)'
#VAR = 'pCO2_(uatm)'
#VAR = 'Phosphate_Concentration_(mmol/m3)'
#VAR = 'Temperature_(degC)'
#VAR = 'pH_tot'
#VAR = 'Omega_Aragonite_(--)'
#VAR = 'Salinity_(psu)'
#VAR = 'Total_Alkalinity_(umol/kg)'
#VAR = 'Inorganic_Carbon_(umol/kg)'
REGION = 'NL'
SECTION = 'SESPB'
SEASON = 'fall'
YEAR = 2015# plotting parameters
CLIM = True # Clim or no clim?
#v = 10
#v_anom=10
CMAP = cmo.cm.seismic
ZMAX = 500
# File to load
my_file = '/home/cyrf0006/github/AZMP-NL/datasets/carbonates/AZMP_carbon_data.csv'
# For colorbar:
vv = vl.variable_parameters(VAR)
num_levels = vv[0]
vmin = vv[1]
vmax = vv[2]
midpoint = vv[3]
colors = vv[4]
ticks = vv[5]
axis_label = vv[6]
extent = vv[7]
v_anom = np.linspace(vv[8], vv[9], vv[10])

# Text for figure:
section_EN = ['BI', 'MB', 'SI', 'WB', 'BB', 'S27', 'FC', 'SEGB', 'SWSPB']
season_EN = ['spring', 'summer', 'fall']
VAR_text = VAR.split('(')[0][0:-1] 
if VAR == 'Oxygen_Saturation_(%)':
    title = 'Oxygen Saturation for section ' + SECTION + ' - ' + SEASON
elif VAR == 'pCO2_(uatm)':
    title = r'$p$$\rm CO_2$ for section ' + SECTION + ' - ' + SEASON
elif VAR == 'pH_Total_(total_scale)':
    title = 'pH Total for section ' + SECTION + ' - ' + SEASON + ' ' + str(YEAR)
elif VAR == 'Omega_Aragonite_(unitless)':
    title = 'Aragonite saturation state for section ' + SECTION + ' - ' + SEASON + ' ' + str(YEAR)
else:
    title=VAR_text
    

# adjust the colorbar
levels = np.linspace(vmin, vmax, num_levels)
midp = np.mean(np.c_[levels[:-1], levels[1:]], axis=1)
vals = np.interp(midp, [vmin, midpoint, vmax], [0, 0.5, 1])
colors = colors(vals)
colors=np.concatenate([[colors[0,:]], colors, [colors[-1,:]]],0)
cmap, norm = from_levels_and_colors(levels, colors, extend=extent)

# Get the data
df = pd.read_csv(my_file)
# set index
df.index = pd.to_datetime(df.Timestamp)
df.drop(columns='Timestamp', inplace=True)

# If sections are trans-region, this is a bit more complicated:
if SECTION == 'xGSL':
    df = df.loc[(df.Station_Name.str.contains('TIDM7')) | (df.Station_Name.str.contains('TCEN')) | (df.Station_Name.str.contains('CH6')) | (df.Station_Name.str.contains('TBB3')) | (df.Station_Name.str.contains('CH2')) | (df.Station_Name.str.contains('CMO1/CH1')) | (df.Station_Name.str.contains('IF12'))]

elif SECTION == 'LC':
    df = df.loc[(df.Station_Name.str.contains('14ML')) | (df.Station_Name.str.contains('CMO3/CH12')) | (df.Station_Name.str.contains('IF37')) | (df.Station_Name.str.contains('TSI3')) | (df.Station_Name.str.contains('TASO3')) | (df.Station_Name.str.contains('IF34')) | (df.Station_Name.str.contains('IF32')) | (df.Station_Name.str.contains('TCEN3')) | (df.Station_Name.str.contains('IF27')) | (df.Station_Name.str.contains('IF28')) | (df.Station_Name.str.contains('CSL_04')) | (df.Station_Name.str.contains('STAB'))]

else: # Extract Region + Section
    df = df[df.Region==REGION]
    df = df.loc[(df.Station_Name.str.contains(SECTION))]

# Extract season
#df = df[df.index.year == YEAR]
if SEASON == 'spring':
    df = df[df.index.month <= 5]
elif SEASON == 'summer':
    df = df[(df.index.month>=6) & (df.index.month<=8)]
elif SEASON == 'fall':
    df = df[df.index.month >= 9]

# Build climatology
df_list = []
for i in df.index.year.unique():
    df_tmp = df[df.index.year == i]
    # Extract variable
    df_tmp = df_tmp[['Depth_(dbar)', 'Station_Name', VAR]]
    df_tmp = df_tmp.pivot(index='Depth_(dbar)', columns='Station_Name') #<--- this is cool!
    df_tmp = df_tmp[VAR]
    # So I re-define a constant vertical axis every 5m.
    depth_range = np.arange(2.5, 2000, 5) # range to look for data
    reg_depth = (depth_range[1:] + depth_range[:-1]) / 2 # mid point of the range
    df_tmp = df_tmp.groupby(pd.cut(df_tmp.index, depth_range)).mean() # <--- This is cool!
    df_tmp.index = reg_depth # replace range by mean depth
    # interpolate vertically and horisontally where possible
    df_tmp.interpolate(axis=0, limit_area='inside', inplace=True)
    df_tmp.interpolate(axis=1, limit_area='inside', inplace=True)
    # Drop depth where there is no values associated to a certain depth (because my axis is every 5m...)
    #df_tmp = df_tmp.dropna(how='all')
    df_list.append(df_tmp)
    #del df_tmp
    

# Create multi-index
df_all = pd.concat(df_list, keys=df.index.year.unique(), sort=True)

# extract year current year
df_year = df_all.xs((YEAR), level=('Timestamp'))

# compute climatology
df_clim = df_all.groupby(level=1).apply(lambda x: x.mean())

# vertically fill NaNs
df_year.interpolate(axis=0, limit_area='inside', inplace=True)
df_clim.interpolate(axis=0, limit_area='inside', inplace=True)
# horizontally fill NaNs
df_year.interpolate(axis=1, limit_area='inside', inplace=True)
df_clim.interpolate(axis=1, limit_area='inside', inplace=True)

# calculate anomaly
df_anom = df_year-df_clim

## ---- Load station lat/lon ---- ##
df_stn = pd.read_excel('/home/cyrf0006/github/AZMP-NL/data/STANDARD_SECTIONS.xlsx')
df_stn = df_stn.drop(['LONG'], axis=1)
df_stn = df_stn.rename(columns={'LONG.1': 'LON'})
df_stn = df_stn.dropna()
if SECTION == 'LC':
    df_stn = df_stn.loc[df_stn.SECTION.str.contains('LAURENTIAN CHANNEL')]
else:        
    df_stn = df_stn[df_stn.STATION.str.contains(SECTION)]
df_stn = df_stn.drop(['SECTION'], axis=1)
df_stn = df_stn.reset_index(drop=True)
# remove hyphen in stn names (normally should keep it everywhere, but simpler here)
#df_stn['STATION'] = df_stn.STATION.str.replace('-', '')

# Sort dataframe according to stations
if (SECTION == 'LC') | (SECTION == 'xGSL'):
    df_year = df_year[df_stn.STATION.to_list()]
    df_clim = df_clim[df_stn.STATION.to_list()]
    df_anom = df_anom[df_stn.STATION.to_list()]

## ---- Compute distance vector ---- ##
distance = np.full((df_clim.keys().shape), np.nan)
lat0 = df_stn[df_stn.index==0]['LAT']
lon0 = df_stn[df_stn.index==0]['LON']
for i, stn in enumerate(df_clim.keys().values):
    #stn = stn.replace('_','-') # replace in case underscore is used
    lat_stn = df_stn[df_stn.STATION==str(stn)]['LAT']
    lon_stn = df_stn[df_stn.STATION==str(stn)]['LON']      
    distance[i] = azst.haversine(lon0, lat0, lon_stn, lat_stn)
XLIM =distance.max()

## ---- Retrieve bathymetry using function ---- ##
if REGION == 'NL':
    bathymetry = azst.section_bathymetry(SECTION)
else:
    bathymetry = []

## ---- plot Figure ---- ##
#XLIM = df_section_itp.index[-1][1]
fig = plt.figure()
# ax1
ax = plt.subplot2grid((1, 1), (0, 0))
if CLIM:
    c = plt.contourf(distance, df_year.index, df_clim, levels, cmap=cmap, extend='both')
else:
    c = plt.contourf(distance, df_year.index, df_year, levels, cmap=cmap, extend='both')

for i in distance:
    plt.plot([i, i], [0, ZMAX], '--k', alpha=.5)
ax.set_ylim([0, ZMAX])
ax.set_xlim([0, XLIM])
#plt.clabel(c_sig1, inline=1, fontsize=10, colors='gray', fmt='%1.1f')
ax.set_ylabel('Depth (m)')
ax.set_xlabel('Distance (km)')
ax.invert_yaxis()
if REGION == 'NL':
    Bgon = plt.Polygon(bathymetry,color=np.multiply([1,.9333,.6667],.4), alpha=1, zorder=10)
    ax.add_patch(Bgon)
cb = plt.colorbar(c)
cb.ax.tick_params(labelsize=8)
cb.set_label(axis_label, fontsize=12, fontweight='normal')
#ax.xaxis.label.set_visible(False)
#ax.tick_params(labelbottom='off')
ax.set_title(title)


fig.set_size_inches(w=8,h=4)
if CLIM:
    fig_name = 'ms_sections_' + VAR_text + '_' + SECTION + '_' + SEASON + '_clim.png' 
else:
    fig_name = 'ms_sections_' + VAR_text + '_' + SECTION + '_' + SEASON + '_' + str(YEAR) + '.png' 
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)



# montage ms_sections_Temperature_BB_fall_clim.png  ms_sections_Temperature_LC_fall_clim.png -tile 1x2 -geometry +10+10  -background white  Figure09a.png

# montage  ms_sections_Oxygen_Saturation_SESPB_fall_clim.png ms_sections_pCO2_SESPB_fall_clim.png -tile 1x2 -geometry +10+10  -background white  Figure09.png
