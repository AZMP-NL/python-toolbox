# -*- coding: utf-8 -*-
"""

To convert cabonate measurements from all regions to CO2sys

Frederic.Cyr@dfo-mpo.gc.ca
April 2020

CO2sys usage:
CO2dict = CO2SYS(PAR1, PAR2, PAR1TYPE, PAR2TYPE, SAL, TEMPIN, TEMPOUT, PRESIN, PRESOUT,
    SI, PO4, pHSCALEIN, K1K2CONSTANTS, KSO4CONSTANTS, NH3=0.0, H2S=0.0, KFCONSTANT=1)
    
where: (* From https://github.com/mvdh7/PyCO2SYS)

Required inputs

The required inputs are identical to the MATLAB version:

    PAR1 - frst known carbonate system parameter value.
    PAR2 - second known carbonate system parameter value.
    PAR1TYPE - integer identifying which parameters PAR1 are.
    PAR2TYPE - integer identifying which parameters PAR2 are.

The possible known carbonate system parameters are 1: total alkalinity in μmol·kg−1, 2: dissolved inorganic carbon in μmol·kg−1, 3: pH (dimensionless), 4: dissolved CO2 partial pressure in μatm, 5: dissolved CO2 fugacity in μatm, and 6: carbonate ion concentration in μmol·kg−1.

Here and throughout the inputs and outputs, "kg−1" refers to the total mass of seawater (solvent + solutes), not just the mass of H2O.

    SAL - practical salinity.
    TEMPIN - temperature of input carbonate system parameters.
    TEMPOUT - temperature at which to calculate outputs.
    PRESIN - pressure of input carbonate system parameters.
    PRESOUT - pressure at which to calculate outputs.

All temperatures are in °C and pressures are in dbar. Pressure is within the water column as typically measured by a CTD sensor, i.e. not including atmospheric pressure. The 'input' conditions could represent conditions in the laboratory during a measurement, while the 'output' conditions could be those observed in situ during sample collection.

    SI - total silicate concentration.
    PO4 - total phosphate concentration.

Nutrient concentrations are all in μmol·kg−1.

    pHSCALEIN - pH scale(s) that pH values in PAR1 or PAR2 are on.

The options are 1: Total scale, 2: Seawater scale, 3: Free scale, and 4: NBS scale, as defined by ZW01.

    K1K2CONSTANTS - which set of constants to use for carbonic acid dissociation.

The options are integers from 1 to 15 inclusive. From the original MATLAB documentation:

%   1 = Roy, 1993                                         T:    0-45  S:  5-45. Total scale. Artificial seawater.
%   2 = Goyet & Poisson                                   T:   -1-40  S: 10-50. Seaw. scale. Artificial seawater.
%   3 = HANSSON              refit BY DICKSON AND MILLERO T:    2-35  S: 20-40. Seaw. scale. Artificial seawater.
%   4 = MEHRBACH             refit BY DICKSON AND MILLERO T:    2-35  S: 20-40. Seaw. scale. Artificial seawater.
%   5 = HANSSON and MEHRBACH refit BY DICKSON AND MILLERO T:    2-35  S: 20-40. Seaw. scale. Artificial seawater.
%   6 = GEOSECS (i.e., original Mehrbach)                 T:    2-35  S: 19-43. NBS scale.   Real seawater.
%   7 = Peng    (i.e., original Mehrbach but without XXX) T:    2-35  S: 19-43. NBS scale.   Real seawater.
%   8 = Millero, 1979, FOR PURE WATER ONLY (i.e., Sal=0)  T:    0-50  S:     0.
%   9 = Cai and Wang, 1998                                T:    2-35  S:  0-49. NBS scale.   Real and artificial seawater.
%  10 = Lueker et al, 2000                                T:    2-35  S: 19-43. Total scale. Real seawater.
%  11 = Mojica Prieto and Millero, 2002.                  T:    0-45  S:  5-42. Seaw. scale. Real seawater
%  12 = Millero et al, 2002                               T: -1.6-35  S: 34-37. Seaw. scale. Field measurements.
%  13 = Millero et al, 2006                               T:    0-50  S:  1-50. Seaw. scale. Real seawater.
%  14 = Millero        2010                               T:    0-50  S:  1-50. Seaw. scale. Real seawater.
%  15 = Waters, Millero, & Woosley 2014                   T:    0-50  S:  1-50. Seaw. scale. Real seawater.

    KSO4CONSTANTS - which sets of constants to use for bisulfate dissociation and borate:chlorinity ratio.

The options are integers from 1 to 4 inclusive. From the original MATLAB documentation:

%  1 = KSO4 of Dickson 1990a   & TB of Uppstrom 1974  (PREFERRED)
%  2 = KSO4 of Khoo et al 1977 & TB of Uppstrom 1974
%  3 = KSO4 of Dickson 1990a   & TB of Lee et al. 2010
%  4 = KSO4 of Khoo et al 1977 & TB of Lee et al. 2010
    
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
from PyCO2SYS import CO2SYS 
#import seawater as swx
import gsw

#my_region=''
my_year=2019
my_season='spring'
my_depth='bottom'
my_variable='pH'

#### Load Carbonate data into a Pandas DataFrame ######
if my_year != 2019:
    from PyCO2SYS.meta import version
    print('Use PyCO2SYS v{}'.format(version))
    #### Load Carbonate data into a Pandas DataFrame ######
    df = pd.read_excel('/home/cyrf0006/github/AZMP-NL/datasets/carbonates/AZMP_OA_CO2stats.xlsx')
    df = df.set_index('timestamp', drop=False)

else:
    # 1. NL
    df_NL = pd.read_excel('/home/cyrf0006/github/AZMP-NL/datasets/carbonates/AZMP_NL_plot.xlsx')
    df_NL = df_NL.set_index('timestamp', drop=False)
    # CO2sys
    CO2dict_NL = CO2SYS(df_NL.TA, df_NL.TIC, 1, 2, df_NL.salinity, 20, df_NL.temperature, 0, df_NL.depth, 0, 0, 1, 1, 1)
    co2sys_NL = pd.DataFrame.from_dict(CO2dict_NL)
    co2sys_NL.index = df_NL.index
    # try to by-pass Gary's calculation
    df_NL['pH'] = co2sys_NL['pHoutTOTAL']
    df_NL['Omega_C'] = co2sys_NL['OmegaCAout']
    df_NL['Omega_A'] = co2sys_NL['OmegaARout']
    df_NL['pH'] = co2sys_NL['pHoutTOTAL']
    # Calculate in situ o2
    SA = gsw.SA_from_SP(df_NL.salinity, df_NL.depth , df_NL.latitude, df_NL.longitude)
    CT = gsw.CT_from_t(SA, df_NL.temperature, df_NL.depth)
    df_NL['satO2'] = gsw.O2sol(SA, CT, df_NL.depth , df_NL.latitude, df_NL.longitude)
    df_NL['O2_ml_l'] = df_NL['satO2_perc']/100*df_NL['satO2']/44.6
    del SA, CT
    
    # 2. Gulf
    df_IML = pd.read_excel('/home/cyrf0006/github/AZMP-NL/datasets/carbonates/IMLJune2019(Gibb).xlsx', header=1, parse_dates = {'timestamp' : [6, 7]})
    df_IML = df_IML.drop(0)
    df_IML.index = pd.to_datetime(df_IML.timestamp)
    df_IML.drop(columns='timestamp', inplace=True)
    df_IML = df_IML.rename(columns={'Unnamed: 64' : 'pH'})
    df_IML = df_IML.rename(columns={'Unnamed: 65' : 'Omega_C'})
    df_IML = df_IML.rename(columns={'Unnamed: 66' : 'Omega_A'})
    df_IML = df_IML.rename(columns={'Unnamed: 67' : 'strate'})
    # Drop empty values
    df_IML.dropna(subset=['pH  labo'], inplace=True)
    # CO2sys
    CO2dict_IML = CO2SYS(df_IML.ALKW_01, df_IML.pH, 1, 3, df_IML.PSAL, 20, df_IML.TE90, 0, df_IML.PRES, 0, 0, 1, 1, 1)
    co2sys_IML = pd.DataFrame.from_dict(CO2dict_IML)
    co2sys_IML.index = df_IML.index
    #df_IML['Omega_C'] = co2sys_IML['OmegaCAout'] # here we use the one derive by Michel
    #df_IML['Omega_A'] = co2sys_IML['OmegaARout']
    df_IML['TIC'] = co2sys_IML['TCO2']
    # Calculate in situ o2
    df_IML = df_IML.astype({'OXY_02':'float'})
    SA = gsw.SA_from_SP(df_IML.PSAL, df_IML.PRES , df_IML.Latitude, df_IML.Longitude)
    CT = gsw.CT_from_t(SA, df_IML.TE90, df_IML.PRES)
    df_IML['satO2'] = gsw.O2sol(SA, CT, df_IML.PRES , df_IML.Latitude, df_IML.Longitude)
    df_IML['satO2_perc'] = df_IML['OXY_02']*44.66/df_IML['satO2']*100  
    del SA, CT


    # 3. Maritimes
    df_MAR = pd.read_excel('/home/cyrf0006/github/AZMP-NL/datasets/carbonates/2019_Spring_AZMP_Mar_TIC-TA.xlsx', parse_dates = {'timestamp' : [7, 8]})
    df_MAR.index = pd.to_datetime(df_MAR.timestamp)
    df_MAR.drop(columns='timestamp', inplace=True)
    # Drop empty values
    df_MAR.dropna(subset=['TA (umol/kg)'], inplace=True)
    # Flagg data
    #df_MAR.dropna()
    #df_MAR.drop(df_MAR.loc[df_MAR['TIC flag']>=3].index, inplace=True)
    # CO2sys
    CO2dict_MAR = CO2SYS(df_MAR['TA (umol/kg)'], df_MAR['TIC (umol/kg)'], 1, 2, df_MAR.Sal00, 20, df_MAR.T068C, 0, df_MAR.PrDM, 0, 0, 1, 1, 1)
    co2sys_MAR = pd.DataFrame.from_dict(CO2dict_MAR)
    co2sys_MAR.index = df_MAR.index
    df_MAR['pH'] = co2sys_MAR['pHoutTOTAL']
    df_MAR['Omega_C'] = co2sys_MAR['OmegaCAout']
    df_MAR['Omega_A'] = co2sys_MAR['OmegaARout']

    # -- add oxygen -- #
    df_o2 = pd.read_csv('WinklerOxygen_COR2019001_Final.dat', header=9, sep=',')
    # rename sample no.
    sampleNo = df_o2.Sample
    sampleNo = sampleNo.map(lambda x: x.replace('_1', ''))
    sampleNo = sampleNo.map(lambda x: x.replace('_2', ''))
    df_o2.Sample = sampleNo
    df_o2 = df_o2.groupby('Sample').mean() 
    df_o2 = df_o2['O2_Concentration(ml/l)']
    df_o2 = df_o2.reset_index()
    df_o2 = df_o2.rename(columns={'Sample' : 'sample_id'})
    df_o2 = df_o2.astype({'sample_id':'int'})
    # merge O2 into carbonate data
    MAR_index = df_MAR.index
    df_MAR = df_MAR.merge(df_o2, on='sample_id', how='left')
    df_MAR.index = MAR_index 
    # Calculate in situ o2
    SA = gsw.SA_from_SP(df_MAR.Sal00, df_MAR.PrDM , df_MAR.latitude, df_MAR.longitude)
    CT = gsw.CT_from_t(SA, df_MAR.T068C, df_MAR.PrDM)
    df_MAR['satO2'] = gsw.O2sol(SA, CT, df_MAR.PrDM , df_MAR.latitude, df_MAR.longitude)
    df_MAR['satO2_perc'] = df_MAR['O2_Concentration(ml/l)']*44.66/df_MAR['satO2']*100  # to umol/l
    del SA, CT
 
    # 4. merge relevant parameters
    df_IML = df_IML.rename(columns={'Station' : 'StationID'})
    df_IML = df_IML.rename(columns={'Latitude' : 'latitude'})
    df_IML = df_IML.rename(columns={'Longitude' : 'longitude'})
    df_IML = df_IML.rename(columns={'PRES' : 'depth'})
    df_IML = df_IML.rename(columns={'TE90' : 'temperature'})
    df_IML = df_IML.rename(columns={'PSAL' : 'salinity'})
    df_IML = df_IML.rename(columns={'ALKW_01' : 'TA'})
    df_IML = df_IML.rename(columns={'OXY_02' : 'O2_ml_l'})
    df_MAR = df_MAR.rename(columns={'Station' : 'StationID'})
    df_MAR = df_MAR.rename(columns={'PrDM' : 'depth'})
    df_MAR = df_MAR.rename(columns={'T068C' : 'temperature'})
    df_MAR = df_MAR.rename(columns={'Sal00' : 'salinity'})
    df_MAR = df_MAR.rename(columns={'TA (umol/kg)' : 'TA'})
    df_MAR = df_MAR.rename(columns={'TIC (umol/kg)' : 'TIC'})
    df_MAR = df_MAR.rename(columns={'O2_Concentration(ml/l)' : 'O2_ml_l'})

    df_MAR['Region'] ='MAR'
    df_IML['Region'] = 'Gulf' 
    
    df = pd.concat([
    df_NL[['StationID', 'latitude', 'longitude', 'depth', 'temperature', 'salinity','satO2_perc', 'O2_ml_l', 'TA', 'TIC', 'pH', 'Omega_A', 'Omega_C', 'Region']],
    df_IML[['StationID', 'latitude', 'longitude', 'depth', 'temperature', 'salinity','satO2_perc', 'O2_ml_l', 'TA', 'TIC', 'pH', 'Omega_A', 'Omega_C', 'Region']],
    df_MAR[['StationID', 'latitude', 'longitude', 'depth', 'temperature', 'salinity','satO2_perc', 'O2_ml_l', 'TA', 'TIC', 'pH', 'Omega_A', 'Omega_C', 'Region']]
    ], sort=True)


print (str(my_year)+' '+my_season+' '+my_depth+' '+my_variable)

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

df = df[df.index.year==my_year]
#df = df.loc[my_year] ###locates year
#df=df[df.Region.str.contains(my_region)] ###locates region
    
#######locates season##########  <------- check this for other season
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
    print('All seasons selected')


if df[my_variable].isna().values.all():# or df.size == 0:
    print ('!!! no data for this season !!!')
    sys.exit()

df.dropna(subset=[my_variable, 'pH', 'Omega_A', 'StationID', 'depth'], axis=0, inplace=True)
df = df.reset_index(drop=True)
# Set depth as float (had some problem with some data set)
df = df.astype({'depth':'float', my_variable:'float'})  
#######locates depth -- either surface, bottom, or within a range of depths############
if my_depth == 'surface':
    df = df.loc[df.groupby('StationID')['depth'].idxmin()] #group by station then pull "min or max depth"
    df = df.loc[df.depth <20] #take all depths >10m (for bottom) to eliminate lone surface samples, all depths <20m (for surface) to eliminate lone deep sample
if my_depth == 'bottom':
    df = df.loc[df.groupby('StationID')['depth'].idxmax()] #group by station then pull "min or max depth"
    df = df.loc[df.depth >10] #take all depths >10m (for bottom) to eliminate lone surface samples
    if my_year == 2019: # Some stations to be flagged (deep bottle not bottom)
        df.drop(df[df.StationID=='GULD_04'].index, inplace=True)
        df.drop(df[df.StationID=='HL_07'].index, inplace=True)
        df.drop(df[df.StationID=='LL_08'].index, inplace=True)

if my_depth == 'range':
    df = df.loc[(df.depth >=135) & (df.depth <=165)]
    df = df.loc[df.groupby('StationID')['depth'].idxmax()] #group by station then pull "min or max depth"


############create bathymetry################
dataFile = '/home/cyrf0006/data/GEBCO/GEBCO_2014_1D.nc'
###map boundaries###
lonLims = [-72, -41.5]  
latLims = [40.5, 58.4]
#latLims = [40.5, 55.1]

##### this is the dataset that will be plotted ##########
lat_data = np.array(df.latitude)
lon_data = np.array(df.longitude)
data = np.array(df[my_variable])
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
if my_variable == 'Omega_A':
    ax.text(-69, 56, str(my_year), horizontalalignment='left', color='white',
            fontsize=12, fontweight='bold', transform=ccrs.PlateCarree())

# Save figure
fig.savefig('AZMP_OA_'+str(my_year)+'_'+my_season+'_'+my_variable+'_'+my_depth+'.png', format='png', dpi=300, bbox_inches='tight')
#lt.close('all') 


## gm convert -trim AZMP_OA_2019_spring_Omega_A_bottom.png  AZMP_OA_2019_spring_Omega_A_bottom.png
## gm convert -trim AZMP_OA_2019_spring_satO2_perc_bottom.png AZMP_OA_2019_spring_satO2_perc_bottom.png
## gm convert -trim AZMP_OA_2019_spring_pH_bottom.png AZMP_OA_2019_spring_pH_bottom.png
## montage AZMP_OA_2019_spring_Omega_A_bottom.png  AZMP_OA_2019_spring_pH_bottom.png AZMP_OA_2019_spring_satO2_perc_bottom.png  -tile 1x3 -geometry +10+10  -background white AZMP_OA_2019_spring.png


## gm convert -trim AZMP_OA_2018_spring_Omega_A_bottom.png  AZMP_OA_2018_spring_Omega_A_bottom.png
## gm convert -trim AZMP_OA_2018_spring_satO2_perc_bottom.png AZMP_OA_2018_spring_satO2_perc_bottom.png
## gm convert -trim AZMP_OA_2018_spring_pH_bottom.png AZMP_OA_2018_spring_pH_bottom.png
## montage AZMP_OA_2018_spring_Omega_A_bottom.png  AZMP_OA_2018_spring_pH_bottom.png AZMP_OA_2018_spring_satO2_perc_bottom.png  -tile 1x3 -geometry +10+10  -background white AZMP_OA_2018_spring.png

## montage AZMP_OA_2018_spring.png  AZMP_OA_2019_spring.png -tile 2x1 -geometry +10+10  -background white AZMP_OA_2018-2019_spring.png


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
##         ('Station', '@StationID'),
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

## # Another one in map-like display
# pH
Viridis7.reverse()
Spectral7.reverse()
mapper = linear_cmap(field_name='pH', palette=Spectral7 ,low=7.5 ,high=8.2)
output_file('OA_2019_pH_map.html')
source = ColumnDataSource(df)
p = figure(plot_width=900, plot_height=600)
p.scatter('longitude', 'latitude', source=source, size=8,
          line_color=mapper,color=mapper)
p.title.text = 'Carbonate parameters in the Atlantic Zone - Spring 2019, bottom'
p.xaxis.axis_label = 'Longitude'
p.yaxis.axis_label = 'Latitude'
p.add_tools(HoverTool(
    tooltips=[
        ('pH', '@pH'),
        ('Saturation Aragonite', '@Omega_A'),
        ('Saturation Calcite', '@Omega_C'),
        ('Total Alkalinity', '@TA'),
        ('Inorg. Carbon', '@TIC'),
        ('O2 sat (%)', '@satO2_perc'),
        ('O2 dissolved (ml/l)', '@O2_ml_l'),
        ('Region', '@Region'),
        ('Station', '@StationID'),
        ('Depth (m)', '@depth'),
        ('Temperature', '@temperature'),
        ('Salinity', '@salinity'),
        ('Latitude', '@latitude'),
        ('Longitude', '@longitude')      
    ],

    mode='mouse'
))
color_bar = ColorBar(color_mapper=mapper['transform'], width=8,  location=(0,0),title="pH")
p.add_layout(color_bar, 'right')
show(p)

# 2. Omega_A
Spectral8.reverse()
mapper = linear_cmap(field_name='Omega_A', palette=Spectral8 ,low=0 ,high=2)
output_file('OA_2019_omegaA_map.html')
source = ColumnDataSource(df)
p = figure(plot_width=900, plot_height=600)
p.scatter('longitude', 'latitude', source=source, size=8,
          line_color=mapper,color=mapper)
p.title.text = 'Carbonate parameters in the Atlantic Zone - Spring 2019, bottom'
p.xaxis.axis_label = 'Longitude'
p.yaxis.axis_label = 'Latitude'
p.add_tools(HoverTool(
    tooltips=[
        ('pH', '@pH'),
        ('Saturation Aragonite', '@Omega_A'),
        ('Saturation Calcite', '@Omega_C'),
        ('Total Alkalinity', '@TA'),
        ('Inorg. Carbon', '@TIC'),
        ('O2 sat (%)', '@satO2_perc'),
        ('O2 dissolved (ml/l)', '@O2_ml_l'),
        ('Region', '@Region'),
        ('Station', '@StationID'),
        ('Depth (m)', '@depth'),
        ('Temperature', '@temperature'),
        ('Salinity', '@salinity'),
        ('Latitude', '@latitude'),
        ('Longitude', '@longitude')      
    ],

    mode='mouse'
))
color_bar = ColorBar(color_mapper=mapper['transform'], width=8,  location=(0,0),title="Omega Ar")
p.add_layout(color_bar, 'right')
show(p)

# 2. Omega_C
mapper = linear_cmap(field_name='Omega_C', palette=Spectral8 ,low=0 ,high=2)
output_file('OA_2019_omegaC_map.html')
source = ColumnDataSource(df)
p = figure(plot_width=900, plot_height=600)
p.scatter('longitude', 'latitude', source=source, size=8,
          line_color=mapper,color=mapper)
p.title.text = 'Carbonate parameters in the Atlantic Zone - Spring 2019, bottom'
p.xaxis.axis_label = 'Longitude'
p.yaxis.axis_label = 'Latitude'
p.add_tools(HoverTool(
    tooltips=[
        ('pH', '@pH'),
        ('Saturation Aragonite', '@Omega_A'),
        ('Saturation Calcite', '@Omega_C'),
        ('Total Alkalinity', '@TA'),
        ('Inorg. Carbon', '@TIC'),
        ('O2 sat (%)', '@satO2_perc'),
        ('O2 dissolved (ml/l)', '@O2_ml_l'),
        ('Region', '@Region'),
        ('Station', '@StationID'),
        ('Depth (m)', '@depth'),
        ('Temperature', '@temperature'),
        ('Salinity', '@salinity'),
        ('Latitude', '@latitude'),
        ('Longitude', '@longitude')      
    ],

    mode='mouse'
))
color_bar = ColorBar(color_mapper=mapper['transform'], width=8,  location=(0,0),title="Omega Ca")
p.add_layout(color_bar, 'right')
show(p)


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

