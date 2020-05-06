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

my_season = 'all'
my_depth = 'all'
df = pd.read_excel('/home/cyrf0006/github/AZMP-NL/datasets/carbonates/AZMP_OA_CO2stats.xlsx')
df = df.set_index('timestamp', drop=False)

    
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


#df.dropna(subset=[my_variable, 'pH', 'Omega_A', 'StationID', 'depth'], axis=0, inplace=True)
# Set depth as float (had some problem with some data set)
df = df.astype({'depth':'float'})  
# locates depth -- either surface, bottom, or within a range of depths
if my_depth == 'surface':
    df = df.loc[df.groupby('StationID')['depth'].idxmin()] #group by station and get min
    df = df.loc[df.depth <10] #take all depths <10m (for surface) 
elif my_depth == 'bottom':
    df = df.loc[df.groupby('StationID')['depth'].idxmax()] #group by station and get max
    df = df.loc[df.depth >10] #take all depths >10m (for bottom)
    if my_year == 2019: # Some stations to be flagged (deep bottle not bottom)
        df.drop(df[df.StationID=='GULD_04'].index, inplace=True)
        df.drop(df[df.StationID=='HL_07'].index, inplace=True)
        df.drop(df[df.StationID=='LL_08'].index, inplace=True)
elif my_depth == 'range':
    df = df.loc[(df.depth >=135) & (df.depth <=165)]
    df = df.loc[df.groupby('StationID')['depth'].idxmax()]
else:
    print('no depth range specified, take all')

df = df.drop(columns='timestamp') 
df = df.drop(columns='Unnamed: 0')
df.reset_index(inplace=True)      
df.timestamp = df.timestamp.dt.strftime(date_format = '%Y-%m-%d')   


## ## ---- Test with Bokeh plot ---- ##
from bokeh.plotting import figure, output_file, show
from bokeh.models import ColumnDataSource
from bokeh.models.tools import HoverTool
from bokeh.models import DatetimeTickFormatter
from bokeh.transform import factor_cmap, factor_mark
# For colorbar
from bokeh.models import ColorBar
from bokeh.palettes import Spectral8, Viridis7, Spectral7, Viridis10
from bokeh.transform import linear_cmap


mapper = linear_cmap(field_name='depth', palette=np.flip(Viridis10), low=0, high=300)
output_file('OA_T-S.html')
source = ColumnDataSource(df)
p = figure(plot_width=900, plot_height=600)
p.scatter('salinity', 'temperature', source=source, size=8, line_color=mapper,color=mapper)
p.title.text = 'Carbonate parameters in the Atlantic Zone - 2014-2019'
p.xaxis.axis_label = 'Salinity'
p.yaxis.axis_label = 'Temperature'
p.add_tools(HoverTool(
    tooltips=[
        ('Temperature', '@temperature'),
        ('Salinity', '@salinity'),
        ('Depth (m)', '@depth'),
        ('pH', '@pH'),
        ('Saturation Aragonite', '@Omega_A'),
        ('Saturation Calcite', '@Omega_C'),
        ('Total Alkalinity', '@TA'),
        ('Inorg. Carbon', '@TIC'),
        ('O2 sat (%)', '@satO2_perc'),
        ('O2 dissolved (ml/l)', '@O2_ml_l'),
        ('Region', '@Region'),
        ('Station', '@StationID'),
        ('Sample date', '@timestamp'),
        ('Latitude', '@latitude'),
        ('Longitude', '@longitude')      
    ],

    mode='mouse'
))
color_bar = ColorBar(color_mapper=mapper['transform'], width=8,  location=(0,0),title="T-S diagram")
p.add_layout(color_bar, 'right')
show(p)
