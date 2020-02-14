"""Some tools to plot sections from AZMP bottle files

** This script whould be moved in ~/github/AZMP once stable

Check azmp_utils.masterfile_section_to_multiindex.py

References
----------

Atlantic Zone Monitoring Program @NAFC:
https://azmp-nl.github.io/


-> run in /home/cyrf0006/AZMP/dig_project/btl_work
"""

__author__ = 'Frederic.Cyr@dfo-mpo.gc.ca'
__version__ = '0.1'

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import azmp_sections_tools as azst
from matplotlib.patches import Polygon
from scipy.io import loadmat
import cmocean
import os

## ----  Some info that should eventually put as function inputs ---- ##
SECTION = 'SESPB'
SECTION_bathy = 'SESPB'
VAR = 'PO4'
YEAR = 2014
SEASON = 'fall'
#levels = np.linspace(0,1,11) # f_pw
#levels = np.linspace(0,16,9) # SIO
#levels = np.linspace(70,100,11) # O2
#levels = np.linspace(0,18,10) # NO3
v = 20
v_sig = [25, 26, 27, 27.5]
v_anom = 10
CMAP = cmocean.cm.thermal

ZMAX = 400

# derived parameters
if VAR == 'temperature':
    v = np.arange(-2,11,1)
    v_anom = np.linspace(-3.5, 3.5, 15)
    v_anom = np.delete(v_anom, np.where(v_anom==0)) 
    CMAP = cmocean.cm.thermal
elif VAR == 'salinity':
    v = np.arange(29,36,.5)
    v_anom = np.linspace(-1.5, 1.5, 16)
    CMAP = cmocean.cm.haline
elif VAR == 'oxygen':
    v = np.arange(5.5,8.5,.2)
    v_anom = np.linspace(-1, 1, 11)
    CMAP = cmocean.cm.thermal
elif VAR == 'satO2_perc':
    v = np.arange(70, 100, 2)
    v_anom = np.linspace(-10, 10, 11)
    CMAP = cmocean.cm.thermal
elif VAR == 'PO4':
    v = np.arange(0, 2, .1)
    v_anom = np.linspace(-1, 1, 11)
    CMAP = cmocean.cm.thermal
elif VAR == 'NO3':
    v = np.arange(0, 16, 1)
    v_anom = np.linspace(-5, 5, 11)
    CMAP = cmocean.cm.thermal
elif VAR == 'SIO':
    v = np.arange(0, 20, 1)
    v_anom = np.linspace(-5, 5, 11)
    CMAP = cmocean.cm.thermal
else:
    v = 10
    v_anom = 10
    CMAP = cmocean.cm.thermal

## ---- Retrieve bottle data ---- ##
SECTION_FILE = 'bottle_data_multiIndex_' + SECTION + '.pkl'
df = pd.read_pickle(SECTION_FILE)

## ---- Compute climatology ---- ##
df_clim = df.xs((VAR, SEASON),level=('variable', 'season'))
df_clim = df_clim.groupby(level=0).apply(lambda x: x.mean())
# Compute also sigma-t
df_sigt_clim = df.xs(('sigmat', SEASON),level=('variable', 'season'))
df_sigt_clim = df_sigt_clim.groupby(level=0).apply(lambda x: x.mean())

## ---- Compute current year & anomaly ---- ##
if YEAR in df.index.get_level_values(1).unique():
    df_year = df.xs((YEAR, VAR, SEASON),level=('year', 'variable', 'season'))
    df_year = df_year.groupby(level=0).apply(lambda x: x.mean())
    df_sigt_year = df.xs((YEAR, 'sigmat', SEASON),level=('year', 'variable', 'season'))
    df_sigt_year = df_sigt_year.groupby(level=0).apply(lambda x: x.mean())

    df_anom = df_year-df_clim

    # Drop NaNs
    df_year = df_year.dropna(axis=0,how='all')
    df_sigt_year = df_sigt_year.dropna(axis=0,how='all')
    df_clim = df_clim.dropna(axis=0,how='all')
    df_sigt_clim = df_sigt_clim.dropna(axis=0,how='all')
    df_anom = df_anom.dropna(axis=0,how='all')

    if df_year.size == 0:
        print(' !!! Empty section [return None] !!!')
        
    ## ---- Load station lat/lon ---- ##
    df_stn = pd.read_excel('/home/cyrf0006/github/AZMP-NL/data/STANDARD_SECTIONS.xlsx')
    df_stn = df_stn.drop(['SECTION', 'LONG'], axis=1)
    df_stn = df_stn.rename(columns={'LONG.1': 'LON'})
    df_stn = df_stn.dropna()
    df_stn = df_stn[df_stn.STATION.str.contains(SECTION)]
    df_stn = df_stn.reset_index(drop=True)
    # remove hyphen in stn names (normally should keep it everywhere, but simpler here)
    df_stn['STATION'] = df_stn.STATION.str.replace('-', '')

    ## ---- Compute distance vector ---- ##
    distance = np.full((df_clim.index.shape), np.nan)
    lat0 = df_stn[df_stn.index==0]['LAT']
    lon0 = df_stn[df_stn.index==0]['LON']
    for i, stn in enumerate(df_clim.index):
        lat_stn = df_stn[df_stn.STATION==stn]['LAT']
        lon_stn = df_stn[df_stn.STATION==stn]['LON']      
        distance[i] = azst.haversine(lon0, lat0, lon_stn, lat_stn)
        
    XLIM =distance.max()

    # Distance for current year
    distance_year = np.full((df_year.index.shape), np.nan)
    for i, stn in enumerate(df_year.index):
        lat_stn = df_stn[df_stn.STATION==stn]['LAT']
        lon_stn = df_stn[df_stn.STATION==stn]['LON']      
        distance_year[i] = azst.haversine(lon0, lat0, lon_stn, lat_stn)

    distance_sigt = np.full((df_sigt_year.index.shape), np.nan)
    for i, stn in enumerate(df_sigt_year.index):
        lat_stn = df_stn[df_stn.STATION==stn]['LAT']
        lon_stn = df_stn[df_stn.STATION==stn]['LON']      
        distance_sigt[i] = azst.haversine(lon0, lat0, lon_stn, lat_stn)

    ## ---- Retrieve bathymetry using function ---- ##
    bathymetry = azst.section_bathymetry(SECTION_bathy)

    ## ---- plot Figure ---- ##
    #XLIM = df_section_itp.index[-1][1]
    fig = plt.figure()
    # ax1
    ax = plt.subplot2grid((3, 1), (0, 0))
    c = plt.contourf(distance_year, df_year.columns, df_year.T, v, cmap=CMAP, extend='both')
    c_sig1 = plt.contour(distance_sigt, df_sigt_year.columns, df_sigt_year.T, v_sig, colors='gray', linewidths=1)
    for i in distance_year:
        plt.plot([i, i], [0, ZMAX], '--k', alpha=.5)
    ax.set_ylim([0, ZMAX])
    ax.set_xlim([0,  XLIM])
    plt.clabel(c_sig1, inline=1, fontsize=10, colors='gray', fmt='%1.1f')
    ax.set_ylabel('Depth (m)', fontWeight = 'bold')
    ax.invert_yaxis()
    Bgon = plt.Polygon(bathymetry,color=np.multiply([1,.9333,.6667],.4), alpha=1, zorder=10)
    ax.add_patch(Bgon)
    plt.colorbar(c)
    ax.xaxis.label.set_visible(False)
    ax.tick_params(labelbottom='off')
    ax.set_title(VAR + ' for section ' + SECTION + ' - ' + SEASON + ' ' + str(YEAR))

    # ax2
    ax2 = plt.subplot2grid((3, 1), (1, 0))
    c = plt.contourf(distance, df_clim.columns, df_clim.T, v, cmap=CMAP, extend='both')
    c_sig2 = plt.contour(distance, df_sigt_clim.columns, df_sigt_clim.T, v_sig, colors='gray', linewidths=1)
    for i in distance:
        plt.plot([i, i], [0, ZMAX], '--k', alpha=.5)
    ax2.set_ylim([0, ZMAX])
    ax2.set_xlim([0,  XLIM])
    plt.clabel(c_sig2, inline=1, fontsize=10, colors='gray', fmt='%1.1f')
    ax2.set_ylabel('Depth (m)', fontWeight = 'bold')
    ax2.invert_yaxis()
    Bgon = plt.Polygon(bathymetry,color=np.multiply([1,.9333,.6667],.4), alpha=1, zorder=10)
    ax2.add_patch(Bgon)
    plt.colorbar(c)
    ax2.xaxis.label.set_visible(False)
    ax2.tick_params(labelbottom='off')
    ax2.set_title('1999-' + str(df.index.get_level_values(1).max()) + ' climatology')

    # ax3
    ax3 = plt.subplot2grid((3, 1), (2, 0))
    c = plt.contourf(distance_year, df_anom.columns, df_anom.T, v_anom, cmap=cmocean.cm.balance, extend='both')
    ax3.set_ylim([0, ZMAX])
    ax3.set_xlim([0,  XLIM])
    ax3.set_ylabel('Depth (m)', fontWeight = 'bold')
    ax3.set_xlabel('Distance (km)', fontWeight = 'bold')
    ax3.invert_yaxis()
    Bgon = plt.Polygon(bathymetry,color=np.multiply([1,.9333,.6667],.4), alpha=1, zorder=10)
    ax3.add_patch(Bgon)
    plt.colorbar(c)
    ax3.set_title(r'Anomaly')

    fig.set_size_inches(w=8,h=12)
    fig_name = 'btl_' + VAR + '_' + SECTION + '_' + SEASON + '_' + str(YEAR) + '.png' 
    fig.savefig(fig_name, dpi=200)
    os.system('convert -trim ' + fig_name + ' ' + fig_name)
else:
    print(' --> No data for this year')

#montage temperature_BB_summer_2018.png salinity_BB_summer_2018.png  -tile 2x1 -geometry +10+10  -background white BB_summer_2018.png 
#montage temperature_SI_summer_2018.png salinity_SI_summer_2018.png  -tile 2x1 -geometry +10+10  -background white SI_summer_2018.png 
#montage temperature_FC_summer_2018.png salinity_FC_summer_2018.png  -tile 2x1 -geometry +10+10  -background white FC_summer_2018.png 
#montage temperature_WB_summer_2018.png salinity_WB_summer_2018.png  -tile 2x1 -geometry +10+10  -background white WB_summer_2018.png 

# Save section in csv
df_year['distance'] = distance_year
df_year.set_index('distance', append=True, inplace=True)
df_clim['distance'] = distance
df_clim.set_index('distance', append=True, inplace=True)
df_clim.T.to_csv('test.csv', float_format='%.3f')  
