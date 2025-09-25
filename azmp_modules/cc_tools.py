"""Some utilitary tools for project on Carbonate Chemistry (aka cc_)

Contains following functions:
- 
-

** idea: add function 'plot_nafo_division(dict)'

----------

Atlantic Zone Monitoring Program @NAFC:
https://azmp-nl.github.io/

"""

__author__ = 'Frederic.Cyr@dfo-mpo.gc.ca'
__version__ = '0.1'

import h5py
import os
import sys
import netCDF4
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr
#from mpl_toolkits.basemap import Basemap
from scipy.interpolate import griddata
from scipy.interpolate import interp1d  # to remove NaNs in profiles
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
#from shapely.ops import cascaded_union
from area import area # external fns to compute surface area
from seawater import extras as swx
import gsw

## from matplotlib.colors import Normalize
from matplotlib.colors import from_levels_and_colors
## import numpy as np
import h5py
import cmocean
import cmocean.cm as cmo
import cartopy. crs as ccrs
import cartopy.feature as cpf
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

# AZMP modules
import cc_variable_list2 as vl
import cc_variable_list_NL as vlNL
import azmp_sections_tools as azst

def HUD2014_to_multiindex(section, z_vec):
    """
    This function reads the HUD2014 Excel sheet and pickle a Multi-index DataFrame organized as:


    Depth                            0    2    4        6        8        10   \
    station year variable    season                                             
    BB01    1999 temperature spring  NaN  NaN  NaN      NaN      NaN      NaN   
                             summer  NaN  NaN  NaN   9.6315   9.6265  9.59314   
                             fall    NaN  NaN  NaN  3.23428  3.23473  3.23507   
                 salinity    spring  NaN  NaN  NaN      NaN      NaN      NaN   
                             summer  NaN  NaN  NaN  34.1785  34.1793  34.1823   
                                       [...]                    

    INPUT:
    - section = section name (ex: 'BB', 'SI', etc.)
    - z_vec = binned depth vector (see example)

    OUTPUT:
    - in addition of pickling the Multi-index, it returns it for verification 
                                       
    usage ex:
    import numpy as np
    import azmp_utils as azu
    azu.masterfile_section_to_multiindex('BB', np.arange(0,350, 2))

    Multi-index manipulation examples:
    df_BB = pd.read_pickle('bottle_data_multiIndex_BB.pkl')
    # 1. Single station average vertical profile
    A = df_BB.xs(('BB01', 'NO3'),level=('station', 'variable'))
    A.groupby(level=0).apply(lambda x: x.mean()).mean()

    # 2. Single year section
    B = df_BB.xs((2016, 'NO3'),level=('year', 'variable'))
    B.groupby(level=0).apply(lambda x: x.mean())

    # 3. 1999-2016 section climato
    C = df_BB.xs(('NO3'),level=('variable'))
    C.groupby(level=0).apply(lambda x: x.mean())

    """

    ## ---- List of variable to export ---- ##
    varname = pd.Series(['temperature', 'salinity', 'sigmat', 'oxygen', 'PO4', 'SIO', 'NO3', 'chlor', 'fluor', 'satO2_perc', 'NPratio', 'f_pw'])
    varname.name='variable'

    ## ----  Load data ---- ##
    df = pd.read_excel('/home/cyrf0006/github/AZMP-NL/data/HUDSON2014_surveys_with_ancillary_data_31Aug2015.xlsx')

    ## ---- Remove space in column names ---- ##
    cols = df.columns
    cols = cols.map(lambda x: x.replace(' ', '_') if isinstance(x, (str, unicode)) else x)
    df.columns = cols

    

    
    # Set date as index
    df = df.set_index('sas_date')
    # Drop other time-related columns
    df = df.drop(['Day', 'Month', 'Year'], axis=1)
    # Keep only targeted section
    df = df[df.section == section]
    # Other derived variables
    df['satO2'] = swx.satO2(df['salinity'], df['temp'])
    df['satO2_perc'] = df['oxygen']/df['satO2']*100
    df['NPratio'] = df['NO3']/df['PO4']
    df['f_pw'] = (df.NO3 - (17.499*df.PO4 - 3.072)) / ((12.368*df.PO4 - 10.549) - (17.499*df.PO4 - 3.072))
    # rename temperature
    df = df.rename(columns={'temp':'temperature'})
    
    sname_unique = pd.Series(df.sname.unique())
    sname_unique.name='station'

    df_list_station = []
    for i, stn in enumerate(sname_unique):

        df_sname = df[df.sname==stn]
        years_unique = df_sname.index.year.unique()
        years_unique.name='year'

        df_list_year = []
        for j, year in enumerate(years_unique):

            df_year = df_sname[df_sname.index.year == year]

            # Select only seasons
            df_spring = df_year[(df_year.index.month>=4) & (df_year.index.month<=6)]
            df_summer = df_year[(df_year.index.month>=7) & (df_year.index.month<=9)]
            df_fall = df_year[(df_year.index.month>=10) & (df_year.index.month<=12)]

            df_list_var = []
            for k, var in enumerate(varname):

                df_season_clean = pd.DataFrame(index=['spring', 'summer', 'fall'], columns=z_vec)
                df_season_clean.index.name='season'
                df_season_clean.columns.name='Depth'

                # Spring
                var_itp = np.full((z_vec.shape), np.nan)
                series_var = df_spring[var]
                if series_var.size>1: # <---- Here I end up ignoring some data if only one sample per profile...
                    series_z = df_spring.depth
                    idx_good = np.argwhere((~np.isnan(series_var)))            
                    interp = interp1d(series_z.values, series_var.values)  
                    idx_interp = np.where((z_vec>=series_z.min()) & (z_vec<=series_z.max()))
                    var_itp[idx_interp] = interp(z_vec[idx_interp]) # interpolate only where possible (1st to last good idx)
                    var_itp_series = var_itp
                    df_season_clean.loc['spring'] = var_itp

                # Summer
                var_itp = np.full((z_vec.shape), np.nan)
                series_var = df_summer[var]
                if series_var.size>1:
                    series_z = df_summer.depth
                    idx_good = np.argwhere((~np.isnan(series_var)))            
                    interp = interp1d(series_z.values, series_var.values)  
                    idx_interp = np.where((z_vec>=series_z.min()) & (z_vec<=series_z.max()))
                    var_itp[idx_interp] = interp(z_vec[idx_interp]) # interpolate only where possible (1st to last good idx)
                    var_itp_series = var_itp
                    df_season_clean.loc['summer'] = var_itp

                # Fall
                var_itp = np.full((z_vec.shape), np.nan)
                series_var = df_fall[var]
                if series_var.size>1:
                    series_z = df_fall.depth
                    idx_good = np.argwhere((~np.isnan(series_var)))            
                    interp = interp1d(series_z.values, series_var.values)  
                    idx_interp = np.where((z_vec>=series_z.min()) & (z_vec<=series_z.max()))
                    var_itp[idx_interp] = interp(z_vec[idx_interp]) # interpolate only where possible (1st to last good idx)
                    var_itp_series = var_itp
                    df_season_clean.loc['fall'] = var_itp


                df_list_var.append(df_season_clean)


            df_list_year.append(pd.concat(df_list_var,keys=varname))


        df_list_station.append(pd.concat(df_list_year,keys=years_unique))


    section_mindex = pd.concat(df_list_station,keys=sname_unique)  

    section_mindex.to_pickle('bottle_data_multiIndex_' + section +'.pkl')

    return section_mindex


def seasonal_map(VARIABLE, YEAR, SEASON, DEPTH):

    '''
    This function was created to generate seasonal maps for the Carbonate dataset paper.
    This script can also be used for AZMP reporting (usually labeled azmp_cc_zonal_YYYY.py).

    usage ex:
    import cc_tools as cc
    cc.seasonal_map('Temperature_(degC)', 2017, 'summer', 'surface')

    variable options are:
    VARIABLE = 'Temperature_(degC)'
    VARIABLE = 'Salinity_(psu)'
    VARIABLE = 'Oxygen_Saturation_(%)'
    VARIABLE = 'Total_Alkalinity_(umol/kg)'
    VARIABLE = 'Inorganic_Carbon_(umol/kg)'
    VARIABLE = 'pCO2_(uatm)',
    VARIABLE = 'pH_Total_(total_scale)'

    See another usage example in a loop in:
    azmp_zonal_cc_2017_ms.py

    Frederic.Cyr@dfo-mpo.gc.ca
    September 2022
    '''


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

    # Figure name (to avoid parenthesis, etc.)
    if VARIABLE == 'Temperature_(degC)':
        FIG_VAR = 'Temperature'
        FIG_LABEL = 'T'
    elif VARIABLE == 'Salinity_(psu)':
        FIG_VAR = 'Salinity'
        FIG_LABEL = 'S'
    elif VARIABLE == 'Total_Alkalinity_(umol/kg)':
        FIG_VAR = 'TA'
        FIG_LABEL = 'TA'
    elif VARIABLE == 'Inorganic_Carbon_(umol/kg)':
        FIG_VAR = 'DIC'
        FIG_LABEL = 'DIC'
    elif VARIABLE == 'pCO2_(uatm)':
        FIG_VAR = 'pCO2'
        FIG_LABEL = r'$pCO_2$'
    elif VARIABLE == 'Omega_Aragonite_(unitless)':
        FIG_VAR = 'OmegaA'
        FIG_LABEL = r'$\Omega_{\rm arg}$'
    elif VARIABLE == 'Omega_Calcite_(unitless)':
        FIG_VAR = 'OmegaC'
        FIG_LABEL = r'$\Omega_{\rm cal}$'
    elif VARIABLE == 'pH_Total_(total_scale)':
        FIG_VAR = 'pH'
        FIG_LABEL = 'pH'
    elif VARIABLE == 'Oxygen_Saturation_(%)':
        FIG_VAR = 'DO_perc'    
        FIG_LABEL = r'$O_2$'
    elif VARIABLE == 'Dissolved_Oxygen_(mL/L)':
        FIG_VAR = 'DO'
        FIG_LABEL = r'$O_2$'

    
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
        df = df[(df.index.month>=7) & (df.index.month<=9)]
    elif SEASON == 'fall':
        dfng = df[~df.Region.str.contains('MAR')]
        dfng = dfng.assign(x=dfng.index.strftime('%m-%d')).query("'10-15' <= x <= '12-31'").drop('x',axis=1)
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
    # locate depth -- either surface, bottom, or within a range of depths
    if DEPTH == 'surface':
        #df = df.loc[df.groupby('Station_Name')['Depth_(dbar)'].idxmin()] #group by station then pull "min or max depth"
        df = df.loc[df.groupby('Longitude_(degEast)')['Depth_(dbar)'].idxmin()] #group by station then pull "min or max depth"
        df = df.loc[df['Depth_(dbar)'] <20] #take all depths >10m (for bottom) to eliminate lone surface samples, all depths <20m (for surface) to eliminate lone deep sample
    if DEPTH == 'bottom':
        #df = df.loc[df.groupby('Station_Name')['Depth_(dbar)'].idxmax()] #group by station then pull "min or max depth"
        df = df.loc[df.groupby('Longitude_(degEast)')['Depth_(dbar)'].idxmax()] #group by station then pull "min or max depth"
        df = df.loc[df['Depth_(dbar)'] >10] #take all depths >10m (for bottom) to eliminate lone surface samples
        if YEAR == 2019: # Some stations to be flagged (deep bottle not bottom)
            df.drop(df[df.Station_Name=='GULD_04'].index, inplace=True)
            df.drop(df[df.Station_Name=='HL_07'].index, inplace=True)
            df.drop(df[df.Station_Name=='LL_08'].index, inplace=True)

    if DEPTH == 'range':
        df = df.loc[(df['Depth_(dbar)'] >=135) & (df.depth <=165)]
        df = df.loc[df.groupby('Station_Name')['Depth_(dbar)'].idxmax()] #group by station then pull "min or max depth"


    ## Load bathymetry
    dataFile = '/home/cyrf0006/data/GEBCO/GEBCO_2014_1D.nc'
    # map boundaries
    lonLims = [-72, -41.5]  
    latLims = [40.5, 58.4]
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


    ## draw map with bathymetry
    print('--- Now plot ---')
    fig = plt.figure(figsize=(7,5))
    ax = fig.add_subplot(111, projection=ccrs.Mercator())
    ax.set_extent([-72, -41.5, 40.5, 58.4], crs=ccrs.PlateCarree())
    ax.add_feature(cpf.NaturalEarthFeature('physical', 'coastline', '50m', edgecolor='k', alpha=0.7, linewidth=0.6, facecolor='black'), zorder=1)#cpf.COLORS['land']))
    m=ax.gridlines(linewidth=0.5, color='black', draw_labels=True, alpha=0.5)
    m.right_labels=False
    m.top_labels=False
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
    SEASON_TEXT = ax.text(-72, 56, SEASON + ' ' +  str(YEAR), horizontalalignment='left', verticalalignment = 'top', color='white', fontsize=12, fontweight='bold', transform=ccrs.PlateCarree())
    DEPTH_TEXT = ax.text(-72, 56, DEPTH + ' ' + FIG_LABEL, horizontalalignment='left', verticalalignment = 'bottom', color='white', fontsize=12, fontweight='bold', transform=ccrs.PlateCarree())

    # add letter ID (for Fgiures 6-8 in manuscript):
    if SEASON == 'fall':
        if DEPTH == 'surface':
            if ((VARIABLE=='Temperature_(degC)') | (VARIABLE=='Total_Alkalinity_(umol/kg)') | (VARIABLE=='pH_Total_(total_scale)')):
                ID_TEXT = ax.text(-72, 57, 'a) ', horizontalalignment='right', verticalalignment = 'bottom', color='black', fontsize=12, fontweight='bold', transform=ccrs.PlateCarree())
            elif ((VARIABLE=='Salinity_(psu)') | (VARIABLE=='Inorganic_Carbon_(umol/kg)') | (VARIABLE=='Omega_Aragonite_(unitless)')):
                  ID_TEXT = ax.text(-72, 57, 'c) ', horizontalalignment='right', verticalalignment = 'bottom', color='black', fontsize=12, fontweight='bold', transform=ccrs.PlateCarree())
            elif ((VARIABLE=='Oxygen_Saturation_(%)') | (VARIABLE=='pCO2_(uatm)') | (VARIABLE=='Omega_Calcite_(unitless)')):
                  ID_TEXT = ax.text(-72, 57, 'e) ', horizontalalignment='right', verticalalignment = 'bottom', color='black', fontsize=12, fontweight='bold', transform=ccrs.PlateCarree())
        elif DEPTH == 'bottom':
            if ((VARIABLE=='Temperature_(degC)') | (VARIABLE=='Total_Alkalinity_(umol/kg)') | (VARIABLE=='pH_Total_(total_scale)')):
                ID_TEXT = ax.text(-72, 57, 'b) ', horizontalalignment='right', verticalalignment = 'bottom', color='black', fontsize=12, fontweight='bold', transform=ccrs.PlateCarree())
            elif ((VARIABLE=='Salinity_(psu)') | (VARIABLE=='Inorganic_Carbon_(umol/kg)') | (VARIABLE=='Omega_Aragonite_(unitless)')):
                  ID_TEXT = ax.text(-72, 57, 'd) ', horizontalalignment='right', verticalalignment = 'bottom', color='black', fontsize=12, fontweight='bold', transform=ccrs.PlateCarree())
            elif ((VARIABLE=='Oxygen_Saturation_(%)') | (VARIABLE=='pCO2_(uatm)') | (VARIABLE=='Omega_Calcite_(unitless)')):
                  ID_TEXT = ax.text(-72, 57, 'f) ', horizontalalignment='right', verticalalignment = 'bottom', color='black', fontsize=12, fontweight='bold', transform=ccrs.PlateCarree())    

         
    # Save figure
    fig_name = 'AZMP_OA_'+str(YEAR)+'_'+SEASON+'_'+FIG_VAR+'_'+DEPTH+'.png'
    fig.savefig(fig_name, dpi=300)
    os.system('convert -trim ' + fig_name + ' ' + fig_name)


    ## Save in French
    if SEASON == 'summer':
        SEASON_TEXT.set_text('été ' +  str(YEAR))
    elif SEASON == 'fall':
         SEASON_TEXT.set_text('automne ' +  str(YEAR))
    elif SEASON == 'spring':
         SEASON_TEXT.set_text('printemps ' +  str(YEAR))
         
    if DEPTH == 'bottom':
        DEPTH_TEXT.set_text('fond')

    fig_name = 'AZMP_OA_'+str(YEAR)+'_'+SEASON+'_'+FIG_VAR+'_'+DEPTH+'_FR.png'
    fig.savefig(fig_name, dpi=300)
    os.system('convert -trim ' + fig_name + ' ' + fig_name)


def seasonal_map_NL(VARIABLE, YEAR, SEASON, DEPTH):

    """

    To produce AZMP BGC ResDoc figures

    Frederic.Cyr@dfo-mpo.gc.ca
    Feb 2023

    """

    # English/French station names
    stationFile = '/home/cyrf0006/github/AZMP-NL/data/STANDARD_SECTIONS.xlsx'
    if SEASON == 'summer':
        SEC_toplot = ['MB', 'SI', 'BB', 'S27', 'FC']
        SEC_toplot_french = [' BM', ' IS', ' BB', 'S27', 'BF']
    elif SEASON == 'fall':
        SEC_toplot = ['BB', 'FC', 'SEGB', 'SPB']
        SEC_toplot_french = [' BB', 'BF', ' GBSE',' BSP']
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
    index_SPB = df.SECTION[df.SECTION=="SOUTHWEST ST PIERRE BANK"].index.tolist()
    index_S27 = df.SECTION[df.SECTION=="STATION 27"].index.tolist()


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
        df = df[(df.index.month>=6) & (df.index.month<=10)]
    elif SEASON == 'fall':
        df = df[(df.index.month>=10) & (df.index.month<=12)]
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
    elif DEPTH == 'bottom':
        df = df.loc[df.groupby('Station_Name')['Depth_(dbar)'].idxmax()] #group by station then pull "min or max depth"
        df = df.loc[df['Depth_(dbar)'] >10] #take all depths >10m (for bottom) to eliminate lone surface samples

    elif DEPTH == 'range':
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
        latLims = [40.5, 52.6]   

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


    ## ---- plot map ---- ##
    print('--- Now plot ---')
    fig = plt.figure(figsize=(7,5))
    ax = fig.add_subplot(111, projection=ccrs.Mercator())
    ax.set_extent([lonLims[0], lonLims[1], latLims[0], latLims[1]], crs=ccrs.PlateCarree())
    #ax.set_extent([-72, -41.5, 40.5, 55.1], crs=ccrs.PlateCarree())
    ax.add_feature(cpf.NaturalEarthFeature('physical', 'coastline', '50m', edgecolor='k', alpha=0.7, linewidth=0.6, facecolor='black'), zorder=15)#cpf.COLORS['land']))
    m=ax.gridlines(linewidth=0.5, color='black', draw_labels=True, alpha=0.5)
    m.top_labels = False
    m.right_labels = False
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
        season_text = ax.text(-42, 57.8, '('+SEASON+')', horizontalalignment='right', color='black', verticalalignment='top', zorder=200,
                fontsize=12, fontweight='bold', transform=ccrs.PlateCarree())    
    else:
         ax.text(-42, 52+.05, str(YEAR), horizontalalignment='right', color='black', zorder=200,
                fontsize=12, fontweight='bold', transform=ccrs.PlateCarree())   
         season_text = ax.text(-42, 52, '('+SEASON+')', horizontalalignment='right', color='black', verticalalignment='top', zorder=200,
                fontsize=12, fontweight='bold', transform=ccrs.PlateCarree())   

    ## plot AZMP section names
    ax_text = SEC_toplot.copy()
    for i, sec in enumerate(SEC_toplot):
        sec = SEC_toplot[i]
        sec_index = eval('index_' + sec)
        if sec == 'FC':
            ax_text[i] = ax.text(stationLon[sec_index][-4], stationLat[sec_index][-4], sec, horizontalalignment='left', verticalalignment='top', color='black', zorder=200, fontsize=12, fontweight='bold', transform=ccrs.PlateCarree())
        elif sec == 'S27':
            ax_text[i] = ax.text(stationLon[sec_index][-8], stationLat[sec_index][-8], sec, horizontalalignment='left', color='black', zorder=200, fontsize=12, fontweight='bold', transform=ccrs.PlateCarree())
        elif sec == 'SPB':
            ax_text[i] = ax.text(stationLon[sec_index][4], stationLat[sec_index][4], ' '+sec, horizontalalignment='left', verticalalignment='top', color='black', zorder=200, fontsize=12, fontweight='bold', transform=ccrs.PlateCarree())
        else:
            ax_text[i] = ax.text(stationLon[sec_index][-1], stationLat[sec_index][-1], ' '+sec, horizontalalignment="left", color="black", zorder=200, fontsize=12, fontweight="bold", transform=ccrs.PlateCarree())

    # Save figure
    fig_name = 'NL_OA_'+str(YEAR)+'_'+SEASON+'_'+FIG_VAR+'_'+DEPTH+'.png'
    fig.savefig(fig_name, dpi=300)
    os.system('convert -trim ' + fig_name + ' ' + fig_name)


    ## Save in French
    if SEASON == 'summer':
        season_text.set_text("(été)")
    elif SEASON == 'fall':
        season_text.set_text("(automne)")

    for i, sec in enumerate(SEC_toplot_french):
        ax_text[i].set_text(SEC_toplot_french[i])

    fig_name = 'NL_OA_'+str(YEAR)+'_'+SEASON+'_'+FIG_VAR+'_'+DEPTH+'_FR.png'
    fig.savefig(fig_name, dpi=300)
    os.system('convert -trim ' + fig_name + ' ' + fig_name)


def section_plot(VAR, YEAR, SEASON, SECTION, ZMAX, REGION='NL', FRENCH=False):

    """

    To produce section plots.
    This is an adaptation of the script 
    cs_section_plot.py

    Possible variables:
    'Omega_Aragonite_(unitless)'
    'pH_Total_(total_scale)'
    'Oxygen_Saturation_(%)'
    'Dissolved_Oxygen_(mL/L)'
    'Phosphate_Concentration_(mmol/m3)'
    'Temperature_(degC)'
    'pH_tot'
    'Omega_Aragonite_(--)'
    'Salinity_(psu)'
    'Total_Alkalinity_(umol/kg)'
    'Inorganic_Carbon_(umol/kg)'
    'AOU'
    ...

    Originally made for work on OA15 paper
    Frederic.Cyr@dfo-mpo.gc.ca
    Sept. 2023

    """

    # File to load
    my_file = '/home/cyrf0006/github/AZMP-NL/datasets/carbonates/AZMP_carbon_data.csv'
    # For colorbar:
    if REGION == 'NL':
        import cc_variable_list_NL as vl
    else:
        import cc_variable_list2 as vl
    # Ajust plotting parameters
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
    section_FR = ['IB', 'MB', 'IS', 'WB', 'BB', 'S27', 'BF', 'GBSE', 'SWSPB']
    SECTION_FR = section_FR[section_EN.index(SECTION)]
    season_EN = ['spring', 'summer', 'fall']
    season_FR = ['printemps', 'été', 'automne']
    SEASON_FR = season_FR[season_EN.index(SEASON)]
    VAR_text = VAR.split('(')[0][0:-1] 
    if VAR == 'Oxygen_Saturation_(%)':
        title = 'Oxygen Saturation for section ' + SECTION + ' - ' + SEASON + ' ' + str(YEAR)
        title_FR = 'Saturation en oxygène pour la section ' + SECTION_FR + ' - ' + SEASON_FR + ' ' + str(YEAR)
    elif VAR == 'pH_Total_(total_scale)':
        title = 'pH Total for section ' + SECTION + ' - ' + SEASON + ' ' + str(YEAR)
        title_FR = 'pH Total pour la section ' + SECTION_FR + ' - ' + SEASON_FR + ' ' + str(YEAR)
    elif VAR == 'Omega_Aragonite_(unitless)':
        title = 'Aragonite saturation state for section ' + SECTION + ' - ' + SEASON + ' ' + str(YEAR)
        title_FR = 'Saturation en aragonite pour la section ' + SECTION_FR + ' - ' + SEASON_FR + ' ' + str(YEAR)
    elif VAR == 'Inorganic_Carbon_(umol/kg)':
        title = 'Inorganic Carbon for section' + SECTION + ' - ' + SEASON + ' ' + str(YEAR)
        title_FR = 'Carbone inorganique pour la section ' + SECTION_FR + ' - ' + SEASON_FR + ' ' + str(YEAR)
    elif VAR == 'Total_Alkalinity_(umol/kg)':
        title = 'Total Alkalinity for section' + SECTION + ' - ' + SEASON + ' ' + str(YEAR)
        title_FR = 'Alkalinité totale pour la section ' + SECTION_FR + ' - ' + SEASON_FR + ' ' + str(YEAR)
    elif VAR == 'AOU':
        VAR_text = 'AOU'
        title = 'AOU for section ' + SECTION + ' - ' + SEASON + ' ' + str(YEAR)
        title_FR = 'UAO pour la section ' + SECTION_FR + ' - ' + SEASON_FR + ' ' + str(YEAR)
    else:
        title=VAR_text
        title_FR=VAR_text
    
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
    # Extract NL data
    df = df[df.Region==REGION]
    # Extract section
    df = df.loc[(df.Station_Name.str.contains(SECTION))]

    # Extract season
    #df = df[df.index.year == YEAR]
    if SEASON == 'spring':
        df = df[df.index.month <= 5]
    elif SEASON == 'summer':
        df = df[(df.index.month>=6) & (df.index.month<=8)]
    elif SEASON == 'fall':
        df = df[df.index.month >= 9]

    # Calculate Apparent Oxygen Utilisation (AOU) if needed
    if VAR == 'AOU':
        SA = gsw.SA_from_SP(df['Salinity_(psu)'], df['Depth_(dbar)'], df['Longitude_(degEast)'], df['Latitude_(degNorth)'])
        CT = gsw.CT_from_t(SA, df['Temperature_(degC)'], df['Depth_(dbar)'])
        O2sol = gsw.O2sol(SA, CT, df['Depth_(dbar)'],  df['Longitude_(degEast)'], df['Latitude_(degNorth)']) # in umol/kg
        O2sol = O2sol/43.570 # in ml/l 
        df['AOU'] = O2sol - df['Dissolved_Oxygen_(mL/L)'] 

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
        df_list.append(df_tmp)

    # Break if empty (no data at all for this season)
    if len(df_list) == 0:
        print('No data for this season [return]')
        return
        
    # Create multi-index
    df_all = pd.concat(df_list, keys=df.index.year.unique(), sort=True)

    # extract year current year
    if YEAR in df_all.index.levels[0]:
        df_year = df_all.xs((YEAR), level=('Timestamp'))
    else:
        print('No data for this year [return]')
        return
        
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
    df_stn = df_stn.drop(['SECTION', 'LONG'], axis=1)
    df_stn = df_stn.rename(columns={'LONG.1': 'LON'})
    df_stn = df_stn.dropna()
    df_stn = df_stn[df_stn.STATION.str.contains(SECTION)]
    df_stn = df_stn.reset_index(drop=True)

    ## ---- Compute distance vector ---- ##
    distance = np.full((df_clim.keys().shape), np.nan)
    lat0 = df_stn[df_stn.index==0]['LAT']
    lon0 = df_stn[df_stn.index==0]['LON']
    for i, stn in enumerate(df_clim.keys().values):
        stn = stn.replace('_','-') # replace in case underscore is used
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
    fig = plt.figure()
    # ax1
    ax = plt.subplot2grid((3, 1), (0, 0))
    c = plt.contourf(distance, df_year.index, df_year, levels, cmap=cmap, extend='both')
    #c_sig1 = plt.contour(distance_sigt, df_sigt_year.columns, df_sigt_year, v_sig, colors='gray', linewidths=1)
    for i in distance:
        plt.plot([i, i], [0, ZMAX], '--k', alpha=.5)
    ax.set_ylim([0, ZMAX])
    ax.set_xlim([0, XLIM])
    #plt.clabel(c_sig1, inline=1, fontsize=10, colors='gray', fmt='%1.1f')
    ax.set_ylabel('Depth (m)', fontweight = 'bold')
    ax.invert_yaxis()
    if REGION == 'NL':
        Bgon = plt.Polygon(bathymetry,color=np.multiply([1,.9333,.6667],.4), alpha=1, zorder=10)
        ax.add_patch(Bgon)
    cb = plt.colorbar(c)
    cb.ax.tick_params(labelsize=8)
    cb.set_label(axis_label, fontsize=12, fontweight='normal')
    ax.xaxis.label.set_visible(False)
    ax.tick_params(labelbottom='off')
    ax.set_title(title)

    # ax2
    ax2 = plt.subplot2grid((3, 1), (1, 0))
    c = plt.contourf(distance, df_clim.index, df_clim, levels, cmap=cmap, extend='both')
    #c_sig2 = plt.contour(distance, df_sigt_clim.columns, df_sigt_clim, v_sig, colors='gray', linewidths=1)
    for i in distance:
        plt.plot([i, i], [0, ZMAX], '--k', alpha=.5)
    ax2.set_ylim([0, ZMAX])
    ax2.set_xlim([0,  XLIM])
    #plt.clabel(c_sig2, inline=1, fontsize=10, colors='gray', fmt='%1.1f')
    ax2.set_ylabel('Depth (m)', fontweight = 'bold')
    ax2.invert_yaxis()
    if REGION == 'NL':
        Bgon = plt.Polygon(bathymetry,color=np.multiply([1,.9333,.6667],.4), alpha=1, zorder=10)
        ax2.add_patch(Bgon)
    cb = plt.colorbar(c)
    cb.ax.tick_params(labelsize=8)
    cb.set_label(axis_label, fontsize=12, fontweight='normal')
    ax2.xaxis.label.set_visible(False)
    ax2.tick_params(labelbottom='off')
    ax2.set_title('Climatology (2014-2020)')

    # ax3
    ax3 = plt.subplot2grid((3, 1), (2, 0))
    #c = plt.contourf(distance, df_anom.index, df_anom, v_anom, cmap=cmo.cm.seismic, extend='both')
    c = plt.contourf(distance, df_anom.index, df_anom, v_anom, cmap=plt.cm.RdBu_r, extend='both')
    ax3.set_ylim([0, ZMAX])
    ax3.set_xlim([0,  XLIM])
    ax3.set_ylabel('Depth (m)', fontweight = 'bold')
    ax3.set_xlabel('Distance (km)', fontweight = 'bold')
    ax3.invert_yaxis()
    if REGION == 'NL':
        Bgon = plt.Polygon(bathymetry,color=np.multiply([1,.9333,.6667],.4), alpha=1, zorder=10)
        ax3.add_patch(Bgon)
    cb = plt.colorbar(c)
    cb.ax.tick_params(labelsize=8)
    cb.set_label(axis_label, fontsize=12, fontweight='normal')
    ax3.set_title(r'Anomaly')

    fig.set_size_inches(w=8,h=12)
    fig_name = 'btl_' + VAR_text + '_' + SECTION + '_' + SEASON + '_' + str(YEAR) + '.png' 
    fig.savefig(fig_name, dpi=200)
    os.system('convert -trim ' + fig_name + ' ' + fig_name)

    ## French figure
    if FRENCH:
        ax.set_title(title_FR)
        ax2.set_title('Climatologie (2014-2020)')
        ax3.set_title(r'Anomalie')
        ax.set_ylabel('Profondeur (m)', fontweight = 'bold')
        ax2.set_ylabel('Profondeur (m)', fontweight = 'bold')
        ax3.set_ylabel('Profondeur (m)', fontweight = 'bold')
        fig.set_size_inches(w=8,h=12)
        fig_name = 'btl_' + VAR_text + '_' + SECTION + '_' + SEASON + '_' + str(YEAR) + '_FR.png' 
        fig.savefig(fig_name, dpi=200)
        os.system('convert -trim ' + fig_name + ' ' + fig_name)

    # For NL ResDoc
    ## montage btl_Omega_Aragonite_SI_summer_2019.png btl_Omega_Aragonite_SI_summer_2020.png -tile 2x1 -geometry +10+10  -background white NL_OA_Omega_Aragonite_SI_2019-2020_summer.png
    ## montage btl_Omega_Aragonite_BB_summer_2019.png btl_Omega_Aragonite_BB_summer_2020.png -tile 2x1 -geometry +10+10  -background white NL_OA_Omega_Aragonite_BB_2019-2020_summer.png
    ## montage btl_Omega_Aragonite_FC_summer_2019.png btl_Omega_Aragonite_FC_summer_2020.png -tile 2x1 -geometry +10+10  -background white NL_OA_Omega_Aragonite_FC_2019-2020_summer.png
    
    ## montage btl_Oxygen_Saturation_SI_summer_2019.png btl_Oxygen_Saturation_SI_summer_2020.png -tile 2x1 -geometry +10+10  -background white NL_OA_Oxygen_Saturation_SI_2019-2020_summer.png
    ## montage btl_Oxygen_Saturation_BB_summer_2019.png btl_Oxygen_Saturation_BB_summer_2020.png -tile 2x1 -geometry +10+10  -background white NL_OA_Oxygen_Saturation_BB_2019-2020_summer.png
    ## montage btl_Oxygen_Saturation_FC_summer_2019.png btl_Oxygen_Saturation_FC_summer_2020.png -tile 2x1 -geometry +10+10  -background white NL_OA_Oxygen_Saturation_FC_2019-2020_summer.png
    
    ## montage btl_pH_Total_SI_summer_2019.png btl_pH_Total_SI_summer_2020.png -tile 2x1 -geometry +10+10  -background white NL_OA_pH_Total_SI_2019-2020_summer.png
    ## montage btl_pH_Total_BB_summer_2019.png btl_pH_Total_BB_summer_2020.png -tile 2x1 -geometry +10+10  -background white NL_OA_pH_Total_BB_2019-2020_summer.png
    ## montage btl_pH_Total_FC_summer_2019.png btl_pH_Total_FC_summer_2020.png -tile 2x1 -geometry +10+10  -background white NL_OA_pH_Total_FC_2019-2020_summer.png
    
    ## In French:
    ## montage btl_Omega_Aragonite_SI_summer_2019_FR.png btl_Omega_Aragonite_SI_summer_2020_FR.png -tile 2x1 -geometry +10+10  -background white NL_OA_Omega_Aragonite_SI_2019-2020_summer_FR.png
    ## montage btl_Omega_Aragonite_BB_summer_2019_FR.png btl_Omega_Aragonite_BB_summer_2020_FR.png -tile 2x1 -geometry +10+10  -background white NL_OA_Omega_Aragonite_BB_2019-2020_summer_FR.png
    ## montage btl_Omega_Aragonite_FC_summer_2019_FR.png btl_Omega_Aragonite_FC_summer_2020_FR.png -tile 2x1 -geometry +10+10  -background white NL_OA_Omega_Aragonite_FC_2019-2020_summer_FR.png
    
    ## montage btl_Oxygen_Saturation_SI_summer_2019_FR.png btl_Oxygen_Saturation_SI_summer_2020_FR.png -tile 2x1 -geometry +10+10  -background white NL_OA_Oxygen_Saturation_SI_2019-2020_summer_FR.png
    ## montage btl_Oxygen_Saturation_BB_summer_2019_FR.png btl_Oxygen_Saturation_BB_summer_2020_FR.png -tile 2x1 -geometry +10+10  -background white NL_OA_Oxygen_Saturation_BB_2019-2020_summer_FR.png
    ## montage btl_Oxygen_Saturation_FC_summer_2019_FR.png btl_Oxygen_Saturation_FC_summer_2020_FR.png -tile 2x1 -geometry +10+10  -background white NL_OA_Oxygen_Saturation_FC_2019-2020_summer_FR.png
    
    ## montage btl_pH_Total_SI_summer_2019_FR.png btl_pH_Total_SI_summer_2020_FR.png -tile 2x1 -geometry +10+10  -background white NL_OA_pH_Total_SI_2019-2020_summer_FR.png
    ## montage btl_pH_Total_BB_summer_2019_FR.png btl_pH_Total_BB_summer_2020_FR.png -tile 2x1 -geometry +10+10  -background white NL_OA_pH_Total_BB_2019-2020_summer_FR.png
    ## montage btl_pH_Total_FC_summer_2019_FR.png btl_pH_Total_FC_summer_2020_FR.png -tile 2x1 -geometry +10+10  -background white NL_OA_pH_Total_FC_2019-2020_summer_FR.png
