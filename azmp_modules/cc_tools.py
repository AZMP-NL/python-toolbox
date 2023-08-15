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
from shapely.ops import cascaded_union
from area import area # external fns to compute surface area
from seawater import extras as swx

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
sys.path.append('/home/jcoyne/Documents/AZMP-NL_python-toolbox/python-toolbox/cc_tools')
import cc_variable_list as vl

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

    # Only year in revieww
    df = df[df.index.year==YEAR]

    # Season
    if SEASON == 'spring':
        df = df[(df.index.month>=3) & (df.index.month<=6)]
    elif SEASON == 'summer':
        df = df[(df.index.month>=7) & (df.index.month<=9)]
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
    SEASON_TEXT = ax.text(-72, 56, SEASON + ' ' +  str(YEAR), horizontalalignment='left', verticalalignment = 'top', color='white',
        fontsize=12, fontweight='bold', transform=ccrs.PlateCarree())
    DEPTH_TEXT = ax.text(-72, 56, DEPTH, horizontalalignment='left', verticalalignment = 'bottom', color='white',
        fontsize=12, fontweight='bold', transform=ccrs.PlateCarree())

         
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

