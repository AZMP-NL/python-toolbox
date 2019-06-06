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
from mpl_toolkits.basemap import Basemap
from scipy.interpolate import griddata
from scipy.interpolate import interp1d  # to remove NaNs in profiles
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from shapely.ops import cascaded_union
from area import area # external fns to compute surface area
from seawater import extras as swx


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
