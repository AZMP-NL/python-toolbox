"""
TO BE EDITED!!!!!!!

This script reads the master file Excel sheet and pickle a Multi-index DataFrame organized as:

depth                            0    2    4     6     8     10       12   \
Station year variable    season                                             
STN27   2014 temperature spring  NaN  NaN  NaN   NaN   NaN   NaN      NaN   
                         summer  NaN  NaN  NaN   NaN   NaN   NaN      NaN   
                         fall    NaN  NaN  NaN  3.46  3.46  3.46  4.37366   
             salinity    spring  NaN  NaN  NaN   NaN   NaN   NaN      NaN   
                         summer  NaN  NaN  NaN   NaN   NaN   NaN      NaN                            
                                   [...]                    

Multi-index manipulation examples:
# to check indices level names
df_mindex.index.names

df = pd.read_pickle('AZMP_OA_multiIndex.pkl')
# 1. Single station average vertical profile
A = df.xs(('BB-01', 'NO3'),level=('station', 'variable'))
A.groupby(level=0).apply(lambda x: x.mean()).mean()

# 2. Single year section
B = df.xs((2016, 'NO3'),level=('year', 'variable'))
B.groupby(level=0).apply(lambda x: x.mean())

# 3. 1999-2016 section climato
C = df.xs(('NO3'),level=('variable'))
C.groupby(level=0).apply(lambda x: x.mean())

"""
import numpy as np
import pandas as pd
## from scipy.interpolate import griddata
from scipy.interpolate import interp1d  # to remove NaNs in profiles
## from shapely.geometry import Point
## from shapely.geometry.polygon import Polygon
## from shapely.ops import cascaded_union
## from area import area # external fns to compute surface area
## from seawater import extras as swx



## ---- Some parameters ---- ##
infile = '/home/cyrf0006/AZMP/cc/AZMP_OA_CO2stats.xlsx'
z_vec = np.arange(0,350, 5)
varname = pd.Series(['TIC', 'TA', 'pH', 'pCO2', 'Omega_C', 'Omega_A'])
varname.name='variable'

## ----  Load  data ---- ##
df = pd.read_excel(infile)
# Set date as index
df = df.set_index('timestamp')
# Drop other time-related or unnecessary columns
df = df.drop(['timestamp.1', 'Name'], axis=1)
# Keep only targeted section
#df = df[df.section == section]


sname_unique = pd.Series(df['Station-ID'].unique())
sname_unique.name='station'

df_list_station = []
for i, stn in enumerate(sname_unique):

    df_sname = df[df['Station-ID']==stn]
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
            df_season_clean.columns.name='depth'

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


df_mindex = pd.concat(df_list_station,keys=sname_unique)  

#section_mindex.to_pickle('AZMP_OA_multiIndex_' + section +'.pkl')
df_mindex.to_pickle('AZMP_OA_multiIndex.pkl')

