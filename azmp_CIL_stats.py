'''
CIL area and other parameters along hydrographic sections.


** This script must be run into ~/AZMP/state_reports/sections_plots/CIL **

*** This script will also generate .csv files (e.g., 'CIL_area_BB.csv') necessary for the IROC input (replace iroc_CIL_area.py)

**** When the original .pkl/.csv files are generated, it is posisble to only update the recent years with azmp_CIL_stats_update.py

Final update November 2021
Frederic.Cyr@dfo-mpo.gc.ca


'''
import os
import netCDF4
import xarray as xr
#from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import pandas as pd
import numpy as  np
from scipy.interpolate import interp1d  # to remove NaNs in profiles
from scipy.interpolate import griddata
import azmp_sections_tools as azst


## ---- Region parameters ---- ## <-------------------------------Would be nice to pass this in a config file '2017.report'
SECTION = 'FC'
SEASON = 'fall'
YEAR = [1928, 2023] # For IROC
YEAR = [1950, 2023] # For AZMP reporting
YEAR_MIN = YEAR[0]
YEAR_MAX = YEAR[1]
YEAR_CLIM = [1991, 2020]
dlat = 1 # how far from station we search
dlon = 1
dz = 1 # vertical bins
dc = .05 # grid resolution
v = np.arange(-2,13,1)

# Years to flag
flag_SI_summer = np.array([1928, 1932, 1936, 1937, 1939, 1940, 1941, 1966, 1967, 1972, 1981, 1989])
flag_BB_summer = np.array([1968])
flag_WB_summer = np.arange(1950, 1973)
flag_WB_summer = np.append(flag_WB_summer, [2011, 2019])       
flag_WB_summer = flag_WB_summer[flag_WB_summer != 1960]
flag_BB_summer = flag_BB_summer[(flag_BB_summer>=YEAR_MIN) & (flag_BB_summer<=YEAR_MAX)]
flag_WB_summer = flag_WB_summer[(flag_WB_summer>=YEAR_MIN) & (flag_WB_summer<=YEAR_MAX)]
flag_BB_fall = np.array([1993])

# CIL surface (Note that there is a bias because )
def area(vs):
    a = 0
    x0,y0 = vs[0]
    for [x1,y1] in vs[1:]:
        dx = x1-x0
        dy = y1-y0
        a += 0.5*(y0*dx - x0*dy)
        x0 = x1
        y0 = y1
    return a

## ---- Get Stations ---- ## 
df_stn = pd.read_excel('~/github/AZMP-NL/utils/STANDARD_SECTIONS.xlsx')
df_stn = df_stn.drop(['SECTION', 'LONG'], axis=1)
df_stn = df_stn.rename(columns={'LONG.1': 'LON'})
df_stn = df_stn.dropna()
df_stn = df_stn[df_stn.STATION.str.contains(SECTION+'-')]
df_stn = df_stn.reset_index(drop=True)
stn_list = df_stn.STATION.values

# Derive regular grid
latLims = np.array([df_stn.LAT.min() - dlat, df_stn.LAT.max() + dlat])
lonLims = np.array([df_stn.LON.min() - dlon, df_stn.LON.max() + dlon])
lon_reg = np.arange(lonLims[0]+dc/2, lonLims[1]-dc/2, dc)
lat_reg = np.arange(latLims[0]+dc/2, latLims[1]-dc/2, dc)

## --------- Get Bathymetry -------- ####
bathy_file = 'operation_files/' + SECTION + '_bathy.npy'
if os.path.isfile(bathy_file):
    print('Load saved bathymetry!')
    Zitp = np.load(bathy_file)

else:
    print('Get bathy...')
    dataFile = '~/data/GEBCO/GEBCO_2023_sub_ice_topo.nc' # Maybe find a better way to handle this file


    dataset = xr.open_dataset(dataFile)
    dataset = dataset.isel(lon=(dataset.lon>=lonLims[0])*(dataset.lon<=lonLims[1]))
    dataset = dataset.isel(lat=(dataset.lat>=latLims[0])*(dataset.lat<=latLims[1]))

    # Extract latitude and longitude
    lon = dataset.lon.values
    lat = dataset.lat.values

    # Extract the elevation
    Z = dataset.elevation.values

    #Set up the lats and lons
    lon_grid, lat_grid = np.meshgrid(lon_reg,lat_reg)
    lon_grid_bathy, lat_grid_bathy = np.meshgrid(lon,lat)

    #Interpolate
    Zitp = griddata((lon_grid_bathy.flatten(),lat_grid_bathy.flatten()), Z.flatten(), (lon_grid,lat_grid), method='linear')
    np.save(bathy_file, Zitp)
    print(' -> Done!')
    print(' -> Done!')

#QUESTION FOR FRED: In other scripts, we use the stn clim to fill everything, should this be made different
## --------- Get climatology -------- ####
clim_file = 'operation_files/df_temperature_' + SECTION + '_' + SEASON + '_clim.pkl'
df_section_clim = pd.read_pickle(clim_file)

## --------- Import yearly data -------- ####
stn_file = 'operation_files/df_stn_mindex_T_'+SECTION+'_'+SEASON+'.pkl'
df_stn_mindex = pd.read_pickle(stn_file)
stn_man_file = 'operation_files/df_stn_man_mindex_T_'+SECTION+'_'+SEASON+'.pkl'
df_stn_man_mindex = pd.read_pickle(stn_man_file)
itp_file = 'operation_files/df_itp_mindex_T_'+SECTION+'_'+SEASON+'.pkl'
df_itp_mindex = pd.read_pickle(itp_file)

## --------- Loop on years -------- ####
years = np.arange(YEAR[0], YEAR[1]+1)
years_series = pd.Series(years)
years_series.name='year'

#Set up the variables that will be saved
CIL_area = np.full((years.size,2), np.nan)
cil_vol_itp = np.full(years.shape, np.nan)
cil_core_itp = np.full(years.shape, np.nan)
cil_coredepth_itp = np.full(years.shape, np.nan)
section_meanT = np.full(years.shape, np.nan)
section_meanT_shelf = np.full(years.shape, np.nan)
section_meanT_cap = np.full(years.shape, np.nan)

#Cycle through each year
for idx, YEAR in enumerate(years):

    #Isolate for the data of interest
    stn_years = df_stn_mindex.index.get_level_values('year').values
    df_section_stn = df_stn_mindex.iloc[stn_years == YEAR].droplevel('year')
    itp_years = df_itp_mindex.index.get_level_values('year').values
    df_section_itp = df_itp_mindex.iloc[itp_years == YEAR].droplevel('year')


    #QUESTION FOR FRED: Why was climatology fill only done for interpolation
    #QUESTION FOR FRED: Why was climatology fill only done for stations that were initially present?
    ## ---- Fill NaNs with climatology ---- ##
    if df_section_itp.shape[0] != 0:
        df_section_itp_orig = df_section_itp.copy()
        df_section_itp = df_section_itp.T
        df_section_itp = df_section_itp.reindex(columns=df_section_clim.T.columns, fill_value=np.nan)
        stn_clim = np.array(df_section_clim.T)
        year_data = np.array(df_section_itp)
        year_data[pd.isnull(year_data)] = stn_clim[pd.isnull(year_data)]
        df_section_itp.iloc[:,:] = year_data.astype(float)
        df_section_itp = df_section_itp.T
    else:
        df_section_itp_orig = df_section_itp.copy()

    # Compute distance vector for option #1 - exact station
    distance_stn = np.full((df_section_stn.index.shape), np.nan)
    for i, stn in enumerate(df_section_stn.index):
        distance_stn[i] = azst.haversine(
            df_stn.LON[0],
            df_stn.LAT[0],
            df_stn[df_stn.STATION==df_section_stn.index[i]].LON.values,
            df_stn[df_stn.STATION==df_section_stn.index[i]].LAT.values)

    # Compute distance vector for option #2 - interp field
    distance_itp = np.full((df_section_itp.index.shape), np.nan)
    for i, stn in enumerate(df_section_itp.index):
        distance_itp[i] = azst.haversine(
            df_stn.LON[0],
            df_stn.LAT[0],
            df_stn[df_stn.STATION==df_section_itp.index[i]].LON.values,
            df_stn[df_stn.STATION==df_section_itp.index[i]].LAT.values)
    distance_itp_orig = np.full((df_section_itp_orig.index.shape), np.nan)
    for i, stn in enumerate(df_section_itp_orig.index):
        distance_itp_orig[i] = azst.haversine(
            df_stn.LON[0],
            df_stn.LAT[0],
            df_stn[df_stn.STATION==df_section_itp_orig.index[i]].LON.values,
            df_stn[df_stn.STATION==df_section_itp_orig.index[i]].LAT.values)


    ## ---- Plot to check the result ---- ##
    XLIM = azst.haversine(df_stn.LON[0], df_stn.LAT[0],df_stn.iloc[-1].LON,df_stn.iloc[-1].LAT)
    #Station based
    if df_section_stn.index.size > 1:
        fig, ax = plt.subplots()
        c = plt.contourf(distance_stn, df_section_stn.columns, df_section_stn.T, v)
        c_cil_stn = plt.contour(distance_stn, df_section_stn.columns, df_section_stn.T, [0,], colors='k', linewidths=2)
        ax.set_ylim([0, 400])
        ax.set_xlim([0,  XLIM])
        ax.set_ylabel('Depth (m)', fontweight = 'bold')
        ax.set_xlabel('Distance (km)')
        ax.invert_yaxis()
        plt.colorbar(c)
        fig.savefig('CIL_' + SECTION + '_' + str(YEAR) + '_1.png', dpi=150)
        plt.close()
        # CIL area
        cil_vol_stn = 0
        for path in c_cil_stn.get_paths()[:]:
            vs = path.vertices
            if vs.shape[0] != 0:
                cil_vol_stn = cil_vol_stn + np.abs(area(vs))/1000
    else:
        cil_vol_stn = np.nan

    #Interpolation based
    if df_section_itp.index.size > 0:
        fig, ax = plt.subplots()
        c = plt.contourf(distance_itp, df_section_itp.columns, df_section_itp.T, v, alpha=.75)
        cfill = plt.contourf(distance_itp_orig, df_section_itp_orig.columns, df_section_itp_orig.T, v, alpha=.75)
        c_cil_itp = plt.contour(distance_itp, df_section_itp.columns, df_section_itp.T, [0,], colors='k', linewidths=2)
        ax.set_ylim([0, 400])
        ax.set_xlim([0,  XLIM])
        ax.set_ylabel('Depth (m)', fontweight = 'bold')
        ax.set_xlabel('Distance (km)')
        ax.invert_yaxis()
        plt.colorbar(c)
        fig.savefig('CIL_' + SECTION + '_' + str(YEAR) + '_2.png', dpi=150)
        plt.close()

        ## 3 variables
        # A) CIL area
        cil_vol_itp = 0
        for path in c_cil_itp.get_paths()[:]:
            vs = path.vertices
            if vs.shape[0] != 0:
                cil_vol_itp = cil_vol_itp + np.abs(area(vs))/1000

        if cil_vol_itp != 0:
            # B) CIL core + core depth
            narrow_depth = df_section_itp[df_section_itp.columns[(df_section_itp.columns>=25) & (df_section_itp.columns<=250)]]
            cil_core_itp[idx] = np.nanmin(narrow_depth.values)
            z_idx = np.where(narrow_depth==np.nanmin(narrow_depth.values))
            cil_coredepth_itp[idx] = np.array(narrow_depth.columns[z_idx[1].min()])

            # C) Section mean T
            section_meanT[idx] = df_section_itp.mean().mean()
            if SECTION == 'FC':
                section_meanT_shelf[idx] = df_section_itp.iloc[0:20,:].mean().mean() # FC-01 to FC-20
                section_meanT_cap[idx] = df_section_itp.iloc[14:,:].mean().mean() # FC-15 to FC-38
            plt.close('all')
        else:
           cil_core_itp[idx] = np.nan
           cil_coredepth_itp[idx] = np.nan
           section_meanT[idx] = np.nan
           section_meanT_shelf[idx] = np.nan
           section_meanT_cap[idx] = np.nan
    else:
        cil_vol_itp = np.nan
        
    CIL_area[idx,:] = cil_vol_stn, cil_vol_itp
    plt.close('all')
    print(' ->  '+str(YEAR)+' done! ')

# Merge CIL timeseries
# A) For IROC:
df = pd.DataFrame(CIL_area)
df.index = years

# B) for AZMP reporting: 
df_CIL= pd.DataFrame([CIL_area[:,1], cil_core_itp, cil_coredepth_itp]).T
df_CIL.index = years_series
df_CIL.columns = ['vol_itp', 'core_itp', 'core_depth_itp']

# Manually flag some years (check section plots for justification) 
if (SECTION=='BB') & (SEASON=='summer'):
    for i in flag_BB_summer:        
        df.loc[i]=np.nan
        df_CIL.loc[i]=np.nan
elif (SECTION=='WB') & (SEASON=='summer'):
    for i in flag_WB_summer:        
        df.loc[i]=np.nan
        df_CIL.loc[i]=np.nan
elif (SECTION=='SI') & (SEASON=='summer'):
    for i in flag_SI_summer:        
        df.loc[i]=np.nan
        df_CIL.loc[i]=np.nan

# Save IROC datadata in csv        
df.columns = ['station-ID', 'interp_field']
outfile = 'CIL_area_' + SECTION + '.csv'
df.to_csv(outfile)

# Save CIL timeseries
picklename = 'df_CIL_' + SECTION + '_' + SEASON + '.pkl'
df_CIL.to_pickle(picklename)

'''
# Save section average temp
df_sectionT = pd.Series(section_meanT)
df_sectionT.index = years_series
picklename = 'df_' + SECTION + '_meanT_' + SEASON + '.pkl'
df_sectionT.to_pickle(picklename)
if SECTION == 'FC':
    df_sectionT = pd.Series(section_meanT_shelf)
    df_sectionT.index = years_series
    picklename = 'df_' + SECTION + '_meanT_shelf_' + SEASON + '.pkl'
    df_sectionT.to_pickle(picklename)
    df_sectionT = pd.Series(section_meanT_cap)
    df_sectionT.index = years_series
    picklename = 'df_' + SECTION + '_meanT_cap_' + SEASON + '.pkl'
    df_sectionT.to_pickle(picklename)
'''


