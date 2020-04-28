'''
CIL area

** This script must be run into ~/AZMP/state_reports/sections_plots/CIL **


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
SECTION = 'BB'
SEASON = 'summer'
CLIM_YEAR = [1950, 2019]
dlat = 2 # how far from station we search
dlon = 2
dz = 1 # vertical bins
dc = .2 # grid resolution

# Years to flag
flag_BB_summer = [1982]
flag_WB_summer = [1953, 1956, 1959, 1982, 2019]
flag_SI_summer = [1989]

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
df_stn = pd.read_excel('/home/cyrf0006/github/AZMP-NL/data/STANDARD_SECTIONS.xlsx')
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
print('Get bathy...')
dataFile = '/home/cyrf0006/data/GEBCO/GEBCO_2014_1D.nc' # Maybe find a better way to handle this file
lon_grid, lat_grid = np.meshgrid(lon_reg,lat_reg)
# Load data
dataset = netCDF4.Dataset(dataFile)
x = [-179-59.75/60, 179+59.75/60] # to correct bug in 30'' dataset?
y = [-89-59.75/60, 89+59.75/60]
spacing = dataset.variables['spacing']
# Compute Lat/Lon
nx = int((x[-1]-x[0])/spacing[0]) + 1  # num pts in x-dir
ny = int((y[-1]-y[0])/spacing[1]) + 1  # num pts in y-dir
lon = np.linspace(x[0],x[-1],nx)
lat = np.linspace(y[0],y[-1],ny)
# interpolate data on regular grid (temperature grid)
# Reshape data
zz = dataset.variables['z']
Z = zz[:].reshape(ny, nx)
Z = np.flipud(Z) # <------------ important!!!
# Reduce data according to Region params
idx_lon = np.where((lon>=lonLims[0]) & (lon<=lonLims[1]))
idx_lat = np.where((lat>=latLims[0]) & (lat<=latLims[1]))
Z = Z[idx_lat[0][0]:idx_lat[0][-1]+1, idx_lon[0][0]:idx_lon[0][-1]+1]
lon = lon[idx_lon[0]]
lat = lat[idx_lat[0]]
# interpolate data on regular grid (temperature grid)
lon_grid_bathy, lat_grid_bathy = np.meshgrid(lon,lat)
lon_vec_bathy = np.reshape(lon_grid_bathy, lon_grid_bathy.size)
lat_vec_bathy = np.reshape(lat_grid_bathy, lat_grid_bathy.size)
z_vec = np.reshape(Z, Z.size)
Zitp = griddata((lon_vec_bathy, lat_vec_bathy), z_vec, (lon_grid, lat_grid), method='linear')
del zz, Z
print(' -> Done!')



## Loop on climatological years
years = np.arange(CLIM_YEAR[0], CLIM_YEAR[1]+1)
years_series = pd.Series(years)
years_series.name='year'

cil_vol_itp_clim = np.full(years.shape, np.nan)
cil_core_itp_clim = np.full(years.shape, np.nan)
cil_coredepth_itp_clim = np.full(years.shape, np.nan)
section_meanT = np.full(years.shape, np.nan)
section_meanT_shelf = np.full(years.shape, np.nan)
section_meanT_cap = np.full(years.shape, np.nan)
df_itp_temp = []
for idx, YEAR in enumerate(years):
    ## -------- Get CTD data -------- ##
    year_file = '/home/cyrf0006/data/dev_database/netCDF/' + str(YEAR) + '.nc'
    print('Get ' + year_file)
    ds = xr.open_dataset(year_file)

    # Remame problematic datasets
    print('!!Remove MEDBA & MEDTE data!!')
    print('  ---> I Should be improme because I remove good data!!!!')
    ds = ds.where(ds.instrument_ID!='MEDBA', drop=True)
    ds = ds.where(ds.instrument_ID!='MEDTE', drop=True)

    # Select Region
    ds = ds.where((ds.longitude>lonLims[0]) & (ds.longitude<lonLims[1]), drop=True)
    ds = ds.where((ds.latitude>latLims[0]) & (ds.latitude<latLims[1]), drop=True)

    # Select time (save several options here)
    if SEASON == 'summer':
        ds = ds.sel(time=((ds['time.month']>=7)) & ((ds['time.month']<=9)))
    elif SEASON == 'spring':
        ds = ds.sel(time=((ds['time.month']>=3)) & ((ds['time.month']<=5)))
    elif SEASON == 'fall':
        ds = ds.sel(time=((ds['time.month']>=10)) & ((ds['time.month']<=12)))
    else:
        print('!! no season specified, used them all! !!')

    # Extract temperature    
    da_temp = ds['temperature']
    lons = np.array(ds.longitude)
    lats = np.array(ds.latitude)
    bins = np.arange(dz/2.0, 500, dz)
    da_temp = da_temp.groupby_bins('level', bins).mean(dim='level')
    del ds
    
    # 1. Temperature to Pandas Dataframe
    print('Process temperature')
    df = da_temp.to_pandas()
    df.columns = bins[0:-1] #rename columns with 'bins'
    # Remove empty rows.
    idx_empty_rows = df.isnull().all(1).values.nonzero()[0]
    df = df.dropna(axis=0,how='all')
    lons = np.delete(lons,idx_empty_rows)
    lats = np.delete(lats,idx_empty_rows)

    ## --- fill 3D cube --- ##  
    z = df.columns.values
    V_temp = np.full((lat_reg.size, lon_reg.size, z.size), np.nan)

    # Aggregate on regular grid
    for i, xx in enumerate(lon_reg):
        for j, yy in enumerate(lat_reg):    
            idx_coords = np.where((lons>=xx-dc/2) & (lons<xx+dc/2) & (lats>=yy-dc/2) & (lats<yy+dc/2))
            tmp = np.array(df.iloc[idx_coords].mean(axis=0))
            idx_good = np.argwhere((~np.isnan(tmp)) & (tmp<30))
            if np.size(idx_good)==1:
                V_temp[j,i,:] = np.array(df.iloc[idx_coords].mean(axis=0))
            elif np.size(idx_good)>1: # vertical interpolation between pts
                interp = interp1d(np.squeeze(z[idx_good]), np.squeeze(tmp[idx_good]))  # <---- Attention: strange syntax but works
                idx_interp = np.arange(np.int(idx_good[0]),np.int(idx_good[-1]+1))
                V_temp[j,i,idx_interp] = interp(z[idx_interp]) # interpolate only where possible (1st to last good idx)

    # horizontal interpolation at each depth
    lon_grid, lat_grid = np.meshgrid(lon_reg,lat_reg)
    lon_vec = np.reshape(lon_grid, lon_grid.size)
    lat_vec = np.reshape(lat_grid, lat_grid.size)
    for k, zz in enumerate(z):
        # Meshgrid 1D data (after removing NaNs)
        tmp_grid = V_temp[:,:,k]
        tmp_vec = np.reshape(tmp_grid, tmp_grid.size)
        # griddata (after removing nans)
        idx_good = np.argwhere(~np.isnan(tmp_vec))
        if idx_good.size>5: # will ignore depth where no data exist
            LN = np.squeeze(lon_vec[idx_good])
            LT = np.squeeze(lat_vec[idx_good])
            TT = np.squeeze(tmp_vec[idx_good])
            zi = griddata((LN, LT), TT, (lon_grid, lat_grid), method='linear')
            V_temp[:,:,k] = zi
        else:
            continue
    print(' -> Done!')    

    
    # mask using bathymetry (I don't think it is necessary, but make nice figures)
    for i, xx in enumerate(lon_reg):
        for j,yy in enumerate(lat_reg):
            if Zitp[j,i] > -10: # remove shallower than 10m
                V_temp[j,i,:] = np.nan

    
    ## ---- Extract section info (test 2 options) ---- ##
    ## Temperature
    temp_coords = np.array([[x,y] for x in lat_reg for y in lon_reg])
    VT = V_temp.reshape(temp_coords.shape[0],V_temp.shape[2])
    ZZ = Zitp.reshape(temp_coords.shape[0],1)
    section_only = []
    df_section_itp = pd.DataFrame(index=stn_list, columns=z)
    for stn in stn_list:
        #2.  From interpolated field (closest to station) 
        station = df_stn[df_stn.STATION==stn]
        idx_opti = np.argmin(np.sum(np.abs(temp_coords - np.array(list(zip(station.LAT.values,station.LON.values)))), axis=1))
        Tprofile = VT[idx_opti,:]
        # remove data below bottom
        bottom_depth = -ZZ[idx_opti]
        Tprofile[z>=bottom_depth]=np.nan
        # store in dataframe
        df_section_itp.loc[stn] = Tprofile

    ## CIL Calculation
    # Compute distance vector for option #2 - interp field
    distance_itp = np.full((df_section_itp.index.shape), np.nan)
    for i, stn in enumerate(df_section_itp.index):
        distance_itp[i] = azst.haversine(df_stn.LON[0], df_stn.LAT[0], df_stn[df_stn.STATION==df_section_itp.index[i]].LON, df_stn[df_stn.STATION==df_section_itp.index[i]].LAT)

    # Option 2 only
    #if df_section_itp.index.size > 0:
    #if df_section_itp.dropna(how='all').size > 0:
    if df_section_itp.dropna(axis=1, how='all').shape[1] > 1: # e.g., WB 1982, only surface obs are available
        fig, ax = plt.subplots()
        c = plt.contourf(distance_itp, df_section_itp.columns, df_section_itp.T)
        c_cil_itp = plt.contour(distance_itp, df_section_itp.columns, df_section_itp.T, [0,], colors='k', linewidths=2)
        ax.set_ylabel('Depth (m)', fontWeight = 'bold')
        ax.set_xlabel('Distance (km)')
        ax.invert_yaxis()
        plt.title(str(YEAR))
        plt.colorbar(c)
        fig_name = 'temp_section' + SECTION + '_' + SEASON + '_' + '_' + str(YEAR) + '_2.png'
        fig.savefig(fig_name, dpi=150)
        fig.clf()
        plt.close('all')
        # CIL area
        cil_vol_itp = 0
        CIL = c_cil_itp.collections[0]
        for path in CIL.get_paths()[:]:
            vs = path.vertices
            cil_vol_itp = cil_vol_itp + np.abs(area(vs))/1000
        cil_vol_itp_clim[idx] = cil_vol_itp 
        # CIL core
        #cil_core_itp_clim[idx] = np.nanmin(df_section_itp.values)
        narrow_depth = df_section_itp[df_section_itp.columns[(df_section_itp.columns>=25) & (df_section_itp.columns<=250)]]
        cil_core_itp_clim[idx] = np.nanmin(narrow_depth.values)
        z_idx = np.where(narrow_depth==np.nanmin(narrow_depth.values))
        cil_coredepth_itp_clim[idx] = np.array(narrow_depth.columns[z_idx[1].min()])
        # Section mean T
        section_meanT[idx] = df_section_itp.mean().mean()
        if SECTION == 'FC':
            section_meanT_shelf[idx] = df_section_itp.iloc[0:20,:].mean().mean() # FC-01 to FC-20
            section_meanT_cap[idx] = df_section_itp.iloc[14:,:].mean().mean() # FC-15 to FC-38

# Merge CIL timeseries
df_CIL= pd.DataFrame([cil_vol_itp_clim, cil_core_itp_clim, cil_coredepth_itp_clim]).T
df_CIL.index = years_series
df_CIL.columns = ['vol_itp', 'core_itp', 'core_depth_itp']
# Manually flag some years (check section plots for justification) 
if (SECTION=='SI') & (SEASON=='summer'):
    df_CIL.loc[flag_SI_summer]=np.nan
if (SECTION=='BB') & (SEASON=='summer'):
    df_CIL.loc[flag_BB_summer]=np.nan
if (SECTION=='WB') & (SEASON=='summer'):
    df_CIL.loc[flag_WB_summer]=np.nan
# Save CIL timeseries
picklename = 'df_CIL_' + SECTION + '_' + SEASON + '.pkl'
df_CIL.to_pickle(picklename)
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
    
