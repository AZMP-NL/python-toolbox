'''
This script is an extension of azmp_CIL_stats.py

Instead of running the entire time series (e.g. 1950-2020), this script will only append the most recent year to the existing pickle.

** This script must be run into ~/AZMP/state_reports/sections_plots/CIL **

*** This script will also generate .csv files (e.g., 'CIL_area_BB.csv') necessary for the IROC input

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
import sys

print('BE CAREFUL NOT TO OVERWRITE RESULT FROM azmp_CIL_stats.py')

## ---- Region parameters ---- ## <-------------------------------Would be nice to pass this in a config file '2017.report'
SECTION = 'WB'
SEASON = 'summer'
YEAR = 2021 # year to update
dlat = 1 # how far from station we search
dlon = 1
dz = 1 # vertical bins
dc = .05 # grid resolution
v = np.arange(-2,13,1)
FLAG = True # whether or not it is a valid year

if SECTION == 'SI':
    df_CIL = pd.read_pickle('df_CIL_SI_summer.pkl')
    df_sectionT = pd.read_pickle('df_SI_meanT_summer.pkl')
elif SECTION == 'WB':
    df_CIL = pd.read_pickle('df_CIL_WB_summer.pkl')
    df_sectionT = pd.read_pickle('df_WB_meanT_summer.pkl')
elif SECTION == 'BB':
    df_CIL = pd.read_pickle('df_CIL_BB_summer.pkl')    
    df_sectionT = pd.read_pickle('df_BB_meanT_summer.pkl')
elif SECTION == 'FC':
    df_CIL = pd.read_pickle('df_CIL_FC_summer.pkl')
    df_sectionT = pd.read_pickle('df_FC_meanT_summer.pkl')
    df_sectionT_cap = pd.read_pickle('df_FC_meanT_cap_summer.pkl')
    df_sectionT_shelf = pd.read_pickle('df_FC_meanT_shelf_summer.pkl')
else:
    disp('Check section!')


# IF A year to flag, we just put NaNs (no need to calculate anything)
if FLAG == True:
    # Append to existing df_CIL
    df_update = pd.DataFrame([np.nan, np.nan, np.nan]).T
    df_update.index = [YEAR]
    df_update.columns = ['vol_itp', 'core_itp', 'core_depth_itp']
    df_CIL = df_CIL.append(df_update)  
    # Save CIL timeseries
    picklename = 'df_CIL_' + SECTION + '_' + SEASON + '.pkl'
    df_CIL.to_pickle(picklename)

    # Append section average temp
    df_sectionT_update = pd.Series(np.nan)
    df_sectionT_update.index = [YEAR]
    df_sectionT = df_sectionT.append(df_sectionT_update)
    picklename = 'df_' + SECTION + '_meanT_' + SEASON + '.pkl'
    df_sectionT.to_pickle(picklename)
    sys.exit()
    
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
bathy_file = './' + SECTION + '_bathy.npy'
if os.path.isfile(bathy_file):
    print('Load saved bathymetry!')
    Zitp = np.load(bathy_file)

else:
    print('Bathymetry file not found. Check your path or run azmp_CIL_stats.py')
    sys.exit()

## --------- Get climatology -------- ####
clim_file = './df_' + SECTION + '_itp_clim.pkl'
if os.path.isfile(clim_file):
    print('Load saved climatology!')
    df_section_itp_clim = pd.read_pickle('df_' + SECTION + '_itp_clim.pkl')
    df_section_stn_clim = pd.read_pickle('df_' + SECTION + '_stn_clim.pkl')  
else:
    print('Climatological file not found. Check your path or run azmp_CIL_stats.py')
    sys.exit()

    
## ---- Open netDF file and calculate CIL area for a single year ---- ##
year_file = '/home/cyrf0006/data/dev_database/netCDF/' + str(YEAR) + '.nc'
print('Get ' + year_file)
ds = xr.open_dataset(year_file)

# Remove problematic datasets
ds = ds.where(ds.instrument_ID!='MEDBA', drop=True)
ds = ds.where(ds.instrument_ID!='MEDTE', drop=True)

# Select Region
ds = ds.where((ds.longitude>lonLims[0]) & (ds.longitude<lonLims[1]), drop=True)
ds = ds.where((ds.latitude>latLims[0]) & (ds.latitude<latLims[1]), drop=True)

# Select time (save several options here)
if SEASON == 'summer':
    ds = ds.sel(time=((ds['time.month']>=6)) & ((ds['time.month']<=8)))
elif SEASON == 'spring':
    ds = ds.sel(time=((ds['time.month']>=3)) & ((ds['time.month']<=5)))
elif SEASON == 'fall':
    ds = ds.sel(time=((ds['time.month']>=9)) & ((ds['time.month']<=12)))
else:
    print('!! no season specified, used them all! !!')

# Extract temperature    
da_temp = ds['temperature']
lons = np.array(ds.longitude)
lats = np.array(ds.latitude)
bins = np.arange(dz/2.0, 500, dz)
da_temp = da_temp.groupby_bins('level', bins).mean(dim='level')

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
V = np.full((lat_reg.size, lon_reg.size, z.size), np.nan)

# Aggregate on regular grid
for i, xx in enumerate(lon_reg):
    for j, yy in enumerate(lat_reg):    
        idx = np.where((lons>=xx-dc/2) & (lons<xx+dc/2) & (lats>=yy-dc/2) & (lats<yy+dc/2))
        tmp = np.array(df.iloc[idx].mean(axis=0))
        idx_good = np.argwhere((~np.isnan(tmp)) & (tmp<30))
        if np.size(idx_good)==1:
            V[j,i,:] = np.array(df.iloc[idx].mean(axis=0))
        elif np.size(idx_good)>1: # vertical interpolation between pts
            interp = interp1d(np.squeeze(z[idx_good]), np.squeeze(tmp[idx_good]))
            idx_interp = np.arange(np.int(idx_good[0]),np.int(idx_good[-1]+1))
            V[j,i,idx_interp] = interp(z[idx_interp]) # interpolate 

# horizontal interpolation at each depth
lon_grid, lat_grid = np.meshgrid(lon_reg,lat_reg)
lon_vec = np.reshape(lon_grid, lon_grid.size)
lat_vec = np.reshape(lat_grid, lat_grid.size)
for k, zz in enumerate(z):
    # Meshgrid 1D data (after removing NaNs)
    tmp_grid = V[:,:,k]
    tmp_vec = np.reshape(tmp_grid, tmp_grid.size)
    # griddata (after removing nans)
    idx_good = np.argwhere(~np.isnan(tmp_vec))
    if idx_good.size>10: # At least N pts for interpolation
        LN = np.squeeze(lon_vec[idx_good])
        LT = np.squeeze(lat_vec[idx_good])
        TT = np.squeeze(tmp_vec[idx_good])
        zi = griddata((LN, LT), TT, (lon_grid, lat_grid), method='linear')
        V[:,:,k] = zi
    else:
        continue
print(' -> Done!')    

# mask using bathymetry (I don't think it is necessary, but make nice figures)
for i, xx in enumerate(lon_reg):
    for j,yy in enumerate(lat_reg):
        if Zitp[j,i] > -10: # remove shallower than 10m
            V[j,i,:] = np.nan
                
## ---- Extract section info (test 2 options) ---- ##
temp_coords = np.array([[x,y] for x in lat_reg for y in lon_reg])
VV = V.reshape(temp_coords.shape[0],V.shape[2])
ZZ = Zitp.reshape(temp_coords.shape[0],1)
section_only = []
df_section_itp = pd.DataFrame(index=stn_list, columns=z)
for stn in stn_list:
    # 1. Section only (by station name)
    ds_tmp = ds.where(ds.comments == stn, drop=True)  
    section_only.append(ds_tmp)

    #2.  From interpolated field (closest to station)
    station = df_stn[df_stn.STATION==stn]
    idx_opti = np.argmin(np.sum(np.abs(temp_coords - np.array(list(zip(station.LAT,station.LON)))), axis=1))
    Tprofile = VV[idx_opti,:]
    # remove data below bottom
    bottom_depth = -ZZ[idx_opti]
    Tprofile[z>=bottom_depth]=np.nan
    # store in dataframe
    df_section_itp.loc[stn] = Tprofile

# convert option #1 to dataframe    
ds_section = xr.concat(section_only, dim='time')
da = ds_section['temperature']
df_section_stn = da.to_pandas()
df_section_stn.index = ds_section.comments.values
del ds, da

# Compute distance vector for option #1 - exact station
distance_stn = np.full((df_section_stn.index.shape), np.nan)
for i, stn in enumerate(df_section_stn.index):
    distance_stn[i] = azst.haversine(df_stn.LON[0], df_stn.LAT[0], df_stn[df_stn.STATION==df_section_stn.index[i]].LON, df_stn[df_stn.STATION==df_section_stn.index[i]].LAT)

# Compute distance vector for option #2 - interp field
distance_itp = np.full((df_section_itp.index.shape), np.nan)
for i, stn in enumerate(df_section_itp.index):
    distance_itp[i] = azst.haversine(df_stn.LON[0], df_stn.LAT[0], df_stn[df_stn.STATION==df_section_itp.index[i]].LON, df_stn[df_stn.STATION==df_section_itp.index[i]].LAT)

## ---- Fill NaNs with climatology (new from April 2021) ---- ##
print('Fill NaNs with climatology')
df_section_itp_orig = df_section_itp.copy()
df_fill = df_section_itp_clim[df_section_itp_orig.isna()]
df_section_itp[df_section_itp_orig.isna()] = df_section_itp_clim[df_section_itp_orig.isna()]

## ---- Plot to check the result ---- ##
XLIM = azst.haversine(df_stn.LON[0], df_stn.LAT[0],df_stn.iloc[-1].LON,df_stn.iloc[-1].LAT)
# station-based
if df_section_stn.index.size > 0:
    fig, ax = plt.subplots()
    c = plt.contourf(distance_stn, df_section_stn.columns, df_section_stn.T, v)
    c_cil_stn = plt.contour(distance_stn, df_section_stn.columns, df_section_stn.T, [0,], colors='k', linewidths=2)
    ax.set_ylim([0, 400])
    ax.set_xlim([0,  XLIM])
    ax.set_ylabel('Depth (m)', fontWeight = 'bold')
    ax.set_xlabel('Distance (km)')
    ax.invert_yaxis()
    plt.colorbar(c)
    fig.savefig('CIL_' + SECTION + '_' + str(YEAR) + '_1.png', dpi=150)
    plt.close()
    # CIL area
    cil_vol_stn = 0
    CIL = c_cil_stn.collections[0]
    for path in CIL.get_paths()[:]:
        vs = path.vertices
        cil_vol_stn = cil_vol_stn + np.abs(area(vs))/1000
else:
    cil_vol_stn = np.nan

# interpolation-based
if df_section_itp.index.size > 0:
    fig, ax = plt.subplots()
    c = plt.contourf(distance_itp, df_section_itp.columns, df_section_itp.T, v, alpha=.75)
    cfill = plt.contourf(distance_itp, df_section_itp.columns, df_section_itp_orig.T, v, alpha=.75)
    c_cil_itp = plt.contour(distance_itp, df_section_itp.columns, df_section_itp.T, [0,], colors='k', linewidths=2)
    ax.set_ylim([0, 400])
    ax.set_xlim([0,  XLIM])
    ax.set_ylabel('Depth (m)', fontWeight = 'bold')
    ax.set_xlabel('Distance (km)')
    ax.invert_yaxis()
    plt.colorbar(c)
    fig.savefig('CIL_' + SECTION + '_' + str(YEAR) + '_2.png', dpi=150)
    plt.close()

    ## 3 variables
    # A) CIL area
    cil_vol_itp = 0
    CIL = c_cil_itp.collections[0]
    for path in CIL.get_paths()[:]:
        vs = path.vertices
        cil_vol_itp = cil_vol_itp + np.abs(area(vs))/1000

    if cil_vol_itp != 0:
        # B) CIL core + core depth
        narrow_depth = df_section_itp[df_section_itp.columns[(df_section_itp.columns>=25) & (df_section_itp.columns<=250)]]
        cil_core_itp_update = np.nanmin(narrow_depth.values) # update!
        z_idx = np.where(narrow_depth==np.nanmin(narrow_depth.values))
        cil_coredepth_itp_update = np.array(narrow_depth.columns[z_idx[1].min()]) #update!
        # C) Section mean T
        section_meanT_update = df_section_itp.mean().mean() # update!
        if SECTION == 'FC':
            section_meanT_shelf_update = df_section_itp.iloc[0:20,:].mean().mean() # FC-01 to FC-20
            section_meanT_cap_update = df_section_itp.iloc[14:,:].mean().mean() # FC-15 to FC-38
        plt.close('all')
    else:
       cil_core_itp_update = np.nan
       cil_coredepth_itp_update = np.nan
       section_meanT_update = np.nan
       section_meanT_shelf_update = np.nan
       section_meanT_cap_update = np.nan
else:
    cil_vol_itp = np.nan

CIL_area_update = cil_vol_stn, cil_vol_itp


## ---- Update existing time series  ---- ##

# A) For IROC:
# previous values
outfile = 'CIL_area_' + SECTION + '.csv'
df_iroc = pd.read_csv(outfile)
df_iroc.set_index('Unnamed: 0', inplace=True)
del df_iroc.index.name
# update
iroc_update = pd.DataFrame(CIL_area_update).T
iroc_index = iroc_update.index.tolist()
iroc_index[0] = YEAR
iroc_update.index = iroc_index
iroc_update.columns = ['station-ID', 'interp_field']
# merge and save
df_iroc = df_iroc.append(iroc_update)   
df_iroc.to_csv(outfile)

# B) for AZMP reporting: 
df_update = pd.DataFrame([CIL_area_update[1], cil_core_itp_update, cil_coredepth_itp_update]).T
df_update.index = [YEAR]
df_update.columns = ['vol_itp', 'core_itp', 'core_depth_itp']
df_CIL = df_CIL.append(df_update)  
# Save CIL timeseries
picklename = 'df_CIL_' + SECTION + '_' + SEASON + '.pkl'
df_CIL.to_pickle(picklename)

# C) For NAFO and AZMP
# Append section average temp
df_sectionT_update = pd.Series(section_meanT_update)
df_sectionT_update.index = [YEAR]
df_sectionT = df_sectionT.append(df_sectionT_update)
picklename = 'df_' + SECTION + '_meanT_' + SEASON + '.pkl'
df_sectionT.to_pickle(picklename)
# distinguish shelf/cap when FC
if SECTION == 'FC':
    #FC shelf
    df_shelf_update = pd.Series(section_meanT_shelf_update)
    df_shelf_update.index = [YEAR]
    df_sectionT_shelf = df_sectionT_shelf.append(df_shelf_update)
    picklename = 'df_' + SECTION + '_meanT_shelf_' + SEASON + '.pkl'
    df_sectionT_shelf.to_pickle(picklename)
    # Flemish Cap
    df_cap_update = pd.Series(section_meanT_cap_update)
    df_cap_update.index = [YEAR]
    df_sectionT_cap = df_sectionT_cap.append(df_cap_update)
    picklename = 'df_' + SECTION + '_meanT_cap_' + SEASON + '.pkl'
    df_sectionT_cap.to_pickle(picklename)
    
    
