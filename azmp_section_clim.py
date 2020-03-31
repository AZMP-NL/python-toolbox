'''
CIL area

This script is still in progress...

 Example how to read MIndex:
    df_itp_mindex = pd.read_pickle('df_itp_mindex.pkl')
    C = df_itp_mindex.xs(('year'),level=(0))
    C.groupby(level=0).apply(lambda x: x.mean())

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
CLIM_YEAR = [1981, 2010]
dlat = 2 # how far from station we search
dlon = 2
z1 = 2
dz = 5 # vertical bins
dc = .2 # grid resolution

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
print(' -> Done!')



## Loop on climatological years
years = np.arange(CLIM_YEAR[0], CLIM_YEAR[1]+1)
years_series = pd.Series(years)
years_series.name='year'

cil_vol_stn_clim = np.full(years.shape, np.nan)
cil_vol_itp_clim = np.full(years.shape, np.nan)
cil_core_stn_clim = np.full(years.shape, np.nan)
cil_core_itp_clim = np.full(years.shape, np.nan)
df_stn_temp = []
df_itp_temp = []
df_stn_sal = []
df_itp_sal = []
df_stn_sig = []
df_itp_sig = []
for idx, YEAR in enumerate(years):
    ## -------- Get CTD data -------- ##
    year_file = '/home/cyrf0006/data/dev_database/netCDF/' + str(YEAR) + '.nc'
    print('Get ' + year_file)
    ds = xr.open_mfdataset(year_file)

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
        ds = ds.sel(time=((ds['time.month']>=6)) & ((ds['time.month']<=9)))
    elif SEASON == 'spring':
        ds = ds.sel(time=((ds['time.month']>=3)) & ((ds['time.month']<=5)))
    elif SEASON == 'fall':
        ds = ds.sel(time=((ds['time.month']>=10)) & ((ds['time.month']<=12)))
    else:
        print('!! no season specified, used them all! !!')

    # Extract temperature    
    da_temp = ds['temperature']
    da_sal = ds['salinity']
    da_sig = ds['sigma-t']
    lons = np.array(ds.longitude)
    lats = np.array(ds.latitude)
    bins = np.arange(z1, 2000, dz)
    da_temp = da_temp.groupby_bins('level', bins).mean(dim='level')
    da_sal = da_sal.groupby_bins('level', bins).mean(dim='level')
    da_sig = da_sig.groupby_bins('level', bins).mean(dim='level')

    # 1. Temperature to Pandas Dataframe
    print('Process temperature')
    df = da_temp.to_pandas()
    df.columns = bins[0:-1] #rename columns with 'bins'
    # Remove empty columns.
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
            if (np.unique(LT).size == 1) | (np.unique(LN).size == 1): # May happend for sampling along 47N section (cannot grid single latitude)
                zi = griddata((LN, LT), TT, (lon_grid, lat_grid), method='nearest')
            else:
                zi = griddata((LN, LT), TT, (lon_grid, lat_grid), method='linear')
            V_temp[:,:,k] = zi
        else:
            continue
    print(' -> Done!')    

    
    # 2. Salinity to Pandas Dataframe
    print('Process salinity')
    df = da_sal.to_pandas()
    df.columns = bins[0:-1] #rename columns with 'bins'
    # Remove empty columns
    idx_empty_rows = df.isnull().all(1).values.nonzero()[0]
    df = df.dropna(axis=0,how='all')
    if np.size(idx_empty_rows):
        lons = np.delete(lons,idx_empty_rows)
        lats = np.delete(lats,idx_empty_rows)

    ## --- fill 3D cube --- ##  
    z = df.columns.values
    V_sal = np.full((lat_reg.size, lon_reg.size, z.size), np.nan)

    # Aggregate on regular grid
    for i, xx in enumerate(lon_reg):
        for j, yy in enumerate(lat_reg):    
            idx_coords = np.where((lons>=xx-dc/2) & (lons<xx+dc/2) & (lats>=yy-dc/2) & (lats<yy+dc/2))
            tmp = np.array(df.iloc[idx_coords].mean(axis=0))
            idx_good = np.argwhere((~np.isnan(tmp)))
            if np.size(idx_good)==1:
                V_sal[j,i,:] = np.array(df.iloc[idx_coords].mean(axis=0))
            elif np.size(idx_good)>1: # vertical interpolation between pts
                interp = interp1d(np.squeeze(z[idx_good]), np.squeeze(tmp[idx_good]))  # <---- Attention: strange syntax but works
                idx_interp = np.arange(np.int(idx_good[0]),np.int(idx_good[-1]+1))
                V_sal[j,i,idx_interp] = interp(z[idx_interp]) # interpolate only where possible (1st to last good idx)

    # horizontal interpolation at each depth
    lon_grid, lat_grid = np.meshgrid(lon_reg,lat_reg)
    lon_vec = np.reshape(lon_grid, lon_grid.size)
    lat_vec = np.reshape(lat_grid, lat_grid.size)
    for k, zz in enumerate(z):
        # Meshgrid 1D data (after removing NaNs)
        tmp_grid = V_sal[:,:,k]
        tmp_vec = np.reshape(tmp_grid, tmp_grid.size)
        #print 'interpolate depth layer ' + np.str(k) + ' / ' + np.str(z.size) 
        # griddata (after removing nans)
        idx_good = np.argwhere(~np.isnan(tmp_vec))
        if idx_good.size>7: # will ignore depth where no data exist
            LN = np.squeeze(lon_vec[idx_good])
            LT = np.squeeze(lat_vec[idx_good])
            TT = np.squeeze(tmp_vec[idx_good])
            if (np.unique(LT).size == 1) | (np.unique(LN).size == 1): # May happend for sampling along 47N section (cannot grid single latitude)
                zi = griddata((LN, LT), TT, (lon_grid, lat_grid), method='nearest')
            else:
                zi = griddata((LN, LT), TT, (lon_grid, lat_grid), method='linear')
            V_sal[:,:,k] = zi
        else:
            continue
    print(' -> Done!')   

    # 3. Sigma-t to Pandas Dataframe
    print('Process sigma-t')
    df = da_sig.to_pandas()
    df.columns = bins[0:-1] #rename columns with 'bins'
    # Remove empty columns
    idx_empty_rows = df.isnull().all(1).values.nonzero()[0]
    df = df.dropna(axis=0,how='all')
    if np.size(idx_empty_rows):
        lons = np.delete(lons,idx_empty_rows)
        lats = np.delete(lats,idx_empty_rows)

    ## --- fill 3D cube --- ##  
    z = df.columns.values
    V_sig = np.full((lat_reg.size, lon_reg.size, z.size), np.nan)

    # Aggregate on regular grid
    for i, xx in enumerate(lon_reg):
        for j, yy in enumerate(lat_reg):    
            idx_coords = np.where((lons>=xx-dc/2) & (lons<xx+dc/2) & (lats>=yy-dc/2) & (lats<yy+dc/2))
            tmp = np.array(df.iloc[idx_coords].mean(axis=0))
            idx_good = np.argwhere((~np.isnan(tmp)))
            if np.size(idx_good)==1:
                V_sig[j,i,:] = np.array(df.iloc[idx_coords].mean(axis=0))
            elif np.size(idx_good)>1: # vertical interpolation between pts
                interp = interp1d(np.squeeze(z[idx_good]), np.squeeze(tmp[idx_good]))  # <---- Attention: strange syntax but works
                idx_interp = np.arange(np.int(idx_good[0]),np.int(idx_good[-1]+1))
                V_sig[j,i,idx_interp] = interp(z[idx_interp]) # interpolate only where possible (1st to last good idx)

    # horizontal interpolation at each depth
    lon_grid, lat_grid = np.meshgrid(lon_reg,lat_reg)
    lon_vec = np.reshape(lon_grid, lon_grid.size)
    lat_vec = np.reshape(lat_grid, lat_grid.size)
    for k, zz in enumerate(z):
        # Meshgrid 1D data (after removing NaNs)
        tmp_grid = V_sig[:,:,k]
        tmp_vec = np.reshape(tmp_grid, tmp_grid.size)
        #print 'interpolate depth layer ' + np.str(k) + ' / ' + np.str(z.size) 
        # griddata (after removing nans)
        idx_good = np.argwhere(~np.isnan(tmp_vec))
        if idx_good.size>7: # will ignore depth where no data exist
            LN = np.squeeze(lon_vec[idx_good])
            LT = np.squeeze(lat_vec[idx_good])
            TT = np.squeeze(tmp_vec[idx_good])
            if (np.unique(LT).size == 1) | (np.unique(LN).size == 1): # May happend for sampling along 47N section (cannot grid single latitude)
                zi = griddata((LN, LT), TT, (lon_grid, lat_grid), method='nearest')
            else:
                zi = griddata((LN, LT), TT, (lon_grid, lat_grid), method='linear')
            V_sig[:,:,k] = zi
        else:
            continue
    print(' -> Done!')
    
    # mask using bathymetry (I don't think it is necessary, but make nice figures)
    for i, xx in enumerate(lon_reg):
        for j,yy in enumerate(lat_reg):
            if Zitp[j,i] > -10: # remove shallower than 10m
                V_temp[j,i,:] = np.nan
                V_sal[j,i,:] = np.nan
                V_sig[j,i,:] = np.nan

    
    ## ---- Extract section info (test 2 options) ---- ##
    ## Temperature & Salinity
    temp_coords = np.array([[x,y] for x in lat_reg for y in lon_reg])
    VT = V_temp.reshape(temp_coords.shape[0],V_temp.shape[2])
    VS = V_sal.reshape(temp_coords.shape[0],V_temp.shape[2])
    VSi = V_sig.reshape(temp_coords.shape[0],V_temp.shape[2])
    ZZ = Zitp.reshape(temp_coords.shape[0],1)
    section_only = []
    df_section_itp = pd.DataFrame(index=stn_list, columns=z)
    df_section_itp_S = pd.DataFrame(index=stn_list, columns=z)
    df_section_itp_Si = pd.DataFrame(index=stn_list, columns=z)
    for stn in stn_list:
        # 1. Section only (by station name)
        ds_tmp = ds.where(ds.comments == stn, drop=True)  
        section_only.append(ds_tmp)

        #2.  From interpolated field (closest to station) 
        station = df_stn[df_stn.STATION==stn]
        idx_opti = np.argmin(np.sum(np.abs(temp_coords - np.array(list(zip(station.LAT.values,station.LON.values)))), axis=1))
        Tprofile = VT[idx_opti,:]
        Sprofile = VS[idx_opti,:]
        Siprofile = VSi[idx_opti,:]
        # remove data below bottom
        bottom_depth = -ZZ[idx_opti]
        Tprofile[z>=bottom_depth]=np.nan
        Sprofile[z>=bottom_depth]=np.nan
        Siprofile[z>=bottom_depth]=np.nan

        # store in dataframe
        df_section_itp.loc[stn] = Tprofile
        df_section_itp_S.loc[stn] = Sprofile
        df_section_itp_Si.loc[stn] = Siprofile

    # convert option #1 to dataframe    
    ds_section = xr.concat(section_only, dim='time')
    da = ds_section['temperature']
    df_section_stn = da.to_pandas()
    df_section_stn.index = ds_section.comments.values
    da = ds_section['salinity']
    df_section_stn_S = da.to_pandas()
    df_section_stn_S.index = ds_section.comments.values
    da = ds_section['sigma-t']
    df_section_stn_Si = da.to_pandas()
    df_section_stn_Si.index = ds_section.comments.values
    
    # drop indices with only NaNs (option #2)
    df_section_itp = df_section_itp.dropna(axis=0, how='all')
    df_section_itp_S = df_section_itp_S.dropna(axis=0, how='all')
    df_section_itp_Si = df_section_itp_Si.dropna(axis=0, how='all')

    
    ## CIL Calculation
    # Compute distance vector for option #1 - exact station
    distance_stn = np.full((df_section_stn.index.shape), np.nan)
    for i, stn in enumerate(df_section_stn.index):
        distance_stn[i] = azst.haversine(df_stn.LON[0], df_stn.LAT[0], df_stn[df_stn.STATION==df_section_stn.index[i]].LON, df_stn[df_stn.STATION==df_section_stn.index[i]].LAT)    
    # Compute distance vector for option #2 - interp field
    distance_itp = np.full((df_section_itp.index.shape), np.nan)
    for i, stn in enumerate(df_section_itp.index):
        distance_itp[i] = azst.haversine(df_stn.LON[0], df_stn.LAT[0], df_stn[df_stn.STATION==df_section_itp.index[i]].LON, df_stn[df_stn.STATION==df_section_itp.index[i]].LAT)
    if df_section_stn.index.size > 1:
        fig, ax = plt.subplots()
        c = plt.contourf(distance_stn, df_section_stn.columns, df_section_stn.T)
        c_cil_stn = plt.contour(distance_stn, df_section_stn.columns, df_section_stn.T, [0,], colors='k', linewidths=2)
        ax.set_ylabel('Depth (m)', fontWeight = 'bold')
        ax.set_xlabel('Distance (km)')
        ax.invert_yaxis()
        plt.colorbar(c)
        plt.title(str(YEAR))
        fig_name = 'temp_section' + SECTION + '_' + SEASON + '_' + '_' + str(YEAR) + '_1.png'
        fig.savefig(fig_name, dpi=150)
        #plt.clf()
        plt.close('all')
        #fig.clf()
        # CIL area
        cil_vol_stn = 0
        CIL = c_cil_stn.collections[0]
        for path in CIL.get_paths()[:]:
            vs = path.vertices
            cil_vol_stn = cil_vol_stn + np.abs(area(vs))/1000
        cil_vol_stn_clim[idx] = cil_vol_stn 
        # CIL core
        cil_core_stn_clim[idx] = np.nanmin(df_section_stn.values)
        # HERE I SOULD ADD CORE DEPTH!!
        
    if df_section_itp.index.size > 1:
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
        #plt.close()
        #plt.clf()
        #fig.clf()
        plt.close('all')
        # CIL area
        cil_vol_itp = 0
        CIL = c_cil_itp.collections[0]
        for path in CIL.get_paths()[:]:
            vs = path.vertices
            cil_vol_itp = cil_vol_itp + np.abs(area(vs))/1000
        cil_vol_itp_clim[idx] = cil_vol_itp 
        # CIL core
        cil_core_itp_clim[idx] = np.nanmin(df_section_itp.values)

    # Store results - Temperature
    df_section_stn.index.name='station'
    df_section_stn.columns.name='depth'
    df_stn_temp.append(df_section_stn)
    df_section_itp.index.name='station'
    df_section_itp.columns.name='depth'
    df_itp_temp.append(df_section_itp)
    # Store results - Salinity
    df_section_stn_S.index.name='station'
    df_section_stn_S.columns.name='depth'
    df_stn_sal.append(df_section_stn_S)
    df_section_itp_S.index.name='station'
    df_section_itp_S.columns.name='depth'
    df_itp_sal.append(df_section_itp_S)
    # Store results - Sigma-t
    df_section_stn_Si.index.name='station'
    df_section_stn_Si.columns.name='depth'
    df_stn_sig.append(df_section_stn_Si)
    df_section_itp_Si.index.name='station'
    df_section_itp_Si.columns.name='depth'
    df_itp_sig.append(df_section_itp_Si)
    # Store results - CIL            
    ## cil_vol_stn_clim[idx] = cil_vol_stn
    ## cil_vol_itp_clim[idx] = cil_vol_itp
    ## cil_core_stn_clim[idx] = np.nanmin(df_section_stn.values)
    ## cil_core_itp_clim[idx] = np.nanmin(df_section_itp.values)
    
## # Save temperature (not necessary)
df_stn_mindex = pd.concat(df_stn_temp,keys=years_series)
df_itp_mindex = pd.concat(df_itp_temp,keys=years_series)
## df_stn_mindex.to_pickle('df_stn_mindex_temp.pkl')
## df_itp_mindex.to_pickle('df_itp_mindex_temp.pkl')
## # Save Salinity (not necessary)
df_stn_mindex_S = pd.concat(df_stn_sal,keys=years_series)
df_itp_mindex_S = pd.concat(df_itp_sal,keys=years_series)
## df_stn_mindex_S.to_pickle('df_stn_mindex_sal.pkl')
## df_itp_mindex_S.to_pickle('df_itp_mindex_sal.pkl')
## # Save Sigma-t (not necessary)
df_stn_mindex_Si = pd.concat(df_stn_sig,keys=years_series)
df_itp_mindex_Si = pd.concat(df_itp_sig,keys=years_series)
## df_stn_mindex_Si.to_pickle('df_stn_mindex_sigt.pkl')
## df_itp_mindex_Si.to_pickle('df_itp_mindex_sigt.pkl')
# Save Climatology (this is really what is needed) - 'itp' only because no. station not always the same for 'stn'
df_clim =  df_itp_mindex.groupby(level=1).apply(lambda x: x.mean())
picklename = 'df_temperature_' + SECTION + '_' + SEASON + '_clim.pkl'
df_clim.to_pickle(picklename)
df_clim_S =  df_itp_mindex_S.groupby(level=1).apply(lambda x: x.mean())
picklename = 'df_salinity_' + SECTION + '_' + SEASON + '_clim.pkl'
df_clim_S.to_pickle(picklename)
df_clim_Si =  df_itp_mindex_Si.groupby(level=1).apply(lambda x: x.mean())
picklename = 'df_sigma-t_' + SECTION + '_' + SEASON + '_clim.pkl'
df_clim_Si.to_pickle(picklename)
# Save CIL timseries
df_CIL= pd.DataFrame([cil_vol_stn_clim, cil_vol_itp_clim, cil_core_stn_clim, cil_core_itp_clim]).T
df_CIL.index = years_series
df_CIL.columns = ['vol_stn', 'vol_itp', 'core_stn', 'core_itp']
picklename = 'df_CIL_' + SECTION + '_' + SEASON + '.pkl'
df_CIL.to_pickle(picklename)

# plot temperature
fig, ax = plt.subplots()
c = plt.contourf(df_clim.index, df_clim.columns, df_clim.T)
c_cil_itp = plt.contour(df_clim.index, df_clim.columns, df_clim.T, [0,], colors='k', linewidths=2)
#ax.set_ylim([0, 400])
#ax.set_xlim([0,  XLIM])
ax.set_ylabel('Depth (m)', fontWeight = 'bold')
#ax.set_xlabel('Distance (km)')
ax.invert_yaxis()
plt.colorbar(c)
figname = 'temperature' + SECTION + '_' + SEASON + '_clim.png'

# plot salinity
fig.savefig(figname, dpi=150)
fig, ax = plt.subplots()
c = plt.contourf(df_clim_S.index, df_clim_S.columns, df_clim_S.T)
ax.set_ylabel('Depth (m)', fontWeight = 'bold')
ax.invert_yaxis()
plt.colorbar(c)
figname = 'salinity' + SECTION + '_' + SEASON + '_clim.png'
fig.savefig(figname, dpi=150)

# plot sigma-t
fig.savefig(figname, dpi=150)
fig, ax = plt.subplots()
c = plt.contourf(df_clim_Si.index, df_clim_Si.columns, df_clim_Si.T)
ax.set_ylabel('Depth (m)', fontWeight = 'bold')
ax.invert_yaxis()
plt.colorbar(c)
figname = 'sigma-t' + SECTION + '_' + SEASON + '_clim.png'
fig.savefig(figname, dpi=150)
