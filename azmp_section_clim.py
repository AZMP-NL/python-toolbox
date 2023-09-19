'''
CIL area

This script is still in progress...

 Example how to read MIndex:
    df_itp_mindex = pd.read_pickle('df_itp_mindex.pkl')
    C = df_itp_mindex.xs(('year'),level=(0))
    C.groupby(level=0).apply(lambda x: x.mean())

'''
import os
import sys
import netCDF4
import xarray as xr
#from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import pandas as pd
import numpy as  np
from scipy.interpolate import interp1d  # to remove NaNs in profiles
from scipy.interpolate import griddata
sys.path.append('/home/jcoyne/Documents/AZMP-NL_python-toolbox/python-toolbox/azmp_modules')
import azmp_sections_tools as azst


## ---- Region parameters ---- ## <-------------------------------Would be nice to pass this in a config file '2017.report'
SECTION = 'BB'
SEASON = 'summer'
CLIM_YEAR = [1950, 2022]
#CLIM_YEAR = [1991, 2010]
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
df_stn = pd.read_excel('operation_files/STANDARD_SECTIONS.xlsx')
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
    dataFile = '/home/jcoyne/Documents/Datasets/GEBCO_2023/GEBCO_2023_sub_ice_topo.nc' # Maybe find a better way to handle this file

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


## Loop on climatological years
years = np.arange(CLIM_YEAR[0], CLIM_YEAR[1]+1)
years_series = pd.Series(years)
years_series.name='year'

cil_vol_stn_clim = np.full(years.shape, np.nan)
cil_vol_stn_man_clim = np.full(years.shape, np.nan)
cil_vol_itp_clim = np.full(years.shape, np.nan)
cil_core_stn_clim = np.full(years.shape, np.nan)
cil_core_stn_man_clim = np.full(years.shape, np.nan)
cil_core_itp_clim = np.full(years.shape, np.nan)
df_stn_temp = []
df_stn_man_temp = []
df_itp_temp = []
df_stn_sal = []
df_stn_man_sal = []
df_itp_sal = []
for idx, YEAR in enumerate(years):
    ## -------- Get CTD data -------- ##
    year_file = '/home/jcoyne/Documents/CASH/Combined_Data/CASTS_new-vertical_v2/' + str(YEAR) + '.nc'
    print('Get ' + year_file)
    ds = xr.open_dataset(year_file)

    # Select Region
    ds = ds.isel(time = (ds.longitude>lonLims[0]) & (ds.longitude<lonLims[1]))
    ds = ds.isel(time = (ds.latitude>latLims[0]) & (ds.latitude<latLims[1]))
    ds = ds.isel(time = ds.source == 'NAFC-Oceanography')

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
    lons = np.array(ds.longitude)
    lats = np.array(ds.latitude)
    bins = np.arange(z1, 2000, dz)
    da_temp = da_temp.groupby_bins('level', bins).mean(dim='level')
    da_sal = da_sal.groupby_bins('level', bins).mean(dim='level')

    # 1. Temperature to Pandas Dataframe
    print('Process temperature')
    df = da_temp.to_pandas()
    df.columns = bins[0:-1] #rename columns with 'bins'
    # Remove empty columns.
    idx_empty_rows = df.isnull().all(1).values.nonzero()[0]
    df = df.dropna(axis=0,how='all')
    lonsT = np.delete(lons,idx_empty_rows)
    latsT = np.delete(lats,idx_empty_rows)

    ## --- fill 3D cube --- ##  
    z = df.columns.values
    V_temp = np.full((lat_reg.size, lon_reg.size, z.size), np.nan)

    # Aggregate on regular grid
    for i, xx in enumerate(lon_reg):
        for j, yy in enumerate(lat_reg):    
            idx_coords = np.where((lonsT>=xx-dc/2) & (lonsT<xx+dc/2) & (latsT>=yy-dc/2) & (latsT<yy+dc/2))
            tmp = np.array(df.iloc[idx_coords].mean(axis=0))
            idx_good = np.argwhere((~np.isnan(tmp)) & (tmp<30))
            if np.size(idx_good)==1:
                V_temp[j,i,:] = np.array(df.iloc[idx_coords].mean(axis=0))
            elif np.size(idx_good)>1: # vertical interpolation between pts
                interp = interp1d(np.squeeze(z[idx_good]), np.squeeze(tmp[idx_good]))  # <---- Attention: strange syntax but works
                idx_interp = np.arange(int(idx_good[0]),int(idx_good[-1]+1))
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
    lonsS = np.delete(lons,idx_empty_rows)
    latsS = np.delete(lats,idx_empty_rows)

    ## --- fill 3D cube --- ##  
    z = df.columns.values
    V_sal = np.full((lat_reg.size, lon_reg.size, z.size), np.nan)

    # Aggregate on regular grid
    for i, xx in enumerate(lon_reg):
        for j, yy in enumerate(lat_reg):    
            idx_coords = np.where((lonsS>=xx-dc/2) & (lonsS<xx+dc/2) & (latsS>=yy-dc/2) & (latsS<yy+dc/2))
            tmp = np.array(df.iloc[idx_coords].mean(axis=0))
            idx_good = np.argwhere((~np.isnan(tmp)))
            if np.size(idx_good)==1:
                V_sal[j,i,:] = np.array(df.iloc[idx_coords].mean(axis=0))
            elif np.size(idx_good)>1: # vertical interpolation between pts
                interp = interp1d(np.squeeze(z[idx_good]), np.squeeze(tmp[idx_good]))  # <---- Attention: strange syntax but works
                idx_interp = np.arange(int(idx_good[0]),int(idx_good[-1]+1))
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

    # mask using bathymetry (I don't think it is necessary, but make nice figures)
    for i, xx in enumerate(lon_reg):
        for j,yy in enumerate(lat_reg):
            if Zitp[j,i] > -10: # remove shallower than 10m
                V_temp[j,i,:] = np.nan
                V_sal[j,i,:] = np.nan

    ## ---- Extract section info using station_ID (test 2 options) ---- ##
    ## Temperature & Salinity
    temp_coords = np.array([[x,y] for x in lat_reg for y in lon_reg])
    VT = V_temp.reshape(temp_coords.shape[0],V_temp.shape[2])
    VS = V_sal.reshape(temp_coords.shape[0],V_temp.shape[2])
    ZZ = Zitp.reshape(temp_coords.shape[0],1)
    # Initialize
    section_only = []
    df_section_itp = pd.DataFrame(index=stn_list, columns=z)
    df_section_itp_S = pd.DataFrame(index=stn_list, columns=z)
    for stn in stn_list:
        # 1. Section only (by station name)
        ds_tmp = ds.isel(time = ds.station_ID == stn)
        if ds_tmp.time.size > 1:
            ds_tmp_mean = ds_tmp.isel(time=[0])
            ds_tmp_mean.temperature[:] = ds_tmp.mean(dim='time',skipna=True).temperature.values
            ds_tmp_mean.salinity[:] = ds_tmp.mean(dim='time',skipna=True).salinity.values
            section_only.append(ds_tmp_mean)
        elif ds_tmp.time.size == 1:
            section_only.append(ds_tmp)


        #2.  From interpolated field (closest to station)
        station = df_stn[df_stn.STATION==stn]
        idx_opti = np.argmin(np.sum(np.abs(temp_coords - np.array([station.LAT.values,station.LON.values]).T), axis=1))
        # *** There's a problem here since often the optimum is NaN-only profile for station 1...
        Tprofile = VT[idx_opti,:]
        Sprofile = VS[idx_opti,:]
        # remove data below bottom
        bottom_depth = -ZZ[idx_opti]
        Tprofile[z>=bottom_depth]=np.nan
        Sprofile[z>=bottom_depth]=np.nan

        # store in dataframe
        df_section_itp.loc[stn] = Tprofile
        df_section_itp_S.loc[stn] = Sprofile

    # convert option #1 to dataframe    
    if len(section_only) > 0:
        ds_section = xr.concat(section_only, dim='time')
        station_ID = np.array(ds_section.station_ID.values.astype(str))
        if ds_section.time.size > 0:
            ds_section = ds_section.groupby(ds_section.station_ID.astype(str)).mean(dim='time')
        ds_section = ds_section.groupby_bins('level', bins).mean(dim='level')
        da = ds_section['temperature'].T
        df_section_stn = da.to_pandas()
        df_section_stn.columns = bins[0:-1]
        df_section_stn.index = np.unique(station_ID)
        da = ds_section['salinity'].T
        df_section_stn_S = da.to_pandas()
        df_section_stn_S.columns = bins[0:-1]
        df_section_stn_S.index = np.unique(station_ID)
        df_section_stn = df_section_stn.interpolate(axis=1,method='linear').where(df_section_stn.bfill(axis=1).notnull())
        df_section_stn_S = df_section_stn_S.interpolate(axis=1,method='linear').where(df_section_stn_S.bfill(axis=1).notnull())
    else:
        df_section_stn = pd.DataFrame(columns=z)
        df_section_stn_S = pd.DataFrame(columns=z)








    ## ---- Extract section info using station_ID_manual (test 3 options) ---- ##
    ## Temperature & Salinity
    temp_coords = np.array([[x,y] for x in lat_reg for y in lon_reg])
    VT = V_temp.reshape(temp_coords.shape[0],V_temp.shape[2])
    VS = V_sal.reshape(temp_coords.shape[0],V_temp.shape[2])
    ZZ = Zitp.reshape(temp_coords.shape[0],1)
    # Initialize
    section_only = []
    df_section_itp = pd.DataFrame(index=stn_list, columns=z)
    df_section_itp_S = pd.DataFrame(index=stn_list, columns=z)
    for stn in stn_list:
        # 1. Section only (by station name)
        ds_tmp = ds.isel(time = ds.station_ID_manual == 'AZMP_NAFC_'+stn)
        if ds_tmp.time.size > 1:
            ds_tmp_mean = ds_tmp.isel(time=[0])
            ds_tmp_mean.temperature[:] = ds_tmp.mean(dim='time',skipna=True).temperature.values
            ds_tmp_mean.salinity[:] = ds_tmp.mean(dim='time',skipna=True).salinity.values
            section_only.append(ds_tmp_mean)
        elif ds_tmp.time.size == 1:
            section_only.append(ds_tmp)

        #2.  From interpolated field (closest to station)
        station = df_stn[df_stn.STATION==stn]
        idx_opti = np.argmin(np.sum(np.abs(temp_coords - np.array([station.LAT.values,station.LON.values]).T), axis=1))
        # *** There's a problem here since often the optimum is NaN-only profile for station 1...
        Tprofile = VT[idx_opti,:]
        Sprofile = VS[idx_opti,:]
        # remove data below bottom
        bottom_depth = -ZZ[idx_opti]
        Tprofile[z>=bottom_depth]=np.nan
        Sprofile[z>=bottom_depth]=np.nan

        # store in dataframe
        df_section_itp.loc[stn] = Tprofile
        df_section_itp_S.loc[stn] = Sprofile

    # convert option #1 to dataframe
    if len(section_only) > 0:
        ds_section = xr.concat(section_only, dim='time')
        station_ID = np.array(ds_section.station_ID_manual.values.astype(str))
        station_ID = np.array([i.split('_')[-1] for i in station_ID])
        if ds_section.time.size > 0:
            ds_section = ds_section.groupby(ds_section.station_ID_manual.astype(str)).mean(dim='time')
        ds_section = ds_section.groupby_bins('level', bins).mean(dim='level')
        da = ds_section['temperature'].T
        df_section_stn_man = da.to_pandas()
        df_section_stn_man.columns = bins[0:-1]
        df_section_stn_man.index = np.unique(station_ID)
        da = ds_section['salinity'].T
        df_section_stn_man_S = da.to_pandas()
        df_section_stn_man_S.columns = bins[0:-1]
        df_section_stn_man_S.index = np.unique(station_ID)
        df_section_stn_man = df_section_stn_man.interpolate(axis=1,method='linear').where(df_section_stn_man.bfill(axis=1).notnull())
        df_section_stn_man_S = df_section_stn_man_S.interpolate(axis=1,method='linear').where(df_section_stn_man_S.bfill(axis=1).notnull())
    else:
        df_section_stn_man = pd.DataFrame(columns=z)
        df_section_stn_man_S = pd.DataFrame(columns=z)




    # drop indices with only NaNs (option #2)
    df_section_itp = df_section_itp.dropna(axis=0, how='all')
    df_section_itp_S = df_section_itp_S.dropna(axis=0, how='all')
    df_section_itp = df_section_itp.astype(float)
    df_section_itp_S = df_section_itp_S.astype(float)


    ## CIL Calculation
    # Compute distance vector for option #1 - exact station
    distance_stn = np.full((df_section_stn.index.shape), np.nan)
    for i, stn in enumerate(df_section_stn.index):
        distance_stn[i] = azst.haversine(df_stn.LON[0], df_stn.LAT[0], df_stn[df_stn.STATION==df_section_stn.index[i]].LON, df_stn[df_stn.STATION==df_section_stn.index[i]].LAT)    
    
    # Compute distance vector for option #1 - station_manual
    distance_stn_man = np.full((df_section_stn_man.index.shape), np.nan)
    for i, stn in enumerate(df_section_stn_man.index):
        distance_stn_man[i] = azst.haversine(df_stn.LON[0], df_stn.LAT[0], df_stn[df_stn.STATION==df_section_stn_man.index[i]].LON, df_stn[df_stn.STATION==df_section_stn_man.index[i]].LAT)
    
    # Compute distance vector for option #2 - interp field
    distance_itp = np.full((df_section_itp.index.shape), np.nan)
    for i, stn in enumerate(df_section_itp.index):
        distance_itp[i] = azst.haversine(df_stn.LON[0], df_stn.LAT[0], df_stn[df_stn.STATION==df_section_itp.index[i]].LON, df_stn[df_stn.STATION==df_section_itp.index[i]].LAT)

    if df_section_stn.index.size > 1:
        fig, ax = plt.subplots()
        levels = np.arange(-2,10)
        c = plt.contourf(distance_stn, df_section_stn.columns, df_section_stn.T,levels)
        c_cil_stn = plt.contour(distance_stn, df_section_stn.columns, df_section_stn.T, [0,], colors='k', linewidths=2)
        ax.set_ylim([0, 400])
        ax.set_ylabel('Depth (m)', fontweight = 'bold')
        ax.set_xlabel('Distance (km)')
        ax.invert_yaxis()
        plt.xlim([0,800])
        plt.colorbar(c)
        plt.title(str(YEAR))
        fig_name = 'temp_section_' + SECTION + '_' + SEASON + '_' + '_' + str(YEAR) + '_1.png'
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
        #if df_section_stn.index.size/df_stn.index.size >= 0.75:
        cil_vol_stn_clim[idx] = cil_vol_stn 
        # CIL core
        cil_core_stn_clim[idx] = np.nanmin(df_section_stn.values)
        # HERE I SOULD ADD CORE DEPTH!!
        
    if df_section_stn_man.index.size > 1:
        fig, ax = plt.subplots()
        levels = np.arange(-2,10)
        c = plt.contourf(distance_stn_man, df_section_stn_man.columns, df_section_stn_man.T,levels)
        c_cil_stn = plt.contour(distance_stn_man, df_section_stn_man.columns, df_section_stn_man.T, [0,], colors='k', linewidths=2)
        ax.set_ylim([0, 400])
        ax.set_ylabel('Depth (m)', fontweight = 'bold')
        ax.set_xlabel('Distance (km)')
        ax.invert_yaxis()
        plt.xlim([0,800])
        plt.colorbar(c)
        plt.title(str(YEAR))
        fig_name = 'temp_section_' + SECTION + '_' + SEASON + '_' + '_' + str(YEAR) + '_1man.png'
        fig.savefig(fig_name, dpi=150)
        #plt.clf()
        plt.close('all')
        #fig.clf()
        # CIL area
        cil_vol_stn = 0
        CIL = c_cil_stn.collections[0]
        for path in CIL.get_paths()[:]:
            vs = path.vertices
            cil_vol_stn_man = cil_vol_stn + np.abs(area(vs))/1000
        #if df_section_stn.index.size/df_stn.index.size >= 0.75:
        cil_vol_stn_man_clim[idx] = cil_vol_stn_man 
        # CIL core
        cil_core_stn_man_clim[idx] = np.nanmin(df_section_stn_man.values)
        # HERE I SOULD ADD CORE DEPTH!!

    if df_section_itp.index.size > 1:
        fig, ax = plt.subplots()
        levels = np.arange(-2,10)
        c = plt.contourf(distance_itp, df_section_itp.columns, df_section_itp.T,levels)
        c_cil_itp = plt.contour(distance_itp, df_section_itp.columns, df_section_itp.T, [0,], colors='k', linewidths=2)
        ax.set_ylim([0, 400])
        ax.set_ylabel('Depth (m)', fontweight = 'bold')
        ax.set_xlabel('Distance (km)')
        ax.invert_yaxis()
        plt.xlim([0,800])
        plt.title(str(YEAR))
        plt.colorbar(c)
        fig_name = 'temp_section_' + SECTION + '_' + SEASON + '_' + '_' + str(YEAR) + '_2.png'
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

    df_section_stn_man.index.name = 'station'
    df_section_stn_man.columns.name = 'depth'
    df_stn_man_temp.append(df_section_stn_man)

    df_section_itp.index.name='station'
    df_section_itp.columns.name='depth'
    df_itp_temp.append(df_section_itp)
    # Store results - Salinity
    df_section_stn_S.index.name='station'
    df_section_stn_S.columns.name='depth'
    df_stn_sal.append(df_section_stn_S)

    df_section_stn_man_S.index.name = ' station'
    df_section_stn_man_S.columns.name = 'depth'
    df_stn_man_sal.append(df_section_stn_man_S)

    df_section_itp_S.index.name='station'
    df_section_itp_S.columns.name='depth'
    df_itp_sal.append(df_section_itp_S)
    # Store results - CIL            
    ## cil_vol_stn_clim[idx] = cil_vol_stn
    ## cil_vol_itp_clim[idx] = cil_vol_itp
    ## cil_core_stn_clim[idx] = np.nanmin(df_section_stn.values)
    ## cil_core_itp_clim[idx] = np.nanmin(df_section_itp.values)
    del ds, df

    
## # Save temperature (not necessary)
df_stn_mindex = pd.concat(df_stn_temp,keys=years_series)
df_stn_man_mindex = pd.concat(df_stn_man_temp,keys=years_series)
df_itp_mindex = pd.concat(df_itp_temp,keys=years_series)
## df_stn_mindex.to_pickle('df_stn_mindex_temp.pkl')
## df_itp_mindex.to_pickle('df_itp_mindex_temp.pkl')
## # Save Salinity (not necessary)
df_stn_mindex_S = pd.concat(df_stn_sal,keys=years_series)
df_stn_man_mindex_S = pd.concat(df_stn_man_sal,keys=years_series)
df_itp_mindex_S = pd.concat(df_itp_sal,keys=years_series)
## df_stn_mindex_S.to_pickle('df_stn_mindex_sal.pkl')
## df_itp_mindex_S.to_pickle('df_itp_mindex_sal.pkl')
# Save Climatology (this is really what is needed) - 'itp' only because no. station not always the same for 'stn'
df_clim =  df_stn_mindex.groupby(level=1).apply(lambda x: x.mean())
picklename = 'df_temperature_' + SECTION + '_' + SEASON + '_clim.pkl'
df_clim.to_pickle(picklename)
df_clim_S =  df_stn_mindex_S.groupby(level=1).apply(lambda x: x.mean())
picklename = 'df_salinity_' + SECTION + '_' + SEASON + '_clim.pkl'
df_clim_S.to_pickle(picklename)
# Save CIL timseries
df_CIL= pd.DataFrame([cil_vol_stn_clim, cil_vol_stn_man_clim, cil_vol_itp_clim, cil_core_stn_clim, cil_core_stn_man_clim, cil_core_itp_clim]).T
df_CIL.index = years_series
df_CIL.columns = ['vol_stn', 'vol_stn_man', 'vol_itp', 'core_stn', 'core_stn_man', 'core_itp']
picklename = 'df_CIL_' + SECTION + '_' + SEASON + '.pkl'
df_CIL.to_pickle(picklename)





## TESTING AREA - Add climatology to blank areas for each year

#This part adds values to nans for stations that are covered
place1 = []
stn_years = np.unique(df_stn_mindex.index.get_level_values('year').values)
stn_man_years = np.unique(df_stn_man_mindex.index.get_level_values('year').values)
itp_years = np.unique(df_itp_mindex.index.get_level_values('year').values)


for YEAR in years:

    if np.isin(YEAR,stn_years):
        #Fill in the stn output with climatology when not available
        stn_clim = df_stn_mindex.groupby(level=1).apply(lambda x: x.mean()).T
        stn_orig = df_stn_mindex.T[[YEAR]].copy()
        stn_orig = stn_orig.T.droplevel('year').T
        stn_orig = stn_orig.reindex(columns=stn_clim.columns,fill_value=np.nan)


        stn_merged = stn_orig.copy()
        place2 = np.array(stn_clim,dtype=float)
        place3 = np.array(stn_orig,dtype=float)
        place3[np.isnan(place3)] = place2[np.isnan(place3)]
        stn_merged.iloc[:,:] = place3.astype(float)

        place1.append(stn_merged.T)


        '''
        #THIS WILL BE COMMENTED OUT
        #Create a figure showing where the AZMP section has been filled in by climatology
        plt.figure(figsize=(10,6))
        plt.pcolor(
            stn_merged.columns,
            stn_merged.index,
            stn_merged,
            alpha=0.5,vmin=-2,vmax=20,snap=True
            )
        plt.pcolor(
            stn_orig.columns,
            stn_orig.index,
            stn_orig,
            alpha=1,vmin=-2,vmax=20,snap=True
            )
        plt.ylabel('Depth (m)')
        plt.ylim([0,1500])
        plt.xlabel('Distance (km)')
        plt.gca().invert_yaxis()
        plt.xticks(rotation=90)
        plt.title(SECTION+'-'+SEASON+', '+str(YEAR))
        plt.colorbar(label='Temperature ($^\circ$C)')
        plt.tight_layout()
        

        fig, ax = plt.subplots()
        levels = np.arange(-2,10)
        c = plt.contourf(distance_itp, df_section_itp.columns, df_section_itp.T,levels)
        c_cil_itp = plt.contour(distance_itp, df_section_itp.columns, df_section_itp.T, [0,], colors='k', linewidths=2)
        ax.set_ylim([0, 400])
        ax.set_ylabel('Depth (m)', fontweight = 'bold')
        ax.set_xlabel('Distance (km)')
        ax.invert_yaxis()
        plt.xlim([0,800])
        plt.title(str(YEAR))
        plt.colorbar(c)
        fig_name = 'temp_section_' + SECTION + '_' + SEASON + '_' + '_' + str(YEAR) + '_2.png'
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
        '''






        ## CIL Calculation
        # Compute distance vector for option #1 - exact station
        distance_stn = np.full((stn_merged.T.index.shape), np.nan)
        for i, stn in enumerate(stn_merged.T.index):
            distance_stn[i] = azst.haversine(df_stn.LON[0], df_stn.LAT[0], df_stn[df_stn.STATION==stn_merged.T.index[i]].LON, df_stn[df_stn.STATION==stn_merged.T.index[i]].LAT)    

        levels = np.arange(-2,10)
        c = plt.contourf(distance_stn, stn_merged.T.columns, stn_merged,levels)
        c_cil_stn = plt.contour(distance_stn, stn_merged.T.columns, stn_merged, [0,], colors='k', linewidths=2)
        plt.clf()
        #fig.clf()
        # CIL area
        cil_vol_stn = 0
        CIL = c_cil_stn.collections[0]
        for path in CIL.get_paths()[:]:
            vs = path.vertices
            cil_vol_stn = cil_vol_stn + np.abs(area(vs))/1000


        #if df_section_stn.index.size/df_stn.index.size >= 0.75:
        cil_vol_stn_clim[np.where(YEAR == years)[0][0]] = cil_vol_stn 
        # CIL core
        cil_core_stn_clim[np.where(YEAR == years)[0][0]] = np.nanmin(stn_merged.values)


    if np.isin(YEAR,stn_man_years):
        #Same, but for stn_man method
        stn_clim = df_stn_mindex.groupby(level=1).apply(lambda x: x.mean()).T
        stn_orig = df_stn_man_mindex.T[[YEAR]].copy()
        stn_orig = stn_orig.T.droplevel('year').T
        stn_orig = stn_orig.reindex(columns=stn_clim.columns,fill_value=np.nan)


        stn_merged = stn_orig.copy()
        place2 = np.array(stn_clim)
        place3 = np.array(stn_orig)
        place3[pd.isnull(place3)] = place2[pd.isnull(place3)]
        stn_merged.iloc[:,:] = place3.astype(float)

        place1.append(stn_merged.T)


        ## CIL Calculation
        # Compute distance vector for option #1 - exact station
        distance_stn = np.full((stn_merged.T.index.shape), np.nan)
        for i, stn in enumerate(stn_merged.T.index):
            distance_stn[i] = azst.haversine(df_stn.LON[0], df_stn.LAT[0], df_stn[df_stn.STATION==stn_merged.T.index[i]].LON, df_stn[df_stn.STATION==stn_merged.T.index[i]].LAT)    

        levels = np.arange(-2,10)
        c = plt.contourf(distance_stn, stn_merged.T.columns, stn_merged,levels)
        c_cil_itp = plt.contour(distance_stn, stn_merged.T.columns, stn_merged, [0,], colors='k', linewidths=2)
        plt.clf()
        #fig.clf()
        # CIL area
        cil_vol_itp = 0
        CIL = c_cil_itp.collections[0]
        for path in CIL.get_paths()[:]:
            vs = path.vertices
            cil_vol_itp = cil_vol_itp + np.abs(area(vs))/1000


        #if df_section_stn.index.size/df_stn.index.size >= 0.75:
        cil_vol_stn_man_clim[np.where(YEAR == years)[0][0]] = cil_vol_itp 
        # CIL core
        cil_core_stn_man_clim[np.where(YEAR == years)[0][0]] = np.nanmin(stn_merged.values)


    if np.isin(YEAR,itp_years):
        #Same, but for interpolation method
        stn_clim = df_stn_mindex.groupby(level=1).apply(lambda x: x.mean()).T
        stn_orig = df_itp_mindex.T[[YEAR]].copy()
        stn_orig = stn_orig.T.droplevel('year').T
        stn_orig = stn_orig.reindex(columns=stn_clim.columns,fill_value=np.nan)


        stn_merged = stn_orig.copy()
        place2 = np.array(stn_clim)
        place3 = np.array(stn_orig)
        place3[pd.isnull(place3)] = place2[pd.isnull(place3)]
        stn_merged.iloc[:,:] = place3.astype(float)

        place1.append(stn_merged.T)



        ## CIL Calculation
        # Compute distance vector for option #1 - exact station
        distance_stn = np.full((stn_merged.T.index.shape), np.nan)
        for i, stn in enumerate(stn_merged.T.index):
            distance_stn[i] = azst.haversine(df_stn.LON[0], df_stn.LAT[0], df_stn[df_stn.STATION==stn_merged.T.index[i]].LON, df_stn[df_stn.STATION==stn_merged.T.index[i]].LAT)    

        levels = np.arange(-2,10)
        c = plt.contourf(distance_stn, stn_merged.T.columns, stn_merged,levels)
        c_cil_itp = plt.contour(distance_stn, stn_merged.T.columns, stn_merged, [0,], colors='k', linewidths=2)
        plt.clf()
        #fig.clf()
        # CIL area
        cil_vol_itp = 0
        CIL = c_cil_itp.collections[0]
        for path in CIL.get_paths()[:]:
            vs = path.vertices
            cil_vol_itp = cil_vol_itp + np.abs(area(vs))/1000


        #if df_section_stn.index.size/df_stn.index.size >= 0.75:
        cil_vol_itp_clim[np.where(YEAR == years)[0][0]] = cil_vol_itp 
        # CIL core
        cil_core_itp_clim[np.where(YEAR == years)[0][0]] = np.nanmin(stn_merged.values)
    print(YEAR)

#place2 = pd.concat(place1,keys=stn_years)
plt.close('all')

# Save CIL timseries
df_CIL= pd.DataFrame([cil_vol_stn_clim, cil_vol_stn_man_clim, cil_vol_itp_clim, cil_core_stn_clim, cil_core_stn_man_clim, cil_core_itp_clim]).T
df_CIL.index = years_series
df_CIL.columns = ['vol_stn', 'vol_stn_man', 'vol_itp', 'core_stn', 'core_stn_man', 'core_itp']
picklename = 'df_CIL_' + SECTION + '_' + SEASON + '.pkl'
df_CIL.to_pickle(picklename)



















# plot temperature
fig, ax = plt.subplots()
c = plt.contourf(df_clim.index, df_clim.columns, df_clim.T)
c_cil_itp = plt.contour(df_clim.index, df_clim.columns, df_clim.T, [0,], colors='k', linewidths=2)
ax.set_ylim([0, 400])
#ax.set_xlim([0,  XLIM])
ax.set_ylabel('Depth (m)', fontweight = 'bold')
#ax.set_xlabel('Distance (km)')
ax.invert_yaxis()
plt.colorbar(c)
figname = 'temperature_' + SECTION + '_' + SEASON + '_clim.png'

# plot salinity
fig.savefig(figname, dpi=150)
fig, ax = plt.subplots()
c = plt.contourf(df_clim_S.index, df_clim_S.columns, df_clim_S.T)
ax.set_ylim([0, 400])
ax.set_ylabel('Depth (m)', fontweight = 'bold')
ax.invert_yaxis()
plt.colorbar(c)
figname = 'salinity_' + SECTION + '_' + SEASON + '_clim.png'
fig.savefig(figname, dpi=150)

