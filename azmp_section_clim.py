import os
import sys
import netCDF4
import xarray as xr
import warnings
import matplotlib.pyplot as plt
import cmocean as cm
import pandas as pd
import numpy as  np
from scipy.interpolate import interp1d  # to remove NaNs in profiles
from scipy.interpolate import griddata
sys.path.append(os.path.expanduser('~/github/AZMP-NL/python-toolbox/azmp_modules'))
import azmp_sections_tools as azst


'''
Sample code for improving the layout and runtime of section calculations.
Working as of January 2024

Typically, we want FC, BB, SI for spring, summer, fall

Remember to move section plots to AZMP_lines/
Move the .pkl files into operation_files/
'''


## ---- Script Functions ---- ##
#CIL surface (Note that there is a bias because )
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

#3D temperature or salinity
def interpolate_variable_3D(df,lat_reg,lon_reg,lonsV,latsV,dc):
    #Create 3D cube
    z = df.columns.values
    V = np.full((lat_reg.size, lon_reg.size, z.size), np.nan)
    #Cycle through each lat and lon
    for i,xx in enumerate(lon_reg):
        for ii,yy in enumerate(lat_reg):
            idx_coords = np.where((lonsV>=xx-dc/2) & (lonsV<xx+dc/2) & (latsV>=yy-dc/2) & (latsV<yy+dc/2))
            var = np.array(df.iloc[idx_coords].mean(axis=0))
            idx_good = np.argwhere((~np.isnan(var)))
            #Determine if vertical interpolation is needed
            if np.size(idx_good) == 1:
                V[ii,i,:] = np.array(df.iloc[idx_coords].mean(axis=0))
            elif np.size(idx_good) > 1:
                interp = interp1d(np.squeeze(z[idx_good]), np.squeeze(var[idx_good]))
                idx_interp = np.arange(int(idx_good[0]),int(idx_good[-1]+1))
                #interpolate only where possible (1st to last good idx)
                V[ii,i,idx_interp] = interp(z[idx_interp])

    #Next, interpolate horizontally
    lon_grid, lat_grid = np.meshgrid(lon_reg,lat_reg)
    lon_vec = np.reshape(lon_grid, lon_grid.size)
    lat_vec = np.reshape(lat_grid, lat_grid.size)
    for i, xx in enumerate(z):
        #Meshgrid 1D data (after removing NaNs)
        V_grid = V[:,:,i]
        V_vec = np.reshape(V_grid, V_grid.size)
        #Interpolate using griddata (after removing nans)
        idx_good = np.argwhere(~np.isnan(V_vec))
        #We will ignore depth where minimal data exists
        if idx_good.size>5: 
            LN = np.squeeze(lon_vec[idx_good])
            LT = np.squeeze(lat_vec[idx_good])
            VV = np.squeeze(V_vec[idx_good])
            #May happen for sampling along 47N section (cannot linear griddata single latitude)
            if (np.unique(LT).size == 1) | (np.unique(LN).size == 1):
                zi = griddata((LN, LT), VV, (lon_grid, lat_grid), method='nearest')
            else:
                zi = griddata((LN, LT), VV, (lon_grid, lat_grid), method='linear')
            V[:,:,i] = zi
        else:
            continue
    return V,z

#Temperature and salinity based upon station ID or station ID manual
def name_variable(ds,stn_list,station_ID,manual=False):
    #Set up a place to record potential stations
    section_only = []
    for stn in stn_list:
        #Determine where potential stations are
        if manual:
            ds_V = ds.isel(time = station_ID == 'AZMP_NAFC_'+stn)
            ds_V['station_ID'] = (['time'], np.tile('AZMP_NAFC_'+stn, ds_V.time.size))
        else:
            ds_V = ds.isel(time = station_ID == stn)
            ds_V['station_ID'] = (['time'], np.tile(stn, ds_V.time.size))
        if ds_V.time.size > 1:
            #If more than 1 measurement, average (shouldn't happen for station_ID)
            ds_V_mean = ds_V.isel(time=[0])
            ds_V_mean.temperature[:,0] = ds_V.mean(dim='time', skipna=True).temperature.values
            ds_V_mean.salinity[:,0] = ds_V.mean(dim='time', skipna=True).salinity.values
            section_only.append(ds_V_mean)
        elif ds_V.time.size == 1:
            section_only.append(ds_V)
    #Convert to a dataframe
    if len(section_only) > 0:
        ds_section = xr.concat(section_only, dim='time')
        station_ID = ds_section.station_ID.values.astype(str)
        if manual:
            station_ID = np.array([i.split('_')[-1] for i in station_ID])
        da = ds_section['temperature'].T
        df_section_stn_T = da.to_pandas()
        df_section_stn_T.columns = bins[0:-1]
        df_section_stn_T.index = np.unique(station_ID)
        da = ds_section['salinity'].T
        df_section_stn_S = da.to_pandas()
        df_section_stn_S.columns = bins[0:-1]
        df_section_stn_S.index = np.unique(station_ID)
        #Interpolate along depth where possible (between points)
        df_section_stn_T = df_section_stn_T.interpolate(axis=1,method='linear').where(df_section_stn_T.bfill(axis=1).notnull())
        df_section_stn_S = df_section_stn_S.interpolate(axis=1,method='linear').where(df_section_stn_S.bfill(axis=1).notnull())
    else:
        df_section_stn_T = pd.DataFrame(columns=z)
        df_section_stn_S = pd.DataFrame(columns=z)

    return df_section_stn_T,df_section_stn_S

#Temperature plot for each method
def plot_temperature(distance, df, year, section, season, method=''):
    #Create the figure
    fig,ax = plt.subplots()
    levels=np.arange(-2,10)
    #Plot the temperature
    c = plt.contourf(
        distance,
        df.columns,
        df.T,
        levels,cmap=cm.cm.thermal
        )
    #Plot and record CIL
    if np.nanmin(df.T) <= 0:
        c_cil = plt.contour(
            distance,
            df.columns,
            df.T,
            [0,], colors='k', linewidths=2
            )
        cil_vol = 0
        CIL = c_cil.collections[0]
        for path in CIL.get_paths()[:]:
            vs = path.vertices
            cil_vol = cil_vol + np.abs(area(vs))/1000
        cil_core = np.nanmin(df.values)
    else:
        cil_vol = 0
        cil_core = np.nan
    #Other stuff
    ax.set_ylim([0, 400])
    ax.set_ylabel('Depth (m)', fontweight = 'bold')
    ax.set_xlabel('Distance (km)')
    ax.invert_yaxis()
    plt.xlim([0,800])
    plt.colorbar(c)
    plt.title(str(year))
    #Save the figure
    fig_name = 'temp_section_' + section + '_' + season + '_' + '_' + str(year) + '_'+method+'.png'
    fig.savefig(fig_name, dpi=150)
    plt.close('all')
    return cil_vol, cil_core

#Fill in temperature with climatology
def temperature_clim_fill(clim, year_data, df_stn):
    #Merge the climatology and yearly data
    year_merged = year_data.copy()
    stn_clim = np.array(clim)
    year_data = np.array(year_data)
    year_data[pd.isnull(year_data)] = stn_clim[pd.isnull(year_data)]
    year_merged.iloc[:,:] = year_data.astype(float)

    #Compute the new distance vector
    distance_stn = np.full((year_merged.T.index.shape), np.nan)
    for i, stn in enumerate(year_merged.T.index):
        distance_stn[i] = azst.haversine(
            df_stn.LON[0],
            df_stn.LAT[0],
            df_stn[df_stn.STATION==year_merged.T.index[i]].LON,
            df_stn[df_stn.STATION==year_merged.T.index[i]].LAT
            )
    #Calculate CIL
    #Plot the figure
    levels = np.arange(-2,10)
    c = plt.contourf(
        distance_stn,
        year_merged.T.columns,
        year_merged,
        levels
        )
    c_cil = plt.contour(
        distance_stn,
        year_merged.T.columns,
        year_merged,
        [0,],colors='k',linewidths=2)
    plt.clf()
    #CIL area
    cil_vol = 0
    CIL = c_cil.collections[0]
    for path in CIL.get_paths()[:]:
        vs = path.vertices
        cil_vol = cil_vol + np.abs(area(vs))/1000
    cil_core = np.nanmin(year_merged)
    return cil_vol, cil_core

'''
SECTION = 'SI'
SEASON = 'summer'
CLIM_YEAR = [1950, 2023]
dlat = 2 # how far from station we search
dlon = 2
z1 = 2
dz = 5 # vertical bins
dc = .2 # grid resolution
CASTS_path = '~/data/CASTS/'
bath_path = '~/data/GEBCO/GEBCO_2023_sub_ice_topo.nc'
'''
def section_clim(SECTION,SEASON,CLIM_YEAR,dlat,dlon,z1,dz,dc,CASTS_path,bath_path):

    ## ---- Get Stations ---- ## 
    df_stn = pd.read_excel('~/github/AZMP-NL/utils/STANDARD_SECTIONS.xlsx')
    df_stn = df_stn.drop(['SECTION', 'LONG'], axis=1)
    df_stn = df_stn.rename(columns={'LONG.1': 'LON'})
    df_stn = df_stn.dropna()
    df_stn = df_stn[df_stn.STATION.str.contains(SECTION+'-')]
    df_stn = df_stn.reset_index(drop=True)
    stn_list = df_stn.STATION.values

    #Derive regular grid
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
        dataFile = bath_path # Maybe find a better way to handle this file

        dataset = xr.open_dataset(dataFile)
        dataset = dataset.isel(lon=(dataset.lon>=lonLims[0])*(dataset.lon<=lonLims[1]))
        dataset = dataset.isel(lat=(dataset.lat>=latLims[0])*(dataset.lat<=latLims[1]))

        #Extract latitude and longitude
        lon = dataset.lon.values
        lat = dataset.lat.values

        #Extract the elevation
        Z = dataset.elevation.values

        #Set up the lats and lons
        lon_grid, lat_grid = np.meshgrid(lon_reg,lat_reg)
        lon_grid_bathy, lat_grid_bathy = np.meshgrid(lon,lat)

        #Interpolate
        Zitp = griddata((lon_grid_bathy.flatten(),lat_grid_bathy.flatten()), Z.flatten(), (lon_grid,lat_grid), method='linear')
        np.save(bathy_file, Zitp)
        print(' -> Done!')


    #Loop on climatological years
    years = np.arange(CLIM_YEAR[0], CLIM_YEAR[1]+1)
    years_series = pd.Series(years)
    years_series.name='year'

    #Set up the empty arrays and lists
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

    #Cycle through each year
    for idx, YEAR in enumerate(years):

        ## -------- Get CTD data -------- ##
        year_file = CASTS_path + str(YEAR) + '.nc'
        print('Get ' + year_file)
        ds = xr.open_dataset(year_file)

        #Select Region
        ds = ds.isel(time = (ds.longitude>lonLims[0]) & (ds.longitude<lonLims[1]))
        ds = ds.isel(time = (ds.latitude>latLims[0]) & (ds.latitude<latLims[1]))
        ds = ds.isel(time = ds.source == 'NAFC-Oceanography')

        #Take the station ID and station ID manual
        station_ID = ds.station_ID.values.astype(str)
        station_ID_manual = ds.station_ID_manual.values.astype(str)

        #Do the bin average here
        bins = np.arange(z1, 2000, dz)
        warnings.simplefilter(action='ignore', category=FutureWarning)
        ds = ds.groupby_bins('level', bins).mean(dim='level')

        #Select time (save several options here)
        if SEASON == 'summer':
            station_ID = station_ID[(ds['time.month'].values>=6)*(ds['time.month'].values<=9)]
            station_ID_manual = station_ID_manual[(ds['time.month'].values>=6)*(ds['time.month'].values<=9)]
            ds = ds.sel(time=((ds['time.month']>=6)) & ((ds['time.month']<=9)))
        elif SEASON == 'spring':
            station_ID = station_ID[(ds['time.month'].values>=3)*(ds['time.month'].values<=5)]
            station_ID_manual = station_ID_manual[(ds['time.month'].values>=3)*(ds['time.month'].values<=5)]
            ds = ds.sel(time=((ds['time.month']>=3)) & ((ds['time.month']<=5)))
        elif SEASON == 'fall':
            station_ID = station_ID[(ds['time.month'].values>=10)*(ds['time.month'].values<=12)]
            station_ID_manual = station_ID_manual[(ds['time.month'].values>=10)*(ds['time.month'].values<=12)]
            ds = ds.sel(time=((ds['time.month']>=10)) & ((ds['time.month']<=12)))
        else:
            print('!! no season specified, used them all! !!')

        #Extract temperature and salinity
        da_temp = ds['temperature'].T
        da_sal = ds['salinity'].T
        lons = np.array(ds.longitude)[0,:]
        lats = np.array(ds.latitude)[0,:]


        ## -------- Method 1: Interpolation -------- ##
        #Temperature to Pandas Dataframe
        print('Process temperature')
        df = da_temp.to_pandas()
        #Rename columns with 'bins'
        df.columns = bins[0:-1]
        # Remove empty columns.
        idx_empty_rows = df.isnull().all(1).values.nonzero()[0]
        df = df.dropna(axis=0,how='all')
        lonsT = np.delete(lons,idx_empty_rows)
        latsT = np.delete(lats,idx_empty_rows)
        V_temp,z = interpolate_variable_3D(df, lat_reg, lon_reg, lonsT, latsT, dc)
        print(' -> Done!')

        #Salinity to Pandas Dataframe
        print('Process salinity')
        df = da_sal.to_pandas()
        #Rename columns with 'bins'
        df.columns = bins[0:-1]
        # Remove empty columns
        idx_empty_rows = df.isnull().all(1).values.nonzero()[0]
        df = df.dropna(axis=0,how='all')
        lonsS = np.delete(lons,idx_empty_rows)
        latsS = np.delete(lats,idx_empty_rows)
        V_sal,z = interpolate_variable_3D(df, lat_reg, lon_reg, lonsS, latsS, dc)
        print(' -> Done!')

        #Mask using bathymetry (not completely necessary)
        for i, xx in enumerate(lon_reg):
            for ii,yy in enumerate(lat_reg):
                #Remove measurements shallower than 10m
                if Zitp[ii,i] > -10:
                    V_temp[ii,i,:] = np.nan
                    V_sal[ii,i,:] = np.nan

        #Determine the indices closest to the stations of interest
        V_coords = np.array([[x,y] for x in lat_reg for y in lon_reg])
        VT = V_temp.reshape(V_coords.shape[0],V_temp.shape[2])
        VS = V_sal.reshape(V_coords.shape[0],V_sal.shape[2])
        ZZ = Zitp.reshape(V_coords.shape[0],1)
        #Set up the dataframes
        df_section_itp_T = pd.DataFrame(index=stn_list, columns=z)
        df_section_itp_S = pd.DataFrame(index=stn_list, columns=z)
        for stn in stn_list:
            #Determine index closest to station
            station = df_stn[df_stn.STATION==stn]
            idx_opti = np.argmin(np.sum(np.abs(V_coords - np.array([station.LAT.values,station.LON.values]).T), axis=1))
            Tprofile = VT[idx_opti,:]
            Sprofile = VS[idx_opti,:]
            # remove data below bottom
            bottom_depth = -ZZ[idx_opti]
            Tprofile[z>=bottom_depth]=np.nan
            Sprofile[z>=bottom_depth]=np.nan
            # store in dataframe
            df_section_itp_T.loc[stn] = Tprofile
            df_section_itp_S.loc[stn] = Sprofile

        # Drop indices with only NaNs (option #2)
        df_section_itp_T = df_section_itp_T.dropna(axis=0, how='all')
        df_section_itp_S = df_section_itp_S.dropna(axis=0, how='all')
        df_section_itp_T = df_section_itp_T.astype(float)
        df_section_itp_S = df_section_itp_S.astype(float)

        #Compute the distance between stations
        distance_itp = np.full(df_section_itp_T.index.shape, np.nan)
        for i,stn in enumerate(df_section_itp_T.index):
            distance_itp[i] = azst.haversine(
                df_stn.LON[0],
                df_stn.LAT[0],
                df_stn[df_stn.STATION==df_section_itp_T.index[i]].LON,
                df_stn[df_stn.STATION==df_section_itp_T.index[i]].LAT,
                )

        #Plot the CIL and record stats
        if df_section_itp_T.index.size > 1:
            cil_vol,cil_core = plot_temperature(distance_itp, df_section_itp_T, YEAR, SECTION, SEASON, method='itp')
            cil_vol_itp_clim[idx] = cil_vol
            cil_core_itp_clim[idx] = cil_core

        #Record the results
        df_section_itp_T.index.name = 'station'
        df_section_itp_T.columns.name = 'depth'
        df_itp_temp.append(df_section_itp_T)
        df_section_itp_S.index.name = 'station'
        df_section_itp_S.columns.name = 'depth'
        df_itp_sal.append(df_section_itp_S)

        ## -------- Method 2: Station_ID -------- ##
        df_section_stn_T,df_section_stn_S = name_variable(ds, stn_list, station_ID.squeeze(), manual=False)

        #Compute the distance between stations
        distance_stn = np.full(df_section_stn_T.index.shape, np.nan)
        for i,stn in enumerate(df_section_stn_T.index):
            distance_stn[i] = azst.haversine(
                df_stn.LON[0],
                df_stn.LAT[0],
                df_stn[df_stn.STATION==df_section_stn_T.index[i]].LON,
                df_stn[df_stn.STATION==df_section_stn_T.index[i]].LAT,
                )

        #Plot the CIL and record stats
        if df_section_stn_T.index.size > 1:
            cil_vol,cil_core = plot_temperature(distance_stn, df_section_stn_T, YEAR, SECTION, SEASON, method='stn')
            cil_vol_stn_clim[idx] = cil_vol
            cil_core_stn_clim[idx] = cil_core

        #Record the results
        df_section_stn_T.index.name = 'station'
        df_section_stn_T.columns.name = 'depth'
        df_stn_temp.append(df_section_stn_T)
        df_section_stn_S.index.name = 'station'
        df_section_stn_S.columns.name = 'depth'
        df_stn_sal.append(df_section_stn_S)

        ## -------- Method 3: Station_ID_manual -------- ##
        df_section_stn_man_T,df_section_stn_man_S = name_variable(ds, stn_list, station_ID_manual.squeeze(), manual=True)

        #Compute the distance between stations
        distance_stn_man = np.full(df_section_stn_man_T.index.shape, np.nan)
        for i,stn in enumerate(df_section_stn_man_T.index):
            distance_stn_man[i] = azst.haversine(
                df_stn.LON[0],
                df_stn.LAT[0],
                df_stn[df_stn.STATION==df_section_stn_man_T.index[i]].LON,
                df_stn[df_stn.STATION==df_section_stn_man_T.index[i]].LAT,
                )

        #Plot the CIL and record stats
        if df_section_stn_man_T.index.size > 1:
            cil_vol,cil_core = plot_temperature(distance_stn_man, df_section_stn_man_T, YEAR, SECTION, SEASON, method='stnman')
            cil_vol_stn_man_clim[idx] = cil_vol
            cil_core_stn_man_clim[idx] = cil_core

        #Record the results
        df_section_stn_man_T.index.name = 'station'
        df_section_stn_man_T.columns.name = 'depth'
        df_stn_man_temp.append(df_section_stn_man_T)
        df_section_stn_man_S.index.name = 'station'
        df_section_stn_man_S.columns.name = 'depth'
        df_stn_man_sal.append(df_section_stn_man_S)


    #Save temperature (not necessary)
    df_stn_mindex = pd.concat(df_stn_temp,keys=years_series)
    df_stn_man_mindex = pd.concat(df_stn_man_temp,keys=years_series)
    df_itp_mindex = pd.concat(df_itp_temp,keys=years_series)

    #Save Salinity (not necessary)
    df_stn_mindex_S = pd.concat(df_stn_sal,keys=years_series)
    df_stn_man_mindex_S = pd.concat(df_stn_man_sal,keys=years_series)
    df_itp_mindex_S = pd.concat(df_itp_sal,keys=years_series)

    #Save Climatology - currently set to stn
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

    ## -------- Fill temperature with climatology for CIL -------- ##
    #We need to add climatology to blank areas
    #This is done retroactively

    #Determine which years are covered in each method
    stn_years = np.unique(df_stn_mindex.index.get_level_values('year').values)
    stn_man_years = np.unique(df_stn_man_mindex.index.get_level_values('year').values)
    itp_years = np.unique(df_itp_mindex.index.get_level_values('year').values)

    #Cycle through each of the years
    for YEAR in years:

        ## -------- Method 1: Interpolation -------- ##
        if np.isin(YEAR, itp_years):
            #Import the climatology
            stn_clim = df_stn_mindex.groupby(level=1).apply(lambda x: x.mean()).T
            #Import the year of interest data
            year_data = df_itp_mindex.T[[YEAR]].copy()
            year_data = year_data.T.droplevel('year').T
            year_data = year_data.reindex(columns=stn_clim.columns, fill_value=np.nan)

            #Determine the new CIL
            cil_vol,cil_core = temperature_clim_fill(stn_clim, year_data, df_stn)

            #Record the CIL vol
            cil_vol_itp_clim[np.where(YEAR == years)[0][0]] = cil_vol
            #Record the CIL core
            cil_core_itp_clim[np.where(YEAR == years)[0][0]] = cil_core

        ## -------- Method 2: Station_ID -------- ##
        if np.isin(YEAR, stn_years):
            #Import the climatology
            stn_clim = df_stn_mindex.groupby(level=1).apply(lambda x: x.mean()).T
            #Import the year of interest data
            year_data = df_stn_mindex.T[[YEAR]].copy()
            year_data = year_data.T.droplevel('year').T
            year_data = year_data.reindex(columns=stn_clim.columns, fill_value=np.nan)

            #Determine the new CIL
            cil_vol,cil_core = temperature_clim_fill(stn_clim, year_data, df_stn)

            #Record the CIL vol
            cil_vol_stn_clim[np.where(YEAR == years)[0][0]] = cil_vol
            #Record the CIL core
            cil_core_stn_clim[np.where(YEAR == years)[0][0]] = cil_core

        ## -------- Method 3: Station_ID_manual -------- ##
        if np.isin(YEAR, stn_years):
            #Import the climatology
            stn_clim = df_stn_man_mindex.groupby(level=1).apply(lambda x: x.mean()).T
            #Import the year of interest data
            year_data = df_stn_man_mindex.T[[YEAR]].copy()
            year_data = year_data.T.droplevel('year').T
            year_data = year_data.reindex(columns=stn_clim.columns, fill_value=np.nan)

            #Determine the new CIL
            cil_vol,cil_core = temperature_clim_fill(stn_clim, year_data, df_stn)

            #Record the CIL vol
            cil_vol_stn_man_clim[np.where(YEAR == years)[0][0]] = cil_vol
            #Record the CIL core
            cil_core_stn_man_clim[np.where(YEAR == years)[0][0]] = cil_core

        print(' -> '+str(YEAR)+' done! ')

    # Save CIL timseries
    df_CIL= pd.DataFrame([cil_vol_stn_clim, cil_vol_stn_man_clim, cil_vol_itp_clim, cil_core_stn_clim, cil_core_stn_man_clim, cil_core_itp_clim]).T
    df_CIL.index = years_series
    df_CIL.columns = ['vol_stn', 'vol_stn_man', 'vol_itp', 'core_stn', 'core_stn_man', 'core_itp']
    picklename = 'df_CIL_' + SECTION + '_' + SEASON + '.pkl'
    df_CIL.to_pickle(picklename)
    os.system('mv temp_section_*.png AZMP_lines/')