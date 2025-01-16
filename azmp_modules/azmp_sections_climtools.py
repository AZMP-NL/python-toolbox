import os
import sys
import netCDF4
import xarray as xr
import warnings
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import cmocean as cm
import pandas as pd
import skfuzzy as fuzz
import numpy as  np
from shapely.geometry import Point,Polygon
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

    #Interpolate the casts vertically before combining
    df = df.interpolate(axis=1,limit_direction='both', limit_area='inside')

    #Create 3D cube
    z = df.columns.values
    V = np.full((lat_reg.size, lon_reg.size, z.size), np.nan)
    #Cycle through each lat and lon
    for i,xx in enumerate(lon_reg):
        for ii,yy in enumerate(lat_reg):
            idx_coords = np.where((lonsV>=xx-dc/2) & (lonsV<xx+dc/2) & (latsV>=yy-dc/2) & (latsV<yy+dc/2))
            var = np.array(df.iloc[idx_coords].mean(axis=0))
            idx_good = np.argwhere((~np.isnan(var)))
            V[ii,i,:] = var
            '''
            #Determine if vertical interpolation is needed
            if np.size(idx_good) == 1:
                V[ii,i,:] = np.array(df.iloc[idx_coords].mean(axis=0))
            elif np.size(idx_good) > 1:
                interp = interp1d(np.squeeze(z[idx_good]), np.squeeze(var[idx_good]))
                idx_interp = np.arange(int(idx_good[0]),int(idx_good[-1]+1))
                #interpolate only where possible (1st to last good idx)
                V[ii,i,idx_interp] = interp(z[idx_interp])
            '''
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


#3D temperature or salinity interpolation IDW
def IDW_variable_3D(df,df_stn,lonsV,latsV,dc):

    #Interpolate the casts vertically before combining
    df = df.interpolate(axis=1,limit_direction='both', limit_area='inside')

    #Set up the dataframe
    df_section_itp = pd.DataFrame(columns=df_stn.STATION.values,index=df.columns.values)

    #Cycle through each station
    for i,stn in enumerate(df_stn.STATION.values):

        #Determine the distance of all points from the station
        stn_lon = df_stn.LON.values[i]
        stn_lat = df_stn.LAT.values[i]
        distance = (stn_lon - lonsV)**2 + (stn_lat - latsV)**2
        distance = np.sqrt(distance)

        #Isolate the measurements of interest
        lonsV_slice = lonsV[distance <= dc]
        latsV_slice = latsV[distance <= dc]
        df_slice = df.iloc[distance <= dc].values
        distance_slice = distance[distance <= dc]

        #Replace 0 distance values with close-to-zero
        distance_slice[distance_slice == 0.] = 1e-6

        #Isolate for each depth
        interpolated_value = np.full(df_slice.shape[1],np.nan)
        for ii in np.arange(interpolated_value.size):

            #Remove the nans
            temp_slice = df_slice[:,ii]
            distance_slice_nonan = distance_slice[~np.isnan(temp_slice)]
            temp_slice = temp_slice[~np.isnan(temp_slice)]

            #If measurements are present, proceed
            if distance_slice_nonan.size > 0:

                #Determine the power
                A = np.pi*(dc**2)
                r_exp = 1/(2*(temp_slice.size/A)**0.5)
                r_obs = np.mean(distance_slice_nonan)
                R = r_obs/r_exp
                R_min = 1
                R_max = 8
                if R < R_min:
                    mu_R = 0
                elif R > R_max:
                    mu_R = 1
                else:
                    mu_R = fuzz.smf(np.array([0,R]),R_min,R_max)[-1]
                if mu_R <= 0.1:
                    power = 0.1
                elif mu_R >= 0.9:
                    power = 3
                else:
                    slope = (3-0.1)/(0.9-0.1)
                    b = 0.1-slope*0.1
                    power = slope*mu_R+b

                #Find the interpolated value
                numerator = np.nansum(temp_slice / (distance_slice_nonan ** power))
                weights = np.sum(1 / (distance_slice_nonan ** power))
                interpolated_value[ii] = numerator/ weights

        #Fill in the data frame
        df_section_itp[stn] = interpolated_value
    return df_section_itp


#Temperature and salinity based upon station ID or station ID manual
def name_variable(ds,z,bins,stn_list,station_ID,station_ID_strip,manual=False):
    #Set up a place to record potential stations
    section_only = []
    #station_ID_strip = np.array([i.strip().split(' ')[-1] for i in station_ID]).astype(str)
    for stn in stn_list:
        #Determine where potential stations are
        if manual:
            if ds.time.size == 1:
                ds_V = ds.isel(time = [station_ID == 'AZMP_NAFC_'+stn])
            else:
                ds_V = ds.isel(time = station_ID == 'AZMP_NAFC_'+stn)
            ds_V['station_ID'] = (['time'], np.tile('AZMP_NAFC_'+stn, ds_V.time.size))
        else:
            if ds.time.size == 1:
                ds_V = ds.isel(time = [station_ID == stn])
                if ds_V.time.size == 0:
                    ds_V = ds.isel(time = [station_ID_strip == stn])
            else:
                ds_V = ds.isel(time = station_ID == stn)
                if ds_V.time.size == 0:
                    ds_V = ds.isel(time = station_ID_strip == stn)
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
    '''
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
        #CIL area
        cil_vol = 0
        for path in c_cil.allsegs[0]:
            vs = path
            if vs.shape[0] != 0:
                cil_vol = cil_vol + np.abs(area(vs))/1000
        cil_core = np.nanmin(df)
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
    fig_name = 'temp_section_' + section + '_' + season + '_'+str(year) + '_'+method+'.png'
    fig.savefig(fig_name, dpi=150)
    plt.close('all')
    return cil_vol, cil_core
    '''
    #TESTING AREA

    #COMBINE WITH CLIMATOLOGY FIRST
    #ONCE COMBINED, THEN RUN INTERPOLATION
    #ENSURE THAT INTERPOLATION IS DONE AT THE BEGINNING OUTSIDE OF LOOP
    #FROM THERE ONLY INTEPROLATE TO POINTS OUTSIDE OF BATH SHAPE FILE

    #Make sure that all the stations are present
    df = df.T
    df = df.reindex(columns=df_stn.STATION.values, fill_value=np.nan).T

    #Compute the new distance vector
    distance = np.full((df.index.shape), np.nan)
    for i, stn in enumerate(df.index):
        distance[i] = azst.haversine(
            df_stn.LON[0],
            df_stn.LAT[0],
            df_stn[df_stn.STATION==df.index[i]].LON.values,
            df_stn[df_stn.STATION==df.index[i]].LAT.values
            )

    ## ---- Retrieve bathymetry using function ---- ##
    bathymetry = azst.section_bathymetry(section, os.path.expanduser('~/github/AZMP-NL/bathymetry/bottom_profiles/'))

    #Interpolate the data to 1km
    dist_interp = np.arange(int(distance.min())+1,int(distance.max()))
    df_var = np.full((df.columns.size,dist_interp.size),np.nan)
    #Cycle through the depths and interpolate horizontally
    for i in np.arange(df_var.shape[0]):
        #Isolate the measurements of interest
        var_slice = df.iloc[:,i].values
        if np.isnan(var_slice).sum() < var_slice.size:
            #Determine which distances can be interpolated
            bath_shape = Polygon(bathymetry)
            bath_filt = np.full(dist_interp.size, True)
            for ii,value in enumerate(bath_filt):
                if bath_shape.contains(Point(df.columns.values[i],dist_interp[ii])):
                    bath_filt[ii] = False

            df_var[i,bath_filt] = np.interp(
                dist_interp[bath_filt],
                distance,
                var_slice)

    #Define the area of each gridcell and calculate CIL
    grid_area = 1*(dz/1000)
    cil_vol = np.sum(df_var<=0)*grid_area
    if np.nanmin(df) <= 0:
        cil_core = np.nanmin(df)
    else:
        cil_core = 0

    #Plot the figure
    fig,ax = plt.subplots()
    c = ax.pcolormesh(
        dist_interp,
        df.columns.values,
        df_var,
        cmap=cm.cm.thermal,vmin=-2,vmax=10,zorder=1)

    #Plot the CIL
    CIL_hatch = np.ma.masked_greater(df_var, 0)
    ax.pcolor(
        dist_interp,
        df.columns.values,
        CIL_hatch,
        hatch='....',
        alpha=0.,zorder=2
        )

    #Plot the bathymetry
    Bgon = plt.Polygon(bathymetry,color=np.multiply([1,.9333,.6667],.4), alpha=0.8)
    ax.add_patch(Bgon)
    ax.set_ylim([0, 400])
    ax.set_xlim([0, np.max(Bgon.xy[:,0])])
    ax.set_ylabel('Depth (m)', fontweight = 'bold')
    ax.set_xlabel('Distance (km)')
    ax.invert_yaxis()
    plt.colorbar(c)
    plt.title(str(year))
    #Save the figure
    fig_name = 'temp_section_' + section + '_' + season + '_'+str(year) + '_'+method+'.png'
    fig.savefig(fig_name, dpi=150)
    plt.close('all')

    return cil_vol, cil_core


#Fill in temperature with climatology
def temperature_clim_fill(clim, year_data, df_stn, dz, year, section, season, method='',plot=False):
    #Merge the climatology and yearly data
    year_orig = year_data.copy()
    year_merged = year_data.copy()
    stn_clim = np.array(clim)
    year_data = np.array(year_data)
    year_data[pd.isnull(year_data)] = stn_clim[pd.isnull(year_data)]
    year_merged.iloc[:,:] = year_data.astype(float)
    year_merged = year_merged.T

    #Compute the new distance vector
    distance_stn_merged = np.full((year_merged.index.shape), np.nan)
    for i, stn in enumerate(year_merged.index):
        distance_stn_merged[i] = azst.haversine(
            df_stn.LON[0],
            df_stn.LAT[0],
            df_stn[df_stn.STATION==year_merged.index[i]].LON.values,
            df_stn[df_stn.STATION==year_merged.index[i]].LAT.values
            )

    ## ---- Retrieve bathymetry using function ---- ##
    bathymetry = azst.section_bathymetry(section, os.path.expanduser('~/github/AZMP-NL/bathymetry/bottom_profiles/'))

    #Interpolate the data to 1km
    dist_interp = np.arange(int(distance_stn_merged.min())+1,int(distance_stn_merged.max()))
    df_var = np.full((year_merged.columns.size,dist_interp.size),np.nan)
    #Cycle through the depths and interpolate horizontally
    for i in np.arange(df_var.shape[0]):
        #Isolate the measurements of interest
        var_slice = year_merged.iloc[:,i].values
        if np.isnan(var_slice).sum() < var_slice.size:
            df_var[i,:] = np.interp(
                dist_interp,
                distance_stn_merged,
                var_slice,
                right=np.nan,left=np.nan)
    #Define the area of each gridcell
    grid_area = 1*(dz/1000)
    cil_vol = np.sum(df_var<=0)*grid_area
    if np.nanmin(df_var) <= 0:
        cil_core = np.nanmin(df_var)
    else:
        cil_core = 0

    #Plot the figure
    if plot:
        fig,ax = plt.subplots()
        ax.pcolormesh(
            dist_interp,
            year_merged.columns.values,
            df_var,
            cmap=cm.cm.thermal,vmin=-2,vmax=10,alpha=0.5,zorder=1)
        c = ax.pcolormesh(
            distance_stn_merged,
            year_orig.T.columns.values,
            year_orig,
            cmap=cm.cm.thermal,vmin=-2,vmax=10,zorder=1)


        #Plot the CIL
        CIL_hatch = np.ma.masked_greater(df_var, 0)
        ax.pcolor(
            dist_interp,
            year_merged.columns.values,
            CIL_hatch,
            hatch='....',
            alpha=0.,zorder=3
            )

        #Plot the bathymetry
        Bgon = plt.Polygon(bathymetry,color=np.multiply([1,.9333,.6667],.4), alpha=0.8)
        ax.add_patch(Bgon)
        ax.set_ylim([0, 400])
        ax.set_xlim([0, np.max(Bgon.xy[:,0])])
        ax.set_ylabel('Depth (m)', fontweight = 'bold')
        ax.set_xlabel('Distance (km)')
        ax.invert_yaxis()
        plt.colorbar(c)
        plt.title(str(year))
        #Save the figure
        fig_name = 'temp_section_' + section + '_' + season + '_'+str(year) + '_'+method+'.png'
        fig.savefig(fig_name, dpi=150)
        plt.close('all')



    return cil_vol, cil_core, year_merged

'''
SECTION = 'FC'
SEASON = 'summer'
YEARS = [1950,2023]
CLIM_YEAR = [1990, 2021]
dlat = 2 # how far from station we search
dlon = 2
z1 = 2
dz = 5 # vertical bins
dc = .2 # grid resolution
CASTS_path = '~/data/CASTS/'
bath_path = '~/data/GEBCO/GEBCO_2023_sub_ice_topo.nc'
'''
def section_clim(SECTION,SEASON,YEARS,CLIM_YEAR,dlat,dlon,z1,dz,dc,CASTS_path,bath_path):

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
    bathy_file = '~/github/AZMP-NL/bathymetry/' + SECTION + '_bathy.npy'
    if os.path.isfile(os.path.expanduser(bathy_file)):
        print('Load saved bathymetry!')
        Zitp = np.load(os.path.expanduser(bathy_file))

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

    #Loop on years
    years = np.arange(YEARS[0], YEARS[1]+1)
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

    #Check to see if files are available
    file_check = [
    'df_stn_mindex_T_'+SECTION+'_'+SEASON+'.pkl',
    'df_stn_man_mindex_T_'+SECTION+'_'+SEASON+'.pkl',
    'df_itp_mindex_T_'+SECTION+'_'+SEASON+'.pkl',
    'df_stn_mindex_S_'+SECTION+'_'+SEASON+'.pkl',
    'df_stn_man_mindex_S_'+SECTION+'_'+SEASON+'.pkl',
    'df_itp_mindex_S_'+SECTION+'_'+SEASON+'.pkl',
    'df_CIL_'+SECTION+'_'+SEASON+'.pkl',
    ]
    data_available = np.array([os.path.isfile('operation_files/'+i) for i in file_check])
    if np.sum(data_available) == data_available.size:
        print('All temperature and salinity data files available!')
        concat = input('Would you like to use existing files and caculate only missing years[y/n]? ')
        if concat == 'y':
            print(' ->  Calculating missing years only.')

            #Import the data
            df_stn_mindex = pd.read_pickle('operation_files/'+file_check[0])
            df_stn_man_mindex = pd.read_pickle('operation_files/'+file_check[1])
            df_itp_mindex = pd.read_pickle('operation_files/'+file_check[2])
            df_stn_mindex_S = pd.read_pickle('operation_files/'+file_check[3])
            df_stn_man_mindex_S = pd.read_pickle('operation_files/'+file_check[4])
            df_itp_mindex_S = pd.read_pickle('operation_files/'+file_check[5])
            df_CIL = pd.read_pickle('operation_files/'+file_check[6])

            #Determine the years present
            years_present = np.unique(df_itp_mindex.index.get_level_values('year').values)
            years = years[~np.isin(years,years_present)]
            years_series = pd.Series(years)
            years_series.name='year'
            cil_vol_stn_clim = np.full(years.shape, np.nan)
            cil_vol_stn_man_clim = np.full(years.shape, np.nan)
            cil_vol_itp_clim = np.full(years.shape, np.nan)
            cil_core_stn_clim = np.full(years.shape, np.nan)
            cil_core_stn_man_clim = np.full(years.shape, np.nan)
            cil_core_itp_clim = np.full(years.shape, np.nan)


        elif concat == 'n':
            print(' ->  Calculating all years')
    else:
        concat = 'n'


    #Cycle through each year
    for idx, YEAR in enumerate(years):

        ## -------- Get CTD data -------- ##
        year_file = CASTS_path + str(YEAR) + '.nc'
        #print('Get ' + year_file)
        ds = xr.open_dataset(year_file)
        ds.load()

        #Select Region
        ds = ds.isel(time = (ds.longitude>lonLims[0]) & (ds.longitude<lonLims[1]))
        ds = ds.isel(time = (ds.latitude>latLims[0]) & (ds.latitude<latLims[1]))
        #ds = ds.isel(time = ds.source == 'NAFC-Oceanography')

        #Take the station ID and station ID manual
        station_ID = ds.station_ID.values.astype(str)
        station_ID_manual = ds.station_ID_manual.values.astype(str)
        station_ID_strip = np.array([i.strip().split(' ')[-1] for i in station_ID]).astype(str)

        #Do the bin average here
        bins = np.arange(z1, 2000, dz)
        warnings.simplefilter(action='ignore', category=FutureWarning)
        ds = ds.groupby_bins('level', bins).mean(dim='level')

        #Select time (save several options here)
        if SEASON == 'summer':
            station_ID = station_ID[(ds['time.month'].values>=6)*(ds['time.month'].values<=9)]
            station_ID_manual = station_ID_manual[(ds['time.month'].values>=6)*(ds['time.month'].values<=9)]
            station_ID_strip = station_ID_strip[(ds['time.month'].values>=6)*(ds['time.month'].values<=9)]
            ds = ds.sel(time=((ds['time.month']>=6)) & ((ds['time.month']<=9)))
        elif SEASON == 'spring':
            station_ID = station_ID[(ds['time.month'].values>=3)*(ds['time.month'].values<=5)]
            station_ID_manual = station_ID_manual[(ds['time.month'].values>=3)*(ds['time.month'].values<=5)]
            station_ID_strip = station_ID_strip[(ds['time.month'].values>=3)*(ds['time.month'].values<=5)]
            ds = ds.sel(time=((ds['time.month']>=3)) & ((ds['time.month']<=5)))
        elif SEASON == 'fall':
            station_ID = station_ID[(ds['time.month'].values>=10)*(ds['time.month'].values<=12)]
            station_ID_manual = station_ID_manual[(ds['time.month'].values>=10)*(ds['time.month'].values<=12)]
            station_ID_strip = station_ID_strip[(ds['time.month'].values>=10)*(ds['time.month'].values<=12)]
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
        #print('Process temperature')
        df = da_temp.to_pandas()
        #Rename columns with 'bins'
        df.columns = bins[0:-1]
        # Remove empty columns.
        idx_empty_rows = df.isnull().all(1).values.nonzero()[0]
        df = df.dropna(axis=0,how='all')
        lonsT = np.delete(lons,idx_empty_rows)
        latsT = np.delete(lats,idx_empty_rows)
        #V_temp,z = interpolate_variable_3D(df, lat_reg, lon_reg, lonsT, latsT, dc)
        if df.shape[0] != 0:
            df_section_itp_T = IDW_variable_3D(df, df_stn, lonsT, latsT, 0.25)
        else:
            df_section_itp_T = pd.DataFrame(index=stn_list, columns=df.columns.values).T
        #print(' -> Done!')

        #Salinity to Pandas Dataframe
        #print('Process salinity')
        df = da_sal.to_pandas()
        #Rename columns with 'bins'
        df.columns = bins[0:-1]
        # Remove empty columns
        idx_empty_rows = df.isnull().all(1).values.nonzero()[0]
        df = df.dropna(axis=0,how='all')
        lonsS = np.delete(lons,idx_empty_rows)
        latsS = np.delete(lats,idx_empty_rows)
        #V_sal,z = interpolate_variable_3D(df, lat_reg, lon_reg, lonsS, latsS, dc)
        if df.shape[0] != 0:
            df_section_itp_S = IDW_variable_3D(df, df_stn, lonsS, latsS, 1)
        else:
            df_section_itp_S = pd.DataFrame(index=stn_list, columns=df.columns.values).T
        #print(' -> Done!')

        #Mask out the bathymetry
        V_coords = np.array([[x,y] for x in lat_reg for y in lon_reg])
        ZZ = Zitp.reshape(V_coords.shape[0],1)
        for stn in stn_list:

            #Determine the index closest to the stn
            station = df_stn[df_stn.STATION==stn]
            idx_opti = np.argmin(np.sum(np.abs(V_coords - np.array([station.LAT.values,station.LON.values]).T), axis=1))
            bottom_depth = -ZZ[idx_opti]

            #Mask out less than 10m, greater than bottom_depth
            stn_depth = df_section_itp_T.index.values
            stn_temp = df_section_itp_T[stn].values
            stn_temp[~((stn_depth >= 10)*(stn_depth <= bottom_depth))] = np.nan
            df_section_itp_T[stn] = stn_temp
            #Same for salinity
            stn_saln = df_section_itp_S[stn].values
            stn_saln[~((stn_depth >= 10)*(stn_depth <= bottom_depth))] = np.nan
            df_section_itp_S[stn] = stn_saln

        df_section_itp_T = df_section_itp_T.T
        df_section_itp_S = df_section_itp_S.T

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

        #Record the results
        df_section_itp_T.index.name = 'station'
        df_section_itp_T.columns.name = 'depth'
        df_itp_temp.append(df_section_itp_T)
        df_section_itp_S.index.name = 'station'
        df_section_itp_S.columns.name = 'depth'
        #Code is temporary for bad 2024 data at BB summer - should be removed eventually
        if SECTION == 'BB':
            if SEASON == 'summer':
                if YEAR == 2024:
                    df_section_itp_S.iloc[:2,:10] = np.nan
        df_itp_sal.append(df_section_itp_S)

        ## -------- Method 2: Station_ID -------- ##
        z = df.columns.values
        df_section_stn_T,df_section_stn_S = name_variable(ds, z, bins, stn_list, station_ID.squeeze(), station_ID_strip.squeeze(), manual=False)

        #Compute the distance between stations
        distance_stn = np.full(df_section_stn_T.index.shape, np.nan)
        for i,stn in enumerate(df_section_stn_T.index):
            distance_stn[i] = azst.haversine(
                df_stn.LON[0],
                df_stn.LAT[0],
                df_stn[df_stn.STATION==df_section_stn_T.index[i]].LON,
                df_stn[df_stn.STATION==df_section_stn_T.index[i]].LAT,
                )

        #Record the results
        df_section_stn_T.index.name = 'station'
        df_section_stn_T.columns.name = 'depth'
        df_stn_temp.append(df_section_stn_T)
        df_section_stn_S.index.name = 'station'
        df_section_stn_S.columns.name = 'depth'
        #Code is temporary for bad 2024 data at BB summer - should be removed eventually
        if SECTION == 'BB':
            if SEASON == 'summer':
                if YEAR == 2024:
                    df_section_stn_S.iloc[:2,:10] = np.nan
        df_stn_sal.append(df_section_stn_S)

        ## -------- Method 3: Station_ID_manual -------- ##
        df_section_stn_man_T,df_section_stn_man_S = name_variable(ds, z, bins, stn_list, station_ID_manual.squeeze(), station_ID_strip.squeeze(), manual=True)

        #Compute the distance between stations
        distance_stn_man = np.full(df_section_stn_man_T.index.shape, np.nan)
        for i,stn in enumerate(df_section_stn_man_T.index):
            distance_stn_man[i] = azst.haversine(
                df_stn.LON[0],
                df_stn.LAT[0],
                df_stn[df_stn.STATION==df_section_stn_man_T.index[i]].LON,
                df_stn[df_stn.STATION==df_section_stn_man_T.index[i]].LAT,
                )

        #Record the results
        df_section_stn_man_T.index.name = 'station'
        df_section_stn_man_T.columns.name = 'depth'
        df_stn_man_temp.append(df_section_stn_man_T)
        df_section_stn_man_S.index.name = 'station'
        df_section_stn_man_S.columns.name = 'depth'
        #Code is temporary for bad 2024 data at BB summer - should be removed eventually
        if SECTION == 'BB':
            if SEASON == 'summer':
                if YEAR == 2024:
                    df_section_stn_man_S.iloc[:2,:10] = np.nan
        df_stn_man_sal.append(df_section_stn_man_S)




    #Concatenate all temperature measurements
    if concat=='y':
        if np.size(df_stn_temp) != 0:
            place_stn = pd.concat(df_stn_temp,keys=years_series)
            df_stn_mindex = pd.concat([place_stn,df_stn_mindex]).sort_index()
            df_stn_mindex = df_stn_mindex[~df_stn_mindex.index.duplicated(keep='first')]
        if np.size(df_stn_man_temp) != 0:
            place_stn_man = pd.concat(df_stn_man_temp,keys=years_series)
            df_stn_man_mindex = pd.concat([place_stn_man,df_stn_man_mindex]).sort_index()
            df_stn_man_mindex = df_stn_man_mindex[~df_stn_man_mindex.index.duplicated(keep='first')]
        if np.size(df_itp_temp) != 0:
            place_itp = pd.concat(df_itp_temp,keys=years_series)
            df_itp_mindex = pd.concat([place_itp,df_itp_mindex]).sort_index()
            df_itp_mindex = df_itp_mindex[~df_itp_mindex.index.duplicated(keep='first')]
    elif concat=='n':
        df_stn_mindex = pd.concat(df_stn_temp,keys=years_series)
        df_stn_man_mindex = pd.concat(df_stn_man_temp,keys=years_series)
        df_itp_mindex = pd.concat(df_itp_temp,keys=years_series)
    #Save all three versions
    df_stn_mindex.to_pickle('df_stn_mindex_T_'+SECTION+'_'+SEASON+'.pkl')
    df_stn_man_mindex.to_pickle('df_stn_man_mindex_T_'+SECTION+'_'+SEASON+'.pkl')
    df_itp_mindex.to_pickle('df_itp_mindex_T_'+SECTION+'_'+SEASON+'.pkl')

    #Concatenate all salinity measurements
    if concat=='y':
        if np.size(df_stn_sal) != 0:
            place_stn = pd.concat(df_stn_sal,keys=years_series)
            df_stn_mindex_S = pd.concat([place_stn,df_stn_mindex_S]).sort_index()
            df_stn_mindex_S = df_stn_mindex_S[~df_stn_mindex_S.index.duplicated(keep='first')]
        if np.size(df_stn_man_sal) != 0:
            place_stn_man = pd.concat(df_stn_man_sal,keys=years_series)
            df_stn_man_mindex_S = pd.concat([place_stn_man,df_stn_man_mindex_S]).sort_index()
            df_stn_man_mindex_S = df_stn_man_mindex_S[~df_stn_man_mindex_S.index.duplicated(keep='first')]
        if np.size(df_itp_sal) != 0:
            place_itp = pd.concat(df_itp_sal,keys=years_series)
            df_itp_mindex_S = pd.concat([place_itp,df_itp_mindex_S]).sort_index()
            df_itp_mindex_S = df_itp_mindex_S[~df_itp_mindex_S.index.duplicated(keep='first')]
    elif concat=='n':
        df_stn_mindex_S = pd.concat(df_stn_sal,keys=years_series)
        df_stn_man_mindex_S = pd.concat(df_stn_man_sal,keys=years_series)
        df_itp_mindex_S = pd.concat(df_itp_sal,keys=years_series)
    #Save all three versions
    df_stn_mindex_S.to_pickle('df_stn_mindex_S_'+SECTION+'_'+SEASON+'.pkl')
    df_stn_man_mindex_S.to_pickle('df_stn_man_mindex_S_'+SECTION+'_'+SEASON+'.pkl')
    df_itp_mindex_S.to_pickle('df_itp_mindex_S_'+SECTION+'_'+SEASON+'.pkl')

    #Save Climatology - currently set to stn
    #Isolate for climatology years
    stn_years = df_stn_mindex.index.get_level_values('year').values
    clim_filt = (stn_years>=CLIM_YEAR[0])*(stn_years<=CLIM_YEAR[1])
    df_stn_mindex_clim = df_stn_mindex.iloc[clim_filt]
    df_clim =  df_stn_mindex_clim.groupby(level=1).apply(lambda x: x.mean())
    picklename = 'df_temperature_' + SECTION + '_' + SEASON + '_clim.pkl'
    df_clim.to_pickle(picklename)
    stn_years = df_stn_mindex_S.index.get_level_values('year').values
    clim_filt = (stn_years>=CLIM_YEAR[0])*(stn_years<=CLIM_YEAR[1])
    df_stn_mindex_S_clim = df_stn_mindex_S.iloc[clim_filt]
    df_clim_S =  df_stn_mindex_S_clim.groupby(level=1).apply(lambda x: x.mean())
    picklename = 'df_salinity_' + SECTION + '_' + SEASON + '_clim.pkl'
    df_clim_S.to_pickle(picklename)

    ## -------- Fill temperature with climatology for CIL -------- ##
    #We need to add climatology to blank areas
    #This is done retroactively

    #Determine which years are covered in each method
    stn_years = np.unique(df_stn_mindex.index.get_level_values('year').values)
    stn_man_years = np.unique(df_stn_man_mindex.index.get_level_values('year').values)
    itp_years = np.unique(df_itp_mindex.index.get_level_values('year').values)
    years = np.arange(YEARS[0],YEARS[1]+1)
    cil_vol_stn_clim = np.full(years.shape, np.nan)
    cil_vol_stn_man_clim = np.full(years.shape, np.nan)
    cil_vol_itp_clim = np.full(years.shape, np.nan)
    cil_core_stn_clim = np.full(years.shape, np.nan)
    cil_core_stn_man_clim = np.full(years.shape, np.nan)
    cil_core_itp_clim = np.full(years.shape, np.nan)
    df_section_itp_meanT = np.full(years.shape, np.nan)
    df_section_stn_meanT = np.full(years.shape, np.nan)
    df_section_stn_man_meanT = np.full(years.shape, np.nan)
    if SECTION == 'FC':
        df_section_itp_meanT_shelf = np.full(years.shape, np.nan)
        df_section_stn_meanT_shelf = np.full(years.shape, np.nan)
        df_section_stn_man_meanT_shelf = np.full(years.shape, np.nan)
        df_section_itp_meanT_cap = np.full(years.shape, np.nan)
        df_section_stn_meanT_cap = np.full(years.shape, np.nan)
        df_section_stn_man_meanT_cap = np.full(years.shape, np.nan)


    #Cycle through each of the years
    for YEAR in years:

        ## -------- Method 1: Interpolation -------- ##
        if np.isin(YEAR, itp_years):
            #Import the climatology
            stn_clim = df_stn_mindex_clim.groupby(level=1).apply(lambda x: x.mean()).T
            #Import the year of interest data
            year_data = df_itp_mindex.T[[YEAR]].copy()
            year_data = year_data.T.droplevel('year').T
            num_of_stns = year_data.columns.size
            year_data = year_data.reindex(columns=stn_clim.columns, fill_value=np.nan)

            #Determine the new CIL
            cil_vol,cil_core,year_merged = temperature_clim_fill(stn_clim, year_data, df_stn, dz, YEAR, SECTION, SEASON, method='itp')

            #Record the CIL vol
            cil_vol_itp_clim[np.where(YEAR == years)[0][0]] = cil_vol
            #Record the CIL core
            cil_core_itp_clim[np.where(YEAR == years)[0][0]] = cil_core

            #Record the depth-station average
            if num_of_stns/year_merged.index.size >= 0.5:
                df_section_itp_meanT[np.where(YEAR == years)[0][0]] = year_merged.mean(axis=1).mean()
            if SECTION == 'FC':
                df_section_itp_meanT_shelf[np.where(YEAR == years)[0][0]] = year_merged.iloc[:20].mean(axis=1).mean()
                df_section_itp_meanT_cap[np.where(YEAR == years)[0][0]] = year_merged.iloc[14:].mean(axis=1).mean()

        ## -------- Method 2: Station_ID -------- ##
        if np.isin(YEAR, stn_years):
            #Import the climatology
            stn_clim = df_stn_mindex_clim.groupby(level=1).apply(lambda x: x.mean()).T
            #Import the year of interest data
            year_data = df_stn_mindex.T[[YEAR]].copy()
            year_data = year_data.T.droplevel('year').T
            year_data = year_data.reindex(columns=stn_clim.columns, fill_value=np.nan)

            #Determine the new CIL
            cil_vol,cil_core, year_merged = temperature_clim_fill(stn_clim, year_data, df_stn, dz, YEAR, SECTION, SEASON, method='stn')

            #Record the CIL vol
            cil_vol_stn_clim[np.where(YEAR == years)[0][0]] = cil_vol
            #Record the CIL core
            cil_core_stn_clim[np.where(YEAR == years)[0][0]] = cil_core

            #Record the depth-station average
            if num_of_stns/year_merged.index.size >= 0.5:
                df_section_stn_meanT[np.where(YEAR == years)[0][0]] = year_merged.mean(axis=1).mean()
            if SECTION == 'FC':
                df_section_stn_meanT_shelf[np.where(YEAR == years)[0][0]] = year_merged.iloc[:20].mean(axis=1).mean()
                df_section_stn_meanT_cap[np.where(YEAR == years)[0][0]] = year_merged.iloc[14:].mean(axis=1).mean()

        ## -------- Method 3: Station_ID_manual -------- ##
        if np.isin(YEAR, stn_man_years):
            #Import the climatology
            stn_clim = df_stn_mindex_clim.groupby(level=1).apply(lambda x: x.mean()).T
            #Import the year of interest data
            year_data = df_stn_man_mindex.T[[YEAR]].copy()
            year_data = year_data.T.droplevel('year').T
            year_data = year_data.reindex(columns=stn_clim.columns, fill_value=np.nan)

            #Determine the new CIL
            cil_vol,cil_core,year_merged = temperature_clim_fill(stn_clim, year_data, df_stn, dz, YEAR, SECTION, SEASON, method='stnman')

            #Record the CIL vol
            cil_vol_stn_man_clim[np.where(YEAR == years)[0][0]] = cil_vol
            #Record the CIL core
            cil_core_stn_man_clim[np.where(YEAR == years)[0][0]] = cil_core

            #Record the depth-station average
            if num_of_stns/year_merged.index.size >= 0.5:
                df_section_stn_man_meanT[np.where(YEAR == years)[0][0]] = year_merged.mean(axis=1).mean()
            if SECTION == 'FC':
                df_section_stn_man_meanT_shelf[np.where(YEAR == years)[0][0]] = year_merged.iloc[:20].mean(axis=1).mean()
                df_section_stn_man_meanT_cap[np.where(YEAR == years)[0][0]] = year_merged.iloc[14:].mean(axis=1).mean()

        #print(' -> CIL climatology fill: '+str(YEAR)+' done! ')

    # Save CIL timseries
    df_CIL= pd.DataFrame([cil_vol_stn_clim, cil_vol_stn_man_clim, cil_vol_itp_clim, cil_core_stn_clim, cil_core_stn_man_clim, cil_core_itp_clim]).T
    df_CIL.index = years
    df_CIL.index.name = 'year'
    df_CIL.columns = ['vol_stn', 'vol_stn_man', 'vol_itp', 'core_stn', 'core_stn_man', 'core_itp']
    picklename = 'df_CIL_' + SECTION + '_' + SEASON + '_climfill.pkl'
    df_CIL.to_pickle(picklename)

    # Save CIL timeseries for IROC
    df = pd.DataFrame([cil_vol_stn_clim, cil_vol_itp_clim]).T
    df.index = years
    df.index.name = 'year'
    df.columns = ['station-ID', 'interp_field']
    outfile = 'CIL_area_' + SECTION + '.csv'
    df.to_csv(outfile)

    #Save the meanT
    df_meanT= pd.DataFrame([df_section_itp_meanT,df_section_stn_meanT,df_section_stn_man_meanT]).T
    df_meanT.index = years
    df_meanT.index.name = 'year'
    df_meanT.columns = ['itp_meanT','stn_meanT','stn_man_meanT']
    picklename = 'df_' + SECTION + '_meanT_' + SEASON + '.pkl'
    df_meanT.to_pickle(picklename)
    if SECTION == 'FC':
        df_meanT= pd.DataFrame([df_section_itp_meanT_shelf,df_section_stn_meanT_shelf,df_section_stn_man_meanT_shelf]).T
        df_meanT.index = years
        df_meanT.index.name = 'year'
        df_meanT.columns = ['itp_meanT_shelf','stn_meanT_shelf','stn_man_meanT_shelf']
        picklename = 'df_' + SECTION + '_meanT_shelf_' + SEASON + '.pkl'
        df_meanT.to_pickle(picklename)
        df_meanT= pd.DataFrame([df_section_itp_meanT_cap,df_section_stn_meanT_cap,df_section_stn_man_meanT_cap]).T
        df_meanT.index = years
        df_meanT.index.name = 'year'
        df_meanT.columns = ['itp_meanT_cap','stn_meanT_cap','stn_man_meanT_cap']
        picklename = 'df_' + SECTION + '_meanT_cap_' + SEASON + '.pkl'
        df_meanT.to_pickle(picklename)
