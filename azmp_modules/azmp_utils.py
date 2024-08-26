"""Some utilitary tools to for AZMP

Contains following functions:
- get_nafo_divisions()
- get_bottomT_climato(INFILES, LON_REG,  LAT_REG, year_lims=[1981, 2010], season=[], zlims=[10, 1000], dz=5, h5_outputfile=[])
- get_bottomS_climato(INFILES, LON_REG,  LAT_REG, year_lims=[1981, 2010], season=[], zlims=[10, 1000], dz=5, h5_outputfile=[])
- get_bottomT(year_file, season, climato_file):
- get_bottomS(year_file, season, climato_file):
- bottomT_quickplot(h5_outputfile, figure_file=[])
- Tbot_to_GIS_ascii(h5file, ascfile)
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
from netCDF4 import date2num,num2date
import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr
from scipy.interpolate import griddata
from scipy.interpolate import interp1d  # to remove NaNs in profiles
from scipy.interpolate.interpnd import _ndim_coords_from_arrays
from scipy.spatial import cKDTree
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from shapely.ops import cascaded_union
import shapefile 
from area import area # external fns to compute surface area
#from seawater import extras as swx
import datetime
# maps
#os.environ['PROJ_LIB'] = '/home/jcoyne/anaconda3/envs/jcoyne_general/share/proj'
#from mpl_toolkits.basemap import Basemap
import cartopy.crs as ccrs
import cartopy.feature as cfeature


    
def get_nafo_divisions():
    """ Will generate a dict with NAFO division shapes.
    Example to access info:
    In [14]: dict['2J']['lat']
    Out[14]: [55.33, 55.33, 52.25, 52.25]

    ** Needs to be expanded with other divisions **
    (see map_extract_nafo.py for help getting polygons)

    """

    # 0B
    xlon = [-62.466, -57.634, -57.8, -57.818, -58.033, -57.557, -57.352, -57.218, -57.218, -57.218, -59, -65, -66.259]
    xlat = [66.25, 66.25, 64.25, 64.072, 63.582, 62.5, 62, 61, 60.75, 60.2, 61, 61, 61.866]
    div0B = {'lat' : xlat, 'lon' : xlon}

    # 1C
    xlon = [-53.342, -57.634, -57.8, -52.044]
    xlat = [66.25, 66.25, 64.25, 64.25]
    div1C = {'lat' : xlat, 'lon' : xlon}

    # 1D
    xlon = [-52.044, -57.8, -57.818, -58.033, -57.557, -50.037]
    xlat = [64.25, 64.25, 64.072, 63.582, 62.5, 62.5]
    div1D = {'lat' : xlat, 'lon' : xlon}

    # 1E
    xlon = [-50.037, -57.557, -57.352, -57.218, -46.153]
    xlat = [62.5, 62.5, 62, 60.75, 60.75]
    div1E = {'lat' : xlat, 'lon' : xlon}

    # 1F
    xlon = [-46.153, -57.218, -57.2187, -42, -42, -44, -44]
    xlat = [60.75, 60.75, 60.2, 52.25, 59, 59, 60.157]
    div1F = {'lat' : xlat, 'lon' : xlon}    

    # 2G
    xlon = [-64.5, -64.5, -59, -52, -61.9]
    xlat = [60.3, 61, 61, 57.667, 57.667]
    div2G = {'lat' : xlat, 'lon' : xlon}

    # 2H
    xlon = [-61.924, -52.010, -47.531, -60.382]
    xlat = [57.667, 57.667, 55.333, 55.333]
    div2H = {'lat' : xlat, 'lon' : xlon}

    # 2J
    xlon = [-59.778, -47.529, -42, -55.712]
    xlat = [55.33, 55.33, 52.25, 52.25]
    div2J = {'lat' : xlat, 'lon' : xlon}

    # 3K (no land)
    xlon = [-55.425, -55.425, -42, -42, -53.466]
    xlat = [51.583, 52.25, 52.25, 49.25, 49.25]
    div3K = {'lat' : xlat, 'lon' : xlon}

    # 3K-xtended (point on land to close the box)
    xlon = [-55.425, -55.425, -42, -42, -53.466, -57.5]
    xlat = [51.583, 52.25, 52.25, 49.25, 49.25, 49.25]
    div3Kx = {'lat' : xlat, 'lon' : xlon}

    # 3L
    xlon = [-53.466, -46.5, -46.5, -54.5, -54.2]
    xlat = [49.25, 49.25, 46, 46, 46.815]
    div3L = {'lat' : xlat, 'lon' : xlon}

    # 3N
    xlon = [-51, -46.5, -46.5, -50, -51, -51]
    xlat = [46, 46, 39, 39, 39.927, 46]
    div3N = {'lat' : xlat, 'lon' : xlon}

    # 3M
    xlon = [-42, -46.5, -46.5, -42, -42]
    xlat = [49.25, 49.25, 39, 39, 49.25]
    div3M = {'lat' : xlat, 'lon' : xlon}

    # 3O
    xlon = [-54.5, -51, -51, -54.5, -54.5]
    xlat = [46, 46, 39.927, 43.064, 46]
    div3O = {'lat' : xlat, 'lon' : xlon}

    #3Pn
    xlon = [-57.523, -59.308, -59.569, -58.817]
    xlat = [47.631, 47.620, 47.477, 46.843]
    div3Pn = {'lat' : xlat, 'lon' : xlon}
    
    #3Ps
    xlon = [-57.523, -58.82, -54.5, -54.5, -54.0]
    xlat = [47.631, 46.843, 43.064, 46, 47.7]
    div3Ps = {'lat' : xlat, 'lon' : xlon}

    #4R
    xlon = [-59.308, -59.571, -60, -60, -57.114]
    xlat = [47.62, 47.473, 47.883, 49.418, 51.412]
    div4R = {'lat' : xlat, 'lon' : xlon}

    #4S
    xlon = [-60, -64.667, -67.304, -67.304, -57.130, -60]
    xlat = [47.835, 49.417, 49.417, 51.4, 51.4, 49.416]
    div4S = {'lat' : xlat, 'lon' : xlon}

    #4T
    xlon = [-60, -64.667, -67.304, -65.409, -61.384, -60.413]
    xlat = [47.835, 49.416, 49.416, 47.075, 45.583, 47.039]
    div4T = {'lat' : xlat, 'lon' : xlon}

    #4Vn
    xlon = [-60.41, -60., -57.45, -60.]
    xlat = [47.05, 47.83, 45.67, 45.67]
    div4Vn = {'lat' : xlat, 'lon' : xlon}

    #4Vs
    xlon = [-60.18, -60.15, -57.45, -50, -59, -59., -60, -60]
    xlat = [45.78, 45.70, 45.67, 39, 39, 44.17, 44.17, 45.67]
    div4Vs = {'lat' : xlat, 'lon' : xlon}

    #4W
    xlon = [-63.5, -63.33, -63.33, -59, -59, -60, -60]
    xlat = [44.48, 44.33, 39, 39, 44.17, 44.17, 45.67]
    div4W = {'lat' : xlat, 'lon' : xlon}

    #4X
    xlon = [-66.91, -66.91, -67.41, -67.74, -67.3, -66, -65.67, -65.67, -63.33, -63.33, -63.5] 
    xlat = [45.05, 43.84, 43.84, 42.88, 42.33, 42.33, 42, 39, 39, 44.33, 44.48]
    div4X = {'lat' : xlat, 'lon' : xlon}

    #5Y
    xlon = [-71, -66.903, -66.903, -67.407, -67.744, -67.303, -70, -70, -71]
    xlat = [45.059, 45.059, 43.833, 43.833, 42.887, 42.334, 42.334, 41.7, 41.7]
    div5Y = {'lat' : xlat, 'lon' : xlon}

    #5Ze
    xlon = [-70, -70, -66, -65.667, -65.667, -70, -70]
    xlat = [41.967, 42.333, 42.333, 42.0, 39.0, 39.0, 41.667]
    div5Ze = {'lat' : xlat, 'lon' : xlon}   
    
    #5Zw
    xlon = [-71.667, -71.667, -70, -70]
    xlat = [41.346, 39, 39, 41.667]
    div5Zw = {'lat' : xlat, 'lon' : xlon}
               
    dict = {}
    dict['0B'] = div0B
    dict['1C'] = div1C
    dict['1D'] = div1D
    dict['1E'] = div1E
    dict['1F'] = div1F
    dict['2G'] = div2G
    dict['2H'] = div2H
    dict['2J'] = div2J
    dict['3K'] = div3K
    dict['3Kx'] = div3Kx
    dict['3L'] = div3L
    dict['3N'] = div3N
    dict['3M'] = div3M
    dict['3O'] = div3O
    dict['3Ps'] = div3Ps
    dict['3Pn'] = div3Pn
    dict['3O'] = div3O
    dict['4R'] = div4R
    dict['4S'] = div4S
    dict['4T'] = div4T
    dict['4Vn'] = div4Vn
    dict['4Vs'] = div4Vs
    dict['4W'] = div4W
    dict['4X'] = div4X
    dict['5Y'] = div5Y
    dict['5Ze'] = div5Ze
    dict['5Zw'] = div5Zw

    return dict


def build_NLshelf_definition(dataFile = '/home/jcoyne/Documents/Datasets/GEBCO_2023/GEBCO_2023_sub_ice_topo.nc' ):
    """ Will generate a dict with NL shelf 1000m limit shape.

    Example to access info:
    >> dict['lat']

    """

    #REMOVE AFTER
    dataFile = '/home/jcoyne/Documents/Datasets/GEBCO_2023/GEBCO_2023_sub_ice_topo.nc'


    ## ---- Get Bathymetry ---- ####
    dc = .1
    lonLims = [-60, -45] # FC AZMP report region
    latLims = [42, 56]
    lon_reg = np.arange(lonLims[0]+dc/2, lonLims[1]-dc/2, dc)
    lat_reg = np.arange(latLims[0]+dc/2, latLims[1]-dc/2, dc)
    lon_grid, lat_grid = np.meshgrid(lon_reg,lat_reg)
    # Load data
    ##dataset = netCDF4.Dataset(dataFile)
    dataset = xr.open_dataset(dataFile)
    x = [-179-59.75/60, 179+59.75/60] # to correct bug in 30'' dataset?
    y = [-89-59.75/60, 89+59.75/60]
    ##spacing = dataset.variables['spacing']


    # Compute Lat/Lon
    ##nx = int((x[-1]-x[0])/spacing[0]) + 1  # num pts in x-dir
    ##ny = int((y[-1]-y[0])/spacing[1]) + 1  # num pts in y-dir
    ##lon = np.linspace(x[0],x[-1],nx)
    ##lat = np.linspace(y[0],y[-1],ny)
    # interpolate data on regular grid (temperature grid)
    # Reshape data
    ##zz = dataset.variables['z']
    ##zz = dataset.variables['elevation']

    ##Z = zz[:].reshape(ny, nx)
    ##Z = np.flipud(Z) # <------------ important!!!
    

    #Gather lat/lon
    lon = dataset.lon.values
    lat = dataset.lat.values
    nx = lon.size
    ny = lat.size


    # Reduce data according to Region params
    ##idx_lon = np.where((lon>=lonLims[0]) & (lon<=lonLims[1]))
    ##idx_lat = np.where((lat>=latLims[0]) & (lat<=latLims[1]))


    dataset = dataset.isel(lon=(lon>=lonLims[0])*(lon<=lonLims[1]), lat=(lat>=latLims[0])*(lat<=latLims[1]))
    Z = dataset.elevation.values
    lon = lon[(lon>=lonLims[0])*(lon<=lonLims[1])]
    lat = lat[(lat>=latLims[0])*(lat<=latLims[1])]

    #See if there's a way to index instead of interpolate
    #lon[11::24]
    #lat[11::24]
    Zitp = Z[11::24,11::24]
    lon_grid,lat_grid = np.meshgrid(lon[11::24],lat[11::24])


    ##Z = Z[idx_lat[0][0]:idx_lat[0][-1]+1, idx_lon[0][0]:idx_lon[0][-1]+1]
    ##lon = lon[idx_lon[0]]
    ##lat = lat[idx_lat[0]]
    # interpolate data on regular grid (temperature grid)
    ##lon_grid_bathy, lat_grid_bathy = np.meshgrid(lon,lat)
    ##lon_vec_bathy = np.reshape(lon_grid_bathy, lon_grid_bathy.size)
    ##lat_vec_bathy = np.reshape(lat_grid_bathy, lat_grid_bathy.size)
    ##z_vec = np.reshape(Z, Z.size)
    ##Zitp = griddata((lon_vec_bathy, lat_vec_bathy), z_vec, (lon_grid, lat_grid), method='linear')

    # Use matplotlib contour to extract 1000m isobath
    cc = plt.contour(lon_reg, lat_reg, -Zitp, [1000])
    c1000  = cc.allsegs[0][1]
    lon_vec = c1000[:,0]
    lat_vec = c1000[:,1]
    idx_drop = np.where((lon_vec < -56) & (lat_vec < 46))
    lon_vec = np.delete(lon_vec, idx_drop)
    lat_vec = np.delete(lat_vec, idx_drop)
    # add some points in land to close the loop
    lon_vec = np.append(lon_vec, lon_vec[-1])
    lat_vec = np.append(lat_vec, 48)
    lon_vec = np.append(lon_vec, -57)
    lat_vec = np.append(lat_vec, 48)
    lon_vec = np.append(lon_vec, -57)
    lat_vec = np.append(lat_vec, 48)

    contour_mask = np.load('/home/cyrf0006/AZMP/state_reports/bottomT/100m_contour_labrador.npy')
    contour_mask = contour_mask[1:-2, :]
    lon_vec = np.append(lon_vec, np.flipud(contour_mask[:,0]))
    lat_vec = np.append(lat_vec, np.flipud(contour_mask[:,1]))
    lon_vec = np.append(lon_vec, lon_vec[0])
    lat_vec = np.append(lat_vec, lat_vec[0])

    A = np.stack((lon_vec, lat_vec)).T

    return A

def get_NLshelf(infile):
    """ If infile exists, it will load a dict with NLshelf. If not, it will create it.
    usage example"
    import azmp_utils as azu
    NLshelf = azu.get_NLshelf('/home/cyrf0006/github/AZMP-NL/data/NLshelf_definition.npy')

    """

    if os.path.isfile(infile):
        print(infile + ' exist! Reading directly')
        A = np.load(infile)
    else:
        print(infile + ' does not exist! Let''s buil it...')
        A = build_NLshelf_definition()
        np.save(infile, A)
        print ['  -> Done!']
  
    return A

def get_ice_regions():
    """ Will generate a dict with Galbraiths ice region shapes.
    Example to access info:
    In [14]: dict['NFLD']['lat']
    Out[14]: [55.33, 55.33, 52.25, 52.25]

    """

    # NLab
    xlon = [-64.67084639498432, -62.2, -55, -55, -64.67084639498432]
    xlat = [60.5, 56.3333, 56.3333, 60.5, 60.5]
    NLab = {'lat' : xlat, 'lon' : xlon}

    # SLab
    xlon = [-56, -56, -58.225, -60.15, -61.5, -62.2, -40, -40, -56]
    xlat = [52.25, 53.383, 54.3333, 54.3333, 56.0, 56.3333, 56.3333, 52.25, 52.25]
    SLab = {'lat' : xlat, 'lon' : xlon}

    # Nfld
    xlon = [-55, -55, -57.1, -56, -55.4166, -55.4166, -40, -40, -55]
    xlat = [43, 48, 49.7, 51.5, 51.5833, 52.25, 52.25, 43, 43]
    Nfld = {'lat' : xlat, 'lon' : xlon}

               
    dict = {}
    dict['NLab'] = NLab
    dict['SLab'] = SLab
    dict['Nfld'] = Nfld

    return dict


def get_bottomT_climato(
    INFILES,
    lonLims,
    latLims,
    bath_file,
    time_adjust=True,
    year_lims=[1991, 2020],
    season=[],
    h5_outputfile=[]):
    '''
    Generates and returns the bottom temperature climatology.
    This script uses GEBCO data, path needs to be provided.
    If file already exists, it will be read and the output will come from that file.

    INFILES: Path location of bottom temperature .nc files. 
    lonLims (list): Longitude restraints [longitude min, longitude max]
    latLims (list): Latitude restraints [latitude min, latitude max]
    years_lim (list): Climatology years, default is [1991,2020]
    season: The season of interest; spring, summer, fall
    h5_outputfile: Location of output file, type .10.h5
    bath_file: Path location of GEBCO bathymetry data, including file name.

    Examples.
    #fish_hab region: lonLims = [-60, -43], latLims = [39, 56], h5_outputfile='Tbot_climato_fall_0.10.h5'
    #FC AZMP report region: lonLims = [-60, -45], latLims = [42, 56], h5_outputfile='Tbot_climato_fall_0.10.h5'
    #Include 2H in above: lonLims = [-63, -45], latLims = [42, 58], h5_outputfile='Tbot_climato_fall_0.10.h5'
    #NSRF area: lonLims = [-70, -56], latLims = [57, 67], year_lims=[2006, 2021], season='summer', h5_outputfile='Tbot_climato_NSRFx_summer_2006-2021.h5'
    # NSRFxx, down to 53N (for Normandeau's work): lonLims = [-70, -56], latLims = [53, 67]
    '''


    # Because the code is trying to divide by zero
    # (https://stackoverflow.com/questions/14861891/runtimewarning-invalid-value-encountered-in-divide)
    np.seterr(divide='ignore', invalid='ignore')

    ## ---- Check if H5 file exists ---- ##        
    if os.path.isfile(h5_outputfile):
        print(h5_outputfile + ' exist! Reading directly')
        h5f = h5py.File(h5_outputfile,'r')
        ds_temp = h5f['Tbot'][:]
        lons = h5f['lon_reg'][:]
        lats = h5f['lat_reg'][:]
        Zitp = h5f['Zitp'][:]
        h5f.close()

    else:

        #Load in the bottom temperature data
        ds = xr.open_dataset(os.path.expanduser(INFILES))
        #Isolate for the time of interest
        ds = ds.sel(TIME=((ds['TIME.year']>=year_lims[0])*(ds['TIME.year']<=year_lims[1])))
        ds = ds.mean('TIME')

        #Isolate for the region of interest
        ds = ds.sel(X=((ds.LONGITUDE[0,:]>=lonLims[0])*(ds.LONGITUDE[0,:]<=lonLims[1])).values)
        ds = ds.sel(Y=((ds.LATITUDE[:,0]>=latLims[0])*(ds.LATITUDE[:,0]<=latLims[1])).values)

        #Record the temperature, latitude, longitude
        if time_adjust:
            ds_temp = ds.BOTTOM_TEMPERATURE_ADJUSTED.values
        else:
            ds_temp = ds.BOTTOM_TEMPERATURE.values
        lons = ds.LONGITUDE.values
        lats = ds.LATITUDE.values

        ## ---- Region parameters ---- ##
        ds_bath = xr.open_dataset(os.path.expanduser(bath_file))
        ds_bath = ds_bath.isel(lon=(ds_bath.lon>=lonLims[0])*(ds_bath.lon<=lonLims[1]))
        ds_bath = ds_bath.isel(lat=(ds_bath.lat>=latLims[0])*(ds_bath.lat<=latLims[1]))
        Zitp = ds_bath.elevation[::30,::30].values

        #Save the climatology
        h5f = h5py.File(h5_outputfile, 'w')
        h5f.create_dataset('Tbot', data=ds_temp)
        h5f.create_dataset('lon_reg', data=lons)
        h5f.create_dataset('lat_reg', data=lats)
        h5f.create_dataset('Zitp', data=Zitp)
        h5f.close()

    #Fill dict for output
    dict = {}
    dict['Tbot'] = ds_temp
    dict['bathy'] = Zitp
    dict['lon_reg'] = lons
    dict['lat_reg'] = lats
    
    return dict


def get_bottomS_climato(
    INFILES,
    lonLims,
    latLims,
    bath_file,
    time_adjust=True,
    year_lims=[1991, 2020],
    season=[],
    h5_outputfile=[]):
    '''
    Generates and returns the bottom salinity climatology.
    This script uses GEBCO data, path needs to be provided.
    If file already exists, it will be read and the output will come from that file.

    INFILES: Path location of bottom salinity .nc files. 
    lonLims (list): Longitude restraints [longitude min, longitude max]
    latLims (list): Latitude restraints [latitude min, latitude max]
    years_lim (list): Climatology years, default is [1991,2020]
    season: The season of interest; spring, summer, fall
    h5_outputfile: Location of output file, type .10.h5
    bath_file: Path location of GEBCO bathymetry data, including file name.

    Examples.
    #fish_hab region: lonLims = [-60, -43], latLims = [39, 56], h5_outputfile='Tbot_climato_'+season+'_0.10.h5'
    #FC AZMP report region: lonLims = [-60, -45], latLims = [42, 56], h5_outputfile='Tbot_climato_'+season+'_0.10.h5'
    #Include 2H in above: lonLims = [-63, -45], latLims = [42, 58], h5_outputfile='Tbot_climato_'+season+'_0.10.h5'
    #NSRF area: lonLims = [-70, -56], latLims = [57, 67], year_lims=[2006, 2021], season='summer', h5_outputfile='Tbot_climato_NSRFx_'+season+'_2006-2021.h5'
    # NSRFxx, down to 53N (for Normandeau's work): lonLims = [-70, -56], latLims = [53, 67]
    '''

    # Because the code is trying to divide by zero
    # (https://stackoverflow.com/questions/14861891/runtimewarning-invalid-value-encountered-in-divide)
    np.seterr(divide='ignore', invalid='ignore')
    
    ## ---- Check if H5 file exists ---- ##        
    if os.path.isfile(h5_outputfile):
        print(h5_outputfile + ' exist! Reading directly')
        h5f = h5py.File(h5_outputfile,'r')
        ds_saln = h5f['Sbot'][:]
        lons = h5f['lon_reg'][:]
        lats = h5f['lat_reg'][:]
        Zitp = h5f['Zitp'][:]
        h5f.close()

    else:

        #Load in the bottom salinity data
        ds = xr.open_dataset(os.path.expanduser(INFILES))
        #Isolate for the time of interest
        ds = ds.sel(TIME=((ds['TIME.year']>=year_lims[0])*(ds['TIME.year']<=year_lims[1])))
        ds = ds.mean('TIME')

        #Isolate for the region of interest
        ds = ds.sel(X=((ds.LONGITUDE[0,:]>=lonLims[0])*(ds.LONGITUDE[0,:]<=lonLims[1])).values)
        ds = ds.sel(Y=((ds.LATITUDE[:,0]>=latLims[0])*(ds.LATITUDE[:,0]<=latLims[1])).values)

        #Record the temperature, latitude, longitude
        if time_adjust:
            ds_saln = ds.BOTTOM_SALINITY_ADJUSTED.values
        else:
            ds_saln = ds.BOTTOM_SALINITY.values
        lons = ds.LONGITUDE.values
        lats = ds.LATITUDE.values

        ## ---- Region parameters ---- ##
        ds_bath = xr.open_dataset(bath_file)
        ds_bath = ds_bath.isel(lon=(ds_bath.lon>=lonLims[0])*(ds_bath.lon<=lonLims[1]))
        ds_bath = ds_bath.isel(lat=(ds_bath.lat>=latLims[0])*(ds_bath.lat<=latLims[1]))
        Zitp = ds_bath.elevation[::30,::30].values

        #Save the climatology
        h5f = h5py.File(h5_outputfile, 'w')
        h5f.create_dataset('Sbot', data=ds_saln)
        h5f.create_dataset('lon_reg', data=lons)
        h5f.create_dataset('lat_reg', data=lats)
        h5f.create_dataset('Zitp', data=Zitp)
        h5f.close()

    #Fill dict for output
    dict = {}
    dict['Sbot'] = ds_saln
    dict['bathy'] = Zitp
    dict['lon_reg'] = lons
    dict['lat_reg'] = lats
    
    return dict



def get_bottomT(year_file, year, season, climato_file, time_adjust=True, nafo_mask=True, lab_mask=True):
    '''
    Generate and return bottom temperature data corresponding to a certain climatology map.
    Returns a dictionary.
    Returns a file of bottom temperature for the year (.h5), before masking.

    year_file: Path location of bottom temperature year of interest, file name included.
    season: Season of interest.
    climato_file: Location of climatology file.
    nafo_mask (bool): Whether or not regions not of interest are masked.
    lab_mask (bool): Whether or not the Labrador 100m contour is masked.
    '''

    ## ---- Load Climato data ---- ##    
    print('Load ' + climato_file)
    h5f = h5py.File(climato_file, 'r')
    Tbot_climato = h5f['Tbot'][:]
    lon_reg = h5f['lon_reg'][:][0,:]
    lat_reg = h5f['lat_reg'][:][:,0]
    Zitp = h5f['Zitp'][:]
    h5f.close()

    ## ---- Derive some parameters ---- ##    
    lonLims = [lon_reg[0], lon_reg[-1]]
    latLims = [lat_reg[0], lat_reg[-1]]

    ## ---- NAFO divisions ---- ##
    nafo_div = get_nafo_divisions()

    ## ---- Get CABOTS data --- ##
    print('Get ' + year_file)
    ds = xr.open_dataset(year_file)
    #Isolate for the time of interest
    if np.size(year) == 1:
        ds = ds.sel(TIME=ds['TIME.year']==int(year))
        ds = ds.mean('TIME')
    else:
        ds = ds.sel(TIME=np.isin(ds['TIME.year'],year))
    # Selection of a subset region
    ds = ds.sel(X=((ds.LONGITUDE[0,:]>=lonLims[0])*(ds.LONGITUDE[0,:]<=lonLims[1])).values)
    ds = ds.sel(Y=((ds.LATITUDE[:,0]>=latLims[0])*(ds.LATITUDE[:,0]<=latLims[1])).values)
    if time_adjust:
        Tbot = ds.BOTTOM_TEMPERATURE_ADJUSTED.values
    else:
        Tbot = ds.BOTTOM_TEMPERATURE.values

    ## Save data in h5 for further use
    if np.size(year) == 1:
        h5_cube_name = 'operation_files/Tcube_' + season + str(year) + '.h5'
        h5f = h5py.File(h5_cube_name, 'w')
        h5f.create_dataset('temperature', data=Tbot)
        h5f.create_dataset('lon_reg', data=lon_reg)
        h5f.create_dataset('lat_reg', data=lat_reg)
        h5f.create_dataset('Zitp', data=Zitp)
        h5f.close()
    else:
        for i,value in enumerate(year):
            h5_cube_name = 'operation_files/Tcube_' + season + str(value) + '.h5'
            h5f = h5py.File(h5_cube_name, 'w')
            h5f.create_dataset('temperature', data=Tbot[i])
            h5f.create_dataset('lon_reg', data=lon_reg)
            h5f.create_dataset('lat_reg', data=lat_reg)
            h5f.create_dataset('Zitp', data=Zitp)
            h5f.close()

    # Mask data on coastal Labrador
    if lab_mask == True:
        print('Mask coastal labrador')
        contour_mask = np.load('operation_files/100m_contour_labrador.npy')
        polygon_mask = Polygon(contour_mask)
        for i, xx in enumerate(lon_reg):
            for j,yy in enumerate(lat_reg):
                point = Point(lon_reg[i], lat_reg[j])
                if polygon_mask.contains(point): # mask data near Labrador in fall
                    Tbot[j,i] = np.nan 

    if nafo_mask == True:
        # Mask data outside Nafo div.
        print('Mask according to NAFO division for ' + season)
        # Polygons
        polygon3K = Polygon(zip(nafo_div['3K']['lon'], nafo_div['3K']['lat']))
        polygon3L = Polygon(zip(nafo_div['3L']['lon'], nafo_div['3L']['lat']))
        polygon3N = Polygon(zip(nafo_div['3N']['lon'], nafo_div['3N']['lat']))
        polygon3O = Polygon(zip(nafo_div['3O']['lon'], nafo_div['3O']['lat']))
        polygon3Ps = Polygon(zip(nafo_div['3Ps']['lon'], nafo_div['3Ps']['lat']))
        polygon2J = Polygon(zip(nafo_div['2J']['lon'], nafo_div['2J']['lat']))

        if season == 'spring':
            for i, xx in enumerate(lon_reg):
                for j,yy in enumerate(lat_reg):
                    point = Point(lon_reg[i], lat_reg[j])
                    if polygon3L.contains(point) | polygon3N.contains(point) | polygon3O.contains(point) | polygon3Ps.contains(point):
                        pass #nothing to do
                    else:
                        Tbot[:,j,i] = np.nan

        elif season == 'fall':
            for i, xx in enumerate(lon_reg):
                for j,yy in enumerate(lat_reg):
                    point = Point(lon_reg[i], lat_reg[j])
                    if polygon2J.contains(point) | polygon3K.contains(point) | polygon3L.contains(point) | polygon3N.contains(point) | polygon3O.contains(point) | polygon3Ps.contains(point):
                        pass #nothing to do
                    else:
                        pass # No mask in the fall

        else:
            print('no division mask, all data taken')
    else:
            print('no division mask, all data taken')

    print(' -> Done!')    

    # Fill dict for output
    dict = {}
    if np.size(year) == 1:
        dict['Tbot'] = Tbot
        dict['bathy'] = Zitp
        dict['lon_reg'] = lon_reg
        dict['lat_reg'] = lat_reg
    else:
        for i,value in enumerate(year):
            dict[str(value)] = {}
            dict[str(value)]['Tbot'] = Tbot[i]
            dict[str(value)]['bathy'] = Zitp
            dict[str(value)]['lon_reg'] = lon_reg
            dict[str(value)]['lat_reg'] = lat_reg


    return dict

def get_bottomS(year_file, year, season, climato_file, time_adjust=True, nafo_mask=True, lab_mask=True):
    '''
    Generate and return bottom salinity data corresponding to a certain climatology map.
    Returns a dictionary.
    Returns a file of bottom salinity for the year (.h5), before masking.

    year_file: Path location of bottom salinity year of interest, file name included.
    season: Season of interest.
    climato_file: Location of climatology file.
    nafo_mask (bool): Whether or not regions not of interest are masked.
    lab_mask (bool): Whether or not the Labrador 100m contour is masked.
    '''



    """ Generate and returns bottom temperature data corresponding to a certain climatology map
    (previously generated with get_bottomS_climato)
    Function returns:
    - Sbot:  gridded bottom salinity
    - lons, lats: coordinates of good casts used to generate the grid
    *Note: they are not regular coordinates of the grid that can be obtained with get_bottomS_climato
        
    Usage ex (suppose climato file already exist):
    import azmp_utils as azu
    climato_file = 'Sbot_climato_spring_0.25.h5'
    year_file = '/home/cyrf0006/data/dev_database/2017.nc'
    h5f = h5py.File(climato_file, 'r')
    Sbot_climato = h5f['Sbot'][:]
    lon_reg = h5f['lon_reg'][:]
    lat_reg = h5f['lat_reg'][:]
    Zitp = h5f['Zitp'][:]
    h5f.close()
    Sbot_dict = azu.get_bottomS(year_file, 'fall', climato_file)

    Dec. 2018: Added a flag for masking or not.

    """

    #Determine if it's NSRF data
    if 'NSRF' in climato_file:
        NSRF_plot = True
    else:
        NSRF_plot = False

    ## ---- Load Climato data ---- ##    
    print('Load ' + climato_file)
    h5f = h5py.File(climato_file, 'r')
    Sbot_climato = h5f['Sbot'][:]
    lon_reg = h5f['lon_reg'][:][0,:]
    lat_reg = h5f['lat_reg'][:][:,0]
    Zitp = h5f['Zitp'][:]
    h5f.close()

    ## ---- Derive some parameters ---- ##    
    lonLims = [lon_reg[0], lon_reg[-1]]
    latLims = [lat_reg[0], lat_reg[-1]]
    
    ## ---- NAFO divisions ---- ##
    nafo_div = get_nafo_divisions()

    ## Get SFAs data
    myshp = open(os.path.expanduser('~/github/AZMP-NL/utils/SFAs/SFAs_PANOMICS_Fall2020_shp/SFAs_PANOMICS_Fall2020.shp'), 'rb')
    mydbf = open(os.path.expanduser('~/github/AZMP-NL/utils/SFAs/SFAs_PANOMICS_Fall2020_shp/SFAs_PANOMICS_Fall2020.dbf'), 'rb')
    r = shapefile.Reader(shp=myshp, dbf=mydbf, encoding = "ISO8859-1")
    records = r.records()
    shapes = r.shapes()

    # Fill dictionary with shapes
    shrimp_area = {}
    for idx, rec in enumerate(records):
        if rec[1] == 'Eastern Assessment Zone':
            shrimp_area['2'] = np.array(shapes[idx].points)
        elif rec[1] == 'Western Assessment Zone':
            shrimp_area['3'] = np.array(shapes[idx].points)
        else:
            shrimp_area[rec[0]] = np.array(shapes[idx].points)

    ## ---- Get CABOTS data --- ##
    print('Get ' + year_file)
    ds = xr.open_dataset(year_file)
    #Isolate for the time of interest
    if np.size(year) == 1:
        ds = ds.sel(TIME=ds['TIME.year']==int(year))
        ds = ds.mean('TIME')
    else:
        ds = ds.sel(TIME=np.isin(ds['TIME.year'],year))
    # Selection of a subset region
    ds = ds.sel(X=((ds.LONGITUDE[0,:]>=lonLims[0])*(ds.LONGITUDE[0,:]<=lonLims[1])).values)
    ds = ds.sel(Y=((ds.LATITUDE[:,0]>=latLims[0])*(ds.LATITUDE[:,0]<=latLims[1])).values)
    if time_adjust:
        Sbot = ds.BOTTOM_SALINITY_ADJUSTED.values
    else:
        Sbot = ds.BOTTOM_SALINITY.values

    ## Save data in h5 for further use
    if np.size(year) == 1:
        h5_cube_name = 'operation_files/Scube_' + season + str(year) + '.h5'
        h5f = h5py.File(h5_cube_name, 'w')
        h5f.create_dataset('salinity', data=Sbot)
        h5f.create_dataset('lon_reg', data=lon_reg)
        h5f.create_dataset('lat_reg', data=lat_reg)
        h5f.create_dataset('Zitp', data=Zitp)
        h5f.close()
    else:
        for i,value in enumerate(year):
            h5_cube_name = 'operation_files/Scube_' + season + str(value) + '.h5'
            h5f = h5py.File(h5_cube_name, 'w')
            h5f.create_dataset('salinity', data=Sbot[i])
            h5f.create_dataset('lon_reg', data=lon_reg)
            h5f.create_dataset('lat_reg', data=lat_reg)
            h5f.create_dataset('Zitp', data=Zitp)
            h5f.close()

    # Mask data on coastal Labrador
    if lab_mask == True:
        print('Mask coastal labrador')
        contour_mask = np.load('operation_files/100m_contour_labrador.npy')
        polygon_mask = Polygon(contour_mask)
        for i, xx in enumerate(lon_reg):
            for j,yy in enumerate(lat_reg):
                point = Point(lon_reg[i], lat_reg[j])
                if polygon_mask.contains(point): # mask data near Labrador in fall
                    Sbot[j,i] = np.nan 

    # Mask data outside Nafo div.
    if nafo_mask == True:
        print('Mask according to NAFO division for ' + season)
        # Polygons
        polygon3K = Polygon(zip(nafo_div['3K']['lon'], nafo_div['3K']['lat']))
        polygon3L = Polygon(zip(nafo_div['3L']['lon'], nafo_div['3L']['lat']))
        polygon3N = Polygon(zip(nafo_div['3N']['lon'], nafo_div['3N']['lat']))
        polygon3O = Polygon(zip(nafo_div['3O']['lon'], nafo_div['3O']['lat']))
        polygon3Ps = Polygon(zip(nafo_div['3Ps']['lon'], nafo_div['3Ps']['lat']))
        polygon2J = Polygon(zip(nafo_div['2J']['lon'], nafo_div['2J']['lat']))
        sfa2 = Polygon(shrimp_area['2'])
        sfa3 = Polygon(shrimp_area['3'])
        sfa4 = Polygon(shrimp_area['4'])
        sfa5 = Polygon(shrimp_area['5'])
        sfa6 = Polygon(shrimp_area['6'])
        sfa7 = Polygon(shrimp_area['7'])

        if NSRF_plot:
            if season == 'summer':
                for i, xx in enumerate(lon_reg):
                    for j,yy in enumerate(lat_reg):
                        point = Point(lon_reg[i], lat_reg[j])
                        if sfa2.contains(point) | sfa3.contains(point) | sfa4.contains(point):
                            pass #nothing to do but cannot implement negative statement "if not" above
                        else:
                            Sbot[:,j,i] = np.nan
            if season == 'fall':
                for i, xx in enumerate(lon_reg):
                    for j,yy in enumerate(lat_reg):
                        point = Point(lon_reg[i], lat_reg[j])
                        if sfa4.contains(point) | sfa5.contains(point) | sfa6.contains(point) | sfa7.contains(point):
                            pass #nothing to do but cannot implement negative statement "if not" above
                        else:
                            Sbot[:,j,i] = np.nan
        else:
            if season == 'spring':
                for i, xx in enumerate(lon_reg):
                    for j,yy in enumerate(lat_reg):
                        point = Point(lon_reg[i], lat_reg[j])
                        if polygon3L.contains(point) | polygon3N.contains(point) | polygon3O.contains(point) | polygon3Ps.contains(point):
                            pass #nothing to do but cannot implement negative statement "if not" above
                        else:
                            Sbot[:,j,i] = np.nan
            elif season == 'fall':
                for i, xx in enumerate(lon_reg):
                    for j,yy in enumerate(lat_reg):
                        point = Point(lon_reg[i], lat_reg[j])
                        if polygon2J.contains(point) | polygon3K.contains(point) | polygon3L.contains(point) | polygon3N.contains(point) | polygon3O.contains(point) | polygon3Ps.contains(point):
                            pass #nothing to do but cannot implement negative statement "if not" above
                        else:
                            pass
            elif season == 'summer':
                for i, xx in enumerate(lon_reg):
                    for j,yy in enumerate(lat_reg):
                        point = Point(lon_reg[i], lat_reg[j])

            else:
                print('no division mask, all data taken')

        print(' -> Done!')

    # Fill dict for output
    dict = {}
    if np.size(year) == 1:
        dict['Sbot'] = Sbot
        dict['bathy'] = Zitp
        dict['lon_reg'] = lon_reg
        dict['lat_reg'] = lat_reg
    else:
        for i,value in enumerate(year):
            dict[str(value)] = {}
            dict[str(value)]['Sbot'] = Sbot[i]
            dict[str(value)]['bathy'] = Zitp
            dict[str(value)]['lon_reg'] = lon_reg
            dict[str(value)]['lat_reg'] = lat_reg


    return dict

def get_surfT_climato(INFILES, LON_REG,  LAT_REG, year_lims=[1981, 2010], season=[], zlims=[5, 20], dz=5, h5_outputfile=[]):
    """ Generate and returns the climatological surface temperature map from ship observations.

    If the pickled filename exists, the function will by-pass the processing and return only saved climatology.
    
    Usage ex:
    dc = .1
    import numpy as np
    import azmp_utils as azu
    lonLims = [-60, -43] # fish_hab region
    latLims = [39, 56]
    lonLims = [-60, -45] # FC AZMP report region
    latLims = [42, 56]
    lonLims = [-63, -45] # include 2H in above
    latLims = [42, 58]
    lon_reg = np.arange(lonLims[0]+dc/2, lonLims[1]-dc/2, dc)
    lat_reg = np.arange(latLims[0]+dc/2, latLims[1]-dc/2, dc)
    azu.get_surfT_climato('/home/cyrf0006/data/dev_database/netCDF/*.nc', lon_reg, lat_reg, season='fall', h5_outputfile='Tsurf_climato_fall_0.10.h5') 

    OR

    azu.get_surfT_climato('/home/cyrf0006/data/dev_database/netCDF/*.nc', lon_reg, lat_reg, year_lims=[1980, 2019], season='fall', h5_outputfile='Tsurf_climato_fall_0.10_1980_2019.h5')
    
    """

    # Because the code is trying to divide by zero
    # (https://stackoverflow.com/questions/14861891/runtimewarning-invalid-value-encountered-in-divide)
    np.seterr(divide='ignore', invalid='ignore')

    ## ---- Check if H5 file exists ---- ##        
    if os.path.isfile(h5_outputfile):
        print(h5_outputfile + ' exist! Reading directly')
        h5f = h5py.File(h5_outputfile,'r')
        Tsurf = h5f['Tsurf'][:]
        lon_reg = h5f['lon_reg'][:]
        lat_reg = h5f['lat_reg'][:]
        lon_orig = h5f['lon_orig'][:]
        lat_orig = h5f['lat_orig'][:]
        Zitp = h5f['Zitp'][:]
        h5f.close()

    else:

        ## ---- Region parameters ---- ##
        #    add_path('/home/cyrf0006/data/GEBCO/')
        dataFile = '/home/cyrf0006/data/GEBCO/GEBCO_2014_1D.nc' # Maybe find a better way to handle this file
        lonLims = [LON_REG[0], LON_REG[-1]]
        latLims = [LAT_REG[0], LAT_REG[-1]]
        zmin = zlims[0] # do try to compute bottom temp above that depth
        zmax = zlims[1] # do try to compute bottom temp below that depth
        lon_reg = LON_REG
        lat_reg = LAT_REG
        lon_grid, lat_grid = np.meshgrid(lon_reg,lat_reg)
        dc = np.round(np.diff(lon_reg[0:2]), 3)[0]

        ## ---- Bathymetry ---- ####
        print('Load and grid bathymetry')
        # Load data
        dataset = netCDF4.Dataset(dataFile)
        # Extract variables
        #x = dataset.variables['x_range'] 
        #y = dataset.variables['y_range']
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
        del Z, lon, lat, zz, dataset
        print(' -> Done!')

        ## ---- Get CTD data --- ##
        print('Get historical data')
        ds = xr.open_mfdataset(INFILES)
        # Selection of a subset region
        ds = ds.where((ds.longitude>lonLims[0]) & (ds.longitude<lonLims[1]), drop=True)
        ds = ds.where((ds.latitude>latLims[0]) & (ds.latitude<latLims[1]), drop=True)
        # Select time (save several options here)
        if season == 'summer':
            #ds = ds.sel(time=ds['time.season']=='JJA')
            ds = ds.sel(time=((ds['time.month']>=7)) & ((ds['time.month']<=9)))
        elif season == 'spring':
            #ds = ds.sel(time=ds['time.season']=='MAM')
            ds = ds.sel(time=((ds['time.month']>=4)) & ((ds['time.month']<=6)))
        elif season == 'fall':
            #ds = ds.sel(time=ds['time.season']=='SON')
            ds = ds.sel(time=((ds['time.month']>=9)) & ((ds['time.month']<=12)))
        else:
            print('!! no season specified, used them all! !!')

        # Remome problematic datasets
        print('!!Remove MEDBA data!!')
        print('  ---> I Should be improme because I remove good data!!!!')
        ds = ds.where(ds.instrument_ID!='MEDBA', drop=True)
        ds = ds.where(ds.instrument_ID!='MEDTE', drop=True)
 
            
        # Time period for climatology
        ds = ds.sel(time=ds['time.year']>=year_lims[0])
        ds = ds.sel(time=ds['time.year']<=year_lims[1])
        # Vertical limitation (for surface; or CIL, for example)
        ds = ds.sel(level=ds['level']<=zmax*2) #double the search for zmax (restrict later)
        # Vertical binning (on dataArray; more appropriate here
        da_temp = ds['temperature']
        lons = np.array(ds.longitude)
        lats = np.array(ds.latitude)
        bins = np.arange(dz/2.0, ds.level.max(), dz)
        da_temp = da_temp.groupby_bins('level', bins).mean(dim='level')
        #To Pandas Dataframe
        df_temp = da_temp.to_pandas()
        df_temp.columns = bins[0:-1] #rename columns with 'bins'
        idx_empty_rows = df_temp.isnull().all(1).values.nonzero()[0]
        df_temp = df_temp.dropna(axis=0,how='all')
        lons = np.delete(lons,idx_empty_rows)
        lats = np.delete(lats,idx_empty_rows)
        del ds, da_temp
        print(' -> Done!')        

        ## --- fill 3D cube --- ##  
        print('Fill regular cube')
        z = df_temp.columns.values
        V = np.full((lat_reg.size, lon_reg.size, z.size), np.nan)

        # Aggregate on regular grid
        for i, xx in enumerate(lon_reg):
            for j, yy in enumerate(lat_reg):    
                idx = np.where((lons>=xx-dc/2) & (lons<xx+dc/2) & (lats>=yy-dc/2) & (lats<yy+dc/2))
                tmp = np.array(df_temp.iloc[idx].mean(axis=0))
                idx_good = np.argwhere((~np.isnan(tmp)) & (tmp<30))
                if np.size(idx_good)==1:
                    V[j,i,:] = np.array(df_temp.iloc[idx].mean(axis=0))
                elif np.size(idx_good)>1: # vertical interpolation between pts
                    interp = interp1d(np.squeeze(z[idx_good]), np.squeeze(tmp[idx_good]))
                    idx_interp = np.arange(np.int(idx_good[0]),np.int(idx_good[-1]+1))
                    V[j,i,idx_interp] = interp(z[idx_interp]) # interpolate only where possible (1st to last good idx)


        # horiozntal interpolation at each depth
        lon_grid, lat_grid = np.meshgrid(lon_reg,lat_reg)
        lon_vec = np.reshape(lon_grid, lon_grid.size)
        lat_vec = np.reshape(lat_grid, lat_grid.size)
        for k, zz in enumerate(z):
            # Meshgrid 1D data (after removing NaNs)
            tmp_grid = V[:,:,k]
            tmp_vec = np.reshape(tmp_grid, tmp_grid.size)
            #print 'interpolate depth layer ' + np.str(k) + ' / ' + np.str(z.size) 
            # griddata (after removing nans)
            idx_good = np.argwhere(~np.isnan(tmp_vec))
            if idx_good.size>10: # will ignore depth where no data exist
                LN = np.squeeze(lon_vec[idx_good])
                LT = np.squeeze(lat_vec[idx_good])
                TT = np.squeeze(tmp_vec[idx_good])
                zi = griddata((LN, LT), TT, (lon_grid, lat_grid), method='linear')
                V[:,:,k] = zi
            else:
                continue
        print(' -> Done!')    

        # mask using bathymetry
        if zmin>10: # this is CIL
             for i, xx in enumerate(lon_reg):
                for j,yy in enumerate(lat_reg):
                    if Zitp[j,i] > -np.array([zmin, zmax]).mean(): # remove shallower than depth range ave
                        V[j,i,:] = np.nan           
        else: # this is SST            
            for i, xx in enumerate(lon_reg):
                for j,yy in enumerate(lat_reg):
                    if Zitp[j,i] > -10: # remove shallower than 10m
                        V[j,i,:] = np.nan

        # average surface obs
        print('Restricting depth range: ' + str(zmin) + '-' + str(zmax)+'m')
        idx = np.where([(z>=zmin) & (z<=zmax)])[1]
        Tsurf = np.nanmean(V[:,:,idx], 2)

        # Save data for further use
        if np.size(h5_outputfile):
            h5f = h5py.File(h5_outputfile, 'w')
            h5f.create_dataset('Tsurf', data=Tsurf)
            h5f.create_dataset('lon_reg', data=lon_reg)
            h5f.create_dataset('lat_reg', data=lat_reg)
            h5f.create_dataset('lon_orig', data=lons)
            h5f.create_dataset('lat_orig', data=lats)
            h5f.create_dataset('Zitp', data=Zitp)
            h5f.create_dataset('z', data=z)
            h5f.close()

    # Fill dict for output
    dict = {}
    dict['Tsurf'] = Tsurf
    dict['bathy'] = Zitp
    dict['lon_reg'] = lon_reg
    dict['lat_reg'] = lat_reg
    dict['lon_orig'] = lons
    dict['lat_orig'] = lats
    
    return dict

def get_cilT_climato(INFILES, LON_REG,  LAT_REG, year_lims=[1981, 2010], season=[], zlims=[5, 500], dz=5, h5_outputfile=[]):
    """ Generate and returns the climatological surface temperature map from ship observations.

    If the pickled filename exists, the function will by-pass the processing and return only saved climatology.
    
    Usage ex:
    dc = .1
    import numpy as np
    import azmp_utils as azu
    lonLims = [-60, -43] # fish_hab region
    latLims = [39, 56]
    lonLims = [-60, -45] # FC AZMP report region
    latLims = [42, 56]
    lonLims = [-63, -45] # include 2H in above
    latLims = [42, 58]
    lon_reg = np.arange(lonLims[0]+dc/2, lonLims[1]-dc/2, dc)
    lat_reg = np.arange(latLims[0]+dc/2, latLims[1]-dc/2, dc)
    azu.get_surfT_climato('/home/cyrf0006/data/dev_database/netCDF/*.nc', lon_reg, lat_reg, season='fall', h5_outputfile='Tsurf_climato_fall_0.10.h5') 

    OR

    azu.get_cil_climato('/home/cyrf0006/data/dev_database/netCDF/*.nc', lon_reg, lat_reg, year_lims=[1980, 2019], season='fall', h5_outputfile='Tsurf_climato_fall_0.10_1980_2019.h5')
    
    """

    # Because the code is trying to divide by zero
    # (https://stackoverflow.com/questions/14861891/runtimewarning-invalid-value-encountered-in-divide)
    np.seterr(divide='ignore', invalid='ignore')

    ## ---- Check if H5 file exists ---- ##        
    if os.path.isfile(h5_outputfile):
        print(h5_outputfile + ' exist! Reading directly')
        h5f = h5py.File(h5_outputfile,'r')
        Tsurf = h5f['Tsurf'][:]
        lon_reg = h5f['lon_reg'][:]
        lat_reg = h5f['lat_reg'][:]
        lon_orig = h5f['lon_orig'][:]
        lat_orig = h5f['lat_orig'][:]
        Zitp = h5f['Zitp'][:]
        h5f.close()

    else:

        ## ---- Region parameters ---- ##
        #    add_path('/home/cyrf0006/data/GEBCO/')
        dataFile = '/home/cyrf0006/data/GEBCO/GEBCO_2014_1D.nc' # Maybe find a better way to handle this file
        lonLims = [LON_REG[0], LON_REG[-1]]
        latLims = [LAT_REG[0], LAT_REG[-1]]
        zmin = zlims[0] # do try to compute bottom temp above that depth
        zmax = zlims[1] # do try to compute bottom temp below that depth
        lon_reg = LON_REG
        lat_reg = LAT_REG
        lon_grid, lat_grid = np.meshgrid(lon_reg,lat_reg)
        dc = np.round(np.diff(lon_reg[0:2]), 3)[0]

        ## ---- Bathymetry ---- ####
        print('Load and grid bathymetry')
        # Load data
        dataset = netCDF4.Dataset(dataFile)
        # Extract variables
        #x = dataset.variables['x_range'] 
        #y = dataset.variables['y_range']
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
        del Z, lon, lat, zz, dataset
        print(' -> Done!')

        ## ---- Get CTD data --- ##
        print('Get historical data')
        ds = xr.open_mfdataset(INFILES)
        # Selection of a subset region
        ds = ds.where((ds.longitude>lonLims[0]) & (ds.longitude<lonLims[1]), drop=True)
        ds = ds.where((ds.latitude>latLims[0]) & (ds.latitude<latLims[1]), drop=True)
        # Select time (save several options here)
        if season == 'summer':
            #ds = ds.sel(time=ds['time.season']=='JJA')
            ds = ds.sel(time=((ds['time.month']>=7)) & ((ds['time.month']<=9)))
        elif season == 'spring':
            #ds = ds.sel(time=ds['time.season']=='MAM')
            ds = ds.sel(time=((ds['time.month']>=4)) & ((ds['time.month']<=6)))
        elif season == 'fall':
            #ds = ds.sel(time=ds['time.season']=='SON')
            ds = ds.sel(time=((ds['time.month']>=9)) & ((ds['time.month']<=12)))
        else:
            print('!! no season specified, used them all! !!')

        # Remome problematic datasets
        print('!!Remove MEDBA data!!')
        print('  ---> I Should be improme because I remove good data!!!!')
        ds = ds.where(ds.instrument_ID!='MEDBA', drop=True)
        ds = ds.where(ds.instrument_ID!='MEDTE', drop=True)
 
            
        # Time period for climatology
        ds = ds.sel(time=ds['time.year']>=year_lims[0])
        ds = ds.sel(time=ds['time.year']<=year_lims[1])
        # Vertical limitation (for surface; or CIL, for example)
        ds = ds.sel(level=ds['level']<=zmax) #double the search for zmax (restrict later)
        # Vertical binning (on dataArray; more appropriate here
        da_temp = ds['temperature']
        lons = np.array(ds.longitude)
        lats = np.array(ds.latitude)
        bins = np.arange(dz/2.0, ds.level.max(), dz)
        da_temp = da_temp.groupby_bins('level', bins).mean(dim='level')
        #To Pandas Dataframe
        df_temp = da_temp.to_pandas()
        df_temp.columns = bins[0:-1] #rename columns with 'bins'
        idx_empty_rows = df_temp.isnull().all(1).values.nonzero()[0]
        df_temp = df_temp.dropna(axis=0,how='all')
        lons = np.delete(lons,idx_empty_rows)
        lats = np.delete(lats,idx_empty_rows)
        del ds, da_temp
        print(' -> Done!')        

        ## --- fill 3D cube --- ##  
        print('Fill regular cube')
        z = df_temp.columns.values
        V = np.full((lat_reg.size, lon_reg.size, z.size), np.nan)

        # Aggregate on regular grid
        for i, xx in enumerate(lon_reg):
            for j, yy in enumerate(lat_reg):    
                idx = np.where((lons>=xx-dc/2) & (lons<xx+dc/2) & (lats>=yy-dc/2) & (lats<yy+dc/2))
                tmp = np.array(df_temp.iloc[idx].mean(axis=0))
                idx_good = np.argwhere((~np.isnan(tmp)) & (tmp<30))
                if np.size(idx_good)==1:
                    V[j,i,:] = np.array(df_temp.iloc[idx].mean(axis=0))
                elif np.size(idx_good)>1: # vertical interpolation between pts
                    interp = interp1d(np.squeeze(z[idx_good]), np.squeeze(tmp[idx_good]))
                    idx_interp = np.arange(np.int(idx_good[0]),np.int(idx_good[-1]+1))
                    V[j,i,idx_interp] = interp(z[idx_interp]) # interpolate only where possible (1st to last good idx)


        # horiozntal interpolation at each depth
        lon_grid, lat_grid = np.meshgrid(lon_reg,lat_reg)
        lon_vec = np.reshape(lon_grid, lon_grid.size)
        lat_vec = np.reshape(lat_grid, lat_grid.size)
        for k, zz in enumerate(z):
            # Meshgrid 1D data (after removing NaNs)
            tmp_grid = V[:,:,k]
            tmp_vec = np.reshape(tmp_grid, tmp_grid.size)
            #print 'interpolate depth layer ' + np.str(k) + ' / ' + np.str(z.size) 
            # griddata (after removing nans)
            idx_good = np.argwhere(~np.isnan(tmp_vec))
            if idx_good.size>10: # will ignore depth where no data exist
                LN = np.squeeze(lon_vec[idx_good])
                LT = np.squeeze(lat_vec[idx_good])
                TT = np.squeeze(tmp_vec[idx_good])
                zi = griddata((LN, LT), TT, (lon_grid, lat_grid), method='linear')
                V[:,:,k] = zi
            else:
                continue
        print(' -> Done!')    

        # mask using bathymetry
        for i, xx in enumerate(lon_reg):
            for j,yy in enumerate(lat_reg):
                
                if Zitp[j,i] > -10: # remove shallower than 10m
                    V[j,i,:] = np.nan
                    
        # Find CIL and average
        for i, xx in enumerate(lon_reg):
            for j,yy in enumerate(lat_reg):
                for k, zz in enumerate(z):
                    if V[j,i,k] > 0: # remove all temp >0
                        V[j,i,k] = np.nan
                                        
        Tcil = np.nanmean(V, 2)

        # Save data for further use
        if np.size(h5_outputfile):
            h5f = h5py.File(h5_outputfile, 'w')
            h5f.create_dataset('Tcil', data=Tcil)
            h5f.create_dataset('lon_reg', data=lon_reg)
            h5f.create_dataset('lat_reg', data=lat_reg)
            h5f.create_dataset('lon_orig', data=lons)
            h5f.create_dataset('lat_orig', data=lats)
            h5f.create_dataset('Zitp', data=Zitp)
            h5f.create_dataset('z', data=z)
            h5f.close()

    # Fill dict for output
    dict = {}
    dict['Tcil'] = Tcil
    dict['bathy'] = Zitp
    dict['lon_reg'] = lon_reg
    dict['lat_reg'] = lat_reg
    dict['lon_orig'] = lons
    dict['lat_orig'] = lats
    
    return dict


def bottomT_quickplot(h5_outputfile, figure_file=[]):
    """ Using h5 file created by get_bottomT_climato, this function plots the bottom temperature using basemap.
    
    Usage ex:
    import azmp_utils as azu
    azu.bottomT_quickplot('Tbot_climato_fall.h5')
    
    """

    ## ---- Load data ---- ##    
    h5f = h5py.File(h5_outputfile,'r')
    try:
        Tbot = h5f['Tbot'][:]
    except:
        Tbot = h5f['Tsurf'][:]

    lon_reg = h5f['lon_reg'][:]
    lat_reg = h5f['lat_reg'][:]
    lons = h5f['lon_orig'][:]
    lats = h5f['lat_orig'][:]
    #Zitp = h5f['Zitp'][:]
    h5f.close()

    ## ---- Derive some parameters ---- ##    
    lon_0 = np.round(np.mean(lon_reg))
    lat_0 = np.round(np.mean(lat_reg))
    lonLims = [lon_reg[0], lon_reg[-1]]
    latLims = [lat_reg[0], lat_reg[-1]]
        
    ## ---- Plot map ---- ##
    fig, ax = plt.subplots(nrows=1, ncols=1)
    m = Basemap(ax=ax, projection='merc',lon_0=lon_0,lat_0=lat_0, llcrnrlon=lonLims[0],llcrnrlat=latLims[0],urcrnrlon=lonLims[1],urcrnrlat=latLims[1], resolution='h')
    levels = np.linspace(0, 15, 16)
    xi, yi = m(*np.meshgrid(lon_reg, lat_reg))
    c = m.contourf(xi, yi, Tbot, levels, cmap=plt.cm.RdBu_r, extend='both')
    #c = m.pcolor(xi, yi, Tbot, levels, cmap=plt.cm.RdBu_r, extend='both')
    #    x,y = m(*np.meshgrid(lon,lat))
    #    cc = m.contour(x, y, -Z, [100, 500, 1000, 4000], colors='grey');
    lon_casts, lat_casts = m(lons, lats)
    m.plot(lon_casts, lat_casts, '.', color='k')
    m.fillcontinents(color='tan');
    
    m.drawparallels([40, 45, 50, 55, 60], labels=[1,0,0,0], fontsize=12, fontweight='normal');
    m.drawmeridians([-60, -55, -50, -45], labels=[0,0,0,1], fontsize=12, fontweight='normal');
    
    cax = plt.axes([0.85,0.15,0.04,0.7])
    cb = plt.colorbar(c, cax=cax, ticks=levels)
    cb.set_label(r'$\rm T(^{\circ}C)$', fontsize=12, fontweight='normal')
    #plt.subplots_adjust(left=.07, bottom=.07, right=.93, top=.9, wspace=.2, hspace=.2)
    

    #### ---- Save Figure ---- ####
    if np.size(figure_file):
        fig.set_size_inches(w=6, h=9)
        fig.set_dpi(200)
        fig.savefig(figure_file)
        print [figure_file + ' saved!']

    plt.show()
        
        
def Tbot_to_GIS_ascii(h5file, ascfile):
    """ Read bottom temperature H5 file and export it to ascii readagle for GIS
    
    Usage ex:
    import azmp_utils as azu
    azu.Tbot_to_GIS_ascii('Tbot_climato_fall.h5', 'bottom_temp.asc'):
   
    """    

    ## ---- Load data ---- ##    
    h5f = h5py.File(h5file,'r')
    try:
        Tbot = h5f['Tbot'][:]
    except:
        try:
            Tbot = h5f['Tsurf'][:]
        except:
            Tbot = h5f['Sbot'][:]
            
    lon_reg = h5f['lon_reg'][:]
    lat_reg = h5f['lat_reg'][:]
    h5f.close()
    
    #### ---- Save CSV data ---- ####
    # find resolution
    dc = np.round(np.diff(lat_reg)[0], 2)
    # Replace NaNs by -9999
    idx_nan = np.where(np.isnan(Tbot))
    Tbot[idx_nan] = -9999
    # define header
    Tbot_flip = np.flipud(Tbot)
    header = '{0:^1s} {1:^1s}\n{2:^1s} {3:^1s}\n{4:^1s} {5:^1s}\n{6:^1s} {7:^1s}\n{8:^1s} {9:^1s}\n{10:^1s} {11:^1s}'.format('NCOLS', np.str(lon_reg.size), 'NROWS', np.str(lat_reg.size), 'XLLCORNER', np.str(lon_reg.min()), 'YLLCORNER', np.str(lat_reg.min()), 'CELLSIZE', np.str(dc), 'NODATA_VALUE', '-9999')
    np.savetxt(ascfile, Tbot_flip, delimiter=" ", header=header, fmt='%5.2f', comments='')

    
def polygon_temperature_stats(dict, shape, nsrf=False, var='temperature'):
    """ to compute some stats about temperature in nafo sub-division
    (e.g., to build ResDoc Scorecards)
    
    Usage ex (for ResDocs):
    import azmp_utils as azu
    from shapely.geometry.polygon import Polygon
    from shapely.ops import cascaded_union
    Tdict = azu.get_bottomT(year_file='/home/cyrf0006/data/dev_database/2017.nc', season='spring', climato_file='Tbot_climato_spring_0.10.h5')
    nafo_div = azu.get_nafo_divisions()
    polygon3L = Polygon(zip(nafo_div['3L']['lon'], nafo_div['3L']['lat']))
    polygon3N = Polygon(zip(nafo_div['3N']['lon'], nafo_div['3N']['lat']))
    polygon3O = Polygon(zip(nafo_div['3O']['lon'], nafo_div['3O']['lat']))
    shape_3LNO = [polygon3L, polygon3N, polygon3O]
    new_shape = cascaded_union(shape_3LNO)
    dict = azu.polygon_temperature_stats(Tdict, new_shape)

    But the stats can be computed on any polygon (cf. shrimp's SFA aras)

    Modified Feb. 2023: Added var option for salinity
    """

    # Output from
    map = {}
    map_org = {}
    for year in dict:
        if var == 'temperature':
            map[year] = dict[year]['Tbot']
            map_org[year] = dict[year]['Tbot_orig']
        elif 'salinity':
            map[year] = dict[year]['Sbot']
            map_org[year] = dict[year]['Sbot_orig']
    bathy = dict[year]['bathy']
    lon_reg = dict[year]['lon_reg']
    lat_reg = dict[year]['lat_reg']

    # derive mean pixel area
    obj = {'type':'Polygon','coordinates':[[[lon_reg[0],lat_reg[0]],[lon_reg[0],lat_reg[-1]],[lon_reg[-1],lat_reg[-1]],[lon_reg[-1],lat_reg[0]],[lon_reg[0],lat_reg[0]]]]}
    pixel_area = area(obj)/1e6/map[year].size

    #Isolate the bathymetry (remove land pixels)
    bath_mask = bathy.astype(float)*-1
    bath_mask[bath_mask <= 10] = np.nan
    bath_mask[bath_mask > 1000] = np.nan

    # select data in polygon
    for i, xx in enumerate(lon_reg):
        for j,yy in enumerate(lat_reg):
            point = Point(lon_reg[i], lat_reg[j])
            if shape.contains(point) == False:
                bath_mask[j,i] = np.nan

    #Cycle through each of the years
    data_vec = {}
    for year in dict:
        data_vec[year] = map_org[year][~np.isnan(bath_mask)]
    bathy_vec_org = bathy[~np.isnan(bath_mask)]

    #Cycle through each year
    dict_together = {}
    for year in map:

        #Determine if any data is present
        if np.size(data_vec[year]) != 0:

            # total area of the polygon (save in pkl)
            total_polygon_area = np.size(data_vec[year])*pixel_area

            # remove nans
            bathy_vec = bathy_vec_org[~np.isnan(data_vec[year])]
            data_vec[year] = data_vec[year][~np.isnan(data_vec[year])]
            
            # mean temperature all polygon
            Tmean = data_vec[year].mean()

            # mean temperature at depth shallower than 100m, 200m, 300m
            Tmean100 = data_vec[year][bathy_vec>=-100].mean()
            Tmean200 = data_vec[year][bathy_vec>=-200].mean()
            Tmean300 = data_vec[year][bathy_vec>=-300].mean()
            
            # area with temperature < 0
            area_colder_0deg = data_vec[year][data_vec[year]<=0].size*pixel_area
            # area with temperature < 1
            area_colder_1deg = data_vec[year][data_vec[year]<=1].size*pixel_area        
            # area with temperature > 2
            area_warmer_2deg = data_vec[year][data_vec[year]>=2].size*pixel_area
            # area with temperature > 2 & < 4 (shrimp habitat)
            area_shrimp = data_vec[year][(data_vec[year]>=2) & (data_vec[year]<=4)].size*pixel_area
            # area with temperature < 2 (crab habitat)
            area_colder_2deg = data_vec[year][data_vec[year]<=2].size*pixel_area
            # % area with temperature < 2 (crab habitat)
            area_colder_2deg_perc = area_colder_2deg/(data_vec[year].size*pixel_area)*100.0

            # Pandalus Borealis habitat
            Pbor = data_vec[year][(bathy_vec>=-460) & (bathy_vec<=-180) &  (data_vec[year]>=-.2) &  (data_vec[year]<=4.7)].size*pixel_area
            Pbor_perc = Pbor/(data_vec[year].size*pixel_area)*100.0
            # Pandalus Montagui habitat
            Pmon = data_vec[year][(bathy_vec>=-600) & (bathy_vec<=-110) &  (data_vec[year]>=-1) &  (data_vec[year]<=3.7)].size*pixel_area
            Pmon_perc = Pmon/(data_vec[year].size*pixel_area)*100.0

            # Measure of the successfulness of the sampling
            sampled_area= data_vec[year].size*pixel_area

            #Determine the percentage of missing area
            percent_coverage = np.sum(~np.isnan(map_org[year]*bath_mask))/np.sum(~np.isnan(bath_mask))
            percent_coverage = percent_coverage*100


            # Fill dict for output
            dict = {}
            dict['Tmean'] = Tmean # temperature and salinity
            dict['Tmean_sha100'] = Tmean100
            dict['Tmean_sha200'] = Tmean200
            dict['Tmean_sha300'] = Tmean300
            dict['percent_coverage'] = percent_coverage

            if var == 'temperature': # temperature only
                dict['area_colder0'] = area_colder_0deg # <--- now in km2. They are divisded by 1000 in scorecard.
                dict['area_colder1'] = area_colder_1deg 
                dict['area_warmer2'] = area_warmer_2deg
                dict['area_shrimp'] = area_shrimp
                dict['area_colder2'] = area_colder_2deg 
                dict['area_colder2_perc'] = area_colder_2deg_perc 
                dict['area_Pborealis'] = Pbor
                dict['area_Pborealis_perc'] = Pbor_perc
                dict['area_Pmontagui'] = Pmon
                dict['area_Pmontagui_perc'] = Pmon_perc
                dict['sampled_area'] = sampled_area
                dict['total_area'] = total_polygon_area

            
            if nsrf:
                # Area of NSRF seafloor with conditions within a certain depth and temperature range (project with Wojciech)
                Pbor_eaz = data_vec[year][(bathy_vec>=-590) & (bathy_vec<=-180) &  (data_vec[year]>=-.4) &  (data_vec[year]<=4.7)].size*pixel_area
                Pbor_waz = data_vec[year][(bathy_vec>=-520) & (bathy_vec<=-210) &  (data_vec[year]>=-.7) &  (data_vec[year]<=4.0)].size*pixel_area
                Pbor_sfa4 = data_vec[year][(bathy_vec>=-590) & (bathy_vec<=-180) &  (data_vec[year]>=-.7) &  (data_vec[year]<=4.7)].size*pixel_area
                Pmon_eaz = data_vec[year][(bathy_vec>=-600) & (bathy_vec<=-120) &  (data_vec[year]>=-.5) &  (data_vec[year]<=3.7)].size*pixel_area
                Pmon_waz = data_vec[year][(bathy_vec>=-530) & (bathy_vec<=-110) &  (data_vec[year]>=-1.2) &  (data_vec[year]<=2.8)].size*pixel_area
                Pmon_sfa4 = data_vec[year][(bathy_vec>=-590) & (bathy_vec<=-140) &  (data_vec[year]>=-0.9) &  (data_vec[year]<=4.0)].size*pixel_area
                # % of good pixels
                Pbor_eaz_perc = Pbor_eaz/(data_vec[year].size*pixel_area)*100.0
                Pbor_waz_perc = Pbor_waz/(data_vec[year].size*pixel_area)*100.0
                Pbor_sfa4_perc = Pbor_sfa4/(data_vec[year].size*pixel_area)*100.0
                Pmon_eaz_perc = Pmon_eaz/(data_vec[year].size*pixel_area)*100.0
                Pmon_waz_perc = Pmon_waz/(data_vec[year].size*pixel_area)*100.0
                Pmon_sfa4_perc = Pmon_sfa4/(data_vec[year].size*pixel_area)*100.0
                
                # Fill dict for output
                dict['Pbor_eaz_habitat'] = Pbor_eaz
                dict['Pbor_waz_habitat'] = Pbor_waz
                dict['Pbor_sfa4_habitat'] = Pbor_sfa4
                dict['Pmon_eaz_habitat'] = Pmon_eaz
                dict['Pmon_waz_habitat'] = Pmon_waz
                dict['Pmon_sfa4_habitat'] = Pmon_sfa4
                dict['Pbor_eaz_perc'] = Pbor_eaz_perc
                dict['Pbor_waz_perc'] = Pbor_waz_perc
                dict['Pbor_sfa4_perc'] = Pbor_sfa4_perc
                dict['Pmon_eaz_perc'] = Pmon_eaz_perc
                dict['Pmon_waz_perc'] = Pmon_waz_perc
                dict['Pmon_sfa4_perc'] = Pmon_sfa4_perc

        else:
            # Fill dict for output
            dict={}
            dict['Tmean'] = np.nan
            dict['Tmean_sha100'] = np.nan
            dict['Tmean_sha200'] = np.nan
            dict['Tmean_sha300'] = np.nan
            dict['area_colder0'] = np.nan
            dict['area_colder1'] = np.nan 
            dict['area_warmer2'] = np.nan
            dict['area_shrimp'] = np.nan
            dict['area_colder2'] = np.nan 
            dict['area_colder2_perc'] = np.nan 
            dict['area_Pborealis'] = np.nan
            dict['area_Pborealis_perc'] = np.nan
            dict['area_Pmontagui'] = np.nan
            dict['area_Pmontagui_perc'] = np.nan
            dict['sampled_area'] = np.nan
            dict['total_area'] = np.nan
        dict_together[year] = dict

    return dict_together

#def polygon_salinity_stats(dict, shape):
#    """ Almost the same as previous, but for salinity
#    This time, it was no designed for ResDoc, but rather to give average salinity over SFA areas.
#    And unlike previous, only the average is given.
#
#    Dec. 2018
#    """

#    # Output from 
#    map = dict['Sbot']
#    bathy = dict['bathy']
#    lon_reg = dict['lon_reg']
#    lat_reg = dict['lat_reg']

#    # derive mean pixel area
#    obj = {'type':'Polygon','coordinates':[[[lon_reg[0],lat_reg[0]],[lon_reg[0],lat_reg[-1]],[lon_reg[-1],lat_reg[-1]],[lon_reg[-1],lat_reg[0]],[lon_reg[0],lat_reg[0]]]]}
#    pixel_area = area(obj)/1e6/map.size
    
#    # select data in polygon
#    data_vec = []
#    bathy_vec = []
#    for i, xx in enumerate(lon_reg):
#        for j,yy in enumerate(lat_reg):
#            point = Point(lon_reg[i], lat_reg[j])            
#            if shape.contains(point):
#                data_vec = np.append(data_vec, map[j,i])
#                bathy_vec = np.append(bathy_vec, bathy[j,i])
#            else:                
#                pass

#    # remove nans            
#    bathy_vec = bathy_vec[~np.isnan(data_vec)]
#    data_vec = data_vec[~np.isnan(data_vec)]
    
#    # mean temperature all polygon
#    Smean = data_vec.mean()

#    # mean temperature at depth shallower than 100m, 200m, 300m
#    Smean100 = data_vec[bathy_vec>=-100].mean()
#    Smean200 = data_vec[bathy_vec>=-200].mean()
#    Smean300 = data_vec[bathy_vec>=-300].mean()
#    
#    # Fill dict for output
#    dict = {}
#    dict['Smean'] = Smean
#    dict['Smean_sha100'] = Smean100
#    dict['Smean_sha200'] = Smean200
#    dict['Smean_sha300'] = Smean300

#    return dict


def masterfile_section_to_multiindex(section, z_vec):
    """
    This function reads the master file Excel sheet and pickle a Multi-index DataFrame organized as:


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

    ## ----  Load biochemical data ---- ##
    df = pd.read_excel('/home/cyrf0006/github/AZMP-NL/data/AZMP_Nutrients_1999_2016_Good_Flags_Only.xlsx')
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



def h5cubes_2nc(h5_fileT, h5_fileS, ncfile): 
    '''
    Function to convert .h5 TCube and Scube into a netCDF file
    Adapted from script h52nc_cube.py.

    Option to add T and S, or just T.

    Originally created for a collaboration with S. Rousseau.
    
    February 2023
    Frederic.Cyr@dfo-mpo.gc.ca
    '''

    h5f = h5py.File(h5_fileT, 'r')
    T = h5f['temperature'][:]
    lon_reg = h5f['lon_reg'][:]
    lat_reg = h5f['lat_reg'][:]
    z = h5f['z'][:]
    h5f.close()

    if h5_fileS:
        h5f = h5py.File(h5_fileS, 'r')
        S = h5f['salinity'][:]
        h5f.close()    

    ## Save Climatology in NetCDF (for sharing purposes)
    print('creating ' + ncfile)
    # create dimension
    ds = netCDF4.Dataset(ncfile, 'w', format='NETCDF4')
    lat_dim = ds.createDimension('lat', len(lat_reg)) # latitude axis
    lon_dim = ds.createDimension('lon', len(lon_reg)) # longitude axis
    depth_dim = ds.createDimension('depth', len(z))
    # title
    ds.title='3D temperature and salinity [optionnal]'
    # create variable
    lat = ds.createVariable('lat', np.float32, ('lat',))
    lat.units = 'degrees_north'
    lat.long_name = 'latitude'
    lon = ds.createVariable('lon', np.float32, ('lon',))
    lon.units = 'degrees_east'
    lon.long_name = 'longitude'
    depth = ds.createVariable('depth', np.float64, ('depth',))
    depth.units = 'm'
    depth.long_name = 'depth'
    temp = ds.createVariable('temp',np.float64,('lat','lon', 'depth')) # note: unlimited dimension is leftmost
    temp.units = 'degC'
    temp.standard_name = 'temperature'
    print(temp)
    print("-- Some pre-defined attributes for variable temp:")
    print("temp.dimensions:", temp.dimensions)
    print("temp.shape:", temp.shape)
    print("temp.dtype:", temp.dtype)
    print("temp.ndim:", temp.ndim)

    # writing data
    lat[:] = lat_reg
    lon[:] = lon_reg
    depth[:] = z
    temp[:,:,:] = T # Appends data along unlimited dimension

    print("-- Wrote data, temp.shape is now ", temp.shape)
    print("-- Min/Max values:", temp[:,:,:].min(), temp[:,:,:].max())

    # Optional salinity
    if h5_fileS:
         sal = ds.createVariable('sal',np.float64,('lat','lon', 'depth')) 
         sal.units = 'psu'
         sal.standard_name = 'salinity'
         sal[:,:,:] = S # Appends data along unlimited dimension        
    
    print(ds)
    # close the Dataset.
    ds.close(); print('Dataset is closed!')

    return None

