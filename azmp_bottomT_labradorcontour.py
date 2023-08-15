
import matplotlib.pyplot as plt
import matplotlib
matplotlib.interactive(True)
import numpy as np
import xarray as xr
from scipy.ndimage import gaussian_filter

'''

The purpose of this script is to create a bathymetry contour for Labrador.
This will create the file 100m_contour_labrador.npy

'''


#For testing purposes remove after
test = np.load('100m_contour_labrador.npy')

#Define the path and load in bathymetry data
bathy_path = '/home/jcoyne/Documents/Datasets/GEBCO_2023/'
ds_bathy = xr.open_dataset(bathy_path + 'GEBCO_2023_sub_ice_topo.nc')

#Process the bathymetry
lonLims = [-65,-55]
latLims = [53,61]
ds_bathy = ds_bathy.isel(lon=(ds_bathy.lon>=lonLims[0])*(ds_bathy.lon<=lonLims[1]))
ds_bathy = ds_bathy.isel(lat=(ds_bathy.lat>=latLims[0])*(ds_bathy.lat<=latLims[1]))

# Extract latitude and longitude
lon = ds_bathy.lon.values
lat = ds_bathy.lat.values

# Extract the elevation
Z = ds_bathy.elevation.values

#Labrador bathymetry is complex, apply smoothing
Z_filtered = gaussian_filter(Z, sigma=5)

#Extract the coordinate of the 100m depth contour
coords = plt.contour(lon,lat,Z_filtered,levels=np.array([-100,0]))
coords = coords.allsegs[0][0]
plt.close()

#Place a coordinate point btw top and bottom to make better shape
coords = np.concatenate((coords, np.array([[lonLims[0],latLims[0]]])))

#Save the file
np.save('100m_contour_labrador.npy', coords)
