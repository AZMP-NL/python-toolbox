import numpy as np
import netCDF4
import h5py
import xarray as xr
import pandas as pd
import datetime
from netCDF4 import date2num,num2date

h5_file = 'Tcube_spring2022.h5'
h5f = h5py.File(h5_file, 'r')
T = h5f['T'][:]
lon_reg = h5f['lon_reg'][:]
lat_reg = h5f['lat_reg'][:]
z = h5f['z'][:]
h5f.close()

## Save Climatology in NetCDF (for sharing purposes)
ncfile = h5_file.split('.')[0] + '.nc'
print('creating ' + ncfile)
# create dimension
ds = netCDF4.Dataset(ncfile, 'w', format='NETCDF4')
lat_dim = ds.createDimension('lat', len(lat_reg)) # latitude axis
lon_dim = ds.createDimension('lon', len(lon_reg)) # longitude axis
depth_dim = ds.createDimension('depth', len(z))
# title
ds.title='3D temperature'
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
data_arr = T 
temp[:,:,:] = data_arr # Appends data along unlimited dimension
print("-- Wrote data, temp.shape is now ", temp.shape)
print("-- Min/Max values:", temp[:,:,:].min(), temp[:,:,:].max())

print(ds)
# close the Dataset.
ds.close(); print('Dataset is closed!')
