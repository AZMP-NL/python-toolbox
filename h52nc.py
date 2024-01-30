import numpy as np
import netCDF4
import h5py
import xarray as xr
import pandas as pd
import datetime
from netCDF4 import date2num,num2date

#h5_file = 'Tbot_climato_NSRFxx_summer_2015-2020.h5'
#h5_file = 'Tbot_climato_NSRFxx_summer_2006-2010.h5'
h5_file = 'Tbot_climato_NSRFxx_summer_2011-2015.h5'
#h5_file = 'Tbot_climato_NSRFxx_summer_2016-2021.h5'
#h5_file = 'Tbot_climato_NSRFxx_summer_2016-2020.h5'
#h5_file = 'Tbot_climato_spring_0.10_2015-2020.h5'
#h5_file = 'Tbot_climato_spring_0.10_2006-2010.h5'
#h5_file = 'Tbot_climato_spring_0.10_2011-2015.h5'
#h5_file = 'Tbot_climato_spring_0.10_2016-2021.h5'
#h5_file = 'Tbot_climato_spring_0.10_2016-2020.h5'
#h5_file = 'Tbot_climato_fall_0.10_2006-2010.h5'
#h5_file = 'Tbot_climato_fall_0.10_2011-2015.h5'
#h5_file = 'Tbot_climato_fall_0.10_2016-2021.h5'
#h5_file = 'Tbot_climato_fall_0.10_2016-2020.h5'

h5f = h5py.File(h5_file, 'r')
Tbot_climato = h5f['Tbot'][:]
lon_reg = h5f['lon_reg'][:]
lat_reg = h5f['lat_reg'][:]
Zitp = h5f['Zitp'][:]
h5f.close()

## Save Climatology in NetCDF (for sharing purposes)
ncfile = h5_file.strip('h5') + 'nc'
print('creating ' + ncfile)
# create dimension
ds = netCDF4.Dataset(ncfile, 'w', format='NETCDF4')
lat_dim = ds.createDimension('lat', len(lat_reg)) # latitude axis
lon_dim = ds.createDimension('lon', len(lon_reg)) # longitude axis
time_dim = ds.createDimension('time', None)
# title
ds.title='NSRF - extended - bottom temperature climatology'
ds.author='Frederic.Cyr@dfo-mpo.gc.ca'
# create variable
lat = ds.createVariable('lat', np.float32, ('lat',))
lat.units = 'degrees_north'
lat.long_name = 'latitude'
lon = ds.createVariable('lon', np.float32, ('lon',))
lon.units = 'degrees_east'
lon.long_name = 'longitude'
time = ds.createVariable('time', np.float64, ('time',))
time.units = 'hours since 1970-01-01'
time.long_name = 'time'
temp = ds.createVariable('temp',np.float64,('time','lat','lon')) # note: unlimited dimension is leftmost
temp.units = 'degC' # degrees Kelvin
temp.standard_name = 'bottom_temperature' # this is a CF standard name
print(temp)

print("-- Some pre-defined attributes for variable temp:")
print("temp.dimensions:", temp.dimensions)
print("temp.shape:", temp.shape)
print("temp.dtype:", temp.dtype)
print("temp.ndim:", temp.ndim)
# writing data
ntimes = 1
lat[:] = lat_reg
lon[:] = lon_reg
data_arr = Tbot_climato 
temp[0,:,:] = data_arr # Appends data along unlimited dimension
print("-- Wrote data, temp.shape is now ", temp.shape)
print("-- Min/Max values:", temp[:,:,:].min(), temp[:,:,:].max())

dates = datetime.datetime(2000,1,1,0)
times = date2num(dates, time.units)
time[:] = times
print(ds)
# close the Dataset.
ds.close(); print('Dataset is closed!')
