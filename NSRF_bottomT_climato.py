'''
Bottom temperature maps for NSRF surveys
Script originally developped for 2019 Northern Shrimp RAP


To generate bottom climato:
import numpy as np
import azmp_utils as azu
dc = .10
lonLims = [-66, -56] # Lab Sea
latLims = [57, 67]
lon_reg = np.arange(lonLims[0]+dc/2, lonLims[1]-dc/2, dc)
lat_reg = np.arange(latLims[0]+dc/2, latLims[1]-dc/2, dc)
Tbot_dict = azu.get_bottomT_climato('/home/cyrf0006/data/dev_database/*.nc', lon_reg, lat_reg, year_lims=[2006, 2018], season='summer', h5_outputfile='Tbot_climato_NSRF_summer_2006-2018.h5')

# Extend to include Hudson Strait
import numpy as np
import azmp_utils as azu
dc = .10
lonLims = [-70, -56] # Lab Sea
latLims = [57, 67]
lon_reg = np.arange(lonLims[0]+dc/2, lonLims[1]-dc/2, dc)
lat_reg = np.arange(latLims[0]+dc/2, latLims[1]-dc/2, dc)
Tbot_dict = azu.get_bottomT_climato('/home/cyrf0006/data/dev_database/netCDF/*.nc', lon_reg, lat_reg, year_lims=[2006, 2021], season='summer', h5_outputfile='Tbot_climato_NSRFx_summer_2006-2021.h5')


Frederic.Cyr@dfo-mpo.gc.ca
February 2018

'''

import os
import netCDF4
import h5py
import xarray as xr
os.environ['PROJ_LIB'] = '/home/cyrf0006/anaconda3/share/proj'
from mpl_toolkits.basemap import Basemap
import numpy as  np
import matplotlib.pyplot as plt
import openpyxl, pprint
import pandas as pd
from scipy.interpolate import griddata
from scipy.interpolate import interp1d  # to remove NaNs in profiles
from scipy.interpolate import RegularGridInterpolator as rgi
import azmp_utils as azu
from scipy.ndimage.filters import gaussian_filter
import scipy.ndimage
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from matplotlib.patches import Rectangle
from matplotlib.patches import Polygon as PP
import cmocean

## ---- Parameters ---- ## 
outfile = 'NSRF_bottom_temp_climato_summer_2006-2018.png'
climato_file = 'Tbot_climato_NSRF_summer_2006-2018.h5'

## ---- Load Climato data ---- ##    
print('Load ' + climato_file)
h5f = h5py.File(climato_file, 'r')
Tbot_climato = h5f['Tbot'][:]
lon_reg = h5f['lon_reg'][:]
lat_reg = h5f['lat_reg'][:]
lon_orig = h5f['lon_orig'][:]
lat_orig = h5f['lat_orig'][:]
Zitp = h5f['Zitp'][:]
h5f.close()

## ---- Derive some parameters ---- ##    
lon_0 = np.round(np.mean(lon_reg))
lat_0 = np.round(np.mean(lat_reg))
lonLims = [lon_reg[0], lon_reg[-1]]
latLims = [lat_reg[0], lat_reg[-1]]
lon_grid, lat_grid = np.meshgrid(lon_reg,lat_reg)
dc = np.diff(lon_reg[0:2])

## ---- NAFO divisions ---- ##
nafo_div = azu.get_nafo_divisions()

## ---- LabSea core areas ---- ##
area1_file = '/home/cyrf0006/AZMP/state_reports/bottomT/Areas_Lab_core_1_modif.csv'
area2_file = '/home/cyrf0006/AZMP/state_reports/bottomT/Areas_Lab_core_2_modif.csv'
A1 = pd.read_csv(area1_file, delimiter=r",")
A2 = pd.read_csv(area2_file, delimiter=r",")

## ---- Plot Climato ---- ##
fig, ax = plt.subplots(nrows=1, ncols=1)
m = Basemap(ax=ax, projection='merc',lon_0=lon_0,lat_0=lat_0, llcrnrlon=lonLims[0],llcrnrlat=latLims[0],urcrnrlon=lonLims[1],urcrnrlat=latLims[1], resolution= 'i')
#levels = np.linspace(-2, 6, 9)
levels = np.linspace(-2, 6, 17)
xi, yi = m(*np.meshgrid(lon_reg, lat_reg))
c = m.contourf(xi, yi, Tbot_climato, levels, cmap=plt.cm.RdBu_r, extend='both')
cc = m.contour(xi, yi, -Zitp, [100, 500, 1000, 4000], colors='grey');
#c = m.pcolormesh(xi,yi,Tbot_climato, clim=[-2,6], cmap=plt.cm.RdBu_r)
#c = m.pcolormesh(xi,yi,Tbot_climato, clim=[-2,6], cmap=cmocean.cm.thermal)
c.set_clim(-2,6)
plt.clabel(cc, inline=1, fontsize=10, fmt='%d')
m.fillcontinents(color='tan');
m.drawparallels([54, 55, 56, 57, 58, 59], labels=[1,0,0,0], fontsize=12, fontweight='normal');
m.drawmeridians([-62, -60, -58, -56], labels=[0,0,0,1], fontsize=12, fontweight='normal');
# plot casts
x, y = m(lon_orig, lat_orig)
m.scatter(x,y, s=50, marker='.',color='k')
# plot areas
#x, y = m(A1.long_dd.values, A1.lat_dd.values)
#m.plot(x,y, '-k', linewidth=3)
#x, y = m(-63, 57.1)
#plt.text(x, y, '  Area 1', fontsize=12, fontweight='bold')
#x, y = m(A2.long_dd.values, A2.lat_dd.values)
#m.plot(x,y, '-k', linewidth=3)
#x, y = m(-60, 54.5)
#plt.text(x, y, ' Area 2', fontsize=12, fontweight='bold')
m.fillcontinents(color='tan');

cax = fig.add_axes([0.16, 0.05, 0.7, 0.025])
cb = plt.colorbar(c, cax=cax, orientation='horizontal')
cb.set_label(r'$\rm T(^{\circ}C)$', fontsize=12, fontweight='normal')
# Save Figure
fig.set_size_inches(w=6, h=9)
fig.set_dpi(300)
fig.savefig(outfile)
os.system('convert -trim ' + outfile + ' ' + outfile)

