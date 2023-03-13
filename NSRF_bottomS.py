'''
Bottom salinity maps for Shrimp RAP

To generate bottom climato:
import numpy as np
import azmp_utils as azu
dc = .10
lonLims = [-66, -56] # Lab Sea
latLims = [57, 67]
lon_reg = np.arange(lonLims[0]+dc/2, lonLims[1]-dc/2, dc)
lat_reg = np.arange(latLims[0]+dc/2, latLims[1]-dc/2, dc)
Sbot_dict = azu.get_bottomS_climato('/home/cyrf0006/data/dev_database/*.nc', lon_reg, lat_reg, year_lims=[2006, 2018], season='summer', h5_outputfile='Sbot_climato_NSRF_summer_2006-2018.h5')

* see: /home/cyrf0006/AZMP/state_reports/bottomT

'''

import os
import netCDF4
import h5py
import xarray as xr
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

def draw_screen_poly( lats, lons, m):
    x, y = m( lons, lats )
    xy = zip(x,y)
    poly = PP( xy, facecolor=[.8, .8, .8])
    plt.gca().add_patch(poly)

## ---- Region parameters ---- ## <-------------------------------Would be nice to pass this in a config file '2017.report'
dataFile = '/home/cyrf0006/data/GEBCO/GRIDONE_1D.nc'
lon_0 = -50
lat_0 = 50
#lonLims = [-60, -44] # FishHab region
#latLims = [39, 56]
proj = 'merc'
zmax = 1000 # do try to compute bottom sal below that depth
zmin = 10
dz = 5 # vertical bins

season = 'summer'
year = 2022
climato_file = 'Sbot_climato_NSRF_summer_2006-2018.h5'
year_file = '/home/cyrf0006/data/dev_database/netCDF/' + year + '.nc'


## ---- Load Climato data ---- ##    
print('Load ' + climato_file)
h5f = h5py.File(climato_file, 'r')
Sbot_climato = h5f['Sbot'][:]
lon_reg = h5f['lon_reg'][:]
lat_reg = h5f['lat_reg'][:]
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

## ---- Get CTD data --- ##
print('Get ' + year_file)
ds = xr.open_mfdataset(year_file)
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
    ds = ds.sel(time=((ds['time.month']>=10)) & ((ds['time.month']<=12)))
else:
    print('!! no season specified, used them all! !!')

# Remome problematic datasets
print('!!Remove MEDBA data!!')
print('  ---> I Should be improme because I remove good data!!!!')
ds = ds.where(ds.instrument_ID!='MEDBA', drop=True)

    
# Vertical binning (on dataArray; more appropriate here
da_sal = ds['salinity']
lons = np.array(ds.longitude)
lats = np.array(ds.latitude)
#bins = np.arange(dz/2.0, ds.level.max(), dz)
bins = np.arange(dz/2.0, 1000, dz)
da_sal = da_sal.groupby_bins('level', bins).mean(dim='level')
#To Pandas Dataframe
df_sal = da_sal.to_pandas()
df_sal.columns = bins[0:-1] #rename columns with 'bins'
# Remove empty columns
idx_empty_rows = df_sal.isnull().all(1).nonzero()[0]
df_sal = df_sal.dropna(axis=0,how='all')
lons = np.delete(lons,idx_empty_rows)
lats = np.delete(lats,idx_empty_rows)
print(' -> Done!')


## --- fill 3D cube --- ##  
print('Fill regular cube')
z = df_sal.columns.values
V = np.full((lat_reg.size, lon_reg.size, z.size), np.nan)

# Aggregate on regular grid
for i, xx in enumerate(lon_reg):
    for j, yy in enumerate(lat_reg):    
        idx = np.where((lons>=xx-dc/2) & (lons<xx+dc/2) & (lats>=yy-dc/2) & (lats<yy+dc/2))
        tmp = np.array(df_sal.iloc[idx].mean(axis=0))
        idx_good = np.argwhere(~np.isnan(tmp))
        if np.size(idx_good)==1:
            V[j,i,:] = np.array(df_sal.iloc[idx].mean(axis=0))
        elif np.size(idx_good)>1: # vertical interpolation between pts
            interp = interp1d(np.squeeze(z[idx_good]), np.squeeze(tmp[idx_good]))  # <---------- Pay attention here, this is a bit unusual, but seems to work!
            idx_interp = np.arange(np.int(idx_good[0]),np.int(idx_good[-1]+1))
            V[j,i,idx_interp] = interp(z[idx_interp]) # interpolate only where possible (1st to last good idx)

                   
# horizontal interpolation at each depth
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
    if idx_good.size>4: # will ignore depth where no data exist
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

# getting bottom salinity
print('Getting bottom Sal.')    
Sbot = np.full([lat_reg.size,lon_reg.size], np.nan) 
for i, xx in enumerate(lon_reg):
    for j,yy in enumerate(lat_reg):
        bottom_depth = -Zitp[j,i] # minus to turn positive
        sal_vec = V[j,i,:]
        ## idx_no_good = np.argwhere(sal_vec>30)
        ## if idx_no_good.size:
        ##     sal_vec[idx_no_good] = np.nan
        idx_good = np.squeeze(np.where(~np.isnan(sal_vec)))
        if idx_good.size:
            idx_closest = np.argmin(np.abs(bottom_depth-z[idx_good]))
        else:
            continue

        if np.abs([idx_closest] - bottom_depth) <= 20:
            Sbot[j,i] = sal_vec[idx_good[idx_closest]]
        elif np.abs(z[idx_closest] - bottom_depth) <= 50:
            #print('used data located [30,50]m from bottom')
            Sbot[j,i] = sal_vec[idx_good[idx_closest]]
            
print(' -> Done!')    

# Mask data outside Nafo div.
print('Mask according to NAFO division for ' + season)
# Polygons
polygon3K = Polygon(zip(nafo_div['3K']['lon'], nafo_div['3K']['lat']))
polygon3L = Polygon(zip(nafo_div['3L']['lon'], nafo_div['3L']['lat']))
polygon3N = Polygon(zip(nafo_div['3N']['lon'], nafo_div['3N']['lat']))
polygon3O = Polygon(zip(nafo_div['3O']['lon'], nafo_div['3O']['lat']))
polygon3Ps = Polygon(zip(nafo_div['3Ps']['lon'], nafo_div['3Ps']['lat']))
polygon2J = Polygon(zip(nafo_div['2J']['lon'], nafo_div['2J']['lat']))

# Contour of data to mask
contour_mask = np.load('100m_contour_labrador.npy')
polygon_mask = Polygon(contour_mask)


if season == 'spring':
    for i, xx in enumerate(lon_reg):
        for j,yy in enumerate(lat_reg):
            point = Point(lon_reg[i], lat_reg[j])
            #if (~polygon3L.contains(point)) & (~polygon3N.contains(point)) & (~polygon3O.contains(point)) & (~polygon3Ps.contains(point)):
            if polygon3L.contains(point) | polygon3N.contains(point) | polygon3O.contains(point) | polygon3Ps.contains(point):
                pass #nothing to do but cannot implement negative statement "if not" above
            else:
                Sbot[j,i] = np.nan
            
                
elif season == 'fall':
    for i, xx in enumerate(lon_reg):
        for j,yy in enumerate(lat_reg):
            point = Point(lon_reg[i], lat_reg[j])
            if polygon2J.contains(point) | polygon3K.contains(point) | polygon3L.contains(point) | polygon3N.contains(point) | polygon3O.contains(point):
                pass #nothing to do but cannot implement negative statement "if not" above
            else:
                Sbot[j,i] = np.nan ### <--------------------- Do mask the fall / OR / 

            if polygon_mask.contains(point): # mask data near Labrador in fall
                Sbot[j,i] = np.nan 
else:
    print('no division mask, all data taken')
            
print(' -> Done!')    

# Salinity anomaly:
anom = Sbot-Sbot_climato


## ---- Plot Anomaly ---- ##
fig, ax = plt.subplots(nrows=1, ncols=1)
m = Basemap(ax=ax, projection='merc',lon_0=lon_0,lat_0=lat_0, llcrnrlon=lonLims[0],llcrnrlat=latLims[0],urcrnrlon=lonLims[1],urcrnrlat=latLims[1], resolution= 'i')
levels = np.array([-1, -.8, -.6, -.4, -.2, .2, .4, .6, .8, 1])
xi, yi = m(*np.meshgrid(lon_reg, lat_reg))
c = m.contourf(xi, yi, anom, levels, cmap=plt.cm.RdBu_r, extend='both')
cc = m.contour(xi, yi, -Zitp, [100, 500, 1000, 4000], colors='grey');
plt.clabel(cc, inline=1, fontsize=10, fmt='%d')
if season=='fall':
    plt.title('Fall Bottom Salinity ' + year + ' Anomaly')
elif season=='spring':
    plt.title('Spring Bottom Salinity ' + year + ' Anomaly')
else:
    plt.title('Bottom Salinity ' + year + '  Anomaly')
m.fillcontinents(color='tan');
m.drawparallels(np.arange(58, 68, 2), labels=[1,0,0,0], fontsize=12, fontweight='normal');
m.drawmeridians(np.arange(-66, -54, 2), labels=[0,0,0,1], fontsize=12, fontweight='normal');
cax = fig.add_axes([0.16, 0.05, 0.7, 0.025])
cb = plt.colorbar(c, cax=cax, orientation='horizontal')
cb.set_label(r'S', fontsize=12, fontweight='normal')
div_toplot = ['2J', '3K', '3L', '3N', '3O', '3Ps']
for div in div_toplot:
    div_lon, div_lat = m(nafo_div[div]['lon'], nafo_div[div]['lat'])
    m.plot(div_lon, div_lat, 'k', linewidth=2)
    ax.text(np.mean(div_lon), np.mean(div_lat), div, fontsize=12, color='black', fontweight='bold')    
# Save Figure
fig.set_size_inches(w=6, h=9)
fig.set_dpi(300)
outfile = 'NSRF_bottom_sal_anomaly_' + season + '_' + year + '.png'
fig.savefig(outfile)
os.system('convert -trim ' + outfile + ' ' + outfile)

## ---- Plot Salinity ---- ##
fig, ax = plt.subplots(nrows=1, ncols=1)
m = Basemap(ax=ax, projection='merc',lon_0=lon_0,lat_0=lat_0, llcrnrlon=lonLims[0],llcrnrlat=latLims[0],urcrnrlon=lonLims[1],urcrnrlat=latLims[1], resolution= 'i')
#levels = np.linspace(-2, 6, 9)
levels = np.linspace(30, 36, 13)
xi, yi = m(*np.meshgrid(lon_reg, lat_reg))
c = m.contourf(xi, yi, Sbot, levels, cmap=plt.cm.RdBu_r, extend='both')
cc = m.contour(xi, yi, -Zitp, [100, 500, 1000, 4000], colors='grey');
plt.clabel(cc, inline=1, fontsize=10, fmt='%d')
if season=='fall':
    plt.title('Fall Bottom Salinity ' + year)
elif season=='spring':
    plt.title('Spring Bottom Salinity ' + year)
else:
    plt.title('Bottom Salinity ' + year)
m.fillcontinents(color='tan');
m.drawparallels(np.arange(58, 68, 2), labels=[1,0,0,0], fontsize=12, fontweight='normal');
m.drawmeridians(np.arange(-66, -54, 2), labels=[0,0,0,1], fontsize=12, fontweight='normal');
x, y = m(lons, lats)
m.scatter(x,y, s=50, marker='.',color='k')
cax = fig.add_axes([0.16, 0.05, 0.7, 0.025])
cb = plt.colorbar(c, cax=cax, orientation='horizontal')
cb.set_label(r'S', fontsize=12, fontweight='normal')
div_toplot = ['2J', '3K', '3L', '3N', '3O', '3Ps']
for div in div_toplot:
    div_lon, div_lat = m(nafo_div[div]['lon'], nafo_div[div]['lat'])
    m.plot(div_lon, div_lat, 'k', linewidth=2)
    ax.text(np.mean(div_lon), np.mean(div_lat), div, fontsize=12, color='black', fontweight='bold')
# Save Figure
fig.set_size_inches(w=6, h=9)
fig.set_dpi(200)
outfile = 'NSRF_bottom_sal_' + season + '_' + year + '.png'
fig.savefig(outfile)
os.system('convert -trim ' + outfile + ' ' + outfile)

## ---- Plot Climato ---- ##
fig, ax = plt.subplots(nrows=1, ncols=1)
m = Basemap(ax=ax, projection='merc',lon_0=lon_0,lat_0=lat_0, llcrnrlon=lonLims[0],llcrnrlat=latLims[0],urcrnrlon=lonLims[1],urcrnrlat=latLims[1], resolution= 'i')
levels = np.linspace(30, 36, 13)
xi, yi = m(*np.meshgrid(lon_reg, lat_reg))
c = m.contourf(xi, yi, Sbot_climato, levels, cmap=plt.cm.RdBu_r, extend='both')
cc = m.contour(xi, yi, -Zitp, [100, 500, 1000, 4000], colors='grey');
plt.clabel(cc, inline=1, fontsize=10, fmt='%d')
if season=='fall':
    plt.title('Fall Bottom Salinity Climatology')
elif season=='spring':
    plt.title('Spring Bottom Salinity Climatology')
else:
    plt.title('Bottom Salinity Climatology (2006-2018)')
m.fillcontinents(color='tan');
m.drawparallels(np.arange(58, 68, 2), labels=[1,0,0,0], fontsize=12, fontweight='normal');
m.drawmeridians(np.arange(-66, -54, 2), labels=[0,0,0,1], fontsize=12, fontweight='normal');
cax = fig.add_axes([0.16, 0.05, 0.7, 0.025])
cb = plt.colorbar(c, cax=cax, orientation='horizontal')
cb.set_label(r'S', fontsize=12, fontweight='normal')
div_toplot = ['2J', '3K', '3L', '3N', '3O', '3Ps']
for div in div_toplot:
    div_lon, div_lat = m(nafo_div[div]['lon'], nafo_div[div]['lat'])
    m.plot(div_lon, div_lat, 'k', linewidth=2)
    ax.text(np.mean(div_lon), np.mean(div_lat), div, fontsize=12, color='black', fontweight='bold')    
# Save Figure
fig.set_size_inches(w=6, h=9)
fig.set_dpi(300)
outfile = 'NSRF_bottom_sal_climato_' + season + '_' + year + '.png'
fig.savefig(outfile)
os.system('convert -trim ' + outfile + ' ' + outfile)


# Convert to a subplot
os.system('montage NSRF_bottom_sal_climato_' + season + '_' + year + '.png NSRF_bottom_sal_' + season + '_' + year + '.png NSRF_bottom_sal_anomaly_' + season + '_' + year + '.png  -tile 3x1 -geometry +10+10  -background white  NSRF_bottomS_' + season + year + '.png') 
os.system('rm NSRF_bottom_sal_climato_' + season + '_' + year + '.png NSRF_bottom_sal_' + season + '_' + year + '.png NSRF_bottom_sal_anomaly_' + season + '_' + year + '.png')
