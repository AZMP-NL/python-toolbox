'''
To generate bottom climato for SA4 only:

import numpy as np
import azmp_utils as azu
dc = .1
lonLims = [-68, -55] # NAFO 4
latLims = [38, 48]
lon_reg = np.arange(lonLims[0]+dc/2, lonLims[1]-dc/2, dc)
lat_reg = np.arange(latLims[0]+dc/2, latLims[1]-dc/2, dc)
azu.get_bottomT_climato('/home/cyrf0006/data/dev_database/*.nc', lon_reg, lat_reg, season='summer', h5_outputfile='Tbot_climato_SA4_summer_0.10.h5') 

'''

## mport netCDF4
import h5py
## import xarray as xr
from mpl_toolkits.basemap import Basemap
import numpy as  np
import matplotlib.pyplot as plt
## import openpyxl, pprint
import pandas as pd
from scipy.interpolate import griddata
from scipy.interpolate import interp1d  # to remove NaNs in profiles
from scipy.interpolate import RegularGridInterpolator as rgi
import azmp_utils as azu
from shapely.geometry.polygon import Polygon
from shapely.ops import cascaded_union


## ---- preamble ---- ##
#years = np.arange(1950, 2019)
years = np.arange(1980, 2019)
lon_0 = -50
lat_0 = 50
proj = 'merc'
plot = False # to plot or not to plot...
season = 'summer'

# load climato
if season == 'fall':
    climato_file = 'Tbot_climato_SA4_fall_0.10.h5'
elif season == 'spring':
    climato_file = 'Tbot_climato_SA4_spring_0.10.h5'
elif season == 'summer':
    climato_file = 'Tbot_climato_SA4_summer_0.10.h5'
h5f = h5py.File(climato_file, 'r')
Tbot_climato = h5f['Tbot'][:]
lon_reg = h5f['lon_reg'][:]
lat_reg = h5f['lat_reg'][:]
Zitp = h5f['Zitp'][:]
h5f.close()

# Derive some map parameters
lon_0 = np.round(np.mean(lon_reg))
lat_0 = np.round(np.mean(lat_reg))
lonLims = [lon_reg[0], lon_reg[-1]]
latLims = [lat_reg[0], lat_reg[-1]]

# NAFO divisions
nafo_div = azu.get_nafo_divisions()
shape_4R = Polygon(zip(nafo_div['4R']['lon'], nafo_div['4R']['lat']))
polygon4Vs = Polygon(zip(nafo_div['4Vs']['lon'], nafo_div['4Vs']['lat']))
polygon4Vn = Polygon(zip(nafo_div['4Vn']['lon'], nafo_div['4Vn']['lat']))
polygon4W = Polygon(zip(nafo_div['4W']['lon'], nafo_div['4W']['lat']))
polygon4X = Polygon(zip(nafo_div['4X']['lon'], nafo_div['4X']['lat']))
shape = [polygon4Vs.buffer(0), polygon4Vn.buffer(0), polygon4W.buffer(0), polygon4X.buffer(0)]
shape_4VWX = cascaded_union(shape)

dict_stats_4R = {}
dict_stats_4VWX = {}


# Loop on years
df_list = []
for year in years:
    print ' ---- ' + np.str(year) + ' ---- '
    year_file = '/home/cyrf0006/data/dev_database/' + np.str(year) + '.nc'
    Tdict = azu.get_bottomT(year_file, season, climato_file)    
    Tbot = Tdict['Tbot']
    lons = Tdict['lons']
    lats = Tdict['lats']
    anom = Tbot-Tbot_climato

    
    # NAFO division stats    
    dict_stats_4R[np.str(year)] = azu.polygon_temperature_stats(Tdict, shape_4R)
    dict_stats_4VWX[np.str(year)] = azu.polygon_temperature_stats(Tdict, shape_4VWX)

    # Append bottom temperature for multi-index export
    df = pd.DataFrame(index=lat_reg, columns=lon_reg)
    df.index.name='latitude'
    df.columns.name='longitude'
    df[:] = Tbot
    df_list.append(df)

    
    if plot:
        # 1.1 - Plot Anomaly
        fig, ax = plt.subplots(nrows=1, ncols=1)
        m = Basemap(ax=ax, projection='merc',lon_0=lon_0,lat_0=lat_0, llcrnrlon=lonLims[0],llcrnrlat=latLims[0],urcrnrlon=lonLims[1],urcrnrlat=latLims[1], resolution= 'i')
        levels = np.linspace(-3.5, 3.5, 8)
        xi, yi = m(*np.meshgrid(lon_reg, lat_reg))
        c = m.contourf(xi, yi, anom, levels, cmap=plt.cm.RdBu_r, extend='both')
        cc = m.contour(xi, yi, -Zitp, [100, 500, 1000, 4000], colors='grey');
        plt.clabel(cc, inline=1, fontsize=10, fmt='%d')
        if season=='fall':
            plt.title('Fall Bottom Temperature Anomaly')
        elif season=='spring':
            plt.title('Spring Bottom Temperature Anomaly')
        else:
            plt.title('Bottom Temperature Anomaly')
        m.fillcontinents(color='tan');
        m.drawparallels([40, 45, 50, 55, 60], labels=[1,0,0,0], fontsize=12, fontweight='normal');
        m.drawmeridians([-60, -55, -50, -45], labels=[0,0,0,1], fontsize=12, fontweight='normal');
        cax = plt.axes([0.85,0.15,0.04,0.7], facecolor='grey')
        cb = plt.colorbar(c, cax=cax)
        cb.set_label(r'$\rm T(^{\circ}C)$', fontsize=12, fontweight='normal')
        div_toplot = ['4R', '4Vn', '4Vs', '4W', '4X']
        for div in div_toplot:
            div_lon, div_lat = m(nafo_div[div]['lon'], nafo_div[div]['lat'])
            m.plot(div_lon, div_lat, 'k', linewidth=2)
            ax.text(np.mean(div_lon), np.mean(div_lat), div, fontsize=12, color='black', fontweight='bold')    
        # Save Figure
        fig.set_size_inches(w=7, h=8)
        fig.set_dpi(200)
        outfile = 'bottom_temp_anomaly_SA4_' + season + '_' + np.str(year) + '.png'
        fig.savefig(outfile)
        plt.close()
        fig.close()
        
        # 1.2 - Plot Temperature
        fig, ax = plt.subplots(nrows=1, ncols=1)
        m = Basemap(ax=ax, projection='merc',lon_0=lon_0,lat_0=lat_0, llcrnrlon=lonLims[0],llcrnrlat=latLims[0],urcrnrlon=lonLims[1],urcrnrlat=latLims[1], resolution= 'i')
        levels = np.linspace(-2, 6, 9)
        xi, yi = m(*np.meshgrid(lon_reg, lat_reg))
        c = m.contourf(xi, yi, Tbot, levels, cmap=plt.cm.RdBu_r, extend='both')
        cc = m.contour(xi, yi, -Zitp, [100, 500, 1000, 4000], colors='grey');
        plt.clabel(cc, inline=1, fontsize=10, fmt='%d')
        if season=='fall':
            plt.title('Fall Bottom Temperature')
        elif season=='spring':
            plt.title('Spring Bottom Temperature')
        else:
            plt.title('Bottom Temperature')
        m.fillcontinents(color='tan');
        m.drawparallels([40, 45, 50, 55, 60], labels=[1,0,0,0], fontsize=12, fontweight='normal');
        m.drawmeridians([-60, -55, -50, -45], labels=[0,0,0,1], fontsize=12, fontweight='normal');
        x, y = m(lons, lats)
        m.scatter(x,y, s=50, marker='.',color='k')
        cax = plt.axes([0.85,0.15,0.04,0.7], facecolor='grey')
        cb = plt.colorbar(c, cax=cax)
        cb.set_label(r'$\rm T(^{\circ}C)$', fontsize=12, fontweight='normal')
        div_toplot = ['4R', '4Vn', '4Vs', '4W', '4X']
        for div in div_toplot:
            div_lon, div_lat = m(nafo_div[div]['lon'], nafo_div[div]['lat'])
            m.plot(div_lon, div_lat, 'k', linewidth=2)
            ax.text(np.mean(div_lon), np.mean(div_lat), div, fontsize=12, color='black', fontweight='bold')
        # Save Figure
        fig.set_size_inches(w=7, h=8)
        fig.set_dpi(200)
        outfile = 'bottom_temp_SA4_' + season + '_' + np.str(year) + '.png'
        fig.savefig(outfile)

df_4R = pd.DataFrame.from_dict(dict_stats_4R, orient='index')
df_4VWX = pd.DataFrame.from_dict(dict_stats_4VWX, orient='index')
outname = 'stats_4R_' + season + '.pkl'
df_4R.to_pickle(outname)
outname = 'stats_4VWX_' + season + '.pkl'
df_4VWX.to_pickle(outname)

