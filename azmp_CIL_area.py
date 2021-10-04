'''
CIL area

This script is still in progress...


'''
import os
import netCDF4
import xarray as xr
os.environ['PROJ_LIB'] = '/home/cyrf0006/anaconda3/share/proj'
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import pandas as pd
import numpy as  np
from scipy.interpolate import interp1d  # to remove NaNs in profiles
from scipy.interpolate import griddata
import azmp_sections_tools as azst


## ---- Region parameters ---- ## <-------------------------------Would be nice to pass this in a config file '2017.report'
SECTION = 'FC'
SEASON = 'summer'
YEAR = 2019
dlat = 2 # how far from station we search
dlon = 2
dz = 1 # vertical bins
dc = .1 # grid resolution

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
df_stn = pd.read_excel('/home/cyrf0006/github/AZMP-NL/data/STANDARD_SECTIONS.xlsx')
df_stn = df_stn.drop(['SECTION', 'LONG'], axis=1)
df_stn = df_stn.rename(columns={'LONG.1': 'LON'})
df_stn = df_stn.dropna()
df_stn = df_stn[df_stn.STATION.str.contains(SECTION)]
df_stn = df_stn.reset_index(drop=True)
stn_list = df_stn.STATION.values

# Derive regular grid
latLims = np.array([df_stn.LAT.min() - dlat, df_stn.LAT.max() + dlat])
lonLims = np.array([df_stn.LON.min() - dlon, df_stn.LON.max() + dlon])
lon_reg = np.arange(lonLims[0]+dc/2, lonLims[1]-dc/2, dc)
lat_reg = np.arange(latLims[0]+dc/2, latLims[1]-dc/2, dc)

## --------- Get Bathymetry -------- ####
print('Get bathy...')
dataFile = '/home/cyrf0006/data/GEBCO/GEBCO_2014_1D.nc' # Maybe find a better way to handle this file
lon_grid, lat_grid = np.meshgrid(lon_reg,lat_reg)
# Load data
dataset = netCDF4.Dataset(dataFile)
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
print(' -> Done!')


## -------- Get CTD data -------- ##
year_file = '/home/cyrf0006/data/dev_database/netCDF/' + str(YEAR) + '.nc'
print('Get ' + year_file)
ds = xr.open_mfdataset(year_file)

# Remame problematic datasets
print('!!Remove MEDBA & MEDTE data!!')
print('  ---> I Should be improme because I remove good data!!!!')
ds = ds.where(ds.instrument_ID!='MEDBA', drop=True)
ds = ds.where(ds.instrument_ID!='MEDTE', drop=True)

# Select Region
ds = ds.where((ds.longitude>lonLims[0]) & (ds.longitude<lonLims[1]), drop=True)
ds = ds.where((ds.latitude>latLims[0]) & (ds.latitude<latLims[1]), drop=True)

# Select time (save several options here)
if SEASON == 'summer':
    ds = ds.sel(time=((ds['time.month']>=7)) & ((ds['time.month']<=9)))
elif SEASON == 'spring':
    ds = ds.sel(time=((ds['time.month']>=3)) & ((ds['time.month']<=5)))
elif SEASON == 'fall':
    ds = ds.sel(time=((ds['time.month']>=10)) & ((ds['time.month']<=12)))
else:
    print('!! no season specified, used them all! !!')

# Extract temperature    
da_temp = ds['temperature']
lons = np.array(ds.longitude)
lats = np.array(ds.latitude)
bins = np.arange(dz/2.0, 500, dz)
da_temp = da_temp.groupby_bins('level', bins).mean(dim='level')
#To Pandas Dataframe
df = da_temp.to_pandas()
df.columns = bins[0:-1] #rename columns with 'bins'
# Remove empty columns
idx_empty_rows = df.isnull().all(1).nonzero()[0]
df = df.dropna(axis=0,how='all')
lons = np.delete(lons,idx_empty_rows)
lats = np.delete(lats,idx_empty_rows)
print(' -> Done!')


## --- fill 3D cube --- ##  
print('Fill regular cube')
z = df.columns.values
V = np.full((lat_reg.size, lon_reg.size, z.size), np.nan)

# Aggregate on regular grid
for i, xx in enumerate(lon_reg):
    for j, yy in enumerate(lat_reg):    
        idx = np.where((lons>=xx-dc/2) & (lons<xx+dc/2) & (lats>=yy-dc/2) & (lats<yy+dc/2))
        tmp = np.array(df.iloc[idx].mean(axis=0))
        idx_good = np.argwhere((~np.isnan(tmp)) & (tmp<30))
        if np.size(idx_good)==1:
            V[j,i,:] = np.array(df.iloc[idx].mean(axis=0))
        elif np.size(idx_good)>1: # vertical interpolation between pts
            interp = interp1d(np.squeeze(z[idx_good]), np.squeeze(tmp[idx_good]))  # <---- Attention: strange syntax but works
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
    if idx_good.size: # will ignore depth where no data exist
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

            
## ---- Check plot (uncomment to check the region) ---- ##
fig, ax = plt.subplots(nrows=1, ncols=1)
m = Basemap(ax=ax, projection='merc',lon_0=-50,lat_0=50, llcrnrlon=-60,llcrnrlat=45,urcrnrlon=-40,urcrnrlat=55, resolution= 'i')
levels = np.linspace(-2, 6, 17)
xi, yi = m(*np.meshgrid(lon_reg, lat_reg))
c = m.contourf(xi, yi, np.squeeze(V[:,:,10]), levels, cmap=plt.cm.RdBu, extend='both')
cc = m.contour(xi, yi, -Zitp, [100, 500, 1000, 4000], colors='grey');
plt.clabel(cc, inline=1, fontsize=10, fmt='%d')
m.fillcontinents(color='tan');
m.drawparallels([40, 45, 50, 55, 60], labels=[0,0,0,0], fontsize=12, fontweight='normal');
m.drawmeridians([-60, -55, -50, -45], labels=[0,0,0,1], fontsize=12, fontweight='normal');
cax = fig.add_axes([0.16, 0.05, 0.7, 0.025])
cb = plt.colorbar(c, cax=cax, orientation='horizontal')
cb.set_label(r'$\rm T(^{\circ}C)$', fontsize=12, fontweight='normal')
x, y = m(lons, lats)
m.scatter(x,y, s=50, marker='.',color='k')
plt.show()



## ---- Extract section info (test 2 options) ---- ##
temp_coords = np.array([[x,y] for x in lat_reg for y in lon_reg])
VV = V.reshape(temp_coords.shape[0],V.shape[2])
ZZ = Zitp.reshape(temp_coords.shape[0],1)
section_only = []
df_section_itp = pd.DataFrame(index=stn_list, columns=z)
for stn in stn_list:
    # 1. Section only (by station name)
    ds_tmp = ds.where(ds.comments == stn, drop=True)  
    section_only.append(ds_tmp)

    #2.  From interpolated field (closest to station)
    station = df_stn[df_stn.STATION==stn]
    idx_opti = np.argmin(np.sum(np.abs(temp_coords - np.array(list(zip(station.LAT,station.LON)))), axis=1))
    Tprofile = VV[idx_opti,:]
    # remove data below bottom
    bottom_depth = -ZZ[idx_opti]
    Tprofile[z>=bottom_depth]=np.nan
    # store in dataframe
    df_section_itp.loc[stn] = Tprofile

# convert option #1 to dataframe    
ds_section = xr.concat(section_only, dim='time')
da = ds_section['temperature']
df_section_stn = da.to_pandas()
df_section_stn.index = ds_section.comments.values

# Compute distance vector for option #1 - exact station
distance_stn = np.full((df_section_stn.index.shape), np.nan)
for i, stn in enumerate(df_section_stn.index):
    distance_stn[i] = azst.haversine(df_stn.LON[0], df_stn.LAT[0], df_stn[df_stn.STATION==df_section_stn.index[i]].LON, df_stn[df_stn.STATION==df_section_stn.index[i]].LAT)
# Compute distance vector for option #2 - interp field
distance_itp = np.full((df_section_itp.index.shape), np.nan)
for i, stn in enumerate(df_section_itp.index):
    distance_itp[i] = azst.haversine(df_stn.LON[0], df_stn.LAT[0], df_stn[df_stn.STATION==df_section_itp.index[i]].LON, df_stn[df_stn.STATION==df_section_itp.index[i]].LAT)

## ---- Plot to check the result ---- ##
v = np.arange(-2,13,1)
## ---- Plot to check the result ---- ##
XLIM = azst.haversine(df_stn.LON[0], df_stn.LAT[0],df_stn.iloc[-1].LON,df_stn.iloc[-1].LAT)
if df_section_stn.index.size > 0:
    fig, ax = plt.subplots()
    c = plt.contourf(distance_stn, df_section_stn.columns, df_section_stn.T, v)
    c_cil_stn = plt.contour(distance_stn, df_section_stn.columns, df_section_stn.T, [0,], colors='k', linewidths=2)
    ax.set_ylim([0, 400])
    ax.set_xlim([0,  XLIM])
    ax.set_ylabel('Depth (m)', fontWeight = 'bold')
    ax.set_xlabel('Distance (km)')
    ax.invert_yaxis()
    plt.colorbar(c)
    fig.savefig('method1.png', dpi=150)
    plt.close()
    # CIL area
    cil_vol_stn = 0
    CIL = c_cil_stn.collections[0]
    for path in CIL.get_paths()[:]:
        vs = path.vertices
        cil_vol_stn = cil_vol_stn + np.abs(area(vs))/1000
else:
    cil_vol_stn = np.nan

if df_section_itp.index.size > 0:
    fig, ax = plt.subplots()
    c = plt.contourf(distance_itp, df_section_itp.columns, df_section_itp.T, v)
    c_cil_itp = plt.contour(distance_itp, df_section_itp.columns, df_section_itp.T, [0,], colors='k', linewidths=2)
    ax.set_ylim([0, 400])
    ax.set_xlim([0,  XLIM])
    ax.set_ylabel('Depth (m)', fontWeight = 'bold')
    ax.set_xlabel('Distance (km)')
    ax.invert_yaxis()
    plt.colorbar(c)
    fig.savefig('method2.png', dpi=150)
    plt.close()
    # CIL area
    cil_vol_itp = 0
    CIL = c_cil_itp.collections[0]
    for path in CIL.get_paths()[:]:
        vs = path.vertices
        cil_vol_itp = cil_vol_itp + np.abs(area(vs))/1000
else:
    cil_vol_itp = np.nan

##         fig, ax = plt.subplots()
## c = plt.contourf(distance_stn, df_section_stn.columns, df_section_stn.T, v)
## c_cil_stn = plt.contour(distance_stn, df_section_stn.columns, df_section_stn.T, [0,], colors='k', linewidths=2)
## ax.set_ylim([0, 400])
## ax.set_xlim([0,  distance_stn[-1]])
## ax.set_ylabel('Depth (m)', fontWeight = 'bold')
## ax.set_xlabel('Distance (km)')
## ax.invert_yaxis()
## plt.colorbar(c)
## fig.savefig('method1.png', dpi=150)

## fig, ax = plt.subplots()
## c = plt.contourf(distance_itp, df_section_itp.columns, df_section_itp.T, v)
## c_cil_itp = plt.contour(distance_itp, df_section_itp.columns, df_section_itp.T, [0,], colors='k', linewidths=2)
## ax.set_ylim([0, 400])
## ax.set_xlim([0,  distance_stn[-1]])
## ax.set_ylabel('Depth (m)', fontWeight = 'bold')
## ax.set_xlabel('Distance (km)')
## ax.invert_yaxis()
## plt.colorbar(c)
## fig.savefig('method2.png', dpi=150)


# CIL surface (Note that there is a bias because )
## def area(vs):
##     a = 0
##     x0,y0 = vs[0]
##     for [x1,y1] in vs[1:]:
##         dx = x1-x0
##         dy = y1-y0
##         a += 0.5*(y0*dx - x0*dy)
##         x0 = x1
##         y0 = y1
##     return a

## cil_vol_stn = 0`
## CIL = c_cil_stn.collections[0]
## for path in CIL.get_paths()[:]:
##     vs = path.vertices
##     cil_vol_stn = cil_vol_stn + np.abs(area(vs))/1000
    
## cil_vol_itp = 0
## CIL = c_cil_itp.collections[0]
## for path in CIL.get_paths()[:]:
##     vs = path.vertices
##     cil_vol_itp = cil_vol_itp + np.abs(area(vs))/1000

print(cil_vol_stn, cil_vol_itp)
