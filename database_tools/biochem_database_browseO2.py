# A first test to read Excel nutrient file and export to Pandas.

# Check in:
#  /home/cyrf0006/research/AZMP_database/biochem

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import datetime
import gsw
from seawater import extras as swx
import matplotlib.dates as mdates
from scipy.interpolate import griddata

# Adjust fontsize/weight
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 18}
plt.rc('font', **font)

## ----  Load data from Excel sheet ---- ##
df = pd.read_excel('~/github/AZMP-NL/data/AZMP_Nutrients_1999_2016_Good_Flags_Only.xlsx')

## ---- Some cleaning ---- ##
# Set date as index
df = df.set_index('sas_date')

# Drop other time-related columns
df = df.drop(['Day', 'Month', 'Year'], axis=1)

# Compute Saturation O2 and add to dataframe
df['satO2'] = swx.satO2(df['salinity'], df['temp'])
df['satO2%'] = df['oxygen']/df['satO2']*100

# Keep BB line only
df_BB = df[df['section']=='BB']

# Keep only some variables
df_BB_O2 = df_BB[['depth','oxygen', 'Longitude', 'satO2']]

# Drop NaNs
#df_BB_O2 = df_BB_O2.dropna(axis=0, how='any')

# depth range
idx = ((df_BB_O2['depth']>0) & (df_BB_O2['depth']<=125))
df_BB_O2_surf = df_BB_O2[idx]

idx = ((df_BB_O2['depth']>=200) & (df_BB_O2['depth']<=350))
df_BB_O2_btm = df_BB_O2[idx]

# longitude range
idx = ((df_BB_O2_btm['Longitude']<-50))
df_BB_O2_btm = df_BB_O2_btm[idx]
idx = ((df_BB_O2_surf['Longitude']<-50))
df_BB_O2_surf = df_BB_O2_surf[idx]

# resample
df_season_surf = df_BB_O2_surf.resample('4M').mean()
df_season_btm = df_BB_O2_btm.resample('4M').mean()


## ----  Build climatology ---- ##
df_O2clim = df_BB.groupby(['depth', 'sname'])['satO2%'].mean() # <---------- THIS IS SICK!!
df_coords = df_BB.groupby(['sname'])['Latitude', 'Longitude'].mean()

# Compute along-transect distance
from math import radians, cos, sin, asin, sqrt
from math import radians, cos, sin, asin, sqrt
def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    km = 6367 * c
    return km

distance = np.zeros(df_coords['Latitude'].shape)
for i in range(df_coords['Latitude'].size):
    distance[i] = haversine(df_coords['Longitude'][0], df_coords['Latitude'][0], df_coords['Longitude'][i], df_coords['Latitude'][i])

# Load bathymetry
bb_bathy = '/home/cyrf0006/github/AZMP-NL/bathymetry/bottom_profiles/bbline.txt'
from numpy import loadtxt
bathy = loadtxt(bb_bathy, delimiter=",", unpack=False)
bathy_x = bathy[:,0]/1000.0
bathy_y = np.abs(bathy[:,1])

bathy_x_close = np.append(bathy_x, np.array([bathy_x[-1], bathy_x[0], bathy_x[0]]))
bathy_y_close = np.append(bathy_y ,np.array([bathy_y.max(), bathy_y.max(), bathy_y[0]]))
bathymetry = zip(bathy_x_close, bathy_y_close)

# Check maximum depth
cast_depth = []
for i in distance:
    min_idx = np.argmin(np.abs(i-bathy_x))
    cast_depth = np.append(cast_depth, bathy_y[min_idx])

   
# Digitize depth
#Pbin = np.arange(2.5, 500, 5)
Pbin = np.array([5,10,20,30,40,50,75,100,150,200,250,500])
O2u = df_O2clim.unstack(level=1) # rearange in matrix (depth vs station)
O2u = O2u[O2u.index<Pbin.max()] # Drop deeper samples
O2list = [];
for stn in O2u.columns:
    O2_vec = O2u[stn].dropna(axis=0, how='any')

    digitized = np.digitize(O2_vec.index, Pbin) #<- this is awesome!
    O2list.append([O2_vec.values[digitized == i].mean() for i in range(0, len(Pbin))])

# To array + interpolation    
O2array = np.transpose(np.array(O2list))
X, Y = np.meshgrid(distance, Pbin)
yi = np.arange(2.5, 500, 5)
xi = np.arange(0,np.max(distance), 10)
Xi, Yi = np.meshgrid(xi,yi)
# back to vector
xvalues = X.ravel()
yvalues = Y.ravel()
zvalues = O2array.ravel()

x=list(xvalues[np.isnan(zvalues)==False]) # remove nans
y=list(yvalues[np.isnan(zvalues)==False])
z=list(zvalues[np.isnan(zvalues)==False])

O2i = griddata((x,y), z, (Xi, Yi), method='linear')


# Plot
fig, ax = plt.subplots()
#ctf = plt.contourf(distance, Pbin, O2array, 30, cmap=plt.cm.RdBu)
ctf = plt.contourf(Xi, Yi, O2i, 30, cmap=plt.cm.RdBu)
cl = plt.colorbar(ctf, orientation='vertical')
ax.set_ylim([0, 500])
ax.set_ylabel('Depth (m)')
ax.set_xlabel('Distance (km)')
ax.set_title(r'$\rm O_2 sat.(\%)$ - Bonavista Section (1999-2017)')
ax.invert_yaxis()
Bgon = plt.Polygon(bathymetry,color=np.multiply([1,.9333,.6667],.4), alpha=1)
ax.add_patch(Bgon)
for i in range(0,len(distance)):
    plt.plot(np.array([distance[i], distance[i]]), np.array([Pbin[0], cast_depth[i]]), '--k', linewidth=0.1)

plt.plot([1, 225, 225, 1, 1], [1,1,125,125,1], '--k')
plt.plot([1, 225, 225, 1, 1], [200,200,350,350,200], '--k')
#ax.text(350, 25, r'$\rm O_2 sat.(\%)$', horizontalalignment='left', verticalalignment='center', fontsize=16, color='k')

#plt.show()




## import matplotlib.mlab as mlab
## from scipy.interpolate import griddata
## def normalize_x(data):
##     data = data.astype(np.float)
##     return (data - xmin) / (xmax - xmin)

## def normalize_y(data):
##     data = data.astype(np.float)
##     return (data - ymin) / (ymax - ymin)

## lon_array = np.array(df_BB_O2['Longitude'])
## depth_array = np.array(df_BB_O2['depth'])
## O2_array = np.array(df_BB_O2['oxygen'])

## xmin = -53
## xmax = -49
## ymin = 0
## ymax = 500

## xi = np.arange(xmin, xmax, .2)
## yi = np.arange(ymin, ymax, 25)

## x_new, xi_new = normalize_x(lon_array), normalize_x(xi)
## y_new, yi_new = normalize_y(depth_array), normalize_y(yi)

## zi = mlab.griddata(x_new, y_new, O2_array, xi_new, yi_new, interp='linear' )
## O2i = griddata((x_new,y_new), O2_array, (xi_new, yi_new), method='linear')



## ---------- Plot surface and bottom with the avarage and STD ----------- ##

surf = df_season_surf['oxygen']/df_season_surf['satO2']*100
btm = df_season_btm['oxygen']/df_season_btm['satO2']*100
x0 = mdates.date2num(df_season_surf.index[0])
xf = mdates.date2num(df_season_surf.index[-1])
xVec = [x0, xf, xf, x0, x0]
x = 0.5
yVec_surf = [surf.mean()+x*surf.std(), surf.mean()+x*surf.std(), surf.mean()-x*surf.std(), surf.mean()-x*surf.std(), surf.mean()+x*surf.std()]
yVec_btm = [btm.mean()+x*btm.std(), btm.mean()+x*btm.std(), btm.mean()-x*btm.std(), btm.mean()-x*btm.std(), btm.mean()+x*btm.std()]
rect_surf = zip(xVec, yVec_surf)
rect_btm = zip(xVec, yVec_btm)

fig = plt.figure(1)
ax = fig.add_subplot(111)

plt.plot(df_season_surf.index, surf, '.-', color='steelblue')
plt.plot(df_season_btm.index, btm, '.-', color='orange')

Rgon1 = plt.Polygon(rect_surf,color='steelblue', alpha=0.2)
ax.add_patch(Rgon1)
Rgon2 = plt.Polygon(rect_btm,color='orange', alpha=0.2)
ax.add_patch(Rgon2)
plt.plot([x0, xf], [surf.mean(), surf.mean()], '--', color='steelblue', alpha=0.5)
plt.plot([x0, xf], [btm.mean(), btm.mean()], '--', color='orange', alpha=0.5)

plt.plot(df_season_surf.index, surf, '.-', color='steelblue')
plt.plot(df_season_btm.index, btm, '.-', color='orange')

plt.xlabel('Year')
plt.ylabel(r'$\rm O_2 sat (\%)$')
plt.title('BB shelf')
plt.legend(['0-125m', '200-350m'])
plt.grid('on')
plt.show()

## Plot climatological section
fig = plt.figure(3)
ax = fig.add_subplot(111)
plt.contourf(X, Y, O2i)
plt.scatter(lon_array, depth_array, c='k', alpha=0.2, marker='.')
plt.title(r'Climatology $\rm O_2$ Bonavista Section')
plt.xlabel('Longitude')
plt.ylabel('Depth (m)')
plt.colorbar()


keyboard
