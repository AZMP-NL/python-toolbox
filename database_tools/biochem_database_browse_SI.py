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

# Keep SI line only
df_SI = df[df['section']=='SI']

# Keep only one season <------------------------------------------------ ONE SEASON ONLY!!!
#df_SI = df_SI[(df_SI.index.month==6) | (df_SI.index.month==7) |  (df_SI.index.month==8)] 


# longitude range
idx = ((df_SI['Longitude']<-50))
df_SI = df_SI[idx]

# Keep only some variables
df_SI_O2 = df_SI[['depth','oxygen', 'Longitude', 'satO2']]
df_SI_nut = df_SI[['depth','PO4', 'NO3', 'SIO', 'Longitude']]
df_SI_nut.to_pickle('SI_nutrients')

# Drop NaNs
#df_SI_O2 = df_SI_O2.dropna(axis=0, how='any')


# depth range
idx = ((df_SI_O2['depth']>0) & (df_SI_O2['depth']<=125))
df_SI_O2_surf = df_SI_O2[idx]
idx = ((df_SI_nut['depth']>0) & (df_SI_nut['depth']<=125))
df_SI_nut_surf = df_SI_nut[idx]

idx = ((df_SI_O2['depth']>=200) & (df_SI_O2['depth']<=350))
df_SI_O2_btm = df_SI_O2[idx]
idx = ((df_SI_nut['depth']>200) & (df_SI_nut['depth']<=350))
df_SI_nut_btm = df_SI_nut[idx]


# resample
df_season_surf = df_SI_O2_surf.resample('AS').mean()
df_season_btm = df_SI_O2_btm.resample('AS').mean()


## ----  Build climatology ---- ##
df_O2clim = df_SI.groupby(['depth', 'sname'])['satO2%'].mean() # <---------- THIS IS SICK!!
df_NO3clim = df_SI.groupby(['depth', 'sname'])['NO3'].mean()
df_PO4clim = df_SI.groupby(['depth', 'sname'])['PO4'].mean()
df_SIOclim = df_SI.groupby(['depth', 'sname'])['SIO'].mean()
df_coords = df_SI.groupby(['sname'])['Latitude', 'Longitude'].mean()


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
si_bathy = '/home/cyrf0006/github/AZMP-NL/bathymetry/bottom_profiles/siline.txt'
from numpy import loadtxt
bathy = loadtxt(si_bathy, delimiter=",", unpack=False)
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

NO3u = df_NO3clim.unstack(level=1)
NO3u = NO3u[NO3u.index<Pbin.max()]

PO4u = df_PO4clim.unstack(level=1)
PO4u = PO4u[PO4u.index<Pbin.max()]

SIOu = df_SIOclim.unstack(level=1)
SIOu = SIOu[SIOu.index<Pbin.max()]

O2list = [];
NO3list = [];
PO4list = [];
SIOlist = [];
  
for stn in O2u.columns:
    O2_vec = O2u[stn].dropna(axis=0, how='any')
    NO3_vec = NO3u[stn].dropna(axis=0, how='any')
    PO4_vec = PO4u[stn].dropna(axis=0, how='any')
    SIO_vec = SIOu[stn].dropna(axis=0, how='any')

    
    digitized = np.digitize(O2_vec.index, Pbin) #<- this is awesome!
    O2list.append([O2_vec.values[digitized == i].mean() for i in range(0, len(Pbin))])

    digitized = np.digitize(NO3_vec.index, Pbin) #<- this is awesome!
    NO3list.append([NO3_vec.values[digitized == i].mean() for i in range(0, len(Pbin))])
    PO4list.append([PO4_vec.values[digitized == i].mean() for i in range(0, len(Pbin))])
    SIOlist.append([SIO_vec.values[digitized == i].mean() for i in range(0, len(Pbin))])

    
# To array + interpolation    
O2array = np.transpose(np.array(O2list))
NO3array = np.transpose(np.array(NO3list))
PO4array = np.transpose(np.array(PO4list))
SIOarray = np.transpose(np.array(SIOlist))
X, Y = np.meshgrid(distance, Pbin)
yi = np.arange(2.5, 500, 5)
xi = np.arange(0,np.max(distance), 10)
Xi, Yi = np.meshgrid(xi,yi)
# back to vector
xvalues = X.ravel()
yvalues = Y.ravel()

# interp O2
zvalues = O2array.ravel()
x=list(xvalues[np.isnan(zvalues)==False]) # remove nans
y=list(yvalues[np.isnan(zvalues)==False])
z=list(zvalues[np.isnan(zvalues)==False])
O2i = griddata((x,y), z, (Xi, Yi), method='linear')

# interp NO3
zvalues = NO3array.ravel()
x=list(xvalues[np.isnan(zvalues)==False]) # remove nans
y=list(yvalues[np.isnan(zvalues)==False])
z=list(zvalues[np.isnan(zvalues)==False])
NO3i = griddata((x,y), z, (Xi, Yi), method='linear')

# interp PO4
zvalues = PO4array.ravel()
x=list(xvalues[np.isnan(zvalues)==False]) # remove nans
y=list(yvalues[np.isnan(zvalues)==False])
z=list(zvalues[np.isnan(zvalues)==False])
PO4i = griddata((x,y), z, (Xi, Yi), method='linear')

# interp SIO
zvalues = SIOarray.ravel()
x=list(xvalues[np.isnan(zvalues)==False]) # remove nans
y=list(yvalues[np.isnan(zvalues)==False])
z=list(zvalues[np.isnan(zvalues)==False])
SIOi = griddata((x,y), z, (Xi, Yi), method='linear')

### ----  Plot timeseries ---- ###
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

lon_array = np.array(df_SI_O2['Longitude'])
depth_array = np.array(df_SI_O2['depth'])

## ---------- Plot surface and bottom with the avarage and STD ----------- ##

# O2
df_season_surf = df_SI_O2_surf.resample('4M').mean()
df_season_btm = df_SI_O2_btm.resample('4M').mean()
#df_season_surf = df_SI_O2_surf.resample('A').mean()
#df_season_btm = df_SI_O2_btm.resample('A').mean()


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
fig.clf()
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
plt.title('SI shelf')
plt.legend(['0-125m', '200-350m'])
plt.grid('on')
fig.set_size_inches(w=10,h=6)
fig_name = 'SI_O2_timeseries.png'
fig.set_dpi(300)
fig.savefig(fig_name)


# NO3
df_season_surf = df_SI_nut_surf.resample('4M').mean()
df_season_btm = df_SI_nut_btm.resample('4M').mean()
#df_season_surf = df_SI_nut_surf.resample('A').mean()
#df_season_btm = df_SI_nut_btm.resample('A').mean()
surf = df_season_surf['NO3']
btm = df_season_btm['NO3']
x0 = mdates.date2num(df_season_surf.index[0])
xf = mdates.date2num(df_season_surf.index[-1])
xVec = [x0, xf, xf, x0, x0]
x = 0.5
yVec_surf = [surf.mean()+x*surf.std(), surf.mean()+x*surf.std(), surf.mean()-x*surf.std(), surf.mean()-x*surf.std(), surf.mean()+x*surf.std()]
yVec_btm = [btm.mean()+x*btm.std(), btm.mean()+x*btm.std(), btm.mean()-x*btm.std(), btm.mean()-x*btm.std(), btm.mean()+x*btm.std()]
rect_surf = zip(xVec, yVec_surf)
rect_btm = zip(xVec, yVec_btm)

fig = plt.figure(1)
fig.clf()
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
plt.ylabel(r'$\rm NO_3 (mmol\,Kg^{-1})$')
plt.title('SI shelf')
plt.legend(['0-125m', '200-350m'])
plt.grid('on')
fig.set_size_inches(w=10,h=6)
fig_name = 'SI_NO3_timeseries.png'
fig.set_dpi(300)
fig.savefig(fig_name)

# PO4
surf = df_season_surf['PO4']
btm = df_season_btm['PO4']
x0 = mdates.date2num(df_season_surf.index[0])
xf = mdates.date2num(df_season_surf.index[-1])
xVec = [x0, xf, xf, x0, x0]
x = 0.5
yVec_surf = [surf.mean()+x*surf.std(), surf.mean()+x*surf.std(), surf.mean()-x*surf.std(), surf.mean()-x*surf.std(), surf.mean()+x*surf.std()]
yVec_btm = [btm.mean()+x*btm.std(), btm.mean()+x*btm.std(), btm.mean()-x*btm.std(), btm.mean()-x*btm.std(), btm.mean()+x*btm.std()]
rect_surf = zip(xVec, yVec_surf)
rect_btm = zip(xVec, yVec_btm)

fig = plt.figure(1)
fig.clf()
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
plt.ylabel(r'$\rm PO_4 (mmol\,Kg^{-1})$')
plt.title('SI shelf')
plt.legend(['0-125m', '200-350m'])
plt.grid('on')
fig.set_size_inches(w=10,h=6)
fig_name = 'SI_PO4_timeseries.png'
fig.set_dpi(300)
fig.savefig(fig_name)

# SIO
surf = df_season_surf['SIO']
btm = df_season_btm['SIO']
x0 = mdates.date2num(df_season_surf.index[0])
xf = mdates.date2num(df_season_surf.index[-1])
xVec = [x0, xf, xf, x0, x0]
x = 0.5
yVec_surf = [surf.mean()+x*surf.std(), surf.mean()+x*surf.std(), surf.mean()-x*surf.std(), surf.mean()-x*surf.std(), surf.mean()+x*surf.std()]
yVec_btm = [btm.mean()+x*btm.std(), btm.mean()+x*btm.std(), btm.mean()-x*btm.std(), btm.mean()-x*btm.std(), btm.mean()+x*btm.std()]
rect_surf = zip(xVec, yVec_surf)
rect_btm = zip(xVec, yVec_btm)

fig = plt.figure(1)
fig.clf()
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
plt.ylabel(r'$\rm SiO (mmol\,Kg^{-1})$')
plt.title('SI shelf')
plt.legend(['0-125m', '200-350m'])
plt.grid('on')
fig.set_size_inches(w=10,h=6)
fig_name = 'SI_SIO_timeseries.png'
fig.set_dpi(300)
fig.savefig(fig_name)

# NO3/PO4
surf = df_season_surf['NO3']/df_season_surf['PO4']
btm = df_season_btm['NO3']/df_season_btm['PO4']
x0 = mdates.date2num(df_season_surf.index[0])
xf = mdates.date2num(df_season_surf.index[-1])
xVec = [x0, xf, xf, x0, x0]
x = 0.5
yVec_surf = [surf.mean()+x*surf.std(), surf.mean()+x*surf.std(), surf.mean()-x*surf.std(), surf.mean()-x*surf.std(), surf.mean()+x*surf.std()]
yVec_btm = [btm.mean()+x*btm.std(), btm.mean()+x*btm.std(), btm.mean()-x*btm.std(), btm.mean()-x*btm.std(), btm.mean()+x*btm.std()]
rect_surf = zip(xVec, yVec_surf)
rect_btm = zip(xVec, yVec_btm)

fig = plt.figure(1)
fig.clf()
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
plt.ylabel(r'$\rm [NO_3]/[PO_4]$')
plt.title('SI shelf')
plt.legend(['0-125m', '200-350m'])
plt.ylim([0, 20])
plt.grid('on')
fig.set_size_inches(w=10,h=6)
fig_name = 'SI_NO3onPO4_timeseries.png'
fig.set_dpi(300)
fig.savefig(fig_name)

# NO3/SiO4
surf = df_season_surf['NO3']/df_season_surf['SIO']
btm = df_season_btm['NO3']/df_season_btm['SIO']
x0 = mdates.date2num(df_season_surf.index[0])
xf = mdates.date2num(df_season_surf.index[-1])
xVec = [x0, xf, xf, x0, x0]
x = 0.5
yVec_surf = [surf.mean()+x*surf.std(), surf.mean()+x*surf.std(), surf.mean()-x*surf.std(), surf.mean()-x*surf.std(), surf.mean()+x*surf.std()]
yVec_btm = [btm.mean()+x*btm.std(), btm.mean()+x*btm.std(), btm.mean()-x*btm.std(), btm.mean()-x*btm.std(), btm.mean()+x*btm.std()]
rect_surf = zip(xVec, yVec_surf)
rect_btm = zip(xVec, yVec_btm)

fig = plt.figure(1)
fig.clf()
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
plt.ylabel(r'$\rm [NO_3]/[SiO_4]$')
plt.title('SI shelf')
plt.legend(['0-125m', '200-350m'])
plt.ylim([0.5, 1.8])
plt.grid('on')
fig.set_size_inches(w=10,h=6)
fig_name = 'SI_NO3onSiO4_timeseries.png'
fig.set_dpi(300)
fig.savefig(fig_name)



## Plot climatological section
# O2 sat%
fig, ax = plt.subplots()
#fig.clf()
ctf = plt.contourf(Xi, Yi, O2i, 30, cmap=plt.cm.RdBu)
cl = plt.colorbar(ctf, orientation='vertical')
ax.set_ylim([0, 500])
ax.set_xlim([0, 300])
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
fig.set_size_inches(w=9,h=6)
fig_name = 'SI_O2.png'
fig.set_dpi(300)
fig.savefig(fig_name)


# NO3
fig, ax = plt.subplots()
#fig.clf()
ctf = plt.contourf(Xi, Yi, NO3i, 30, cmap=plt.cm.RdBu_r)
cl = plt.colorbar(ctf, orientation='vertical')
ax.set_ylim([0, 500])
ax.set_xlim([0, 300])
ax.set_ylabel('Depth (m)')
ax.set_xlabel('Distance (km)')
ax.set_title(r'$\rm NO_3 (mmol\,Kg^{-1})$ - Bonavista Section (1999-2017)')
ax.invert_yaxis()
Bgon = plt.Polygon(bathymetry,color=np.multiply([1,.9333,.6667],.4), alpha=1)
ax.add_patch(Bgon)
for i in range(0,len(distance)):
    plt.plot(np.array([distance[i], distance[i]]), np.array([Pbin[0], cast_depth[i]]), '--k', linewidth=0.1)
plt.plot([1, 225, 225, 1, 1], [1,1,125,125,1], '--k')
plt.plot([1, 225, 225, 1, 1], [200,200,350,350,200], '--k')
fig.set_size_inches(w=9,h=6)
fig_name = 'SI_NO3.png'
fig.set_dpi(300)
fig.savefig(fig_name)


# PO4
fig, ax = plt.subplots()
#fig.clf()
ctf = plt.contourf(Xi, Yi, PO4i, 30, cmap=plt.cm.RdBu_r)
cl = plt.colorbar(ctf, orientation='vertical')
ax.set_ylim([0, 500])
ax.set_xlim([0, 300])
ax.set_ylabel('Depth (m)')
ax.set_xlabel('Distance (km)')
ax.set_title(r'$\rm PO_4 (mmol\,Kg^{-1})$ - Bonavista Section (1999-2017)')
ax.invert_yaxis()
Bgon = plt.Polygon(bathymetry,color=np.multiply([1,.9333,.6667],.4), alpha=1)
ax.add_patch(Bgon)
for i in range(0,len(distance)):
    plt.plot(np.array([distance[i], distance[i]]), np.array([Pbin[0], cast_depth[i]]), '--k', linewidth=0.1)
plt.plot([1, 225, 225, 1, 1], [1,1,125,125,1], '--k')
plt.plot([1, 225, 225, 1, 1], [200,200,350,350,200], '--k')
fig.set_size_inches(w=9,h=6)
fig_name = 'SI_PO4.png'
fig.set_dpi(300)
fig.savefig(fig_name)



# SiO
fig, ax = plt.subplots()
#fig.clf()
ctf = plt.contourf(Xi, Yi, SIOi, 30, cmap=plt.cm.RdBu_r)
cl = plt.colorbar(ctf, orientation='vertical')
ax.set_ylim([0, 500])
ax.set_xlim([0, 300])
ax.set_ylabel('Depth (m)')
ax.set_xlabel('Distance (km)')
ax.set_title(r'$\rm SiO_4 (mmol\,Kg^{-1})$ - Bonavista Section (1999-2017)')
ax.invert_yaxis()
Bgon = plt.Polygon(bathymetry,color=np.multiply([1,.9333,.6667],.4), alpha=1)
ax.add_patch(Bgon)
for i in range(0,len(distance)):
    plt.plot(np.array([distance[i], distance[i]]), np.array([Pbin[0], cast_depth[i]]), '--k', linewidth=0.1)
plt.plot([1, 225, 225, 1, 1], [1,1,125,125,1], '--k')
plt.plot([1, 225, 225, 1, 1], [200,200,350,350,200], '--k')
fig.set_size_inches(w=9,h=6)
fig_name = 'SI_SIO.png'
fig.set_dpi(300)
fig.savefig(fig_name)


# NO3/PO4
fig, ax = plt.subplots()
#fig.clf()
ctf = plt.contourf(Xi, Yi, NO3i/PO4i, 30, cmap=plt.cm.RdBu_r)
cl = plt.colorbar(ctf, orientation='vertical')
ax.set_ylim([0, 500])
ax.set_xlim([0, 300])
ax.set_ylabel('Depth (m)')
ax.set_xlabel('Distance (km)')
ax.set_title(r'$\rm [NO_3]/[PO_4]$ - Bonavista Section (1999-2017)')
ax.invert_yaxis()
Bgon = plt.Polygon(bathymetry,color=np.multiply([1,.9333,.6667],.4), alpha=1)
ax.add_patch(Bgon)
for i in range(0,len(distance)):
    plt.plot(np.array([distance[i], distance[i]]), np.array([Pbin[0], cast_depth[i]]), '--k', linewidth=0.1)
plt.plot([1, 225, 225, 1, 1], [1,1,125,125,1], '--k')
plt.plot([1, 225, 225, 1, 1], [200,200,350,350,200], '--k')
fig.set_size_inches(w=9,h=6)
fig_name = 'SI_NO3onPO4.png'
fig.set_dpi(300)
fig.savefig(fig_name)

