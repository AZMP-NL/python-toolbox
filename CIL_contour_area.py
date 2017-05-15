## Play around with CTD cast to plot transects along AZMP lines ##
# for now, run in /home/cyrf0006/research/AZMP/2017Spring after doing:
#  $ grep -l SEGB 39173/*.cnv > SEGB2017.list


import os
from seabird.cnv import fCNV
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
from matplotlib.mlab import griddata
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

## ------------------- Parameters to be edited ----------------- ##
YLIMS = [0, 500]
CHL_MAX = 12 
St27 = [47.550, -52.590]

#filelist = np.genfromtxt('SEGB2017.list', dtype=str)
#bathy = np.array([52, 180, 170, 126, 89, 81, 79, 71, 70, 71, 70, 65, 66, 171, 396, 1466, 2502, 2854, 3125, 3247, 3380])

#filelist = np.genfromtxt('SESPB2017.list', dtype=str)
#filelist = np.genfromtxt('SWSPB2017.list', dtype=str)
filelist = np.genfromtxt('FC2017.list', dtype=str)
bathy = np.flipud(np.array([104,185,155,124,138,100,105,172,82,80,102,134,164,215,550,880,1152,1159,993,346,300,274,247,167,148,148,133,158,272,326,548,1214,44000,4400,4400]))

## ------------------------------------------------------------- ##

LATlist = []
LONlist = []
Tlist = []
Slist = []
Plist = []
Clist = []
SIGlist = [];
Flist = [];
O2list = [];
PHlist = [];

Pbin = np.arange(2.5, 1200, 5)

    
for fname in filelist:
    profile = fCNV(fname)

    LATlist.append(profile.attributes['LATITUDE'])
    LONlist.append(profile.attributes['LONGITUDE'])

    # Must get profile, remove upcast + 5-m bin average
    P = np.array(profile['PRES'])
    T = np.array(profile['TEMP'])
    S = np.array(profile['PSAL'])
    C = np.array(profile['CNDC'])
    SIG = np.array(profile['sigma_t'])
    F = np.array(profile['flECO-AFL'])
    O2 = np.array(profile['oxigen_ml_L'])
    PH= np.array(profile['ph'])
   
    Ibtm = np.argmax(P)    
    digitized = np.digitize(P[0:Ibtm], Pbin) #<- this is awesome!

    Tlist.append([T[0:Ibtm][digitized == i].mean() for i in range(0, len(Pbin))])
    Slist.append([S[0:Ibtm][digitized == i].mean() for i in range(0, len(Pbin))])
    Plist.append([P[0:Ibtm][digitized == i].mean() for i in range(0, len(Pbin))])
    Clist.append([C[0:Ibtm][digitized == i].mean() for i in range(0, len(Pbin))])
    SIGlist.append([SIG[0:Ibtm][digitized == i].mean() for i in range(0, len(Pbin))])
    Flist.append([F[0:Ibtm][digitized == i].mean() for i in range(0, len(Pbin))])
    O2list.append([O2[0:Ibtm][digitized == i].mean() for i in range(0, len(Pbin))])
    PHlist.append([PH[0:Ibtm][digitized == i].mean() for i in range(0, len(Pbin))])


# List2Array
LATarray = np.array(LATlist)
LONarray = np.array(LONlist)

Tarray = np.transpose(np.array(Tlist))
Sarray = np.transpose(np.array(Slist))
Parray = np.transpose(np.array(Plist))
SIGarray = np.transpose(np.array(SIGlist))
Farray = np.transpose(np.array(Flist))
Oarray = np.transpose(np.array(O2list))
PHarray = np.transpose(np.array(PHlist))

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

distance = np.zeros(np.shape(LATarray))
for i in range(len(LATarray)):
    distance[i] = haversine(LONarray[0], LATarray[0], LONarray[i], LATarray[i])
    
## # for bathymetry:
if 'bathy' in locals():
    bathy_x = np.append(distance, [distance[-1], distance[0], distance[0] ])
    bathy_y = np.append(bathy, [np.max(bathy), np.max(bathy), bathy[0]])
    bathymetry = zip(bathy_x, bathy_y)

# Check maximum depth
## cast_depth = []
## for i in distance:
##     min_idx = np.argmin(np.abs(i-distance_bathy))
##     cast_depth = np.append(cast_depth, bathy[min_idx])

cast_depth = bathy
     
## --------  CIL surface from regular grid ---------- ##
dx = 1 #km
dz = Pbin[1]-Pbin[0] #m
xi = np.arange(distance[0],distance[-1],1)
yi = Pbin;

x, y = np.meshgrid(distance, Pbin) 
z = Tarray
x = x.reshape(np.size(z))
y = y.reshape(np.size(z))
z = z.reshape(np.size(z))
x = x[~np.isnan(z)]
y = y[~np.isnan(z)]
z = z[~np.isnan(z)]

# grid the data
zi = griddata(x,y,z,xi,yi, interp='linear')

# mask data in topo
zi = np.array(zi)
bathy_itp = np.interp(xi,distance, bathy)
for i in range(len(bathy_itp)):
    zi[Pbin>bathy_itp[i],i] = np.nan
    
# Compte CIL volume after removing nans
zi_lin = zi.reshape(np.size(zi))
zi_lin = zi_lin[~np.isnan(zi_lin)]
cil_idx = np.where(zi_lin < 0)
cil_vol = np.size(cil_idx)*dz*dx/1000

# Check which side we are going
if haversine(LONarray[0], LATarray[0], St27[1], St27[0]) > haversine(LONarray[-1], LATarray[-1], St27[1], St27[0]):
    distance = np.abs(distance-distance.max())
    xi = np.abs(xi-xi.max())

# for bathymetry (after computing good side):
if 'bathy' in locals():
    bathy_x = np.append(distance, [distance[-1], distance[0], distance[0] ])
    bathy_y = np.append(bathy, [np.max(bathy), np.max(bathy), bathy[0]])
    bathymetry = zip(bathy_x, bathy_y)   
## -------------------------------------------------- ## 
## --------  CIL surface from raw data contours ---------- ##
fig, axes = plt.subplots(nrows=2, ncols=1)


# S0 - Tcontour
plt.axes(axes[0])
ctf = plt.contourf(distance, Pbin, Tarray, 30, cmap=plt.cm.RdBu_r, y_dir='reverse')
c_cil = plt.contour(distance, Pbin, Tarray, [0,], colors='k', linewidths=2)
#ct = plt.contour(distance, Pbin, SIGarray, 10, colors='k', linewidths=0.5)
cl = plt.colorbar(ctf, orientation='vertical')
axes[0].tick_params(labelbottom='off')
axes[0].set_ylim(YLIMS)
axes[0].set_xlim(distance[0], distance[-1])
axes[0].set_ylabel('Depth (m)')
axes[0].invert_yaxis()
axes[0].invert_xaxis()
axes[0].text(5, YLIMS[1]*.90, r'T($^{\circ}$C)', horizontalalignment='left', verticalalignment='center', fontsize=16, color='k')

# print bathymetry
if 'bathy' in locals():
    Bgon = plt.Polygon(bathymetry,color=np.multiply([1,.9333,.6667],.4), alpha=0.8)
    axes[1].add_patch(Bgon)
for i in range(0,len(distance)):
    plt.plot(np.array([distance[i], distance[i]]), np.array([Pbin[0], cast_depth[i]]), '--k', linewidth=0.1)


# CIL surface
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

CIL = c_cil.collections[0]
vs = CIL.get_paths()[0].vertices
cil_vol2 = np.abs(area(vs))/1000.
print cil_vol2



# S1 - Tgrid
plt.axes(axes[1])
CS = plt.contour(xi,yi,zi,[0,],linewidths=2,colors='k')
CS = plt.contourf(xi,yi,zi,30,cmap=plt.cm.RdBu_r, y_dir='reverse')
cl = plt.colorbar(CS, orientation='vertical')
axes[1].set_xlim(distance[0], distance[-1])
axes[1].set_ylim(YLIMS)
axes[1].set_ylabel('Depth (m)')
axes[1].invert_yaxis()
axes[1].invert_xaxis()
axes[1].text(5, YLIMS[1]*.90, r'T($^{\circ}$C)', horizontalalignment='left', verticalalignment='center', fontsize=16, color='k')

# print bathymetry
if 'bathy' in locals():
    Bgon = plt.Polygon(bathymetry,color=np.multiply([1,.9333,.6667],.4), alpha=0.8)
    axes[0].add_patch(Bgon)
for i in range(0,len(distance)):
    plt.plot(np.array([distance[i], distance[i]]), np.array([Pbin[0], cast_depth[i]]), '--k', linewidth=0.1)

print cil_vol

plt.show()
