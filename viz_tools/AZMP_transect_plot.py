## Play around with CTD cast to plot transects along AZMP lines ##
# for now, run in /home/cyrf0006/research/AZMP/2017Spring after doing:
#  $ grep -l SEGB 39173/*.cnv > SEGB2017.list


import os
from seabird.cnv import fCNV
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
from scipy.interpolate import griddata
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

## ------------------- Parameters to be edited ----------------- ##
YLIMS = [0, 200]
CHL_MAX = 5
St27 = [47.550, -52.590]

#filelist = np.genfromtxt('SEGB2017.list', dtype=str)
#bathy = np.array([52, 180, 170, 126, 89, 81, 79, 71, 70, 71, 70, 65, 66, 171, 396, 1466, 2502, 2854, 3125, 3247, 3380])

#filelist = np.genfromtxt('SESPB2017.list', dtype=str)
#filelist = np.genfromtxt('SWSPB2017.list', dtype=str)
#filelist = np.genfromtxt('FC2017.list', dtype=str)
#filelist = np.genfromtxt('BB2017.list', dtype=str)
#filelist = np.genfromtxt('TEL176_FC.list', dtype=str)
#filelist = np.genfromtxt('TEL176_BB.list', dtype=str)
#filelist = np.genfromtxt('TEL176_WB.list', dtype=str)
#filelist = np.genfromtxt('TEL176_SI.list', dtype=str)
#filelist = np.genfromtxt('TEL176_MB.list', dtype=str)
#filelist = np.genfromtxt('TEL176_BI.list', dtype=str)
#filelist = np.genfromtxt('DIS009_SI.list', dtype=str)
#filelist = np.genfromtxt('FD009_FC.list', dtype=str)
#filelist = np.genfromtxt('FD009_BB.list', dtype=str)
filelist = np.genfromtxt('HUD118_SEGB.list', dtype=str)
#filelist = np.genfromtxt('HUD118_SWSPB.list', dtype=str)



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

    vars = profile.keys()
    
     # Must get profile, remove upcast + 5-m bin average
    P = np.array(profile['PRES'])
    T = np.array(profile['TEMP'])
    S = np.array(profile['PSAL'])
    C = np.array(profile['CNDC'])
    SIG = np.array(profile['sigma_t'])
    F = np.array(profile['flECO-AFL'])
    O2 = np.array(profile['oxigen_ml_L'])

    if 'ph' in vars:
        PH = np.array(profile['ph'])
    else:
        PH = np.array(profile['TEMP'])*np.nan
    
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

# Sort arrays according to distance
I = np.argsort(distance)
distance = distance[I]
LATarray = LATarray[I]
LONarray = LONarray[I]
Tarray = Tarray[:,I]
Sarray = Sarray[:,I]
Parray = Parray[:,I]
SIGarray = SIGarray[:,I]
Farray = Farray[:,I]
Oarray = Oarray[:,I]
PHarray = PHarray[:,I]
     
## # for bathymetry:
## if 'bathy' in locals():
##     bathy_x = np.append(distance, [distance[-1], distance[0], distance[0] ])
##     bathy_y = np.append(bathy, [np.max(bathy), np.max(bathy), bathy[0]])
##     bathymetry = zip(bathy_x, bathy_y)

# for GEBCO bathymetry
import coord_list
import netCDF4
import get_GEBCO
interval = 10000.0 #meters
azimuth = coord_list.calculateBearing(LATarray[0], LONarray[0], LATarray[-1], LONarray[-1])
coords = coord_list.main(interval,azimuth,LATarray[0], LONarray[0], LATarray[-1], LONarray[-1])
lat_max = max(l[0] for l in coords)
lat_min = min(l[0] for l in coords)
lon_max = max(l[1] for l in coords)
lon_min = min(l[1] for l in coords)
lat, lon, Z = get_GEBCO.main('/home/cyrf0006/data/GEBCO/GEBCO_08.nc', [-lat_max, -lat_min, lon_min, lon_max])
lat = -lat
Z = -Z
X, Y = np.meshgrid(lon, lat)  # grid X,Y
X = X.reshape(np.size(Z)) #<--- check if a fonction exists for that
Y = Y.reshape(np.size(Z))
Z = Z.reshape(np.size(Z))

bathy = []
distance_bathy = []
print "Computing bathymetry..."
for i in range(0,len(coords)):
    optim = np.abs(np.array(coords[i]) - np.array(zip(Y,X)))
    optim_sum = optim.sum(axis=1)
    min_idx = np.argmin(optim_sum)
    bathy = np.append(bathy, Z[min_idx])
    distance_bathy = np.append(distance_bathy, haversine(coords[0][1], coords[0][0], coords[i][1], coords[i][0])) 
print "done!"

# Check which direction we are going (approaching St.27 or not)
if haversine(LONarray[0], LATarray[0], St27[1], St27[0]) > haversine(LONarray[-1], LATarray[-1], St27[1], St27[0]):
    distance = np.abs(distance-distance.max())
    if 'bathy' in locals():
        distance_bathy = np.abs(distance_bathy - distance_bathy.max())

# Check maximum depth
cast_depth = []
for i in distance:
    min_idx = np.argmin(np.abs(i-distance_bathy))
    cast_depth = np.append(cast_depth, bathy[min_idx])
    
             
# Bathymetry Polygon
bathy_x = np.append(distance_bathy, [distance_bathy[-1], distance_bathy[0], distance_bathy[0] ])
bathy_y = np.append(bathy, [np.max(bathy), np.max(bathy), bathy[0]])
bathymetry = zip(bathy_x, bathy_y)


        
## ---- now plot ---- ##
fig, axes = plt.subplots(nrows=5, ncols=1)

# S0 - T
plt.axes(axes[0])
ctf = plt.contourf(distance, Pbin, Tarray, 30, cmap=plt.cm.RdBu_r, y_dir='reverse')
c_cil = plt.contour(distance, Pbin, Tarray, [0,], colors='k', linewidths=2)
ct = plt.contour(distance, Pbin, SIGarray, 10, colors='k', linewidths=0.5)
cl = plt.colorbar(ctf, orientation='vertical')
axes[0].tick_params(labelbottom='off')
axes[0].set_ylim(YLIMS)
axes[0].set_ylabel('Depth (m)')
axes[0].invert_yaxis()
axes[0].text(5, YLIMS[1]*.90, r'T($^{\circ}$C)', horizontalalignment='left', verticalalignment='center', fontsize=16, color='k')
# print bathymetry
if 'bathy' in locals():
    Bgon = plt.Polygon(bathymetry,color=np.multiply([1,.9333,.6667],.4), alpha=0.8)
    axes[0].add_patch(Bgon)
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
    
# S1 - S
plt.axes(axes[1])
ctf = plt.contourf(distance, Pbin, Sarray, 30, cmap=plt.cm.RdBu_r)
ct = plt.contour(distance, Pbin, SIGarray, 10, colors='k', linewidths=0.5)
cl = plt.colorbar(ctf, orientation='vertical')
axes[1].tick_params(labelbottom='off')
axes[1].set_ylim(YLIMS)
axes[1].set_ylabel('Depth (m)')
axes[1].invert_yaxis()
axes[1].text(5, YLIMS[1]*.90, r'$\rm S_p$', horizontalalignment='left', verticalalignment='center', fontsize=16, color='k')
plt.plot(np.array([distance[i], distance[i]]), np.array([Pbin[0], cast_depth[i]]), '--k', linewidth=0.1)  
if 'bathy' in locals():
    Bgon = plt.Polygon(bathymetry,color=np.multiply([1,.9333,.6667],.4), alpha=0.8)
    axes[1].add_patch(Bgon)
for i in range(0,len(distance)):
        plt.plot(np.array([distance[i], distance[i]]), np.array([Pbin[0], cast_depth[i]]), '--k', linewidth=0.1)  

# S2 - O2
plt.axes(axes[2])
ctf = plt.contourf(distance, Pbin, Oarray, 30, cmap=plt.cm.RdBu)
ct = plt.contour(distance, Pbin, SIGarray, 10, linewidths=0.5, colors='k')
cl = plt.colorbar(ctf, orientation='vertical')
axes[2].tick_params(labelbottom='off')
axes[2].set_ylim(YLIMS)
axes[2].set_ylabel('Depth (m)')
axes[2].invert_yaxis()
axes[2].text(5, YLIMS[1]*.90, r'$\rm O_2$($\rm ml L^{-1}$)', horizontalalignment='left', verticalalignment='center', fontsize=16, color='k')
if 'bathy' in locals():
    Bgon = plt.Polygon(bathymetry,color=np.multiply([1,.9333,.6667],.4), alpha=0.8)
    axes[2].add_patch(Bgon)
for i in range(0,len(distance)):
    plt.plot(np.array([distance[i], distance[i]]), np.array([Pbin[0], cast_depth[i]]), '--k', linewidth=0.1)

# S3 - CHL
plt.axes(axes[3])
cf = plt.contourf(distance, Pbin, Farray, levels=np.arange(0, CHL_MAX, 1), cmap=plt.cm.PuBuGn, extend="max") #extend='both'
cc = plt.contour(distance, Pbin, SIGarray, 10, linewidths=0.5, colors='k')
cf.cmap.set_under('k')
cf.set_clim(0, CHL_MAX)
cb = plt.colorbar(cf)
axes[3].set_ylim(YLIMS)
axes[3].tick_params(labelbottom='off')
axes[3].set_ylabel('Depth (m)')
axes[3].invert_yaxis()
axes[3].text(5, YLIMS[1]*.90, r'chl-a($\rm mg m^{-3}$)', horizontalalignment='left', verticalalignment='center', fontsize=16, color='k')
if 'bathy' in locals():
    Bgon = plt.Polygon(bathymetry,color=np.multiply([1,.9333,.6667],.4), alpha=0.8)
    axes[3].add_patch(Bgon)
for i in range(0,len(distance)):
    plt.plot(np.array([distance[i], distance[i]]), np.array([Pbin[0], cast_depth[i]]), '--k', linewidth=0.1)

# S4 - PH
plt.axes(axes[4])
ctf = plt.contourf(distance, Pbin, PHarray, 30, cmap=plt.cm.RdBu)
ct = plt.contour(distance, Pbin, SIGarray, 10, linewidths=0.5, colors='k')
cl = plt.colorbar(ctf, orientation='vertical')
axes[4].set_ylim(YLIMS)
axes[4].set_ylabel('Depth (m)')
axes[4].set_xlabel('Along-transect distance (km)')
axes[4].invert_yaxis()
axes[4].text(5, YLIMS[1]*.90, 'pH', horizontalalignment='left', verticalalignment='center', fontsize=16, color='k')
if 'bathy' in locals():
    Bgon = plt.Polygon(bathymetry,color=np.multiply([1,.9333,.6667],.4), alpha=0.8)
    axes[4].add_patch(Bgon)
for i in range(0,len(distance)):
    plt.plot(np.array([distance[i], distance[i]]), np.array([Pbin[0], cast_depth[i]]), '--k', linewidth=0.1)

fig.set_size_inches(w=6,h=8)
fig.tight_layout()
fig_name = 'transect_to-be-rename.pdf'
fig.savefig(fig_name)
#os.system('pdfcrop %s %s &> /dev/null &'%(fig_name, fig_name))


