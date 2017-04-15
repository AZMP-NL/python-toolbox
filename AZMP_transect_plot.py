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

filelist = np.genfromtxt('SEGB2017.list', dtype=str)
bathy = np.array([52, 180, 170, 126, 89, 81, 79, 71, 70, 71, 70, 65, 66, 171, 396, 1466, 2502, 2854, 3125, 3247, 3380])

#filelist = np.genfromtxt('SESPB2017.list', dtype=str)
#filelist = np.genfromtxt('SWSPB2017.list', dtype=str)

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
    
    Ibtm = np.argmax(P)    
    digitized = np.digitize(P[0:Ibtm], Pbin) #<- this is awesome!

    Tlist.append([T[0:Ibtm][digitized == i].mean() for i in range(0, len(Pbin))])
    Slist.append([S[0:Ibtm][digitized == i].mean() for i in range(0, len(Pbin))])
    Plist.append([P[0:Ibtm][digitized == i].mean() for i in range(0, len(Pbin))])
    Clist.append([C[0:Ibtm][digitized == i].mean() for i in range(0, len(Pbin))])
    SIGlist.append([SIG[0:Ibtm][digitized == i].mean() for i in range(0, len(Pbin))])
    Flist.append([F[0:Ibtm][digitized == i].mean() for i in range(0, len(Pbin))])
    O2list.append([O2[0:Ibtm][digitized == i].mean() for i in range(0, len(Pbin))])
    

# List2Array
LATarray = np.array(LATlist)
LONarray = np.array(LONlist)
myCasts = np.zeros((len(LATarray)))-5

Tarray = np.transpose(np.array(Tlist))
Sarray = np.transpose(np.array(Slist))
SIGarray = np.transpose(np.array(SIGlist))
Farray = np.transpose(np.array(Flist))
Oarray = np.transpose(np.array(O2list))

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

# for bathymetry:
bathy_x = np.append(distance, [distance[-1], distance[0], distance[0] ])
bathy_y = np.append(bathy, [np.max(bathy), np.max(bathy), bathy[0]])
bathymetry = zip(bathy_x, bathy_y)

## ---- now plot ---- ##

fig, axes = plt.subplots(nrows=4, ncols=1)

# S0 - T
plt.axes(axes[0])
ctf = plt.contourf(distance, Pbin, Tarray, 30, cmap=plt.cm.RdBu_r, y_dir='reverse')
ct = plt.contour(distance, Pbin, SIGarray, 10, colors='k', linewidths=0.5)
#plt.plot(distance, myCasts)
cl = plt.colorbar(ctf, orientation='vertical')
#plt.plot(grid=True)
axes[0].tick_params(labelbottom='off')
axes[0].set_ylim(0,250)
axes[0].set_ylabel('Depth (m)')
axes[0].invert_yaxis()
for i in range(0,len(distance)):
    plt.plot(np.array([distance[i], distance[i]]), np.array([Pbin[0], bathy[i]]), '--k', linewidth=0.1)
axes[0].text(5, 225, r'T($^{\circ}$C)', horizontalalignment='left', verticalalignment='center', fontsize=16, color='k')
# print bathymetry
Bgon = plt.Polygon(bathymetry,color=np.multiply([1,.9333,.6667],.4), alpha=0.8)
axes[0].add_patch(Bgon)

# S1 - S
plt.axes(axes[1])
ctf = plt.contourf(distance, Pbin, Sarray, 30, cmap=plt.cm.RdBu_r)
ct = plt.contour(distance, Pbin, SIGarray, 10, colors='k', linewidths=0.5)
cl = plt.colorbar(ctf, orientation='vertical')
axes[1].tick_params(labelbottom='off')
axes[1].set_ylim(0,250)
axes[1].set_ylabel('Depth (m)')
axes[1].invert_yaxis()
for i in range(0,len(distance)):
    plt.plot(np.array([distance[i], distance[i]]), np.array([Pbin[0], bathy[i]]), '--k', linewidth=0.1)
axes[1].text(5, 225, r'$\rm S_p$', horizontalalignment='left', verticalalignment='center', fontsize=16, color='k')
# print bathymetry
Bgon = plt.Polygon(bathymetry,color=np.multiply([1,.9333,.6667],.4), alpha=0.8)
axes[1].add_patch(Bgon)


# S2 - O2
plt.axes(axes[2])
ctf = plt.contourf(distance, Pbin, Oarray, 30, cmap=plt.cm.RdBu)
ct = plt.contour(distance, Pbin, SIGarray, 10, linewidths=0.5, colors='k')
cl = plt.colorbar(ctf, orientation='vertical')
axes[2].tick_params(labelbottom='off')
axes[2].set_ylim(0,250)
axes[2].set_ylabel('Depth (m)')
axes[2].invert_yaxis()
for i in range(0,len(distance)):
    plt.plot(np.array([distance[i], distance[i]]), np.array([Pbin[0], bathy[i]]), '--k', linewidth=0.1)
axes[2].text(5, 225, r'$\rm O_2$($\rm ml L^{-1}$)', horizontalalignment='left', verticalalignment='center', fontsize=16, color='k')
# print bathymetry
Bgon = plt.Polygon(bathymetry,color=np.multiply([1,.9333,.6667],.4), alpha=0.8)
axes[2].add_patch(Bgon)

# S3 - CHL
plt.axes(axes[3])
ctf = plt.contourf(distance, Pbin, Farray, 30, cmap=plt.cm.PuBuGn)
ct = plt.contour(distance, Pbin, SIGarray, 10, linewidths=0.5, colors='k')
cl = plt.colorbar(ctf, orientation='vertical')
axes[3].set_ylim(0,250)
axes[3].set_ylabel('Depth (m)')
#axes[3].set_xlabel('Latitude ($^{\circ}$N)')
axes[3].set_xlabel('Along-transect distance (km)')
axes[3].invert_yaxis()
#axes[3].invert_xaxis()
for i in range(0,len(distance)):
    plt.plot(np.array([distance[i], distance[i]]), np.array([Pbin[0], bathy[i]]), '--k', linewidth=0.1)
axes[3].text(5, 225, r'chl-a($\rm mg m^{-3}$)', horizontalalignment='left', verticalalignment='center', fontsize=16, color='k')
# print bathymetry
Bgon = plt.Polygon(bathymetry,color=np.multiply([1,.9333,.6667],.4), alpha=0.8)
axes[3].add_patch(Bgon)

    
fig.set_size_inches(w=6,h=8)
fig.tight_layout()
fig_name = 'transect_to-be-rename.pdf'
fig.savefig(fig_name)
#os.system('pdfcrop %s %s &> /dev/null &'%(fig_name, fig_name))




