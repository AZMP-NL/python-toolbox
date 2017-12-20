## Play around with CTD cast to plot transects along AZMP lines ##
# for now, run in /home/cyrf0006/research/AZMP/2017Spring after doing:
#  $ grep -l SEGB 39173/*.cnv > SEGB2017.list


import os
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
from scipy.interpolate import griddata
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from echonix import ek60, echogram, imaging, raw

## ------------------- Parameters to be edited ----------------- ##

filelist = np.genfromtxt('EK60.list', dtype=str)
decim_factorX = 10
decim_factorZ = 1


## ------------------------------------------------------------- ##


#Pbin = np.arange(2.5, 1200, 5)
d = raw.load_raw(filelist[0])
t0 = d[0].dgheader.datetime[0]
d = raw.load_raw(filelist[-1])
tf = d[0].dgheader.datetime[-1]

print(d[0].configurationheader.surveyname)
    

Sv38, r = ek60.raws_to_sv(filelist, 38000)
Sv120, r = ek60.raws_to_sv(filelist, 120000)
Sv200, r = ek60.raws_to_sv(filelist, 200000)

Sv38 = Sv38[::decim_factorZ, ::decim_factorX]
Sv120 = Sv120[::decim_factorZ, ::decim_factorX]
Sv200 = Sv200[::decim_factorZ, ::decim_factorX]

### Note to myself:
# check here how I can extract time and depth vector and build a dataframe to pickle....


df = pd.DataFrame(np.transpose(EPSarray), index=MSStime, columns=Pbin)
df.to_pickle('MSS_S1.pkl')
# df = pd.read_pickle(file_name) t





    
im = imaging.composite(Sv38, Sv120, Sv200, min = -95, max = -50)

echogram.imshow(im, range=r)

plt.xlabel('Sample')
plt.ylabel('Range /m')

plt.title('Test plot EK60')
plt.colorbar()

plt.show()

keyboard


    Sv38 = Sv38[::decim_factorZ, ::decim_factorX]
    Sv120 = Sv120[::decim_factorZ, ::decim_factorX]
    Sv200 = Sv200[::decim_factorZ, ::decim_factorX]

    Sv38_list.append(Sv38)
    
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
