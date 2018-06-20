# A first test to read Excel nutrient file and export to Pandas.

# Check in:
#  /home/cyrf0006/research/AZMP_database/biochem

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import datetime

# Adjust fontsize/weight
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 18}
plt.rc('font', **font)

## ----  Load data from Excel sheet ---- ##
df = pd.read_excel('/home/cyrf0006/github/AZMP-NL/data/AZMP_Nutrients_1999_2016_Good_Flags_Only.xlsx')

## ---- Some cleaning ---- ##
# Set date as index
df = df.set_index('sas_date')

# Drop other time-related columns
df = df.drop(['Day', 'Month', 'Year'], axis=1)

# Replace HUD115 temperature by NaNs
## idx = ((df.index.month==11) & (df.index.year==2015))
## df['temp'][idx] = np.nan
## idx = ((df.index.month==12) & (df.index.year==2015))
## df['temp'][idx] = np.nan   
#OR
#df['temp'][(df['shiptrip']=='HUD115')]=np.nan


# Try to drop Smith Sound data
df = df.drop(df.index[df['section']=='SS'], axis=0)


#plt.hist(df['PO4'], 50)
#plt.show()


## ---- Get numpy array ---- ##
PO4 = df['PO4'].values
NO3 = df['NO3'].values
SIO = df['SIO'].values
T = df['temp'].values
S = df['salinity'].values
SIG = df['sigmat'].values
O2 = df['oxygen'].values
Z = df['depth'].values


## ---- PO4 - NO3 relationship ---- ##
# from Jones et al. (1998)
PO4_fit = np.linspace(0, np.max(PO4), 2)
NO3_fit_AW = 17.499*PO4_fit - 3.072
NO3_fit_PW = 12.368*PO4_fit - 10.549

fig = plt.figure(1)
plt.clf()
plt.scatter(PO4, NO3, c=df.index.year, alpha=0.5, cmap=plt.cm.RdBu_r)
plt.plot(PO4_fit, NO3_fit_AW, 'r--')
plt.plot(PO4_fit, NO3_fit_PW, 'r--')
plt.ylabel(r'NO3 ($\rm mmol m^{-3}$)')
plt.xlabel(r'PO4 ($\rm mmol m^{-3}$)')
plt.xlim([0,2.5])
plt.ylim([0,25])
plt.colorbar()
plt.text(1.5, 25.1, 'AW', fontweight='bold')
plt.text(2.25, 20.1, 'PW', fontweight='bold')
#plt.show()
plt.text(0.05,21, 'N = %d\n' % np.size(np.where((~np.isnan(PO4)) & (~np.isnan(NO3)))), fontsize=20, fontweight='bold')
#plt.show()

fig.set_size_inches(w=8, h=6)
fig.set_dpi(300)
fig.savefig('PO4-NO3_scatter.png')


## ---- PO4 - NO3 / T relationship ---- ##
fig = plt.figure(1)
plt.clf()
plt.scatter(PO4, NO3, c=Z, alpha=0.5)
plt.plot(PO4_fit, NO3_fit_AW, 'r--')
plt.plot(PO4_fit, NO3_fit_PW, 'r--')
plt.ylabel(r'NO3 ($\rm mmol m^{-3}$)')
plt.xlabel(r'PO4 ($\rm mmol m^{-3}$)')
plt.xlim([0,2.5])
plt.ylim([0,25])
plt.colorbar()
plt.text(1.5, 25.1, 'AW', fontweight='bold')
plt.text(2.25, 20.1, 'PW', fontweight='bold')
#plt.show()
plt.text(0.05,21, 'N = %d\n' % np.size(np.where((~np.isnan(PO4)) & (~np.isnan(NO3)))), fontsize=20, fontweight='bold')
#plt.show()

fig.set_size_inches(w=8, h=6)
fig.set_dpi(300)
fig.savefig('PO4-NO3_Z_scatter.png')


## ---- PO4 - NO3 / SIO relationship ---- ##
ratio = PO4/SIO
idx = np.where((np.log10(ratio)>-1) & (np.log10(ratio)<1))  
fig = plt.figure(1)
plt.clf()
plt.scatter(PO4[idx], NO3[idx], c=np.log10(ratio[idx]), alpha=0.5)
plt.plot(PO4_fit, NO3_fit_AW, 'r--')
plt.plot(PO4_fit, NO3_fit_PW, 'r--')
plt.ylabel(r'NO3 ($\rm mmol m^{-3}$)')
plt.xlabel(r'PO4 ($\rm mmol m^{-3}$)')
plt.xlim([0,2.5])
plt.ylim([0,25])
cb = plt.colorbar()
cb.ax.set_ylabel(r'$\rm log_{10}(PO4/SIO)$', fontsize=20, fontweight='bold')
plt.text(1.5, 25.1, 'AW', fontweight='bold')
plt.text(2.25, 20.1, 'PW', fontweight='bold')
#plt.show()
plt.text(0.05,21, 'N = %d\n' % np.size(np.where((~np.isnan(PO4[idx])) & (~np.isnan(NO3[idx])))), fontsize=20, fontweight='bold')
#plt.show()

fig.set_size_inches(w=8, h=6)
fig.set_dpi(300)
fig.savefig('PO4-NO3_ratioPO4-SIO_scatter.png')


## ---- T-S relationship ---- ##
plt.scatter(S, T, c=df.index.year, alpha=0.5)
plt.ylabel('temp')
plt.xlabel('salinity')
plt.colorbar()
plt.ylim([-2,20])
plt.xlim([30,37])
plt.show()


## ---- T-O2 relationship ---- ##
fig = plt.figure(1)
plt.clf()
plt.scatter(T, O2, c=df.index.year, alpha=0.5, cmap=plt.cm.RdBu_r)
plt.ylabel(r'$\rm O_2 (mg L^{-1})$', fontsize=20, fontweight='bold')
plt.xlabel(r'$\rm T (^{\circ}C)$', fontsize=20, fontweight='bold')
cb = plt.colorbar()
#cb.ax.set_ylabel('Years', fontsize=20, fontweight='bold')
plt.ylim([2,12])
plt.xlim([-2,19])
plt.text(10,10, 'N = %d\n' % np.size(np.where((~np.isnan(T)) & (~np.isnan(O2)))), fontsize=20, fontweight='bold')
#plt.show()

fig.set_size_inches(w=8, h=6)
fig.set_dpi(300)
fig.savefig('T-O2_scatter.png')

## ---- T-S Oxygen relationship ---- ##
fig = plt.figure(1)
plt.clf()
plt.scatter(S, T, c=O2, alpha=0.3)
plt.ylabel('T')
plt.xlabel('S')
cb = plt.colorbar()
cb.ax.set_ylabel('Years')
plt.ylim([-2,20])
plt.xlim([30,37])
plt.show()


## ---- S-NO3 relationship ---- ##
plt.scatter(S, NO3, c=df.index.year, alpha=0.3)
plt.ylabel('NO3')
plt.xlabel('S')
plt.colorbar()
#plt.ylim([-2,20])
#plt.xlim([30,37])
plt.show()


## ---- S-NO3 relationship (depth) ---- ##
zmax = 500
zmin = 50
Ztop = Z[(Z>zmin) & (Z<zmax) & (NO3>0) & (S<=33.5)]
Stop = S[(Z>zmin) & (Z<zmax) & (NO3>0) & (S<=33.5)]
NO3top = NO3[(Z>zmin) & (Z<zmax) & (NO3>0) & (S<=33.5)]

#Z[(Z>zmin) & (Z<zmax) & (NO3>0)].size / Z[(Z>zmin) & (Z<zmax) & (NO3>0) & (S<=33.5)].size

plt.scatter(Stop, NO3top, c=Ztop, alpha=0.3)
plt.ylabel('NO3')
plt.xlabel('S')
plt.colorbar()
#plt.ylim([-2,20])
#plt.xlim([30,37])
plt.show()



idx = [i for i,v in enumerate(NO3) if v > 1e-4]
ratio =  PO4[idx]/NO3[idx]
plt.hist(ratio, np.linspace(0,10,1000))
plt.show()
plt.plot(df['NO3'][idx], ratio, '.k')
plt.show()

