# -*- coding: utf-8 -*-

###########CSL changed from SS to GSL in here#########
###########for ease of importing, data files have been altered due to formatting#########

"""
Created on Wed Feb  6 09:43:50 2019

@author: GIBBO
"""
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np
import seaborn as sns
import datetime
from scipy import stats
import seawater as swx
import cmocean
import co2sys
from sklearn.impute import SimpleImputer

###Adjust fontsize/weight
font = {'family' : 'DejaVu Sans',
        'weight' : 'bold',
        'size'   : 12}
plt.rc('font', **font)
pd.set_option('display.max_rows', 500)

variables = ['TripID', 'Region', 'StationID', 'timestamp', 'latitude', 'longitude', 'depth', 'NO3', 'O2', 'TA', 'temp_pH', 'pH_25', 'PO4', 'salinity', 'SiO', 'temperature', 'TIC', 'TAflag', 'pHflag', 'TICflag']

################### Load NL data into a Pandas DataFrame #####################################
dfNL = pd.read_excel('C:\Users\gibbo\Documents\data\AZMP_NL\AZMP_NL_2014_2019.xlsx', delimiter=',')
dfNL['timestamp'] =  pd.to_datetime(dfNL['Date'], format='%m/%d/%Y')
dfNL = dfNL.reset_index(drop=True)

## Remove space in column names ##
cols = dfNL.columns
cols = cols.map(lambda x: x.replace(' ', '_') if isinstance(x, (str, unicode)) else x)
dfNL.columns = cols

###########################rename columns##########################
dfNL = dfNL.rename(columns={'Ship_Trip' : 'TripID'})
dfNL = dfNL.rename(columns={'Latitude_(Decimal_Degrees)' : 'latitude'})
dfNL = dfNL.rename(columns={'Longiude_(Decimal_Degrees)' : 'longitude'})
dfNL = dfNL.rename(columns={'Station_Name' : 'StationID'})
dfNL = dfNL.rename(columns={'Pressure_(dbar)' : 'depth'})
dfNL = dfNL.rename(columns={'Salinity_(psu)' : 'salinity'})
dfNL = dfNL.rename(columns={'Temperature_(C)' : 'temperature'})
dfNL = dfNL.rename(columns={'Dissolved_Oxygen_(mL/L)' : 'O2'})
dfNL = dfNL.rename(columns={'Calibrated_Oxygen_(mL/L)' : 'oxygen'})
dfNL = dfNL.rename(columns={'Phosphate_Concentration_(mmol_m-3)' : 'PO4'})
dfNL = dfNL.rename(columns={'Silicate_Concentration_(mmol_m-3)' : 'SiO'})
dfNL = dfNL.rename(columns={'Nitrate_Concentration_(mmol_m-3)' : 'NO3'})
dfNL = dfNL.rename(columns={'pH_25_' : 'pH_25'})
dfNL = dfNL.rename(columns={'Discrete_Chlorophyll_a_(mg_m-3)' : 'chla'})
dfNL = dfNL.rename(columns={'Total_P_(mmol_m-3)' : 'PO4'})
dfNL = dfNL.rename(columns={'Total_Si_(mmol_m-3)' : 'SiO'})
dfNL = dfNL.rename(columns={'Total_NO2/3_(mmol_m-3)' : 'NO3'})
dfNL = dfNL.rename(columns={'Measured_TA_(micromol/kg)' : 'TA'})
dfNL = dfNL.rename(columns={'Measure_TCO2_(micromol/kg)' : 'TIC'})
dfNL = dfNL.rename(columns={'Density_(kg_m3)' : 'density'})
dfNL = dfNL.rename(columns={'TA_Data_Flag' : 'TAflag'})
dfNL = dfNL.rename(columns={'TCO2_Data_Flag' : 'TICflag'})
dfNL['Region'] = 'NL'
dfNL['temp_pH']= np.NaN
dfNL['pHflag']= np.NaN

dfNL = dfNL.loc[:,variables]

df = dfNL.copy(deep=True)

##############Make one df of all three Regions (datasets) of Original Data to plot (includes temperature and salinity data without associated carbonate data#######
df = df.reset_index(drop=True)
df['satO2'] = swx.satO2(df['salinity'], df['temperature']) 
df['satO2_perc'] = df['O2']/df['satO2']*100  ######calculate percent oxygen saturation
df['AOU'] = df['O2']-df['satO2']
df['NO3']=df['NO3']/1.025 #####convert nutrient data from uM=mmol/m3/umol/L to umol/kgSW (/1.025)
df['SiO']=df['SiO']/1.025
df['PO4']=df['PO4']/1.025
#df.to_excel("C:\Users\gibbo\Documents\data\AZMP_OA\AZMP_OA_Alldata.xlsx")
dforig = df.copy(deep=True) ###all original data, since some is replaced/modified

#############data clean up for CO2sys = estimate TA from TA-S plot, remove empty rows, remove flagged data, replace nutrient nan with zeros################
#############drop flagged data####################
df.drop(df.loc[df['TAflag']>=3].index, inplace=True)
df.drop(df.loc[df['pHflag']>=3].index, inplace=True)
df.drop(df.loc[df['TICflag']>=3].index, inplace=True)

df.dropna(subset=['TA', 'TIC', 'pH_25'], axis=0, thresh=2, inplace=True)

#################### replace NaN in nutrient by zero (does not change anything in Excel CO2SYS)
no_nan = df.SiO[(df.PO4.isna()) | (df.SiO.isna())].size
perc_nan = float(no_nan)/df.shape[0]*100.0
df.PO4 = df.PO4.replace(np.nan, 0)
df.SiO = df.SiO.replace(np.nan, 0)
print(str(no_nan) + ' out of ' + str(df.shape[0]) + ' nutrient data where replaced by zeros (' + str(perc_nan) + ' %)')

##################below is the CO2sys script#############
###############change any parameters here################
df['par1type'] = (1)
df['par2type'] = (2)
df['par3type'] = (3)
df['presin'] = (0)
df['pHscale'] = (1)
df['k1k2c'] = (4)
df['kso4c'] = (1)

####################I separated the data according to sets of carbonate parameters###########
####################for TA and DIC data#######################
df_shelf = df.copy(deep=True) ###keeps all rows with a pH value of zero, left with TA and DIC
par1type = df_shelf['par1type'] # The first parameter supplied is of type "1", which is "alkalinity"
par2type = df_shelf['par2type']  # The first parameter supplied is of type "2", which is "DIC"
par3type = df_shelf['par3type']  # The first parameter supplied is of type "3", which is "pH"
presin   = df_shelf['presin']  # Pressure at input conditions
tempout  = df_shelf['temperature']  # Temperature at output conditions.
presout  = df_shelf['depth']  # Pressure    at output conditions.
pHscale  = df_shelf['pHscale']  # pH scale at which the input pH is reported ("1" means "Total Scale")
k1k2c    = df_shelf['k1k2c']  # Choice of H2CO3 and HCO3- dissociation constants K1 and K2 ("4" means "Mehrbach refit")
kso4c    = df_shelf['kso4c']  # Choice of HSO4- dissociation constants KSO4 ("1" means "Dickson")
alk_s  = df_shelf['TA']
dic_s  = df_shelf['TIC']
sal_s  = df_shelf['salinity']
pH_s = df_shelf['pH_25'] # pH measured in the lab (at ~25C)
temp_s = df_shelf['temp_pH'] # Temperature at which the pH was measured (~25C)
si_s   = df_shelf['SiO']
p_s    = df_shelf['PO4']
out_shelf, niceheaders = co2sys.CO2SYS(alk_s, dic_s, par1type, par2type, sal_s, temp_s, tempout, presin, presout, si_s, p_s, pHscale, k1k2c, kso4c)
dfCO2out_shelf = pd.DataFrame.from_dict(out_shelf, orient='index')
dfCO2out_shelf = dfCO2out_shelf.transpose()
dfCO2out_shelf = dfCO2out_shelf[['TAlk', 'TCO2', 'pHoutTOTAL', 'pCO2out', 'OmegaCAout', 'OmegaARout']].copy() 
dfCO2out_shelf = dfCO2out_shelf.set_index(df_shelf.index)
dfCO2shelf = pd.merge(df_shelf, dfCO2out_shelf, how='outer', left_index=True, right_index=True)

#############combine data and create/cleanup new dataframe for PCA or plotting############
variables = ['TripID', 'Region', 'StationID', 'timestamp', 'latitude', 'longitude', 'depth', 'temperature', 'salinity', 'satO2_perc', 'NO3', 'PO4', 'SiO', 'TCO2', 'TAlk', 'pHoutTOTAL', 'pCO2out', 'OmegaCAout', 'OmegaARout']
dfCO2sys = dfCO2shelf.loc[:,variables]
dfCO2sys = dfCO2sys.rename(columns={'TAlk' : 'TA'})
dfCO2sys = dfCO2sys.rename(columns={'TCO2' : 'TIC'})
dfCO2sys = dfCO2sys.rename(columns={'pHoutTOTAL' : 'pH'})
dfCO2sys = dfCO2sys.rename(columns={'pCO2out' : 'pCO2'})
dfCO2sys = dfCO2sys.rename(columns={'OmegaCAout' : 'Omega_C'})
dfCO2sys = dfCO2sys.rename(columns={'OmegaARout' : 'Omega_A'})
dfCO2sys = dfCO2sys.dropna(subset=['Omega_C'], axis=0)
#dfCO2sys.to_excel("C:\Users\gibbo\Documents\data\AZMP_OA\AZMP_OA_CO2stats.xlsx")

################recombine with original dataset that includes the NaNs for plotting then cleanup##################
var = ['TripID', 'Region', 'StationID', 'timestamp', 'latitude', 'longitude', 'depth', 'temperature', 'salinity', 'satO2_perc', 'NO3', 'PO4', 'SiO']
dforig = dforig.loc[:,var]
dfCO2var = dfCO2sys.loc[:,['TA', 'TIC', 'pH', 'pCO2', 'Omega_A', 'Omega_C']]
dfplot = pd.merge(dforig, dfCO2var, how='left', left_index=True, right_index=True)
dfplot = dfplot.dropna(axis=0, thresh=8)
#dfplot.to_excel("C:\Users\gibbo\Documents\data\AZMP_OA\AZMP_OA_plot.xlsx")



##########below is the original CO2SYS function################################
#import pandas as pd
#import co2sys

#par1type = 1  # The first parameter supplied is of type "1", which is "alkalinity"
#par2type = 2  # The first parameter supplied is of type "2", which is "DIC"
#par3type = 3  # The first parameter supplied is of type "3", which is "pH"
#presin   = 0  # Pressure at input conditions
#tempout  = 15.267  # Temperature at output conditions.
#presout  = 3.498  # Pressure    at output conditions.
#pHscale  = 1  # pH scale at which the input pH is reported ("1" means "Total Scale")
#k1k2c    = 4  # Choice of H2CO3 and HCO3- dissociation constants K1 and K2 ("4" means "Mehrbach refit")
#kso4c    = 1  # Choice of HSO4- dissociation constants KSO4 ("1" means "Dickson")
#alk_s  = 2134.7
#dic_s  = 1948.8
#sal_s  = 31.066
#pH_s = 7.892
#temp_s = 25
#si_s   = 1.005
#p_s    = 0.1945
#
#out, niceheaders = co2sys.CO2SYS(alk_s, pH_s, par1type, par3type, sal_s, temp_s, tempout, presin, presout, si_s, p_s, pHscale, k1k2c, kso4c)
#dfCO2out = pd.DataFrame.from_dict(out, orient='index')
#dfCO2sys = dfCO2out.transpose()
#dfCO2sys = dfCO2sys[['TAlk', 'TCO2', 'pHoutTOTAL', 'pCO2out', 'OmegaCAout', 'OmegaARout']].copy() 
#dfCO2sys = dfCO2sys.dropna(axis=0, how='any')
#dfCO2sys.to_excel("C:\Users\gibbo\Documents\data\AZMP_OA\AZMP_OA__NL_CO2sys.xlsx")
