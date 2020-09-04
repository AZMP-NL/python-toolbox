# -*- coding: utf-8 -*-
"""
Script to QA/QC and process all ocean cabonate data acquired as part of the AZMP in 2014.
The script also:

1. Convert cabonate measurements made from all regions to a comment format

2. Derive extras parameters using CO2sys

@author:
Olivia.Gibb@dfo-mpo.gc.ca
Frederic.Cyr@dfo-mpo.gc.ca

* Note on CO2sys usage:
CO2dict = CO2SYS(PAR1, PAR2, PAR1TYPE, PAR2TYPE, SAL, TEMPIN, TEMPOUT, PRESIN, PRESOUT,
    SI, PO4, pHSCALEIN, K1K2CONSTANTS, KSO4CONSTANTS, NH3=0.0, H2S=0.0, KFCONSTANT=1)
    
where: (* From https://github.com/mvdh7/PyCO2SYS)

Required inputs

The required inputs are identical to the MATLAB version:

    PAR1 - frst known carbonate system parameter value.
    PAR2 - second known carbonate system parameter value.
    PAR1TYPE - integer identifying which parameters PAR1 are.
    PAR2TYPE - integer identifying which parameters PAR2 are.

The possible known carbonate system parameters are 1: total alkalinity in μmol·kg−1, 2: dissolved inorganic carbon in μmol·kg−1, 3: pH (dimensionless), 4: dissolved CO2 partial pressure in μatm, 5: dissolved CO2 fugacity in μatm, and 6: carbonate ion concentration in μmol·kg−1.

Here and throughout the inputs and outputs, "kg−1" refers to the total mass of seawater (solvent + solutes), not just the mass of H2O.

    SAL - practical salinity.
    TEMPIN - temperature of input carbonate system parameters.
    TEMPOUT - temperature at which to calculate outputs.
    PRESIN - pressure of input carbonate system parameters.
    PRESOUT - pressure at which to calculate outputs.

All temperatures are in °C and pressures are in dbar. Pressure is within the water column as typically measured by a CTD sensor, i.e. not including atmospheric pressure. The 'input' conditions could represent conditions in the laboratory during a measurement, while the 'output' conditions could be those observed in situ during sample collection.

    SI - total silicate concentration.
    PO4 - total phosphate concentration.

Nutrient concentrations are all in μmol·kg−1.

    pHSCALEIN - pH scale(s) that pH values in PAR1 or PAR2 are on.

The options are 1: Total scale, 2: Seawater scale, 3: Free scale, and 4: NBS scale, as defined by ZW01.

    K1K2CONSTANTS - which set of constants to use for carbonic acid dissociation.

The options are integers from 1 to 15 inclusive. From the original MATLAB documentation:

%   1 = Roy, 1993                                         T:    0-45  S:  5-45. Total scale. Artificial seawater.
%   2 = Goyet & Poisson                                   T:   -1-40  S: 10-50. Seaw. scale. Artificial seawater.
%   3 = HANSSON              refit BY DICKSON AND MILLERO T:    2-35  S: 20-40. Seaw. scale. Artificial seawater.
%   4 = MEHRBACH             refit BY DICKSON AND MILLERO T:    2-35  S: 20-40. Seaw. scale. Artificial seawater.
%   5 = HANSSON and MEHRBACH refit BY DICKSON AND MILLERO T:    2-35  S: 20-40. Seaw. scale. Artificial seawater.
%   6 = GEOSECS (i.e., original Mehrbach)                 T:    2-35  S: 19-43. NBS scale.   Real seawater.
%   7 = Peng    (i.e., original Mehrbach but without XXX) T:    2-35  S: 19-43. NBS scale.   Real seawater.
%   8 = Millero, 1979, FOR PURE WATER ONLY (i.e., Sal=0)  T:    0-50  S:     0.
%   9 = Cai and Wang, 1998                                T:    2-35  S:  0-49. NBS scale.   Real and artificial seawater.
%  10 = Lueker et al, 2000                                T:    2-35  S: 19-43. Total scale. Real seawater.
%  11 = Mojica Prieto and Millero, 2002.                  T:    0-45  S:  5-42. Seaw. scale. Real seawater
%  12 = Millero et al, 2002                               T: -1.6-35  S: 34-37. Seaw. scale. Field measurements.
%  13 = Millero et al, 2006                               T:    0-50  S:  1-50. Seaw. scale. Real seawater.
%  14 = Millero        2010                               T:    0-50  S:  1-50. Seaw. scale. Real seawater.
%  15 = Waters, Millero, & Woosley 2014                   T:    0-50  S:  1-50. Seaw. scale. Real seawater.

    KSO4CONSTANTS - which sets of constants to use for bisulfate dissociation and borate:chlorinity ratio.

The options are integers from 1 to 4 inclusive. From the original MATLAB documentation:

%  1 = KSO4 of Dickson 1990a   & TB of Uppstrom 1974  (PREFERRED)
%  2 = KSO4 of Khoo et al 1977 & TB of Uppstrom 1974
%  3 = KSO4 of Dickson 1990a   & TB of Lee et al. 2010
%  4 = KSO4 of Khoo et al 1977 & TB of Lee et al. 2010
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import datetime
from scipy import stats
import seawater as swx
from PyCO2SYS import CO2SYS 
from PyCO2SYS.meta import version
    
# Setup paths (user specific, to be adjusted)
#for cyrf0006:
dataset_main_path = '/home/cyrf0006/github/AZMP-NL/datasets/carbonates/'
dataset_path = '/home/cyrf0006/github/AZMP-NL/datasets/carbonates/data_received/'
dataset_path2 = '/home/cyrf0006/github/AZMP-NL/datasets/carbonates/raw_data/'
# for gibbo
#dataset_path = 'C:\Users\gibbo\Documents\data\AZMP_OA\'

###########################################
## ----------- Load all data ----------- ##
###########################################

## 1.1 --- BIO data from BIOChem (GSL dataset are also here, but we just keep 2014)
df_MAR = pd.read_csv(os.path.join(dataset_path, 'Discrete_Data_ 20190305.csv'), delimiter=',')
df_MAR['timestamp'] = df_MAR.year.astype(str) + '-' + df_MAR.month.astype(str) + '-' + df_MAR.day.astype(str) + ' ' + (df_MAR.time.astype(str)).str.zfill(4)
df_MAR['timestamp'] = pd.to_datetime(df_MAR['timestamp'], format='%Y-%m-%d %H%M')
#df_MAR.index = pd.to_datetime(df_MAR.timestamp)

# Rename columns (HERE! update field name or rename later)
df_MAR = df_MAR.rename(columns={'mission_descriptor' : 'TripID'})
df_MAR = df_MAR.rename(columns={'start_depth' : 'depth'})
df_MAR = df_MAR.rename(columns={'collector_station_name' : 'StationID'})
df_MAR = df_MAR.rename(columns={'Alkalinity_umol/kg' : 'TA'})
#df_MAR = df_MAR.rename(columns={'Chl_a_Holm-Hansen_F' : 'chla'})
df_MAR = df_MAR.rename(columns={'NO2NO3_Tech_F' : 'NO3'})
df_MAR = df_MAR.rename(columns={'O2_CTD_mLL' : 'O2CTD'})
df_MAR = df_MAR.rename(columns={'O2_Winkler_Auto' : 'O2'})
df_MAR = df_MAR.rename(columns={'PO4_Tech_F' : 'PO4'})
df_MAR = df_MAR.rename(columns={'SiO4_Tech_F' : 'SiO'})
df_MAR = df_MAR.rename(columns={'Temp_CTD_1968' : 'temperature'})
df_MAR = df_MAR.rename(columns={'Salinity_CTD' : 'salinity'})
df_MAR = df_MAR.rename(columns={'pH_invitro_total_f' : 'pH_25'}) # this is pH measured at 25 degrees -not from CTD, used to calculate pH insitu
df_MAR = df_MAR.rename(columns={'pH_invitro_temp' : 'temp_pH'})  # temperature at which the pH was measured (~25C)

# combine BIO and IML (GSL 2014) columns (they had different headers)
df_MAR.TA.update(df_MAR.ALKW_01)
df_MAR.TIC.update(df_MAR.TICW_01)
#df_MAR.chla.update(df_MAR.Chl_a_Welschmeyer_sF)
df_MAR.NO3.update(df_MAR.NO2NO3_Tech_SF)
df_MAR.NO3.update(df_MAR.NO2NO3_Tech_Fsh)
df_MAR.PO4.update(df_MAR.PO4_Tech_SF)
df_MAR.PO4.update(df_MAR.PO4_Tech_Fsh)
df_MAR.SiO.update(df_MAR.SiO4_Tech_Fsh)
df_MAR.SiO.update(df_MAR.SiO4_Tech_SF)
df_MAR.temperature.update(df_MAR.Temp_CTD_1990)
df_MAR.pH_25.update(df_MAR.pH_invitro_total_p)
df_MAR['Region'] = 'MAR'
df_MAR.loc[(df_MAR['collector_event_id'].str.contains('IML')), 'Region'] = 'GSL'
#change units of GSL data from umol/L to ml/L (/44.661):
df_MAR['O2GSL'] = df_MAR.loc[df_MAR['Region'].str.contains('GSL')]['O2']/44.661
df_MAR.O2.update(df_MAR.O2GSL)

## 1.2 --- BIO flag data
bio_flags = pd.read_excel(os.path.join(dataset_path2, 'DIS_Data_20190305flags.xlsx'))
flag_vars = ['sample_id', 'method', 'data_qc_code']
bio_flags = bio_flags.loc[:,flag_vars]
bio_flags = bio_flags.pivot_table(values='data_qc_code', index='sample_id', columns='method')
bio_flags = bio_flags.rename(columns={'Alkalinity_umol/kg' : 'TAflag'})
bio_flags = bio_flags.rename(columns={'pH_invitro_total_f' : 'pHflag'})
bio_flags = bio_flags.rename(columns={'TIC' : 'TICflag'})
bio_flags.TAflag.update(bio_flags.ALKW_01)
bio_flags.TICflag.update(bio_flags.TICW_01)
bio_flags.pHflag.update(bio_flags.pH_invitro_total_p)
bio_var = ['TAflag', 'TICflag', 'pHflag']
bio_flags = bio_flags.loc[:,bio_var]

# Merge flags
bio_flags.reset_index(inplace=True)
df_MAR['sample_id'] = df_MAR['collector_event_id'] +"_"+ df_MAR['collector_sample_id'].map(str)
df_MAR = df_MAR.merge(bio_flags, on='sample_id', how='left')

# Only select desired columns
variables = ['TripID', 'Region', 'StationID', 'timestamp', 'latitude', 'longitude', 'depth', 'NO3', 'O2', 'O2sat_perc', 'TA', 'temp_pH', 'pH_25', 'PO4', 'salinity', 'SiO', 'temperature', 'TIC', 'TAflag', 'pHflag', 'TICflag']
df_MAR = df_MAR.loc[:,variables]

# Drop GSL data > 2014 (in 2014, data processed at BIO, the rest at IML).
dropGSL = df_MAR[(df_MAR['Region'] == 'GSL') & (df_MAR['timestamp'].dt.year > 2014)].index
df_MAR.drop(dropGSL, inplace=True)
# Attribute CSL to GSL.
#df_MAR.loc[(df_MAR['StationID'].str.contains('CSL', na=False)), 'Region'] = 'GSL'

# rename the IML station ID to match rest of IML dataset using 2014 Dataset
# It's a bit complicated because depth does not match when rounded
df2014f = pd.read_excel(os.path.join(dataset_path2,'AZMP_OA_IML2014f.xlsx'), encoding='utf-8')
GSL_orig = df_MAR[df_MAR.Region=='GSL']
GSL_orig.dropna(subset=['TA', 'TIC', 'pH_25'], axis=0, thresh=2, inplace=True) 
GSL_2014 = df2014f[['Station-ID', 'Latitude', 'Longitude']] 
GSL_2014 = GSL_2014.rename(columns={'Latitude' : 'latitude'})
GSL_2014 = GSL_2014.rename(columns={'Longitude' : 'longitude'})
GSL_2014 = GSL_2014.rename(columns={'Station-ID' : 'StationID_new'})
GSL_orig.latitude = GSL_orig.latitude.round(4)
GSL_orig.longitude = GSL_orig.longitude.round(4)
GSL_2014.latitude = GSL_2014.latitude.round(4)
GSL_2014.longitude = GSL_2014.longitude.round(4)
# Station correspondance
stn_corresp = GSL_orig.merge(GSL_2014, on=['latitude', 'longitude'], how='left')
stn_corresp = stn_corresp.loc[:,['StationID', 'StationID_new']].drop_duplicates() 
# Updated DataFrame
GSL_update = GSL_orig.merge(stn_corresp, on='StationID', how='left')
GSL_update.drop(columns=['StationID'], inplace=True)   
GSL_update = GSL_update.rename(columns={'StationID_new' : 'StationID'})
GSL_update.index = GSL_orig.index
# Replace in df_MAR
df_MAR.iloc[GSL_update.index] = GSL_update

# Remove the Labrador Sea AR7W Line
df_MAR=df_MAR[~df_MAR.StationID.str.contains('L3')]

## 1.3 --- BIO's Spring 2019
# Load data (no nutrients, no oxygen)
MAR2019 = pd.read_excel(os.path.join(dataset_path, '2019_Spring_AZMP_Mar_TIC-TA.xlsx'), parse_dates = {'timestamp' : [7, 8]})
# swap weird flags
MAR2019 = MAR2019.replace('I23', 3)
MAR2019 = MAR2019.replace('I22', 3)


# add oxygen
df_o2 = pd.read_csv(os.path.join(dataset_path, 'WinklerOxygen_COR2019001_Final.dat'), header=9, sep=',')
# rename sample no.
sampleNo = df_o2.Sample
sampleNo = sampleNo.map(lambda x: x.replace('_1', ''))
sampleNo = sampleNo.map(lambda x: x.replace('_2', ''))
df_o2.Sample = sampleNo
df_o2 = df_o2.groupby('Sample').mean() 
df_o2 = df_o2['O2_Concentration(ml/l)']
df_o2 = df_o2.reset_index()
df_o2 = df_o2.rename(columns={'Sample' : 'sample_id'})
df_o2 = df_o2.astype({'sample_id':'int'})
# merge O2 into carbonate data
MAR_index = MAR2019.index
MAR2019 = MAR2019.merge(df_o2, on='sample_id', how='left')
MAR2019.index = MAR_index

# add nutrients
df_nut = pd.read_excel(os.path.join(dataset_path, 'COR2019-001_Spring_AZMP_NUTS_Final_JB.xls'), header=1)
df_nut = df_nut.replace('contaminated', np.nan)
df_nut = df_nut.replace('DL', np.nan)
# rename second rows (duplicates)
df_nut['SAMPLE ID'] = df_nut.groupby(np.arange(len(df_nut)) // 2)['SAMPLE ID'].transform('mean')
# average duplicates
df_nut = df_nut.groupby('SAMPLE ID').mean() 
# add nitrate nitrite and rename columns
df_nut['NO3'] = df_nut['NITRATE'] + df_nut['NITRITE'] 
df_nut = df_nut.rename(columns={'PHOSPHATE' : 'PO4'})
df_nut = df_nut.rename(columns={'SILICATE' : 'SiO'})     
# reset index and rename
df_nut = df_nut.reset_index()
df_nut = df_nut.rename(columns={'SAMPLE ID' : 'sample_id'})
df_nut = df_nut.astype({'sample_id':'int'})
# merge O2 into carbonate data
df_nut = df_nut[['sample_id', 'PO4', 'SiO', 'NO3']]
MAR_index = MAR2019.index
MAR2019 = MAR2019.merge(df_nut, on='sample_id', how='left')
MAR2019.index = MAR_index

# Rename columns
MAR2019 = MAR2019.rename(columns={'Station' : 'StationID'})
MAR2019 = MAR2019.rename(columns={'PrDM' : 'depth'})
MAR2019 = MAR2019.rename(columns={'T068C' : 'temperature'})
MAR2019 = MAR2019.rename(columns={'Sal00' : 'salinity'})
MAR2019 = MAR2019.rename(columns={'TA (umol/kg)' : 'TA'})
MAR2019 = MAR2019.rename(columns={'TA flag' : 'TAflag'})
MAR2019 = MAR2019.rename(columns={'TIC (umol/kg)' : 'TIC'})
MAR2019 = MAR2019.rename(columns={'TIC flag' : 'TICflag'})
MAR2019 = MAR2019.rename(columns={'O2_Concentration(ml/l)' : 'O2'})
MAR2019 = MAR2019.rename(columns={'cruise_number' : 'TripID'})
 
# select columns
MAR2019 = MAR2019.loc[:,variables]
MAR2019['Region'] ='MAR'
# append to df_MAR
df_MAR = df_MAR.append(MAR2019, ignore_index=True)


## 2. --- NL data (main NL file. Olivia may have changed the headers)
# This one is not used (updated at 2.2 below)
df_NL = pd.read_csv(os.path.join(dataset_path2,'AZMP_OA_NL.csv'), delimiter=',', encoding='ISO-8859-1')
df_NL['timestamp'] = df_NL.Date.astype(str) + ' ' + df_NL.GMT.astype(str)
df_NL['timestamp'] =  pd.to_datetime(df_NL['timestamp'], format='%m/%d/%Y %H:%M:%S')
df_NL.index = df_NL.timestamp

# Add more DO data (3245 to 3571 good data) Because it was provided late.
df_NL2018DO = pd.read_excel(os.path.join(dataset_path2, 'NL_2018_satO2_perc.xlsx'))
df_NL2018DO['timestamp'] = df_NL2018DO.Date.astype(str) + ' ' + df_NL2018DO.GMT.astype(str)
df_NL2018DO['timestamp'] =  pd.to_datetime(df_NL2018DO['timestamp'], format='%Y/%m/%d %H:%M:%S')
df_NL2018DO.index = df_NL2018DO.timestamp
df_NL['Dissolved Oxygen (mL/L)'].update(df_NL2018DO['Dissolved Oxygen (mL/L)'])

# Rename columns
df_NL = df_NL.rename(columns={'Ship_Trip' : 'TripID'})
df_NL = df_NL.rename(columns={'Latitude (Decimal Degrees)' : 'latitude'})
df_NL = df_NL.rename(columns={'Longiude (Decimal Degrees)' : 'longitude'})
df_NL = df_NL.rename(columns={'Station Name' : 'StationID'})
df_NL = df_NL.rename(columns={'Pressure (dbar)' : 'depth'})
df_NL = df_NL.rename(columns={'Salinity (psu)' : 'salinity'})
df_NL = df_NL.rename(columns={'Temperature (C)' : 'temperature'})
df_NL = df_NL.rename(columns={'Dissolved Oxygen (mL/L)' : 'O2'})
df_NL = df_NL.rename(columns={'Calibrated Oxygen (mL/L)' : 'oxygen'})
df_NL = df_NL.rename(columns={'Phosphate Concentration (mmol m-3)' : 'PO4'})
df_NL = df_NL.rename(columns={'Silicate Concentration (mmol m-3)' : 'SiO'})
df_NL = df_NL.rename(columns={'Nitrate Concentration (mmol m-3)' : 'NO3'})
df_NL = df_NL.rename(columns={'pH 25 ' : 'pH 25'})
df_NL = df_NL.rename(columns={'Discrete Chlorophyll a (mg m-3)' : 'chla'})
df_NL = df_NL.rename(columns={'Total P (mmol m-3)' : 'PO4'})
df_NL = df_NL.rename(columns={'Total Si (mmol m-3)' : 'SiO'})
df_NL = df_NL.rename(columns={'Total NO2/3 (mmol m-3)' : 'NO3'})
df_NL = df_NL.rename(columns={'Measured TA (micromol/kg)' : 'TA'})
df_NL = df_NL.rename(columns={'Measure TCO2 (micromol/kg)' : 'TIC'})
df_NL = df_NL.rename(columns={'Density (kg m3)' : 'density'})
df_NL = df_NL.rename(columns={'TA Data Flag' : 'TAflag'})
df_NL = df_NL.rename(columns={'TCO2 Data Flag' : 'TICflag'})
df_NL['Region'] = 'NL'
df_NL['temp_pH']= np.NaN
df_NL['pHflag']= np.NaN

# reset index
df_NL = df_NL.reset_index(drop=True)

# Only select subset variables 
df_NL = df_NL.loc[:,variables]

# 2.2 Updated NL data
df_NLu = pd.read_excel(os.path.join(dataset_main_path, 'Ocean_Carbon_NL.xlsx'))
df_NLu['timestamp'] = pd.to_datetime(df_NLu.Date.astype(str)) + pd.to_timedelta(df_NLu.GMT.astype(str)) 
df_NLu.index = df_NLu.timestamp
# Rename columns
df_NLu = df_NLu.rename(columns={'Ship_Trip' : 'TripID'})
df_NLu = df_NLu.rename(columns={'Station Name' : 'StationID'})
df_NLu = df_NLu.rename(columns={'Latitude (Decimal Degrees)' : 'latitude'})
df_NLu = df_NLu.rename(columns={'Longitude (Decimal Degrees)' : 'longitude'})
df_NLu = df_NLu.rename(columns={'Station Name' : 'StationID'})
df_NLu = df_NLu.rename(columns={'Pressure (dbar)' : 'depth'})
df_NLu = df_NLu.rename(columns={'Salinity (psu)' : 'salinity'})
df_NLu = df_NLu.rename(columns={'Temperature (°C)' : 'temperature'})
df_NLu = df_NLu.rename(columns={'Dissolved Oxygen (mL/L)' : 'O2'})
df_NLu = df_NLu.rename(columns={'Calibrated Oxygen (mL/L)' : 'oxygen'})
df_NLu = df_NLu.rename(columns={'Phosphate Concentration (mmol m-3)' : 'PO4'})
df_NLu = df_NLu.rename(columns={'Silicate Concentration (mmol m-3)' : 'SiO'})
df_NLu = df_NLu.rename(columns={'Nitrate Concentration (mmol m-3)' : 'NO3'})
#df_NLu = df_NLu.rename(columns={'pH 25 ' : 'pH 25'})
df_NLu = df_NLu.rename(columns={'Discrete Chlorophyll a (mg m-3)' : 'chla'})
df_NLu = df_NLu.rename(columns={'Total P (mmol m-3)' : 'PO4'})
df_NLu = df_NLu.rename(columns={'Total Si (mmol m-3)' : 'SiO'})
df_NLu = df_NLu.rename(columns={'Total NO2/3 (mmol m-3)' : 'NO3'})
df_NLu = df_NLu.rename(columns={'Measured TA (µmol/kg)' : 'TA'})
df_NLu = df_NLu.rename(columns={'Measured TCO2 (µmol/kg)' : 'TIC'})
df_NLu = df_NLu.rename(columns={'Density (kg m3)' : 'density'})
df_NLu = df_NLu.rename(columns={'TA Data Flag' : 'TAflag'})
df_NLu = df_NLu.rename(columns={'TCO2 Data Flag' : 'TICflag'})
df_NLu['Region'] = 'NL'
df_NLu['temp_pH']= np.NaN
df_NLu['pHflag']= np.NaN

# Only select subset variables 
df_NLu = df_NLu.loc[:,variables]

## 3. --- IML data
# 2014 very different format, thus will be taken from Biochem (see above)
df2015 = pd.read_excel(os.path.join(dataset_path2, 'AZMP_OA_IML2015.xlsx'), encoding='utf-8')
df2015['CTD Station (nom/no)'] = df2015['CTD Station (nom/no)'].replace(7, 'TESL2')
df2015['CTD Station (nom/no)'] = df2015['CTD Station (nom/no)'].replace(16, 'TSI3')
df2015['CTD Station (nom/no)'] = df2015['CTD Station (nom/no)'].replace(23, 'TASO3')
df2015['CTD Station (nom/no)'] = df2015['CTD Station (nom/no)'].replace(28, 'CMO2/CH9')
df2015['CTD Station (nom/no)'] = df2015['CTD Station (nom/no)'].replace(39, 'CH2')
df2015['CTD Station (nom/no)'] = df2015['CTD Station (nom/no)'].replace(51, 'TCEN3')
df2015['CTD Station (nom/no)'] = df2015['CTD Station (nom/no)'].replace(56, 'CSL_04')
df2015['CTD Station (nom/no)'] = df2015['CTD Station (nom/no)'].replace(61, 'TIDM3')
df2015['CTD Station (nom/no)'] = df2015['CTD Station (nom/no)'].replace(67, 'TIDM9')
df2015f = pd.read_excel(os.path.join(dataset_path2,'AZMP_OA_IML2015f.xlsx'), encoding='utf-8')
df2016 = pd.read_excel(os.path.join(dataset_path2,'AZMP_OA_IML2016.xlsx'), encoding='utf-8')
df2016f = pd.read_excel(os.path.join(dataset_path2,'AZMP_OA_IML2016f.xlsx'), encoding='utf-8')
df2017 = pd.read_excel(os.path.join(dataset_path2,'AZMP_OA_IML2017.xlsx'), encoding='utf-8')
df2017f = pd.read_excel(os.path.join(dataset_path2,'AZMP_OA_IML2017f.xlsx'), encoding='utf-8')
df2017s = pd.read_excel(os.path.join(dataset_path2,'AZMP_OA_IML2017s.xlsx'), encoding='utf-8') ###some items in CTD Fichier were changed for an easier 'TripID'
df2018 = pd.read_excel(os.path.join(dataset_path2,'AZMP_OA_IML2018.xlsx'), encoding='utf-8')
df2018s = pd.read_excel(os.path.join(dataset_path2,'AZMP_OA_IML2018s.xlsx'), encoding='utf-8') ###some items in CTD Fichier were changed for an easier 'TripID'
df2018s = df2018s.rename(columns={'CTD Mission (nom)' : 'Mission  (nom)'})
df2018f = pd.read_excel(os.path.join(dataset_path2,'AZMP_OA_IML2018f.xlsx'), encoding='utf-8')

# 2019 data
xls = pd.ExcelFile(os.path.join(dataset_path,'IMLSpring and Fall 2019 (Gibb).xlsx'))
df2019s = pd.read_excel(xls, 'June 2019-2', header=1)
df2019f = pd.read_excel(xls, 'Fall 2019', header=1)
df2019s = df2019s.drop(0)
df2019f = df2019f.drop(0)
# merge 2 sheets
df2019f = df2019f.rename(columns={'PSAL_BS' : 'psal_bs'})
df2019f = df2019f.rename(columns={'oxy_02' : 'OXY_02'})
df2019f = df2019f.rename(columns={'Q_oxy' : 'Q_OXY'})
df2019f = df2019f.rename(columns={'oxy_02.1' : 'OXY_02.1'})
df2019f = df2019f.rename(columns={'Q_oxy.1' : 'Q_OXY.1'})
df2019 = pd.concat([df2019s, df2019f])
# clean un-necessary columns
df2019.Date = pd.to_datetime(df2019.Date)
df2019 = df2019.drop(columns=['Unnamed: 49'])
df2019 = df2019.drop(columns=['Unnamed: 50'])
# Rename columns because 2019 was provided in a different format.
df2019 = df2019.rename(columns={'Unnamed: 0' : 'Mission  (nom)'})
df2019 = df2019.rename(columns={'Fichier' : 'CTD Fichier (nom)'})
df2019 = df2019.rename(columns={'Station' : 'CTD Station (nom/no)'})
df2019 = df2019.rename(columns={'Latitude' : 'CTD Latitude (degres)'})
df2019 = df2019.rename(columns={'Longitude' : 'CTD Longitude (degres)'})
df2019 = df2019.rename(columns={'Date' : 'CTD Date (jj-mmm-yyyy)'})
df2019 = df2019.rename(columns={'Heure' : 'CTD Heure (GMT)'})
df2019 = df2019.rename(columns={'PRES' : 'CTD PRES (dbar)'})
df2019 = df2019.rename(columns={'PRES_SDEV' : 'CTD PRES_SDEV (dbar)'})
df2019 = df2019.rename(columns={'TE90' : 'CTD TE90 (celsius)'})
df2019 = df2019.rename(columns={'TE90_SDEV' : 'CTD TE90_SDEV (celsius)'})
df2019 = df2019.rename(columns={'PSAL' : 'CTD PSAL (psu)'})
df2019 = df2019.rename(columns={'SIGT' : 'CTD SIGT (kg/m**3)'})
df2019 = df2019.rename(columns={'FLOR' : 'CTD FLOR (mg/m**3)'})
df2019 = df2019.rename(columns={'DOXY' : 'CTD DOXY (ml/l)'})
df2019 = df2019.rename(columns={'PHPH' : 'CTD PHPH (NBS scale)'})
df2019 = df2019.rename(columns={'PHPH_T' : 'CTD PHPH (Total scale)'})
df2019 = df2019.rename(columns={'psal_bs' : 'labo psal_bs (PSU)'})
df2019 = df2019.rename(columns={'Q_PSAL' : 'labo Q_PSAL ((none))'})
df2019 = df2019.rename(columns={'OXY_02' : 'labo oxy_02 (ml/l)'})
df2019 = df2019.rename(columns={'Q_OXY' : 'labo Q_OXY ((none))'})
df2019 = df2019.rename(columns={'OXY_02.1' : 'labo oxy_02 (ml/l).1'})
df2019 = df2019.rename(columns={'Q_OXY.1' : 'labo Q_OXY ((none)).1'})
df2019 = df2019.rename(columns={'NO2_03' : 'labo NO2_03 (mmol/m**3)'})
df2019 = df2019.rename(columns={'Q_NO2' : 'labo Q_NO2 ((none))'})
df2019 = df2019.rename(columns={'NO2_03.1' : 'labo NO2_03 (mmol/m**3).1'})
df2019 = df2019.rename(columns={'Q_NO2.1' : 'labo Q_NO2 ((none)).1'})
df2019 = df2019.rename(columns={'NO2_03.2' : 'labo NO2_03 (mmol/m**3).2'})
df2019 = df2019.rename(columns={'Q_NO2.2' : 'labo Q_NO2 ((none)).2'})
df2019 = df2019.rename(columns={'NOX_03' : 'labo NOX_03 (mmol/m**3)'})
df2019 = df2019.rename(columns={'Q_NOX' : 'labo Q_NOX ((none))'})
df2019 = df2019.rename(columns={'NOX_03.1' : 'labo NOX_03 (mmol/m**3).1'})
df2019 = df2019.rename(columns={'Q_NOX.1' : 'labo Q_NOX ((none)).1'})
df2019 = df2019.rename(columns={'NOX_03.2' : 'labo NOX_03 (mmol/m**3).2'})
df2019 = df2019.rename(columns={'Q_NOX.2' : 'labo Q_NOX ((none)).2'})
df2019 = df2019.rename(columns={'PO4_03' : 'labo PO4_03 (mmol/m**3)'})
df2019 = df2019.rename(columns={'Q_PO4' : 'labo Q_PO4 ((none))'})
df2019 = df2019.rename(columns={'PO4_03.1' : 'labo PO4_03 (mmol/m**3).1'})
df2019 = df2019.rename(columns={'Q_PO4.1' : 'labo Q_PO4 ((none)).1'})
df2019 = df2019.rename(columns={'PO4_03.2' : 'labo PO4_03 (mmol/m**3).2'})
df2019 = df2019.rename(columns={'Q_PO4.2' : 'labo Q_PO4 ((none)).2'})
df2019 = df2019.rename(columns={'Si_03' : 'labo Si_03 (mmol/m**3)'})
df2019 = df2019.rename(columns={'Q_Si' : 'labo Q_Si ((none))'})
df2019 = df2019.rename(columns={'Si_03.1' : 'labo Si_03 (mmol/m**3).1'})
df2019 = df2019.rename(columns={'Q_Si.1' : 'labo Q_Si ((none)).1'})
df2019 = df2019.rename(columns={'LBPHT_02' : 'labo LBPHT_02 (Total scale)'})
df2019 = df2019.rename(columns={'pHT_02' : 'labo pHT_02 (Total scale)'})
df2019 = df2019.rename(columns={'pHT_02.1' : 'labo pHT_02 (Total scale).1'})
df2019 = df2019.rename(columns={'ALKW_01' : 'labo_ALKW_01_(umol/kg**1)'}) # only one alk, so same as 'At' below
df2019 = df2019.rename(columns={'TICW_01' : 'labo_TICW_01_(umol/kg**1)'}) # only since fall 2019
df2019 = df2019.rename(columns={'At' : 'Moyenne At (umol/kg)'})
df2019 = df2019.drop(columns=['At.1']) # N
df2019 = df2019.drop(columns=['At.2']) # empty
df2019 = df2019.rename(columns={'At.3' : 'Flag At '})
df2019 = df2019.rename(columns={'LABT_01' : 'labo LABT_01 (deg C)'})
df2019 = df2019.rename(columns={'pH  labo' : 'Moyenne  pH total in situ'})
df2019 = df2019.drop(columns=['pH  labo.1']) # N
df2019 = df2019.drop(columns=['pH  labo.2']) # empty
df2019 = df2019.rename(columns={'pH  labo.3' : 'Flag pH total in situ'})
df2019 = df2019.rename(columns={'PO4_03.3' : 'Mean PO4_03 (mmol/m**3)'})
df2019 = df2019.rename(columns={'PO4_03.4' : 'Flag PO4_03 '})
df2019 = df2019.rename(columns={'Si_03.2' : 'Mean Si_03 (mmol/m**3)'})
df2019 = df2019.rename(columns={'Si_03.3' : 'Flag Si_03 '})
df2019 = df2019.drop(columns=['pH in situ'])
df2019 = df2019.rename(columns={'Unnamed: 65' : ' WCa out '})
df2019 = df2019.rename(columns={'Unnamed: 66' : ' WAr out '})
df2019 = df2019.rename(columns={'Unnamed: 67' : ' Strate '})

# Station Riki
xls_riki = pd.ExcelFile(os.path.join(dataset_path,'StationRimouski (2013-2018)(GIBBB).xlsx'))
riki2013 = pd.read_excel(xls_riki, 'Rimouski 2013', header=1)
riki2014 = pd.read_excel(xls_riki, 'Rimouski 2014', header=1)
riki2015 = pd.read_excel(xls_riki, 'Rimouski 2015', header=1)
riki2017 = pd.read_excel(xls_riki, 'Rimouski 2017', header=1)
riki2018 = pd.read_excel(xls_riki, 'Rimouski 2018', header=1)


df2019s = df2019s.drop(0)
df2019f = df2019f.drop(0)
# merge 2 sheets
df2019f = df2019f.rename(columns={'PSAL_BS' : 'psal_bs'})
df2019f = df2019f.rename(columns={'oxy_02' : 'OXY_02'})
df2019f = df2019f.rename(columns={'Q_oxy' : 'Q_OXY'})
df2019f = df2019f.rename(columns={'oxy_02.1' : 'OXY_02.1'})
df2019f = df2019f.rename(columns={'Q_oxy.1' : 'Q_OXY.1'})
df2019 = pd.concat([df2019s, df2019f])
# clean un-necessary columns

# Merge all 
df_IML = pd.concat([df2019, df2018f, df2018s, df2018, df2017f, df2017s, df2017, df2016f, df2016, df2015, df2015f], axis=0, sort=False)
df_IML = df_IML.reset_index()

# Clean dataframe
cols = df_IML.columns
cols = cols.map(lambda x: x.replace(' ', '_') if isinstance(x, (str)) else x)
df_IML.columns = cols
df_IML = df_IML.replace(r'^\s+$', np.nan, regex=True)
df_IML = df_IML.rename(columns={'CTD_Date_(jj-mmm-yyyy)' : 'Dates'})
df_IML = df_IML.rename(columns={'CTD_Heure_(GMT)' : 'Times'})
df_IML['timestamp'] = df_IML.Dates.astype(str) + ' ' + df_IML.Times.astype(str)
df_IML['timestamp'] =  pd.to_datetime(df_IML['timestamp'], format='%Y/%m/%d %H:%M:%S')
df_IML = df_IML.rename(columns={'CTD_Fichier_(nom)' : 'TripID'})
df_IML['TripID'] = df_IML['TripID'].replace('-','', regex=True)
df_IML.loc[(df_IML['TripID'].str.contains('iml|TEL|PER', na=False)), 'TripID'] = df_IML['TripID'].str[:10]

# Rename columns
df_IML = df_IML.rename(columns={'CTD_Station_(nom/no)' : 'StationID'})
df_IML = df_IML.rename(columns={'CTD_Longitude_(degres)' : 'longitude'})
df_IML = df_IML.rename(columns={'CTD_Latitude_(degres)' : 'latitude'})
df_IML = df_IML.rename(columns={'CTD_PRES_(dbar)' : 'depth'})
df_IML = df_IML.rename(columns={'CTD_TE90_(celsius)' : 'temperature'})
df_IML = df_IML.rename(columns={'CTD_PSAL_(psu)' : 'salinity'})
df_IML = df_IML.rename(columns={'CTD_DOXY_(ml/l)' : 'oxyCTD'})
df_IML = df_IML.rename(columns={'labo_LABT_01_(deg_C)' : 'temp_pH'})
df_IML = df_IML.rename(columns={'Flag_At_' : 'TAflag'})
df_IML = df_IML.rename(columns={'Flag_pH_total_in_situ' : 'pHflag'})

# Rename StationIDs (to be consistent)
df_IML['StationID'] = df_IML['StationID'].replace('CMO3 (CH12)', 'CMO3/CH12')
df_IML['StationID'] = df_IML['StationID'].replace('TESL3 (RIKI)', 'TESL3')
df_IML['StationID'] = df_IML['StationID'].replace('TSI1 (CG)', 'TSI1')
df_IML['StationID'] = df_IML['StationID'].replace('TSI1/Courant_Gaspe', 'TSI1')
df_IML['StationID'] = df_IML['StationID'].replace('TIDM7/6.4', 'TIDM7')
df_IML['StationID'] = df_IML['StationID'].replace('TIDM5(4.5)', 'TIDM5')
df_IML['StationID'] = df_IML['StationID'].replace('TIDM5/4.5', 'TIDM5')
df_IML['StationID'] = df_IML['StationID'].replace('TIDM2(2.2)', 'TIDM2')
df_IML['StationID'] = df_IML['StationID'].replace('CMO1 (CH1)', 'CMO1/CH1')
df_IML['StationID'] = df_IML['StationID'].replace('CMO2 (CH9)', 'CMO2/CH9')
df_IML['StationID'] = df_IML['StationID'].replace('TIDM9 (Shediac)', 'TIDM9')
df_IML['StationID'] = df_IML['StationID'].replace('TIDM9,SHEDIAC', 'TIDM9')
df_IML['StationID'] = df_IML['StationID'].replace('TIDM9/Shediac', 'TIDM9')
df_IML['StationID'] = df_IML['StationID'].replace('TIDM9ShediacValley', 'TIDM9')
df_IML['StationID'] = df_IML['StationID'].replace('CMO1,CH1', 'CMO1/CH1')
df_IML['StationID'] = df_IML['StationID'].replace('TESL3,Riki', 'TESL3')
df_IML['StationID'] = df_IML['StationID'].replace('TESL3/Riki', 'TESL3')
df_IML['StationID'] = df_IML['StationID'].replace('ES2(13ML)', 'ES2/13ML')
df_IML['StationID'] = df_IML['StationID'].replace('If37', 'IF37')

# Create region tag
df_IML['Region'] = 'GSL'

# Calculate mean values of replicate analyses
df_IML['NO3'] = df_IML[['labo_NOX_03_(mmol/m**3)', 'labo_NOX_03_(mmol/m**3).1', 'labo_NOX_03_(mmol/m**3).2', 'labo_NOx_03_(mmol/m**3)', 'labo_NOx_03_(mmol/m**3).1', 'labo_NOx_03_(mmol/m**3).2']].mean(axis=1)
df_IML['PO4'] = df_IML[['labo_PO4_03_(mmol/m**3)', 'labo_PO4_03_(mmol/m**3).1', 'labo_PO4_03_(mmol/m**3).2']].mean(axis=1)
df_IML['SiO'] = df_IML[['labo_Si_03_(mmol/m**3)', 'labo_Si_03_(mmol/m**3).1', 'labo_Si_03_(mmol/m**3).2']].mean(axis=1)
df_IML['TA'] = df_IML[['labo_ALKW_01_(umol/kg**1)', 'labo_ALKW_01_(umol/kg**1).1']].mean(axis=1)
df_IML['TIC'] = df_IML[['labo_TICW_01_(umol/kg**1)']].mean(axis=1)
df_IML['pH_25'] = df_IML[['labo_LBPHT_02_(Total_scale)', 'labo_LBPHT_02_(Total_scale).1']].mean(axis=1)
df_IML['O2'] = df_IML[['labo_oxy_02_(ml/l)', 'labo_oxy_02_(ml/l).1', 'labo_OXY_02_(ml/l)', 'labo_OXY_02_(ml/l).1']].mean(axis=1)

# Fill missing depths (with zbouteille)
idx_missing_depth = df_IML.loc[df_IML.depth.isna()].index
df_IML.loc[idx_missing_depth, 'depth'] = df_IML.loc[idx_missing_depth, 'CTD_zbouteille_(dbar)']

# Only select subset variables 
df_IML = df_IML.loc[:,variables]

## 4. Merged dataset 
df = pd.concat([df_MAR, df_NLu, df_IML], axis=0)
df = df.reset_index(drop=True)
satO2 = swx.satO2(df['salinity'], df['temperature']) 
df['O2sat_perc'] = df['O2']/satO2*100 # calculate oxygen saturation %
#df['AOU'] = df['O2']-df['satO2']
df['NO3']=df['NO3']/1.025 # convert nutrient data from uM=mmol/m3/umol/L to umol/kgSW
df['SiO']=df['SiO']/1.025
df['PO4']=df['PO4']/1.025
dforig = df.copy(deep=True) # all original data, since some is replaced/modified
# ** This one could be exported for dataset paper.

# **** Note that out of total 22422 samples, 14346 have pH, TA and TIC NaNs...
# I think we should get rid of them now (to same dforig above).
df.drop(df.loc[(df.pH_25.isna()) & (df.TA.isna()) & (df.TIC.isna())].index, inplace=True)


###########################################
## ------ Prepare data for CO2sys ------ ##
###########################################
# (estimate TA from TA-S plot, remove empty rows, remove flagged data, replace nutrient nan with zeros)

# HERE!!!!!n problems with flag: df.TAflag.unique()  -> array([0.0, nan, 'I23', 'I22', 1, 3, 8, 7.0], dtype=object)
# try:
#    df[df.TAflag=='I23'].values    
# make sure not str in columns
df = df.astype({'TAflag':'float', 'pHflag':'float', 'TICflag':'float'})

# Drop flagged data
df.drop(df.loc[df['TAflag']>=3].index, inplace=True)
df.drop(df.loc[df['pHflag']>=3].index, inplace=True)
df.drop(df.loc[df['TICflag']>=3].index, inplace=True)

# Calculate missing TA values based on TA-S relationship
dflinreg = df.copy(deep=True)
dflinreg = dflinreg[dflinreg['Region'].str.contains("GSL")]
dflinreg = dflinreg[(dflinreg.timestamp.dt.month>=3) & (dflinreg.timestamp.dt.month<=6)]
dflinreg.dropna(subset=['TA', 'salinity'], inplace=True)
var_x = 'salinity'
var_y = 'TA'
# plot
linreg = sns.lmplot(var_x, var_y, dflinreg, legend=False, fit_reg=True, ci=95, scatter_kws={'s':20, 'alpha':1}, line_kws={'lw':1}, hue='Region', palette='plasma', height=4, aspect=1.1)
lines = sns.regplot(x=var_x, y=var_y, data=dflinreg, scatter=False, ci=95, color='k', line_kws={'lw': 1})
linreg.set_xlabels(r'Salinity')
linreg.set_ylabels(r'TA $(\rm \mu $mol/kg)')
axes = linreg.axes
axes[0,0].set_xlim(17, 38)
axes[0,0].set_ylim(1700, 2500)
sns.set_style('whitegrid')
sns.despine(top=False, right=False, left=False, bottom=False)   
plt.plot([-2,22], [400,400,], 'black', linewidth=2, linestyle='dashed')
# Regression coeff
x=dflinreg[var_x].values
y=dflinreg[var_y].values
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
T = stats.t.ppf(1-0.025, len(dflinreg.salinity)-2)
ci = std_err*T
print("linear regression line y=",slope, "x +",intercept)
print("r squared: ", r_value**2)
print("p_value: ", p_value)
print("standard error: ", std_err)
print("95% confidence interval: ", ci)
print('')

# Replace (*note that it's okay to do df.salinity since associated to index)
df.loc[((df.timestamp.dt.year==2015) & (df.timestamp.dt.month==6)) & (df['Region'].str.contains('GSL')), 'TA'] = ((slope * df.salinity)+ intercept)


###########################################
## ------------ Apply CO2sys ----------- ##
###########################################
print('Use PyCO2SYS v{}'.format(version))

# Keep only indices with at least 2 parameters
df.dropna(subset=['TA', 'TIC', 'pH_25'], axis=0, thresh=2, inplace=True)

# replace NaN in nutrient by zero (CO2sys default)
no_nan = df.SiO[(df.PO4.isna()) | (df.SiO.isna())].size
perc_nan = float(no_nan)/df.shape[0]*100.0
df.PO4 = df.PO4.replace(np.nan, 0)
df.SiO = df.SiO.replace(np.nan, 0)
print(str(no_nan) + ' out of ' + str(df.shape[0]) + ' nutrient data where replaced by zeros (' + str(perc_nan) + ' %)')

# Isolate index with 3 parameters
df_tmp = df.copy()
df_TA_TIC_pH = df_tmp.dropna(subset=['TA', 'TIC', 'pH_25'], thresh=3, axis=0)
df_tmp.drop(df_TA_TIC_pH.index, inplace=True)
# Isolate index with TA/TIC
df_TA_TIC = df_tmp.dropna(subset=['TA', 'TIC'], thresh=2, axis=0)
df_tmp.drop(df_TA_TIC.index, inplace=True)
# Isolate index with TA/pH
df_TA_pH = df_tmp.dropna(subset=['TA', 'pH_25'], thresh=2, axis=0)
df_tmp.drop(df_TA_pH.index, inplace=True)
# Isolate index with TIC/pH (only one!? TA flaged?)
df_TIC_pH = df_tmp.dropna(subset=['TIC', 'pH_25'], thresh=2, axis=0)
df_tmp.drop(df_TIC_pH.index, inplace=True)
if df_tmp.size:
    print('Might have problem, there are left overs for CO2sys!!')

# apply CO2sys (TA/TIC/pH)
CO2dict = CO2SYS(df_TA_TIC_pH.TA, df_TA_TIC_pH.pH_25, 1, 3, df_TA_TIC_pH.salinity, 20, df_TA_TIC_pH.temperature, 0, df_TA_TIC_pH.depth, df_TA_TIC_pH.SiO, df_TA_TIC_pH.PO4, 1, 4, 1)
co2sys_TA_TIC_pH = pd.DataFrame.from_dict(CO2dict)
co2sys_TA_TIC_pH.index = df_TA_TIC_pH.index
del CO2dict    
# apply CO2sys (TA/TIC)
CO2dict = CO2SYS(df_TA_TIC.TA, df_TA_TIC.TIC, 1, 2, df_TA_TIC.salinity, 20, df_TA_TIC.temperature, 0, df_TA_TIC.depth, df_TA_TIC.SiO, df_TA_TIC.PO4, 1, 4, 1)
co2sys_TA_TIC = pd.DataFrame.from_dict(CO2dict)
co2sys_TA_TIC.index = df_TA_TIC.index
del CO2dict
# apply CO2sys (TA/pH)
CO2dict = CO2SYS(df_TA_pH.TA, df_TA_pH.pH_25, 1, 3, df_TA_pH.salinity, 20, df_TA_pH.temperature, 0, df_TA_pH.depth, df_TA_pH.SiO, df_TA_pH.PO4, 1, 4, 1)
co2sys_TA_pH = pd.DataFrame.from_dict(CO2dict)
co2sys_TA_pH.index = df_TA_pH.index
del CO2dict
# apply CO2sys (TIC/pH)
CO2dict = CO2SYS(df_TIC_pH.TIC, df_TIC_pH.pH_25, 2, 3, df_TIC_pH.salinity, 20, df_TIC_pH.temperature, 0, df_TIC_pH.depth, df_TIC_pH.SiO, df_TIC_pH.PO4, 1, 4, 1)
co2sys_TIC_pH = pd.DataFrame.from_dict(CO2dict)
co2sys_TIC_pH.index = df_TIC_pH.index
del CO2dict

# This is the merged CO2sys results (#index same as df)
df_co2sys = pd.concat([co2sys_TA_TIC_pH, co2sys_TA_TIC, co2sys_TA_pH, co2sys_TIC_pH], axis=0)

# Add new columns in main DataFrame
df['pH_tot'] = df_co2sys['pHoutTOTAL']
df['TAc'] = df_co2sys['TAlk']
df['TICc'] = df_co2sys['TCO2']
df['Omega_C'] = df_co2sys['OmegaCAout']
df['Omega_A'] = df_co2sys['OmegaARout']
df['pCO2'] = df_co2sys['pCO2out']

################################################
## -------- Final Cleaning and Saving ------- ##
################################################

# convert trip and stations to string
df['TripID'] = df['TripID'].astype(str)
df['StationID'] = df['StationID'].astype(str)

# Update some Trip Names
df['TripID'] = df['TripID'].replace('20114', 'HUD114')
df['TripID'] = df['TripID'].replace('20115', 'HUD115')
df['TripID'] = df['TripID'].replace('20116', 'HUD116')
df['TripID'] = df['TripID'].replace('39144', 'TEL144')
df['TripID'] = df['TripID'].replace('39148', 'TEL148')
df['TripID'] = df['TripID'].replace('39157', 'TEL157')
df['TripID'] = df['TripID'].replace('39159', 'TEL159')
df['TripID'] = df['TripID'].replace('39160', 'TEL160')
df['TripID'] = df['TripID'].replace('39173', 'TEL173')
df['TripID'] = df['TripID'].replace('39176', 'TEL176')
#df['TripID'] = df['TripID'].replace('TEL2017177', 'TEL177')
#df['TripID'] = df['TripID'].replace('TEL2018196', 'TEL196')
df['TripID'] = df['TripID'].replace('15009', 'DIS009')
df['TripID'] = df['TripID'].replace('JC001', 'COO001')

# Update some Trip Names
df['StationID'] = df['StationID'].replace('TESL3(IML4)RIKI', 'TESL3')
df['StationID'] = df['StationID'].replace('CM03 (CH12)', 'CMO3/CH12')
df['StationID'] = df['StationID'].replace('CMO3', 'CMO3/CH12')
df['StationID'] = df['StationID'].replace('CMO2', 'CMO2/CH9')
df['StationID'] = df['StationID'].replace('CMO1/CH1', 'CMO1/CH1')
df['StationID'] = df['StationID'].replace('CMO1', 'CMO1/CH1')
df['StationID'] = df['StationID'].replace('TIDM7 (6.2)', 'TIDM7')
df['StationID'] = df['StationID'].replace('TIDM9 Shediac Valley', 'TIDM9')
df['StationID'] = df['StationID'].replace('TIDM9 (Shédiac)', 'TIDM9')
df['StationID'] = df['StationID'].replace('Tesl-2', 'TESL2')
df['StationID'] = df['StationID'].replace('Tesl-3  (Riki-1)', 'TESL3')
df['StationID'] = df['StationID'].replace('Tesl-6', 'TESL6')
df['StationID'] = df['StationID'].replace('nan', 'unknown')


# Reorder columns
df = df[[
    'timestamp', 'Region', 'TripID', 'StationID',
    'latitude', 'longitude', 'depth',
    'temperature', 'salinity', 'O2',
    'NO3', 'PO4', 'SiO',
    'TAc', 'TICc', 'pH_tot', 'Omega_A', 'Omega_C', 'pCO2', 'O2sat_perc',
    'pH_25', 'temp_pH', 'TA', 'TIC']]

# HERE!!! Check names      
df = df.rename(columns={'TripID' : 'Trip_Name'})
df = df.rename(columns={'StationID' : 'Station_Name'})
df = df.rename(columns={'timestamp' : 'Timestamp'})
df = df.rename(columns={'latitude' : 'Latitude_(degNorth)'})
df = df.rename(columns={'longitude' : 'Longitude_(degEast)'})
df = df.rename(columns={'depth' : 'Depth_(dbar)'})
df = df.rename(columns={'temperature' : 'Temperature_(degC)'})
df = df.rename(columns={'salinity' : 'Salinity_(psu)'})
df = df.rename(columns={'NO3' : 'Nitrate_Concentration_(mmol/m3)'})
df = df.rename(columns={'O2' : 'Dissolved_Oxygen_(mL/L)'})
df = df.rename(columns={'O2sat_perc' : 'Oxygen_Saturation_(%)'})
df = df.rename(columns={'TA' : 'Total_Alkalinity_Measured_(umol/kg)'})
df = df.rename(columns={'temp_pH' : 'pH_lab_temp_(degC)'})
df = df.rename(columns={'pH_25' : 'pH_lab_(--)'})
df = df.rename(columns={'PO4' : 'Phosphate_Concentration_(mmol/m3)'})
df = df.rename(columns={'SiO' : 'Silicate_Concentration_(mmol/m3)'})
df = df.rename(columns={'TIC' : 'Inorganic_Carbon_measured_(umol/kg)'})
df = df.rename(columns={'ph_tot' : 'pH_Total_(total_scale)'})
df = df.rename(columns={'TAc' : 'Total_Alkalinity_(umol/kg)'})
df = df.rename(columns={'TICc' : 'Inorganic_Carbon_(umol/kg)'})
df = df.rename(columns={'Omega_C' : 'Omega_Calcite_(--)'})
df = df.rename(columns={'Omega_A' : 'Omega_Aragonite_(--)'})
df = df.rename(columns={'pCO2' : 'pCO2_(uatm)'})

# Save final dataset
df.to_csv(os.path.join(dataset_main_path, 'AZMP_carbon_data.csv'), float_format='%.4f', index=False)

