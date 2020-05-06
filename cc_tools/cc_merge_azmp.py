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
df_MAR['Region'] = 'SS'
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
variables = ['TripID', 'Region', 'StationID', 'timestamp', 'latitude', 'longitude', 'depth', 'NO3', 'O2', 'TA', 'temp_pH', 'pH_25', 'PO4', 'salinity', 'SiO', 'temperature', 'TIC', 'TAflag', 'pHflag', 'TICflag']
df_MAR = df_MAR.loc[:,variables]

# Drop GSL data > 2014 (HERE! why?)
dropGSL = df_MAR[(df_MAR['Region'] == 'GSL') & (df_MAR['timestamp'].dt.year > 2014)].index
df_MAR.drop(dropGSL, inplace=True)
df_MAR.loc[(df_MAR['StationID'].str.contains('CSL', na=False)), 'Region'] = 'GSL'

# Remove the Labrador Sea AR7W Line
df_MAR=df_MAR[~df_MAR.StationID.str.contains('L3')]

# rename the IML station ID to match rest of IML dataset
df_MAR['StationID'] = df_MAR['StationID'].replace('0421', 'TIDM9')
df_MAR['StationID'] = df_MAR['StationID'].replace('0431', 'TIDM8')
df_MAR['StationID'] = df_MAR['StationID'].replace('0441', 'TIDM7')
df_MAR['StationID'] = df_MAR['StationID'].replace('0451', 'TIDM6')
df_MAR['StationID'] = df_MAR['StationID'].replace('0461', 'TIDM5')
df_MAR['StationID'] = df_MAR['StationID'].replace('0471', 'TIDM4')
df_MAR['StationID'] = df_MAR['StationID'].replace('0481', 'TIDM3')
df_MAR['StationID'] = df_MAR['StationID'].replace('0491', 'TIDM2')
df_MAR['StationID'] = df_MAR['StationID'].replace('0501', 'TIDM1')
df_MAR['StationID'] = df_MAR['StationID'].replace('0541', 'TDC1')
df_MAR['StationID'] = df_MAR['StationID'].replace('0551', 'TDC2')
df_MAR['StationID'] = df_MAR['StationID'].replace('0561', 'TDC3')
df_MAR['StationID'] = df_MAR['StationID'].replace('0571', 'TDC4')
df_MAR['StationID'] = df_MAR['StationID'].replace('0581', 'TDC5')
df_MAR['StationID'] = df_MAR['StationID'].replace('0591', 'TDC6')
df_MAR['StationID'] = df_MAR['StationID'].replace('0621', 'IF28')
df_MAR['StationID'] = df_MAR['StationID'].replace('0631', 'IF27')
df_MAR['StationID'] = df_MAR['StationID'].replace('0641', 'TCEN1')
df_MAR['StationID'] = df_MAR['StationID'].replace('0661', 'TCEN3')
df_MAR['StationID'] = df_MAR['StationID'].replace('0681', 'TCEN5')
df_MAR['StationID'] = df_MAR['StationID'].replace('0701', 'CH6')
df_MAR['StationID'] = df_MAR['StationID'].replace('0721', 'TBB1')
df_MAR['StationID'] = df_MAR['StationID'].replace('0731', 'TBB2')
df_MAR['StationID'] = df_MAR['StationID'].replace('0741', 'TBB3')
df_MAR['StationID'] = df_MAR['StationID'].replace('0761', 'TBB5')
df_MAR['StationID'] = df_MAR['StationID'].replace('0781', 'TBB7')
df_MAR['StationID'] = df_MAR['StationID'].replace('0821', 'IF10')
df_MAR['StationID'] = df_MAR['StationID'].replace('0841', 'IF12')
df_MAR['StationID'] = df_MAR['StationID'].replace('0861', 'CMO1/CH1')
df_MAR['StationID'] = df_MAR['StationID'].replace('0881', 'CH2')
df_MAR['StationID'] = df_MAR['StationID'].replace('0911', 'CMO2/CH9')
df_MAR['StationID'] = df_MAR['StationID'].replace('0921', 'IF1')
df_MAR['StationID'] = df_MAR['StationID'].replace('0931', 'TASO5')
df_MAR['StationID'] = df_MAR['StationID'].replace('0951', 'TASO3')
df_MAR['StationID'] = df_MAR['StationID'].replace('0971', 'TASO1')
df_MAR['StationID'] = df_MAR['StationID'].replace('0981', 'TSI1')
df_MAR['StationID'] = df_MAR['StationID'].replace('1001', 'TSI3')
df_MAR['StationID'] = df_MAR['StationID'].replace('1021', 'TSI5')
df_MAR['StationID'] = df_MAR['StationID'].replace('1031', 'TSI6')
df_MAR['StationID'] = df_MAR['StationID'].replace('1061', 'TESL2')
df_MAR['StationID'] = df_MAR['StationID'].replace('1071', 'TESL3/IML4')
df_MAR['StationID'] = df_MAR['StationID'].replace('1111', 'TESL6')
df_MAR['StationID'] = df_MAR['StationID'].replace('1121', 'TESL7')
df_MAR['StationID'] = df_MAR['StationID'].replace('1151', 'CM03/CH12')
df_MAR['StationID'] = df_MAR['StationID'].replace('123A1', '14ML')
df_MAR['StationID'] = df_MAR['StationID'].replace('1401', 'IF37')
df_MAR['StationID'] = df_MAR['StationID'].replace('1421', 'IF32')
df_MAR['StationID'] = df_MAR['StationID'].replace('1441', 'IF34')

## 2. --- NL data (main NL file. Olivia may have changed the headers)
# This is basically the ocean_carbon_2014-2018
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
df_NL = df_NL.reset_index(drop=True)

# Remove space in column names (MAY NOT BE NECESSARY IF WE RENAME)
cols = df_NL.columns
cols = cols.map(lambda x: x.replace(' ', '_') if isinstance(x, (str)) else x)
df_NL.columns = cols

# Rename columns
df_NL = df_NL.rename(columns={'Ship_Trip' : 'TripID'})
df_NL = df_NL.rename(columns={'Latitude_(Decimal_Degrees)' : 'latitude'})
df_NL = df_NL.rename(columns={'Longiude_(Decimal_Degrees)' : 'longitude'})
df_NL = df_NL.rename(columns={'Station_Name' : 'StationID'})
df_NL = df_NL.rename(columns={'Pressure_(dbar)' : 'depth'})
df_NL = df_NL.rename(columns={'Salinity_(psu)' : 'salinity'})
df_NL = df_NL.rename(columns={'Temperature_(C)' : 'temperature'})
df_NL = df_NL.rename(columns={'Dissolved_Oxygen_(mL/L)' : 'O2'})
df_NL = df_NL.rename(columns={'Calibrated_Oxygen_(mL/L)' : 'oxygen'})
df_NL = df_NL.rename(columns={'Phosphate_Concentration_(mmol_m-3)' : 'PO4'})
df_NL = df_NL.rename(columns={'Silicate_Concentration_(mmol_m-3)' : 'SiO'})
df_NL = df_NL.rename(columns={'Nitrate_Concentration_(mmol_m-3)' : 'NO3'})
df_NL = df_NL.rename(columns={'pH_25_' : 'pH_25'})
df_NL = df_NL.rename(columns={'Discrete_Chlorophyll_a_(mg_m-3)' : 'chla'})
df_NL = df_NL.rename(columns={'Total_P_(mmol_m-3)' : 'PO4'})
df_NL = df_NL.rename(columns={'Total_Si_(mmol_m-3)' : 'SiO'})
df_NL = df_NL.rename(columns={'Total_NO2/3_(mmol_m-3)' : 'NO3'})
df_NL = df_NL.rename(columns={'Measured_TA_(micromol/kg)' : 'TA'})
df_NL = df_NL.rename(columns={'Measure_TCO2_(micromol/kg)' : 'TIC'})
df_NL = df_NL.rename(columns={'Density_(kg_m3)' : 'density'})
df_NL = df_NL.rename(columns={'TA_Data_Flag' : 'TAflag'})
df_NL = df_NL.rename(columns={'TCO2_Data_Flag' : 'TICflag'})
df_NL['Region'] = 'NL'
df_NL['temp_pH']= np.NaN
df_NL['pHflag']= np.NaN

# Only select subset variables 
df_NL = df_NL.loc[:,variables]


# 2.2 Updated NL data
df_NLu = pd.read_excel(os.path.join(dataset_main_path, 'Ocean_Carbon.xlsx'))
df_NLu['timestamp'] = pd.to_datetime(df_NLu.Date.astype(str)) + pd.to_timedelta(df_NLu.GMT.astype(str)) 
df_NLu.index = df_NLu.timestamp


# HERE!!!!!


## 3. --- IML data
df2018 = pd.read_excel(os.path.join(dataset_path2,'AZMP_OA_IML2018.xlsx'), encoding='utf-8')
df2017 = pd.read_excel(os.path.join(dataset_path2,'AZMP_OA_IML2017.xlsx'), encoding='utf-8')
df2016 = pd.read_excel(os.path.join(dataset_path2,'AZMP_OA_IML2016.xlsx'), encoding='utf-8')
df2015 = pd.read_excel(os.path.join(dataset_path2,'AZMP_OA_IML2015.xlsx'), encoding='utf-8')
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
df2018f = pd.read_excel(os.path.join(dataset_path2,'AZMP_OA_IML2018f.xlsx'), encoding='utf-8')
df2017f = pd.read_excel(os.path.join(dataset_path2,'AZMP_OA_IML2017f.xlsx'), encoding='utf-8')
df2016f = pd.read_excel(os.path.join(dataset_path2,'AZMP_OA_IML2016f.xlsx'), encoding='utf-8')
df2018s = pd.read_excel(os.path.join(dataset_path2,'AZMP_OA_IML2018s.xlsx'), encoding='utf-8') ###some items in CTD Fichier were changed for an easier 'TripID'
df2017s = pd.read_excel(os.path.join(dataset_path2,'AZMP_OA_IML2017s.xlsx'), encoding='utf-8') ###some items in CTD Fichier were changed for an easier 'TripID'
df_IML = pd.concat([df2018f, df2018s, df2018, df2017f, df2017s, df2017, df2016f, df2016, df2015, df2015f], axis=0, sort=False)
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
df_IML['pH_25'] = df_IML[['labo_LBPHT_02_(Total_scale)', 'labo_LBPHT_02_(Total_scale).1']].mean(axis=1)
df_IML['TIC']= np.NaN
df_IML['O2'] = df_IML[['labo_oxy_02_(ml/l)', 'labo_oxy_02_(ml/l).1', 'labo_OXY_02_(ml/l)', 'labo_OXY_02_(ml/l).1']].mean(axis=1)

# Only select subset variables 
df_IML = df_IML.loc[:,variables]

## 4. Merged dataset 
df = pd.concat([df_MAR, df_NL, df_IML], axis=0)
df = df.reset_index(drop=True)
satO2 = swx.satO2(df['salinity'], df['temperature']) 
df['satO2_perc'] = df['O2']/satO2*100 # calculate oxygen saturation %
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

# Save final dataset
df.to_csv(os.path.join(dataset_main_path, 'AZMP_carbon_data.csv', float_format='%.4f')
