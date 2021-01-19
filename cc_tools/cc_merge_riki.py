# -*- coding: utf-8 -*-
"""
This script is a complement to cc_merge_azmp.py to merge only Station Rimouski data.


NOTE CONCERNING ocean acidification DATA TYPES			
LBPHT_01: pH mesuré au laboratoire : spectrophotometric determination on freshly collected sample (SOP 6b in Dickson et al. 2007)	

LBPHT_02: pH mesuré au laboratoire : spectrophotometric determination (SOP 6b in Dickson et al. 2007) BUT on preserved sample			

pHT_01: pH mesuré au laboratoire : spectrophotometric determination on freshly collected sample (SOP 6b in Dickson et al. 2007) et converti au pH in situ selon Lewis & Wallace 1998			
pHT_02: pH mesuré au laboratoire : spectrophotometric determination (SOP 6b in Dickson et al. 2007) BUT on preserved sample et converti au pH in situ selon Lewis & Wallace 1998			
TICW_01: Determination of total dissolved inorganic carbon in seawater (SOP 2 in Dickson et al. 2007)			

ALKW_01: Determination of total alkalinity in seawater using an open-cell titration (SOP 3b in Dickson et al. 2007)

Columns to keep after CO2sys:
    'timestamp', 'Region', 'TripID', 'StationID',
    'latitude', 'longitude', 'depth',
    'temperature', 'salinity', 'O2',
    'NO3', 'PO4', 'SiO',
    'TAc', 'TICc', 'pH_tot', 'Omega_A', 'Omega_C', 'pCO2', 'O2sat_perc',
    'pH_25', 'temp_pH', 'TA', 'TIC']]
    
@author:
Olivia.Gibb@dfo-mpo.gc.ca
Frederic.Cyr@dfo-mpo.gc.ca

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

# path:
dataset_path = '/home/cyrf0006/github/AZMP-NL/datasets/carbonates/data_received/'

## 3.3 --- Station Riki
# load all
xls_riki = pd.ExcelFile(os.path.join(dataset_path,'StationRimouski (2013-2018)(GIBBB).xlsx'))
riki2013 = pd.read_excel(xls_riki, 'Rimouski 2013', header=1)
riki2014 = pd.read_excel(xls_riki, 'Rimouski 2014', header=1)
riki2015 = pd.read_excel(xls_riki, 'Rimouski 2015', header=1)
riki2016 = pd.read_excel(xls_riki, 'Rimouski 2016', header=1)
riki2017 = pd.read_excel(xls_riki, 'Rimouski 2017', header=1)
riki2018 = pd.read_excel(xls_riki, 'Rimouski 2018', header=1)
# drop header lines
riki2013 = riki2013.drop(0)
riki2014 = riki2014.drop([0,1])
riki2015 = riki2015.drop(0)
riki2016 = riki2016.drop(0)
riki2017 = riki2017.drop(0)
riki2018 = riki2018.drop(0)

## Cleaning the mess...
# lower all column names
riki2013.columns = map(str.lower, riki2013.columns)
riki2014.columns = map(str.lower, riki2014.columns)
riki2015.columns = map(str.lower, riki2015.columns)
riki2016.columns = map(str.lower, riki2016.columns)
riki2017.columns = map(str.lower, riki2017.columns)
riki2018.columns = map(str.lower, riki2018.columns)

## Columns to keep: 'timestamp', 'Region', 'TripID', 'StationID', 'latitude', 'longitude', 'depth',
## 'temperature', 'salinity', 'O2', 'NO3', 'PO4', 'SiO', 'TAc', 'TICc', 'pH_tot', 'Omega_A', 'Omega_C',
## 'pCO2', 'O2sat_perc', 'pH_25', 'temp_pH', 'TA', 'TIC'
 
# Timestamp (parse Date / Heure)  - Date is a mess, but can be solved as is:
# 2014 has a different format than other years!!!
riki2013['timestamp'] = pd.to_datetime(riki2013.date).astype(str) + ' ' + riki2013.heure.astype(str)
riki2013['timestamp'] =  pd.to_datetime(riki2013['timestamp'], format='%Y-%m-%d %H:%M:%S')
riki2014['timestamp'] = pd.to_datetime(riki2014.date).astype(str) + ' ' + riki2014.heure.astype(str)
riki2014['timestamp'] =  pd.to_datetime(riki2014['timestamp'], format='%Y-%m-%d %H:%M:%S')
riki2015['timestamp'] = pd.to_datetime(riki2015.date).astype(str) + ' ' + riki2015.heure.astype(str)
riki2015['timestamp'] =  pd.to_datetime(riki2015['timestamp'], format='%Y-%m-%d %H:%M:%S')
riki2016['timestamp'] = pd.to_datetime(riki2016.date).astype(str) + ' ' + riki2016.heure.astype(str)
riki2016['timestamp'] =  pd.to_datetime(riki2016['timestamp'], format='%Y-%m-%d %H:%M:%S')
riki2017['timestamp'] = pd.to_datetime(riki2017.date).astype(str) + ' ' + riki2017.heure.astype(str)
riki2017['timestamp'] =  pd.to_datetime(riki2017['timestamp'], format='%Y-%m-%d %H:%M:%S')
riki2018['timestamp'] = pd.to_datetime(riki2018.date).astype(str) + ' ' + riki2018.heure.astype(str)
riki2018['timestamp'] =  pd.to_datetime(riki2018['timestamp'], format='%Y-%m-%d %H:%M:%S')

# Region (all GSL, to be added later)

# TripID
riki2013 = riki2013.rename(columns={'mission' : 'TripID'})
riki2014 = riki2014.rename(columns={'mission' : 'TripID'})
riki2015 = riki2015.rename(columns={'mission' : 'TripID'})
riki2016 = riki2016.rename(columns={'unnamed: 0' : 'TripID'})
riki2017 = riki2017.rename(columns={'unnamed: 0' : 'TripID'})
riki2018 = riki2018.rename(columns={'unnamed: 0' : 'TripID'})

# 'latitude' / 'longitude' ->  uniform (to be rename later)
# depth  -> 'pres' uniform (to be rename later)
# temperature -> 'te90' uniform (to be rename later)
# salinity -> 'psal' uniform (to be rename later)
## O2 -> oxy_02 (_zero-2) area available 2015-2018 (to be rename later)

##  NO3 - Deal with duplicates and flags
# 2013
riki2013['q_no3+2'].replace(np.nan, 1, inplace=True)
riki2013.loc[riki2013['q_no3+2']!=1, 'no3+2_mean'] = np.nan
riki2013['NO3'] = riki2013['no3+2_mean']
# 2014
riki2014['q_no3+2'].replace(np.nan, 1, inplace=True)
riki2014.loc[riki2014['q_no3+2']!=1, 'no3+2_mean'] = np.nan
riki2014['NO3'] = riki2014['no3+2_mean']
# 2015
riki2015['q_no3+2'].replace(np.nan, 1, inplace=True)
riki2015.loc[riki2015['q_no3+2']!=1, 'no3+2_mean'] = np.nan
riki2015['NO3'] = riki2015['no3+2_mean']
# 2016
riki2016.loc[riki2016['q_nox']!=1, 'nox_03'] = np.nan
riki2016.loc[riki2016['q_nox.1']!=1, 'nox_03.1'] = np.nan
riki2016.loc[riki2016['q_nox.2']!=1, 'nox_03.2'] = np.nan
riki2016['NO3'] = riki2016[['nox_03', 'nox_03.1', 'nox_03.2']].mean(axis=1)
# 2017
riki2017.loc[riki2017['q_nox']!=1, 'nox_03'] = np.nan
riki2017.loc[riki2017['q_nox.1']!=1, 'nox_03.1'] = np.nan
riki2017.loc[riki2017['q_nox.2']!=1, 'nox_03.2'] = np.nan
riki2017['NO3'] = riki2017[['nox_03', 'nox_03.1', 'nox_03.2']].mean(axis=1)
# 2018
riki2018.loc[riki2018['q_nox']!=1, 'nox_03'] = np.nan
riki2018.loc[riki2018['q_nox.1']!=1, 'nox_03.1'] = np.nan
riki2018.loc[riki2018['q_nox.2']!=1, 'nox_03.2'] = np.nan
riki2018['NO3'] = riki2018[['nox_03', 'nox_03.1', 'nox_03.2']].mean(axis=1)

##  PO4 - Deal with duplicates and flags
# 2013
riki2013['q_po4'].replace(np.nan, 1, inplace=True)
riki2013.loc[riki2013['q_po4']!=1, 'po4_mean'] = np.nan
riki2013['PO4'] = riki2013['po4_mean']
# 2014
riki2014['q_po4'].replace(np.nan, 1, inplace=True)
riki2014.loc[riki2014['q_po4']!=1, 'po4_mean'] = np.nan
riki2014['PO4'] = riki2014['po4_mean']
# 2015
riki2015['q_po4'].replace(np.nan, 1, inplace=True)
riki2015.loc[riki2015['q_po4']!=1, 'po4_mean'] = np.nan
riki2015['PO4'] = riki2015['po4_mean']
# 2016
riki2016.loc[riki2016['q_po4']!=1, 'po4_03'] = np.nan
riki2016.loc[riki2016['q_po4.1']!=1, 'po4_03.1'] = np.nan
riki2016.loc[riki2016['q_po4.2']!=1, 'po4_03.2'] = np.nan
riki2016['PO4'] = riki2016[['po4_03', 'po4_03.1', 'po4_03.2']].mean(axis=1)
# 2017
riki2017.loc[riki2017['q_po4']!=1, 'po4_03'] = np.nan
riki2017.loc[riki2017['q_po4.1']!=1, 'po4_03.1'] = np.nan
riki2017.loc[riki2017['q_po4.2']!=1, 'po4_03.2'] = np.nan
riki2017['PO4'] = riki2017[['po4_03', 'po4_03.1', 'po4_03.2']].mean(axis=1)
# 2018
riki2018.loc[riki2018['q_po4']!=1, 'po4_03'] = np.nan
riki2018.loc[riki2018['q_po4.1']!=1, 'po4_03.1'] = np.nan
riki2018.loc[riki2018['q_po4.2']!=1, 'po4_03.2'] = np.nan
riki2018['PO4'] = riki2018[['po4_03', 'po4_03.1', 'po4_03.2']].mean(axis=1)

##  SiO - Deal with duplicates and flags
# 2013
riki2013['q_si'].replace(np.nan, 1, inplace=True)
riki2013.loc[riki2013['q_si']!=1, 'si_mean'] = np.nan
riki2013['SiO'] = riki2013['si_mean']
# 2014
riki2014['q_si'].replace(np.nan, 1, inplace=True)
riki2014.loc[riki2014['q_si']!=1, 'si_mean'] = np.nan
riki2014['SiO'] = riki2014['si_mean']
# 2015
riki2015['q_si'].replace(np.nan, 1, inplace=True)
riki2015.loc[riki2015['q_si']!=1, 'si_mean'] = np.nan
riki2015['SiO'] = riki2015['si_mean']
# 2016
riki2016.loc[riki2016['q_si']!=1, 'si_03'] = np.nan
riki2016.loc[riki2016['q_si.1']!=1, 'si_03.1'] = np.nan
riki2016.loc[riki2016['q_si.2']!=1, 'si_03.2'] = np.nan
riki2016['SiO'] = riki2016[['si_03', 'si_03.1', 'si_03.2']].mean(axis=1)
# 2017
riki2017.loc[riki2017['q_si']!=1, 'si_03'] = np.nan
riki2017.loc[riki2017['q_si.1']!=1, 'si_03.1'] = np.nan
riki2017.loc[riki2017['q_si.2']!=1, 'si_03.2'] = np.nan
riki2017['SiO'] = riki2017[['si_03', 'si_03.1', 'si_03.2']].mean(axis=1)
# 2018
riki2018.loc[riki2018['q_si']!=1, 'si_03'] = np.nan
riki2018.loc[riki2018['q_si.1']!=1, 'si_03.1'] = np.nan
riki2018.loc[riki2018['q_si.2']!=1, 'si_03.2'] = np.nan
riki2018['SiO'] = riki2018[['si_03', 'si_03.1', 'si_03.2']].mean(axis=1)

# TA (need to deal with flag later)
riki2013 = riki2013.rename(columns={'at' : 'TA'})
riki2014 = riki2014.rename(columns={'at' : 'TA'})
riki2015 = riki2015.rename(columns={'mean at' : 'TA'})
riki2016 = riki2016.rename(columns={'at' : 'TA'})
riki2016 = riki2016.rename(columns={'at.1' : 'TAflag'})
riki2017 = riki2017.rename(columns={'at' : 'TA'})
riki2017 = riki2017.rename(columns={'at.1' : 'TAflag'})
riki2018 = riki2018.rename(columns={'at' : 'TA'})
riki2018 = riki2018.rename(columns={'at.1' : 'TAflag'})

# pH Lab (25 degC)
riki2013['pH_25'] = riki2013[['lbpht_01', 'lbpht_01.1', 'lbpht_01.2']].mean(axis=1)
riki2013['pHflag'] = riki2013[['q_lbpht', 'q_lbpht.1', 'q_lbpht.2']].mean(axis=1)
riki2014['pH_25'] = riki2014[['lbpht_01', 'lbpht_01.1']].mean(axis=1)
riki2014['pHflag'] = riki2014[['q_lbpht', 'q_lbpht.1']].mean(axis=1)
# 2015 and 2016 are special because intrument shifted throughout the season
riki2015['pH_25'] = riki2015[['lbpht_01', 'lbpht_01.1','lbpht_02', 'lbpht_02.1']].mean(axis=1)
riki2015['pHflag'] = riki2015[['q_lbpht', 'q_lbpht.1','q_lbpht.2', 'q_lbpht.3']].replace(9, np.nan).mean(axis=1)
riki2016['pH_25'] = riki2016[['lbpht_01', 'lbpht_01.1','lbpht_02', 'lbpht_02.1']].mean(axis=1)
riki2016['pHflag'] = riki2016[['q_lbpht', 'q_lbpht.1','q_lbpht.2', 'q_lbpht.3']].replace(9, np.nan).mean(axis=1)
riki2017['pH_25'] = riki2017[['lbpht_02', 'lbpht_02.1']].mean(axis=1)
riki2017['pHflag'] = riki2017[['q_lbpht', 'q_lbpht.1']].mean(axis=1)
riki2018['pH_25'] = riki2018[['lbpht_02', 'lbpht_02.1']].mean(axis=1)
riki2018['pHflag'] = riki2018['pH_25']*np.nan # no Flag

# CTD pH
## riki2013 = riki2013.rename(columns={'ctd pht_' : 'CTD PHPH (Total scale)'})
## riki2014 = riki2014.rename(columns={'ctd pht_' : 'CTD PHPH (Total scale)'})
## riki2015 = riki2015.rename(columns={'ctd pht_' : 'CTD PHPH (Total scale)'})
## riki2016 = riki2016.rename(columns={'phph_t' : 'CTD PHPH (Total scale)'})
## riki2017 = riki2017.rename(columns={'phph_t' : 'CTD PHPH (Total scale)'})
## riki2018 = riki2018.rename(columns={'phph' : 'CTD PHPH (NBS scale)'})

# merge
riki = pd.concat([riki2013, riki2014, riki2015, riki2016, riki2017, riki2018], sort=False)

# TA and pH flags
riki['TAflag'].replace(np.nan, 1, inplace=True)
riki['pHflag'].replace(np.nan, 1, inplace=True)

# O2
riki = riki.rename(columns={'oxy_02' : 'O2'})
riki = riki.rename(columns={'q_oxy' : 'O2_flag'})
riki['O2_flag'].replace(np.nan, 1, inplace=True)
riki.loc[riki['O2_flag']!=1, 'O2'] = np.nan

# Rename merged product
riki['Region'] = 'GSL'
riki['StationID'] = 'TESL3'
riki = riki.rename(columns={'pres' : 'depth'})
riki = riki.rename(columns={'te90' : 'temperature'})
riki = riki.rename(columns={'psal' : 'salinity'})
riki['labt_01'] = 'temp_pH'

# Only keep desired columns 
variables = ['TripID', 'Region', 'StationID', 'timestamp', 'latitude', 'longitude', 'depth', 'NO3', 'O2','TA', 'temp_pH', 'pH_25', 'PO4', 'salinity', 'SiO', 'temperature', 'TAflag', 'pHflag']
riki = riki.loc[:,variables]

# Make sure all dtypes are okay
riki = riki.astype({'TripID': 'str',
             'Region' : 'str',
             'StationID' : 'str',
             'latitude' : 'float',
             'longitude' : 'float',
             'depth' : 'float',
             'NO3' : 'float',
             'PO4' : 'float',
             'SiO' : 'float',
             'temperature' : 'float',
             'salinity' : 'float',
             'O2' : 'float',
             'TA' : 'float',
             'temp_pH' : 'float',
             'pH_25' : 'float',
             'TAflag' : 'int',
             'pHflag' : 'int'
            })


riki.to_pickle('riki.pkl')  

