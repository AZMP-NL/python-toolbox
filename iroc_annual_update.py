'''

I provide 4 files/timeseries to the IROC:

Newf_Cartwright_Air_Timeseries.csv
Newf_CIL_Area_Timeseries.csv
Newf_SeaIce_Timeseries.csv
Newf_Station27_Annual.csv


'''

import os
import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 

clim_years = [1991, 2020]
year = 2024

## -- 1. Air T Cartwright -- ##
airT_anom = pd.read_pickle('operation_files/airT_anom.pkl')
airT_anom = airT_anom[['Cartwright']]
airT_anom = airT_anom[airT_anom.index>=1935]
airT_mnth = pd.read_pickle('operation_files/airT_monthly.pkl')
df_clim_period = airT_mnth[(airT_mnth.index.year>=clim_years[0]) & (airT_mnth.index.year<=clim_years[1])]
df_monthly_stack = df_clim_period.groupby([(df_clim_period.index.year),(df_clim_period.index.month)]).mean()
df_monthly_clim = df_monthly_stack.groupby(level=1).mean()
df_stack = airT_mnth.groupby([(airT_mnth.index.year),(airT_mnth.index.month)]).mean() 
df_stack_anom = df_stack.sub(df_monthly_clim, level=1)
df_annual_anom = df_stack_anom.groupby(level=0).mean()
df_annual_star =  df_annual_anom + df_monthly_clim.mean() 
#airT_mnth = airT_mnth.groupby(airT_mnth.index.year).mean()
airT_mnth = df_annual_star[['Cartwright']]
airT_mnth = airT_mnth[airT_mnth.index>=1935]
airT = pd.concat([airT_mnth,airT_anom], axis=1)

#Set up the meta-data
with open('iroc_airT_Cartwright.csv','w',newline='') as fd:
    wr = csv.writer(fd, quoting=csv.QUOTE_ALL)
    wr.writerow(['Station Description:','Cartwright Air Temperature'])
    wr.writerow(['Region:','2'])
    wr.writerow(['Latitude:','53.7 N'])
    wr.writerow(['Longitude:','57.0 W'])
    wr.writerow(['Measurement Depth/Range:','Land Surface'])
    wr.writerow(['Sampling Frequency:',''])
    wr.writerow(['Notes:','5 Year average calculated as a running mean'])
    wr.writerow(['Data Source:','Northwest Atlantic Fisheries Centre - Canada'])
    wr.writerow(['Contact Name (email):','Jonathan Coyne (Jonathan.Coyne@dfo-mpo.gc.ca)'])
    wr.writerow(['Location of source datasets:','Northwest Atlantic Fisheries Centre - Canada'])
    wr.writerow(['Data Protection:','Include on ICES website for free distribution.'])
    wr.writerow(['Decimal Year','Temperature °C','Temperature Anomaly Â°C'])

#Append the data to the bottom
airT.to_csv('iroc_airT_Cartwright.csv', mode='a', sep=',', float_format='%.2f', header=False)



## -- 2. CIL data -- ##
# ** need to run iroc_CIL_area.py for 'BB' and 'SI'
# **** HERE, should load .csv file in ~/AZMP/state_reports/sections_plots/CIL since 2021!! (note from 21 Nov. 2021)
df_BB = pd.read_csv('operation_files/CIL_area_BB.csv')
df_SI = pd.read_csv('operation_files/CIL_area_SI.csv')
df_BB.set_index('year', inplace=True)
df_SI.set_index('year', inplace=True)
# select interp field
df_BB = df_BB[['interp_field']]
df_SI = df_SI[['interp_field']]
df_BB.rename(columns={"interp_field": "BB"}, inplace=True)
df_SI.rename(columns={"interp_field": "SI"}, inplace=True)
# Concatenate
df = pd.concat([df_BB, df_SI], axis=1)
df.index.name='Year'
df = df.replace({0:np.nan})
# Some manual tweaking to remove very low values
flag_SI_itp = np.array([1932,1936,1937,1939,1940,1941,1950,1966,1967,1968,1972,1981,1989,2022])
flag_BB_itp = np.array([1966,1967,1968,1993,2022])
df.BB.loc[df.index<1950]=np.nan
for i in flag_SI_itp:
	df.SI.loc[df.index==i] = np.nan
for i in flag_BB_itp:
	df.BB.loc[df.index==i] = np.nan
df = df.dropna(how='all')

#Set up the meta-data
with open('iroc_CIL.csv','w',newline='') as fd:
    wr = csv.writer(fd, quoting=csv.QUOTE_ALL)
    wr.writerow(['Station Description:','CIL Layer Newfoundland and Labrador'])
    wr.writerow(['Region:','2'])
    wr.writerow(['Latitude:','49.5 N'])
    wr.writerow(['Longitude:','50.5 W'])
    wr.writerow(['Area','POLYGON((-47.947 48.733, -52.967 48.733, -52.967 50.332, -47.947 50.332, -47.947 48.733))'])
    wr.writerow(['Measurement Depth/Range:',''])
    wr.writerow(['Sampling Frequency:',''])
    wr.writerow(['Notes:','Area of Cold Intermediate Layer'])
    wr.writerow(['Data Source:','Northwest Atlantic Fisheries Centre - Canada'])
    wr.writerow(['Contact Name (email):','Jonathan Coyne (Jonathan.Coyne@dfo-mpo.gc.ca)'])
    wr.writerow(['Location of source datasets:','Northwest Atlantic Fisheries Centre - Canada'])
    wr.writerow(['Data Protection:','Include on ICES website for free distribution.'])
    wr.writerow(['Decimal Year','Newfoundland CIL Area (km^2)','Labrador CIL Area (km^2)'])

#Append the data to the bottom
df.to_csv('iroc_CIL.csv', mode='a', sep=',', float_format='%.1f',header=False)




'''
## For comparision
df_orig = pd.read_csv('../2017_orig/Newf_CIL_Area_Timeseries.csv', header=12, names=['Year', 'BB', 'SI'])
df_orig.set_index('Year', inplace=True)
df_BB.plot()
df_orig.BB.plot()
plt.title('CIL calculation on BB')
plt.legend(['new calculation', 'previous version'])
plt.grid()

df_SI.plot()
df_orig.SI.plot()
plt.title('CIL calculation on SI')
plt.legend(['new calculation', 'previous version'])
plt.grid()

# To help writing the IROC text:
df_clim = df[(df.index>=clim_years[0]) & (df.index<=clim_years[1])]
df_std_anom = (df - df_clim.mean()) / df_clim.std()
df_std_anom.plot()

## -- 3. Stn 27 data -- ##
#Manually copy-paste iroc_stn27.csv in
# created by azmp_stn27_analysis.py ** check to make sure you have the good climatology!
df_s27 = pd.read_csv('./Newf_Station27_Annual.csv', header=16, names=['Year', 'T', 'Tanom', 'Tstd', 'S', 'Sanom', 'Sstd'])
df_s27.set_index('Year', inplace=True)
df_s27.Tstd.plot()

## -- 4. sea ice data -- ##
# Provided by P.S. Galbraith for 3 regions
df_ice = pd.read_csv('./Newf_SeaIce_Timeseries.csv', header=12, names=['Year', 'NLab', 'SLab', 'NFLD'])
df_ice.set_index('Year', inplace=True)
df_ice_clim = df_ice[(df_ice.index>=clim_years[0]) & (df_ice.index<=clim_years[1])]

df_ice_std_anom = (df_ice - df_ice_clim.mean()) / df_ice_clim.std()
'''