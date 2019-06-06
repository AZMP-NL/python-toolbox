# -*- coding: utf-8 -*-
'''
To generate Colbourne's and STACFIS composite anomalies


Uses this pickled DataFrame:

/home/cyrf0006/AZMP/state_reports/SSTs/SSTs_merged_monthly.pkl
generated by from azmp_sst_scorecards.py


'''

import numpy as  np
import matplotlib.pyplot as plt
import pandas as pd
import os
import unicodedata
from matplotlib.colors import from_levels_and_colors

clim_year = [1981, 2010]
years = [1980, 2018]
width = 0.7


#### ---- LOAD THE DATA ---- ####
# 1. CIL [years: vol_itp, core_itp, core_depth_itp]
df_CIL_SI = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/sections_plots/CIL/df_CIL_SI_summer.pkl')
df_CIL_BB = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/sections_plots/CIL/df_CIL_BB_summer.pkl')
df_CIL_FC = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/sections_plots/CIL/df_CIL_FC_summer.pkl')

# 2. NAO & AO [years: Value]
nao_winter = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/NAO/NAO_winter.pkl')
ao = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/NAO/AO_annual.pkl')

# 3. Air Temperature
df_air = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/airTemp/airT_monthly.pkl')
df_air = df_air.resample('As').mean()
df_air.index = df_air.index.year

# 4. SSTs (problem: NS data missing prior 1997...)
df_sst = pd.read_pickle('/home/cyrf0006/AZMP/annual_meetings/2019/SSTs_merged_monthly.pkl')
df_sst = df_sst.resample('As').mean()
df_sst.index = df_sst.index.year

df_sst_1997 = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/SSTs/SSTs_bometrics_annual.pkl')
df_sst_1997 = df_sst_1997.resample('As').mean()
df_sst_1997.index = df_sst_1997.index.year


# 5. Bottom temperature
# 3LNO - Spring
df_3LNO_spring = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/bottomT/stats_3LNO_spring.pkl')
df_3LNO_spring.index = pd.to_datetime(df_3LNO_spring.index) # update index to datetime
df_3LNO_spring.index = df_3LNO_spring.index.year
df_3LNO_spring = df_3LNO_spring.Tmean
# 3Ps - Spring
df_3Ps_spring = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/bottomT/stats_3Ps_spring.pkl')
df_3Ps_spring.index = pd.to_datetime(df_3Ps_spring.index) # update index to datetime
df_3Ps_spring.index = df_3Ps_spring.index.year
df_3Ps_spring = df_3Ps_spring.Tmean
# 2J - Fall
df_2J_fall = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/bottomT/stats_2J_fall.pkl')
df_2J_fall.index = pd.to_datetime(df_2J_fall.index) # update index to datetime
df_2J_fall.index = df_2J_fall.index.year
df_2J_fall = df_2J_fall.Tmean
# 3K - Fall
df_3K_fall = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/bottomT/stats_3K_fall.pkl')
df_3K_fall.index = pd.to_datetime(df_3K_fall.index) # update index to datetime
df_3K_fall.index = df_3K_fall.index.year
df_3K_fall = df_3K_fall.Tmean
# 3LNO - Fall
df_3LNO_fall = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/bottomT/stats_3LNO_fall.pkl')
df_3LNO_fall.index = pd.to_datetime(df_3LNO_fall.index) # update index to datetime
df_3LNO_fall.index = df_3LNO_fall.index.year
df_3LNO_fall = df_3LNO_fall.Tmean
# 3LNO - Summer
df_3M_summer = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/bottomT/stats_3M_summer.pkl')
df_3M_summer.index = pd.to_datetime(df_3M_summer.index) # update index to datetime
df_3M_summer.index = df_3M_summer.index.year
df_3M_summer = df_3M_summer.Tmean

# 6. S27
df_s27 = pd.read_pickle('/home/cyrf0006/AZMP/S27/S27_temp_annual.pkl')
df_s27_mean = df_s27['mean']

# 7. Section average Temeprature (should eventually add salinity in these dataFrame, see azmp_CIL_stats.py)
df_SI = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/sections_plots/df_SI_meanT_summer.pkl')
df_BB = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/sections_plots/df_BB_meanT_summer.pkl')
df_FC = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/sections_plots/df_FC_meanT_summer.pkl')
df_FC_shelf = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/sections_plots/df_FC_meanT_shelf_summer.pkl')
df_FC_cap = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/sections_plots/df_FC_meanT_cap_summer.pkl')

# 8. Greenland Fylla and Cape Desolation
df_FB4 = pd.read_csv('/home/cyrf0006/data/IROC_timeseries/Greenland_Fylla_0-50_Annual.csv', header=15, index_col='Year')
df_CD3 = pd.read_csv('/home/cyrf0006/data/IROC_timeseries/Greenland_Desolation_2000_Annual.csv', header=14, index_col='Year')
# Keep only temeprature
df_FB4 = df_FB4.iloc[:,0]
df_CD3 = df_CD3.iloc[:,0]



#### ---- STACFIS - 3LNO ---- ####
df_comp_3LNO = pd.concat([df_s27_mean,
                          df_3LNO_spring, df_3LNO_fall,
                          df_SI, df_BB, df_FC_shelf,
                          df_CIL_SI.vol_itp, df_CIL_BB.vol_itp, df_CIL_FC.vol_itp,
                          df_sst.Avalon_Channel, df_sst.Hybernia, df_sst.Flemish_Pass 
                          ], axis=1)

df_3LNO_clim = df_comp_3LNO[(df_comp_3LNO.index>=clim_year[0]) & (df_comp_3LNO.index<=clim_year[1])]
std_anom_3LNO = (df_comp_3LNO-df_3LNO_clim.mean(axis=0))/df_3LNO_clim.std(axis=0)
# revert CIL volume
std_anom_3LNO['vol_itp'] = std_anom_3LNO['vol_itp']*-1
# mean anomaly
composite_3LNO = std_anom_3LNO.mean(axis=1)
composite_3LNO.to_csv('composite_3LNO.csv', float_format='%.2f')

# Plot
fig, ax = plt.subplots(nrows=1, ncols=1)
sign=composite_3LNO>0
composite_3LNO.plot(kind='bar', color=sign.map({True: 'indianred', False: 'steelblue'}), width = width)
plt.ylabel('Standardized Anomaly', weight='bold', fontsize=14)
plt.title('Composite anomaly 3LNO', weight='bold', fontsize=14)
plt.grid()
plt.ylim([-2.5,2.5])
fig.set_size_inches(w=15,h=7)
fig_name = 'composite_3LNO.png'
fig.savefig(fig_name, dpi=200)
os.system('convert -trim composite_3LNO.png composite_3LNO.png')


#### ---- STACFIS - 3M ---- ####
df_comp_3M = pd.concat([df_FC_cap, 
                        df_sst.Flemish_Cap,
                        df_3M_summer
                        ], axis=1)

df_3M_clim = df_comp_3M[(df_comp_3M.index>=clim_year[0]) & (df_comp_3M.index<=clim_year[1])]
std_anom_3M = (df_comp_3M-df_3M_clim.mean(axis=0))/df_3M_clim.std(axis=0)
composite_3M = std_anom_3M.mean(axis=1)
composite_3M.to_csv('composite_3M.csv', float_format='%.2f')

# Plot
fig, ax = plt.subplots(nrows=1, ncols=1)
sign=composite_3M>0
composite_3M.plot(kind='bar', color=sign.map({True: 'indianred', False: 'steelblue'}), width = width)
plt.ylabel('Standardized Anomaly', weight='bold', fontsize=14)
plt.title('Composite anomaly 3M', weight='bold', fontsize=14)
plt.grid()
plt.ylim([-2.5,2.5])
fig.set_size_inches(w=15,h=7)
fig_name = 'composite_3M.png'
fig.savefig(fig_name, dpi=200)
os.system('convert -trim composite_3M.png composite_3M.png')

#### ---- STACFIS - SA1 ---- ####
df_comp_SA1 = pd.concat([df_sst.Central_Labrador_Sea, df_sst.Greenland_Shelf,
                         df_FB4, df_CD3
                          ], axis=1)

df_SA1_clim = df_comp_SA1[(df_comp_SA1.index>=clim_year[0]) & (df_comp_SA1.index<=clim_year[1])]
std_anom_SA1 = (df_comp_SA1-df_SA1_clim.mean(axis=0))/df_SA1_clim.std(axis=0)
composite_SA1 = std_anom_SA1.mean(axis=1)
composite_SA1.to_csv('composite_SA1.csv', float_format='%.2f')

# Plot
fig, ax = plt.subplots(nrows=1, ncols=1)
sign=composite_SA1>0
composite_SA1.plot(kind='bar', color=sign.map({True: 'indianred', False: 'steelblue'}), width = width)
plt.ylabel('Standardized Anomaly', weight='bold', fontsize=14)
plt.title('Composite anomaly SA1', weight='bold', fontsize=14)
plt.grid()
plt.ylim([-2.5,2.5])
fig.set_size_inches(w=15,h=7)
fig_name = 'composite_SA1.png'
fig.savefig(fig_name, dpi=200)
os.system('convert -trim composite_SA1.png composite_SA1.png')



#### ---- STACFIS - SA234 ---- ####
# ** HEre I ignored SSTs, because not available for NS region prior to 1997
#                          df_sst.Hudson_Strait, df_sst.Hamilton_Bank, df_sst['St.Anthony_Basin'], df_sst.Orphan_Knoll 
#                          df_sst.Avalon_Channel, df_sst.Hybernia, df_sst.Flemish_Pass, df_sst.Flemish_Cap, ...

df_comp_3LNO = pd.concat([df_s27_mean,
                          df_3LNO_spring, df_3LNO_fall, df_3M_summer,
                          df_SI, df_BB, df_FC_shelf,
                          df_CIL_SI, df_CIL_BB, df_CIL_FC
                          ], axis=1)

df_3LNO_clim = df_comp_3LNO[(df_comp_3LNO.index>=clim_year[0]) & (df_comp_3LNO.index<=clim_year[1])]
std_anom_3LNO = (df_comp_3LNO-df_3LNO_clim.mean(axis=0))/df_3LNO_clim.std(axis=0)
composite_3LNO = std_anom_3LNO.mean(axis=1)
composite_3LNO.to_csv('composite_3LNO.csv', float_format='%.2f')



