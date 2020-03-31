# -*- coding: utf-8 -*-
'''
To generate AZMP Station 27 analysis scorecards

Uses these pickled DataFrames (generated in azmp_s27.py and azmp_s27_analysis.py):

S27_temperature_monthly.pkl
S27_salinity_monthly.pkl
S27_CIL_summer_stats.pkl
S27_MLD_monthly.pkl
S27_stratif_monthly.pkl


'''

import numpy as  np
import matplotlib.pyplot as plt
import pandas as pd
import os
import unicodedata
from matplotlib.colors import from_levels_and_colors

year_clim = [1981, 2010]
years = [1980, 2019]

def is_number(s):
    #https://www.pythoncentral.io/how-to-check-if-a-string-is-a-number-in-python-including-unicode/
    try:
        float(s)
        return True
    except ValueError:
        pass 
    try:
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass 
    return False


#### ------------- T-S ---------------- ####
# Load pickled data
df_temp = pd.read_pickle('S27_temperature_monthly.pkl')
df_sal = pd.read_pickle('S27_salinity_monthly.pkl')
#df_temp = df_temp[df_temp.index.year>=years[0]]
#df_sal = df_sal[df_sal.index.year>=years[0]]

# Temperature
#df_temp = df_temp.resample('As').mean() # *
T_0_btm = df_temp.mean(axis=1)
T_0_50 = df_temp[df_temp.columns[(df_temp.columns<=50)]].mean(axis=1)
T_0_20 = df_temp[df_temp.columns[(df_temp.columns<=20)]].mean(axis=1)
T_150_btm = df_temp[df_temp.columns[(df_temp.columns>=150)]].mean(axis=1)
# Salinity
#df_sal = df_sal.resample('As').mean() # *
S_0_btm = df_sal.mean(axis=1)
S_0_50 = df_sal[df_sal.columns[(df_sal.columns<=50)]].mean(axis=1)
S_150_btm = df_sal[df_sal.columns[(df_sal.columns>=150)]].mean(axis=1)
# Merge and compute anomalies 
df_T = pd.concat([T_0_btm, T_0_50, T_150_btm], axis=1, keys=['Temp 0-176m', 'Temp 0-50m', 'Temp 150-176m'])
df_S = pd.concat([S_0_btm, S_0_50, S_150_btm], axis=1, keys=['Sal 0-176m', 'Sal 0-50m', r'Sal 150-176m ${~}$'])
df_T_nico = pd.concat([T_0_btm, T_0_20, T_150_btm], axis=1, keys=['Temp 0-176m', 'Temp 0-20m', 'Temp 150-176m'])

# Save csv (for Guoqi)
df_T.to_csv('stn27_monthly_T.csv', float_format='%.3f')
df_T_nico.to_csv('stn27_monthly_T_nico.csv', float_format='%.3f')
df_T_nico.resample('As').mean().to_csv('stn27_annual_T_nico.csv', float_format='%.3f')
df_S.to_csv('stn27_monthly_S.csv', float_format='%.3f')
# **New from 2020 - annual anomaly is mean of monthly anomalies**
# Temperature anomalies
df_T_stack = df_T.groupby([(df_T.index.year),(df_T.index.month)]).mean()
df_T_stack_ave = df_T_stack['Temp 0-176m']
df_T_stack_surf = df_T_stack['Temp 0-50m']
df_T_stack_bot = df_T_stack['Temp 150-176m']
df_T_unstack_ave = df_T_stack_ave.unstack()
df_T_unstack_surf = df_T_stack_surf.unstack()
df_T_unstack_bot = df_T_stack_bot.unstack()
T_clim_period = df_T[(df_T.index.year>=year_clim[0]) & (df_T.index.year<=year_clim[1])]
T_monthly_stack = T_clim_period.groupby([(T_clim_period.index.year),(T_clim_period.index.month)]).mean()
T_monthly_clim = T_monthly_stack.mean(level=1)
T_monthly_std = T_monthly_stack.std(level=1)
Tave_monthly_anom = df_T_unstack_ave - T_monthly_clim['Temp 0-176m']
Tsurf_monthly_anom = df_T_unstack_surf - T_monthly_clim['Temp 0-50m']
Tbot_monthly_anom = df_T_unstack_bot - T_monthly_clim['Temp 150-176m']
Tave_monthly_stdanom = Tave_monthly_anom / T_monthly_std['Temp 0-176m']
Tsurf_monthly_stdanom = Tsurf_monthly_anom / T_monthly_std['Temp 0-50m']
Tbot_monthly_stdanom = Tbot_monthly_anom / T_monthly_std['Temp 150-176m']
Tave_anom = Tave_monthly_anom.mean(axis=1) 
Tsurf_anom = Tsurf_monthly_anom.mean(axis=1) 
Tbot_anom = Tbot_monthly_anom.mean(axis=1) 
Tave_anom.index = pd.to_datetime(Tave_anom.index, format='%Y')
Tsurf_anom.index = pd.to_datetime(Tsurf_anom.index, format='%Y')
Tbot_anom.index = pd.to_datetime(Tbot_anom.index, format='%Y')
Tave_anom_std = Tave_monthly_stdanom.mean(axis=1) 
Tsurf_anom_std = Tsurf_monthly_stdanom.mean(axis=1) 
Tbot_anom_std = Tbot_monthly_stdanom.mean(axis=1) 
Tave_anom_std.index = pd.to_datetime(Tave_anom_std.index, format='%Y')
Tsurf_anom_std.index = pd.to_datetime(Tsurf_anom_std.index, format='%Y')
Tbot_anom_std.index = pd.to_datetime(Tbot_anom_std.index, format='%Y')
# Salinity anomalies
df_S_stack = df_S.groupby([(df_S.index.year),(df_S.index.month)]).mean()
df_S_stack_ave = df_S_stack['Sal 0-176m']
df_S_stack_surf = df_S_stack['Sal 0-50m']
df_S_stack_bot = df_S_stack['Sal 150-176m ${~}$']
df_S_unstack_ave = df_S_stack_ave.unstack()
df_S_unstack_surf = df_S_stack_surf.unstack()
df_S_unstack_bot = df_S_stack_bot.unstack()
S_clim_period = df_S[(df_S.index.year>=year_clim[0]) & (df_S.index.year<=year_clim[1])]
S_monthly_stack = S_clim_period.groupby([(S_clim_period.index.year),(S_clim_period.index.month)]).mean()
S_monthly_clim = S_monthly_stack.mean(level=1)
S_monthly_std = S_monthly_stack.std(level=1)
Save_monthly_anom = df_S_unstack_ave - S_monthly_clim['Sal 0-176m']
Ssurf_monthly_anom = df_S_unstack_surf - S_monthly_clim['Sal 0-50m']
Sbot_monthly_anom = df_S_unstack_bot - S_monthly_clim['Sal 150-176m ${~}$']
Save_monthly_stdanom = Save_monthly_anom / S_monthly_std['Sal 0-176m']
Ssurf_monthly_stdanom = Ssurf_monthly_anom / S_monthly_std['Sal 0-50m']
Sbot_monthly_stdanom = Sbot_monthly_anom / S_monthly_std['Sal 150-176m ${~}$']
Save_anom = Save_monthly_anom.mean(axis=1) 
Ssurf_anom = Ssurf_monthly_anom.mean(axis=1) 
Sbot_anom = Sbot_monthly_anom.mean(axis=1) 
Save_anom.index = pd.to_datetime(Save_anom.index, format='%Y')
Ssurf_anom.index = pd.to_datetime(Ssurf_anom.index, format='%Y')
Sbot_anom.index = pd.to_datetime(Sbot_anom.index, format='%Y')
Save_anom_std = Save_monthly_stdanom.mean(axis=1) 
Ssurf_anom_std = Ssurf_monthly_stdanom.mean(axis=1) 
Sbot_anom_std = Sbot_monthly_stdanom.mean(axis=1) 
Save_anom_std.index = pd.to_datetime(Save_anom_std.index, format='%Y')
Ssurf_anom_std.index = pd.to_datetime(Ssurf_anom_std.index, format='%Y')
Sbot_anom_std.index = pd.to_datetime(Sbot_anom_std.index, format='%Y')
# merge anomalies
T_anom = pd.concat([Tave_anom,Tsurf_anom,Tbot_anom], axis=1, keys=['Temp 0-176m','Temp 0-50m','Temp 150-176m'])
T_anom_std = pd.concat([Tave_anom_std,Tsurf_anom_std,Tbot_anom_std], axis=1, keys=['Temp 0-176m','Temp 0-50m','Temp 150-176m'])
T_anom = T_anom[T_anom.index.year>=1947]
T_anom_std = T_anom_std[T_anom_std.index.year>=1947]
S_anom = pd.concat([Save_anom,Ssurf_anom,Tbot_anom], axis=1, keys=['Sal 0-176m','Sal 0-50m','Sal 150-176m ${~}$'])
S_anom_std = pd.concat([Save_anom_std,Ssurf_anom_std,Tbot_anom_std], axis=1, keys=['Sal 0-176m','Sal 0-50m','Sal 150-176m ${~}$'])
S_anom = S_anom[S_anom.index.year>=1947]
S_anom_std = S_anom_std[S_anom_std.index.year>=1947]
# Save pkl for climate indices (whole timeseries)
T_anom_std.to_pickle('s27_temp_std_anom.pkl')
S_anom_std.to_pickle('s27_sal_std_anom.pkl')
# keep only relevant window
T_anom_std = T_anom_std[T_anom_std.index.year>=years[0]]
S_anom_std = S_anom_std[S_anom_std.index.year>=years[0]]
# annual clims (for inclusion)
T_clim_period_annual = T_clim_period.resample('As').mean()    
S_clim_period_annual = S_clim_period.resample('As').mean()    

# old way
old_df_T = df_T.resample('As').mean() # ** there is a slight different if annual averages are made before
old_df_S = df_S.resample('As').mean() # ** there is a slight different if annual averages are made before
# clim
old_T_clim_period = old_df_T[(old_df_T.index.year>=year_clim[0]) & (old_df_T.index.year<=year_clim[1])]
old_S_clim_period = old_df_S[(old_df_S.index.year>=year_clim[0]) & (old_df_S.index.year<=year_clim[1])]
# Anomaly
old_T_anom = (old_df_T - old_T_clim_period.mean())
old_T_anom_std = old_T_anom / old_T_clim_period.std()
old_S_anom = (old_df_S - old_S_clim_period.mean())
old_S_anom_std = old_S_anom / old_S_clim_period.std()
# Save pkl for climate indices (whole timeseries)
#old_T_anom_std.to_pickle('s27_temp_stn_anom.pkl')
#old_S_anom_std.to_pickle('s27_sal_stn_anom.pkl')
# Reduce series (for scorecards)
old_T_anom_std = old_T_anom_std[old_T_anom_std.index.year>=years[0]]
old_S_anom_std = old_S_anom_std[old_S_anom_std.index.year>=years[0]]

#### ------------- CIL ---------------- ####
# Load pickled data
df_CIL = pd.read_pickle('S27_CIL_summer_stats.pkl')
df_CIL = df_CIL[df_CIL.index.year>=years[0]]
df_CIL.dropna(inplace=True)
# Compute anomalies 
CIL_clim_period = df_CIL[(df_CIL.index.year>=year_clim[0]) & (df_CIL.index.year<=year_clim[1])]
# Anomaly
CIL_anom = (df_CIL - CIL_clim_period.mean())
CIL_anom_std = CIL_anom / CIL_clim_period.std()

#### ------------- MLD ---------------- ####
# Load pickled data
df_MLD = pd.read_pickle('S27_MLD_monthly.pkl')
# flag some data:
df_MLD[df_MLD.index=='2019-03-15']=np.nan
df_MLD = df_MLD[df_MLD.index.year>=years[0]]
MLD_clim_period = df_MLD[(df_MLD.index.year>=year_clim[0]) & (df_MLD.index.year<=year_clim[1])]
# Monthly clim
MLD_monthly_stack = df_MLD.groupby([(df_MLD.index.year),(df_MLD.index.month)]).mean()
MLD_monthly_clim_stack = MLD_clim_period.groupby([(MLD_clim_period.index.year),(MLD_clim_period.index.month)]).mean()
MLD_monthly_clim = MLD_monthly_clim_stack.mean(level=1)
MLD_monthly_std = MLD_monthly_clim_stack.std(level=1)
# Monthly anomaly
MLD_anom = MLD_monthly_stack.sub(MLD_monthly_clim, level=1)
MLD_std_anom = MLD_anom.div(MLD_monthly_std, level=1)
MLD_std_anom.index = pd.to_datetime(MLD_std_anom.index.get_level_values(0).astype(str) + '-' + MLD_std_anom.index.get_level_values(1).astype(str) + '-1', format="%Y-%m-%d")
# dropna
MLD_std_anom.dropna(inplace=True)
# Annual and seasonal anomalies (average of monthly anomalies rather than anomalies of monthly values))
MLD_winter_anom = MLD_std_anom[(MLD_std_anom.index.month>=1) & (MLD_std_anom.index.month<=3)].resample('As').mean()
MLD_spring_anom = MLD_std_anom[(MLD_std_anom.index.month>=4) & (MLD_std_anom.index.month<=6)].resample('As').mean()
MLD_summer_anom = MLD_std_anom[(MLD_std_anom.index.month>=7) & (MLD_std_anom.index.month<=9)].resample('As').mean()
MLD_fall_anom = MLD_std_anom[(MLD_std_anom.index.month>=10) & (MLD_std_anom.index.month<=12)].resample('As').mean()
MLD_annual_anom = MLD_std_anom.resample('As').mean()
# Mean values for climatology period (for last columns)
MLD_winter_clim = MLD_clim_period[(MLD_clim_period.index.month>=1) & (MLD_clim_period.index.month<=3)].resample('As').mean()
MLD_spring_clim = MLD_clim_period[(MLD_clim_period.index.month>=4) & (MLD_clim_period.index.month<=6)].resample('As').mean()
MLD_summer_clim = MLD_clim_period[(MLD_clim_period.index.month>=7) & (MLD_clim_period.index.month<=9)].resample('As').mean()
MLD_fall_clim = MLD_clim_period[(MLD_clim_period.index.month>=10) & (MLD_clim_period.index.month<=12)].resample('As').mean()
MLD_annual_clim = MLD_clim_period.resample('As').mean()
# concat clim and anomaly
MLD_clim = pd.concat([MLD_winter_clim, MLD_spring_clim, MLD_summer_clim, MLD_fall_clim, MLD_annual_clim], axis=1, keys=['MLD winter', 'MLD spring', 'MLD summer', 'MLD fall', 'MLD annual'])
MLD_anom_std = pd.concat([MLD_winter_anom, MLD_spring_anom, MLD_summer_anom, MLD_fall_anom, MLD_annual_anom], axis=1, keys=['MLD winter', 'MLD spring', 'MLD summer', 'MLD fall', 'MLD annual'])
## # Annual and seasonal anomalies (old way)
## MLD_winter = df_MLD[(df_MLD.index.month>=1) & (df_MLD.index.month<=3)].resample('As').mean()
## MLD_spring = df_MLD[(df_MLD.index.month>=4) & (df_MLD.index.month<=6)].resample('As').mean()
## MLD_summer = df_MLD[(df_MLD.index.month>=7) & (df_MLD.index.month<=9)].resample('As').mean()
## MLD_fall = df_MLD[(df_MLD.index.month>=10) & (df_MLD.index.month<=12)].resample('As').mean()
## MLD_annual = df_MLD.resample('As').mean()
## # Clim period
## MLD_winter_clim = MLD_clim_period[(MLD_clim_period.index.month>=1) & (MLD_clim_period.index.month<=3)].resample('As').mean()
## MLD_spring_clim = MLD_clim_period[(MLD_clim_period.index.month>=4) & (MLD_clim_period.index.month<=6)].resample('As').mean()
## MLD_summer_clim = MLD_clim_period[(MLD_clim_period.index.month>=7) & (MLD_clim_period.index.month<=9)].resample('As').mean()
## MLD_fall_clim = MLD_clim_period[(MLD_clim_period.index.month>=10) & (MLD_clim_period.index.month<=12)].resample('As').mean()
## MLD_annual_clim = MLD_clim_period.resample('As').mean()
## # Anomalies
## MLD_std_anom_winter = (MLD_winter - MLD_winter_clim.mean()) / MLD_winter_clim.std()
## MLD_std_anom_spring = (MLD_spring - MLD_spring_clim.mean()) / MLD_spring_clim.std()
## MLD_std_anom_summer = (MLD_summer - MLD_summer_clim.mean()) / MLD_summer_clim.std()
## MLD_std_anom_fall = (MLD_fall - MLD_fall_clim.mean()) / MLD_fall_clim.std()
## MLD_std_anom_annual = (MLD_annual - MLD_annual_clim.mean()) / MLD_annual_clim.std()


#### ------------- Stratification ---------------- ####
# Load pickled data
df_strat = pd.read_pickle('S27_stratif_monthly.pkl')
# flag some data:
df_strat[df_strat.index=='2019-03-15']=np.nan
df_strat = df_strat[df_strat.index.year>=years[0]]
strat_clim_period = df_strat[(df_strat.index.year>=year_clim[0]) & (df_strat.index.year<=year_clim[1])]
# Monthly clim
strat_monthly_stack = df_strat.groupby([(df_strat.index.year),(df_strat.index.month)]).mean()
strat_monthly_clim_stack = strat_clim_period.groupby([(strat_clim_period.index.year),(strat_clim_period.index.month)]).mean()
strat_monthly_clim = strat_monthly_clim_stack.mean(level=1)
strat_monthly_std = strat_monthly_clim_stack.std(level=1)
# Monthly anomaly
strat_anom = strat_monthly_stack.sub(strat_monthly_clim, level=1)
strat_std_anom = strat_anom.div(strat_monthly_std, level=1)
strat_std_anom.index = pd.to_datetime(strat_std_anom.index.get_level_values(0).astype(str) + '-' + strat_std_anom.index.get_level_values(1).astype(str) + '-1', format="%Y-%m-%d")
strat_std_anom.dropna(inplace=True) # dropna
# Annual and seasonal anomalies (average of monthly anomalies rather than anomalies of monthly values))
strat_winter_anom = strat_std_anom[(strat_std_anom.index.month>=1) & (strat_std_anom.index.month<=3)].resample('As').mean()
strat_spring_anom = strat_std_anom[(strat_std_anom.index.month>=4) & (strat_std_anom.index.month<=6)].resample('As').mean()
strat_summer_anom = strat_std_anom[(strat_std_anom.index.month>=7) & (strat_std_anom.index.month<=9)].resample('As').mean()
strat_fall_anom = strat_std_anom[(strat_std_anom.index.month>=10) & (strat_std_anom.index.month<=12)].resample('As').mean()
strat_annual_anom = strat_std_anom.resample('As').mean()
# Mean values for climatology period (for last columns)
strat_winter_clim = strat_clim_period[(strat_clim_period.index.month>=1) & (strat_clim_period.index.month<=3)].resample('As').mean()
strat_spring_clim = strat_clim_period[(strat_clim_period.index.month>=4) & (strat_clim_period.index.month<=6)].resample('As').mean()
strat_summer_clim = strat_clim_period[(strat_clim_period.index.month>=7) & (strat_clim_period.index.month<=9)].resample('As').mean()
strat_fall_clim = strat_clim_period[(strat_clim_period.index.month>=10) & (strat_clim_period.index.month<=12)].resample('As').mean()
strat_annual_clim = strat_clim_period.resample('As').mean()
# concat clim and anomaly
strat_clim = pd.concat([strat_winter_clim, strat_spring_clim, strat_summer_clim, strat_fall_clim, strat_annual_clim], axis=1, keys=['strat winter', 'strat spring', 'strat summer', 'strat fall', 'strat annual'])
strat_anom_std = pd.concat([strat_winter_anom, strat_spring_anom, strat_summer_anom, strat_fall_anom, strat_annual_anom], axis=1, keys=['strat winter', 'strat spring', 'strat summer', 'strat fall', 'strat annual'])

## # Annual and seasonal anomalies
## strat_winter = df_strat[(df_strat.index.month>=1) & (df_strat.index.month<=3)].resample('As').mean()
## strat_spring = df_strat[(df_strat.index.month>=4) & (df_strat.index.month<=6)].resample('As').mean()
## strat_summer = df_strat[(df_strat.index.month>=7) & (df_strat.index.month<=9)].resample('As').mean()
## strat_fall = df_strat[(df_strat.index.month>=10) & (df_strat.index.month<=12)].resample('As').mean()
## strat_annual = df_strat.resample('As').mean()
## # Clim period
## strat_winter_clim = strat_clim_period[(strat_clim_period.index.month>=1) & (strat_clim_period.index.month<=3)].resample('As').mean()
## strat_spring_clim = strat_clim_period[(strat_clim_period.index.month>=4) & (strat_clim_period.index.month<=6)].resample('As').mean()
## strat_summer_clim = strat_clim_period[(strat_clim_period.index.month>=7) & (strat_clim_period.index.month<=9)].resample('As').mean()
## strat_fall_clim = strat_clim_period[(strat_clim_period.index.month>=10) & (strat_clim_period.index.month<=12)].resample('As').mean()
## strat_annual_clim = strat_clim_period.resample('As').mean()
## # Anomalies
## strat_std_anom_winter = (strat_winter - strat_winter_clim.mean()) / strat_winter_clim.std()
## strat_std_anom_spring = (strat_spring - strat_spring_clim.mean()) / strat_spring_clim.std()
## strat_std_anom_summer = (strat_summer - strat_summer_clim.mean()) / strat_summer_clim.std()
## strat_std_anom_fall = (strat_fall - strat_fall_clim.mean()) / strat_fall_clim.std()
## strat_std_anom_annual = (strat_annual - strat_annual_clim.mean()) / strat_annual_clim.std()
## # concat clim and anomaly
## strat_clim = pd.concat([strat_winter_clim, strat_spring_clim, strat_summer_clim, strat_fall_clim, strat_annual_clim], axis=1, keys=['strat winter', 'strat spring', 'strat summer', 'strat fall', 'strat annual'])
## strat_anom_std = pd.concat([strat_std_anom_winter, strat_std_anom_spring, strat_std_anom_summer, strat_std_anom_fall, strat_std_anom_annual], axis=1, keys=['strat winter', 'strat spring', 'strat summer', 'strat fall', 'strat annual'])


#### ------------- Build the scorecard ---------------- ####
# preamble
year_list = T_anom_std.index.year.astype('str')
year_list = [i[2:4] for i in year_list] # 2-digit year
year_list.append(r'$\rm \overline{x}$') # add 2 extra columns
year_list.append(r'sd')
year_list_FR = year_list[:]
year_list_FR[-1] = u'ET'
# Build the colormap
vmin = -3.49
vmax = 3.49
midpoint = 0
levels = np.linspace(vmin, vmax, 15)
midp = np.mean(np.c_[levels[:-1], levels[1:]], axis=1)
colvals = np.interp(midp, [vmin, midpoint, vmax], [-1, 0., 1])
normal = plt.Normalize(-3.49, 3.49)
reds = plt.cm.Reds(np.linspace(0,1, num=7))
blues = plt.cm.Blues_r(np.linspace(0,1, num=7))
whites = [(1,1,1,1)]*2
colors = np.vstack((blues[0:-1,:], whites, reds[1:,:]))
colors = np.concatenate([[colors[0,:]], colors, [colors[-1,:]]], 0)
cmap, norm = from_levels_and_colors(levels, colors, extend='both')
cmap_r, norm_r = from_levels_and_colors(levels, np.flipud(colors), extend='both')
# Common parameters
hcell, wcell = 0.5, 0.7
hpad, wpad = .8, 0 

### 1. Temperature
my_df = T_anom_std.T
my_df['MEAN'] = T_clim_period_annual.mean()
my_df['SD'] =  T_clim_period_annual.std()
#my_df.rename(columns={'MEAN': r'$\rm \overline{x}$', 'SD': r'sd'}, inplace=True)

# Get text values +  cell color
vals = np.around(my_df.values,1)
vals[vals==-0.] = 0.
vals_color = vals.copy()
vals_color[:,-1] = 0 # No color to last two columns (mean and STD)
vals_color[:,-2] = 0


nrows, ncols = my_df.index.size+1, my_df.columns.size
fig=plt.figure(figsize=(ncols*wcell+wpad, nrows*hcell+hpad)) # size+1 because of year's row
ax = fig.add_subplot(111)
ax.axis('off')
#do the table
header = ax.table(cellText=[['']],
                      colLabels=['-- Vertically averaged temperature --'],
                      loc='center'
                      )
header.set_fontsize(13)
the_table=ax.table(cellText=vals, rowLabels=my_df.index, colLabels=year_list,
                    loc='center', cellColours=cmap(norm(vals_color)), cellLoc='center',
                    bbox=[0, 0, 1, 0.5]
                    )
# change font color to white where needed:
the_table.auto_set_font_size(False)
the_table.set_fontsize(13)
table_props = the_table.properties()
table_cells = table_props['child_artists']
last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
#cell_dict = the_table.get_celld()
#for i in np.arange(1,4):
#    cell_dict[(i,)].set_width(0.3)
for key, cell in the_table.get_celld().items():
    cell_text = cell.get_text().get_text() 
    if is_number(cell_text) == False:
        pass
    elif key[0] == 0: #year's row = no color
        pass
    elif key[1] in last_columns:
         cell._text.set_color('darkslategray')
    elif (np.float(cell_text) <= -2) | (np.float(cell_text) >= 2) :
        cell._text.set_color('white')
    elif (cell_text=='nan'):
        cell._set_facecolor('darkgray')
        cell._text.set_color('darkgray')
plt.savefig("scorecards_s27_T.png", dpi=300)
os.system('convert -trim scorecards_s27_T.png scorecards_s27_T.png')

# French table
header = ax.table(cellText=[['']],
                      colLabels=[u'-- Moyenne verticale de température --'],
                      loc='center'
                      )
header.set_fontsize(13)
the_table=ax.table(cellText=vals, rowLabels=my_df.index, colLabels=year_list_FR,
                    loc='center', cellColours=cmap(norm(vals_color)), cellLoc='center',
                    bbox=[0, 0, 1, 0.5]
                    )
# change font color to white where needed:
the_table.auto_set_font_size(False)
the_table.set_fontsize(13)
table_props = the_table.properties()
table_cells = table_props['child_artists']
last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
for key, cell in the_table.get_celld().items():
    cell_text = cell.get_text().get_text()
    if is_number(cell_text) == False:
        pass
    elif key[0] == 0: #year's row = no color
        pass
    elif key[1] in last_columns:
         cell._text.set_color('darkslategray')
    elif (np.float(cell_text) <= -2) | (np.float(cell_text) >= 2) :
        cell._text.set_color('white')
    elif (cell_text=='nan'):
        cell._set_facecolor('darkgray')
        cell._text.set_color('darkgray')
plt.savefig("scorecards_s27_T_FR.png", dpi=300)
os.system('convert -trim scorecards_s27_T_FR.png scorecards_s27_T_FR.png')



### 2. Salinity
my_df = S_anom_std.T
my_df['MEAN'] = S_clim_period_annual.mean()
my_df['SD'] =  S_clim_period_annual.std()
#my_df.rename(columns={'MEAN': r'$\rm \overline{x}$', 'SD': r'sd'}, inplace=True)

# Get text values +  cell color
vals = np.around(my_df.values,1)
vals[vals==-0.] = 0.
vals_color = vals.copy()
vals_color[:,-1] = 0 # No color to last two columns (mean and STD)
vals_color[:,-2] = 0

nrows, ncols = my_df.index.size, my_df.columns.size
fig=plt.figure(figsize=(ncols*wcell+wpad, nrows*hcell+hpad))
ax = fig.add_subplot(111)
ax.axis('off')
#do the table
header = ax.table(cellText=[['']],
                      colLabels=['-- Vertically averaged salinity --'],
                      loc='center'
                      )
header.set_fontsize(13)
the_table=ax.table(cellText=vals, rowLabels=my_df.index, colLabels=None,
                    loc='center', cellColours=cmap(norm(vals_color)), cellLoc='center',
                    bbox=[0, 0, 1, 0.5]
                    )
# change font color to white where needed:
the_table.auto_set_font_size(False)
the_table.set_fontsize(13)
table_props = the_table.properties()
table_cells = table_props['child_artists']
last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
#cell_dict = the_table.get_celld()
#for i in np.arange(1,3):
#    cell_dict[(i,-1)].set_width(0.5)
for key, cell in the_table.get_celld().items():
    cell_text = cell.get_text().get_text() 
    if is_number(cell_text) == False:
        pass
    ## elif key[0] == 0: #year's row = no color
    ##     pass
    elif key[1] in last_columns:
         cell._text.set_color('darkslategray')
    elif (np.float(cell_text) <= -2) | (np.float(cell_text) >= 2) :
        cell._text.set_color('white')
    elif (cell_text=='nan'):
        cell._set_facecolor('darkgray')
        cell._text.set_color('darkgray')
plt.savefig("scorecards_s27_S.png", dpi=300)
os.system('convert -trim scorecards_s27_S.png scorecards_s27_S.png')

# French table
header = ax.table(cellText=[['']],
                      colLabels=[u'-- Moyenne verticale de salinité --'],
                      loc='center'
                      )
header.set_fontsize(13)
the_table=ax.table(cellText=vals, rowLabels=my_df.index, colLabels=None,
                    loc='center', cellColours=cmap(norm(vals_color)), cellLoc='center',
                    bbox=[0, 0, 1, 0.5]
                    )
# change font color to white where needed:
the_table.auto_set_font_size(False)
the_table.set_fontsize(13)
table_props = the_table.properties()
table_cells = table_props['child_artists']
last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
for key, cell in the_table.get_celld().items():
    cell_text = cell.get_text().get_text()
    if is_number(cell_text) == False:
        pass
    ## elif key[0] == 0: #year's row = no color
    ##     pass
    elif key[1] in last_columns:
         cell._text.set_color('darkslategray')
    elif (np.float(cell_text) <= -2) | (np.float(cell_text) >= 2) :
        cell._text.set_color('white')
    elif (cell_text=='nan'):
        cell._set_facecolor('darkgray')
        cell._text.set_color('darkgray')
plt.savefig("scorecards_s27_S_FR.png", dpi=300)
os.system('convert -trim scorecards_s27_S_FR.png scorecards_s27_S_FR.png')


### 3. CIL
my_df = CIL_anom_std.T
my_df['MEAN'] = CIL_clim_period.mean()
my_df['SD'] =  CIL_clim_period.std()

# Get text values +  cell color
vals = np.around(my_df.values,1)
vals[vals==-0.] = 0.
vals_color = vals.copy()
vals_color[:,-1] = 0 # No color to last two columns (mean and STD)
vals_color[:,-2] = 0
# Reverse some CIL values for colormap (thickness & core depth)
vals_color[2,:] = vals_color[2,:]*-1
vals_color[3,:] = vals_color[3,:]*-1

nrows, ncols = my_df.index.size, my_df.columns.size
fig=plt.figure(figsize=(ncols*wcell+wpad, nrows*hcell+hpad))
ax = fig.add_subplot(111)
ax.axis('off')
#do the table
header = ax.table(cellText=[['']],
                      colLabels=['-- Cold intermediate layer (CIL) properties --'],
                      loc='center'
                      )
header.set_fontsize(13)
the_table=ax.table(cellText=vals, rowLabels=my_df.index, colLabels=None,
                    loc='center', cellColours=cmap(norm(vals_color)), cellLoc='center',
                    bbox=[0, 0, 1, 0.5]
                    )
# change font color to white where needed:
the_table.auto_set_font_size(False)
the_table.set_fontsize(13)
table_props = the_table.properties()
table_cells = table_props['child_artists']
last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
for key, cell in the_table.get_celld().items():
    cell_text = cell.get_text().get_text() 
    if is_number(cell_text) == False:
        pass
    ## elif key[0] == 0: #year's row = no color
    ##     pass
    elif key[1] in last_columns:
         cell._text.set_color('darkslategray')
    elif (np.float(cell_text) <= -2) | (np.float(cell_text) >= 2) :
        cell._text.set_color('white')
    elif (cell_text=='nan'):
        cell._set_facecolor('darkgray')
        cell._text.set_color('darkgray')
plt.savefig("scorecards_s27_CIL.png", dpi=300)
os.system('convert -trim scorecards_s27_CIL.png scorecards_s27_CIL.png')

# French table
my_df.index = [u'CIF temp', u'CIF T coeur', u'CIF prof coeur', u'CIF épaisseur']

header = ax.table(cellText=[['']],
                      colLabels=[u'-- Propriétés de la couche intermédiaire froide (CIF) --'],
                      loc='center'
                      )
header.set_fontsize(13)
the_table=ax.table(cellText=vals, rowLabels=my_df.index, colLabels=None,
                    loc='center', cellColours=cmap(norm(vals_color)), cellLoc='center',
                    bbox=[0, 0, 1, 0.5]
                    )
# change font color to white where needed:
the_table.auto_set_font_size(False)
the_table.set_fontsize(13)
table_props = the_table.properties()
table_cells = table_props['child_artists']
last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
for key, cell in the_table.get_celld().items():
    cell_text = cell.get_text().get_text()
    if is_number(cell_text) == False:
        pass
    ## elif key[0] == 0: #year's row = no color
    ##     pass
    elif key[1] in last_columns:
         cell._text.set_color('darkslategray')
    elif (np.float(cell_text) <= -2) | (np.float(cell_text) >= 2) :
        cell._text.set_color('white')
    elif (cell_text=='nan'):
        cell._set_facecolor('darkgray')
        cell._text.set_color('darkgray')
plt.savefig("scorecards_s27_CIL_FR.png", dpi=300)
os.system('convert -trim scorecards_s27_CIL_FR.png scorecards_s27_CIL_FR.png')


### 4. MLD
my_df = MLD_anom_std.T
my_df['MEAN'] = MLD_clim.mean()
my_df['SD'] =  MLD_clim.std()

# Get text values +  cell color
vals = np.around(my_df.values,1)
vals[vals==-0.] = 0.
vals_color = vals.copy()
vals_color[:,-1] = 0 # No color to last two columns (mean and STD)
vals_color[:,-2] = 0

nrows, ncols = my_df.index.size, my_df.columns.size
fig=plt.figure(figsize=(ncols*wcell+wpad, nrows*hcell+hpad))
ax = fig.add_subplot(111)
ax.axis('off')
#do the table
header = ax.table(cellText=[['']],
                      colLabels=['-- Mixed layer depth (MLD) --'],
                      loc='center'
                      )
header.set_fontsize(13)
the_table=ax.table(cellText=vals, rowLabels=my_df.index, colLabels=None,
                    loc='center', cellColours=cmap(norm(vals_color)), cellLoc='center',
                    bbox=[0, 0, 1, 0.5]
                    )
# change font color to white where needed:
the_table.auto_set_font_size(False)
the_table.set_fontsize(13)
table_props = the_table.properties()
table_cells = table_props['child_artists']
last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
for key, cell in the_table.get_celld().items():
    cell_text = cell.get_text().get_text() 
    if is_number(cell_text) == False:
        pass
    ## elif key[0] == 0: #year's row = no color
    ##     pass
    elif key[1] in last_columns:
         cell._text.set_color('darkslategray')
    elif (np.float(cell_text) <= -2) | (np.float(cell_text) >= 2) :
        cell._text.set_color('white')
    elif (cell_text=='nan'):
        cell._set_facecolor('darkgray')
        cell._text.set_color('darkgray')
plt.savefig("scorecards_s27_MLD.png", dpi=300)
os.system('convert -trim scorecards_s27_MLD.png scorecards_s27_MLD.png')

# French table
my_df.index = [u'PCM hiver', u'PCM printemps', u'PCM été', u'PCM automne', u'PCM annuel']

header = ax.table(cellText=[['']],
                      colLabels=[u'-- Profondeur de la couche de mélange (PCM) --'],
                      loc='center'
                      )
header.set_fontsize(13)
the_table=ax.table(cellText=vals, rowLabels=my_df.index, colLabels=None,
                    loc='center', cellColours=cmap(norm(vals_color)), cellLoc='center',
                    bbox=[0, 0, 1, 0.5]
                    )
# change font color to white where needed:
the_table.auto_set_font_size(False)
the_table.set_fontsize(13)
table_props = the_table.properties()
table_cells = table_props['child_artists']
last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
for key, cell in the_table.get_celld().items():
    cell_text = cell.get_text().get_text()
    if is_number(cell_text) == False:
        pass
    ## elif key[0] == 0: #year's row = no color
    ##     pass
    elif key[1] in last_columns:
         cell._text.set_color('darkslategray')
    elif (np.float(cell_text) <= -2) | (np.float(cell_text) >= 2) :
        cell._text.set_color('white')
    elif (cell_text=='nan'):
        cell._set_facecolor('darkgray')
        cell._text.set_color('darkgray')
plt.savefig("scorecards_s27_MLD_FR.png", dpi=300)
os.system('convert -trim scorecards_s27_MLD_FR.png scorecards_s27_MLD_FR.png')


### 5. strat
my_df = strat_anom_std.T
my_df['MEAN'] = strat_clim.mean()
my_df['SD'] =  strat_clim.std()

# Get text values +  cell color
vals = np.around(my_df.values,1)
vals[vals==-0.] = 0.
vals_color = vals.copy()
vals_color[:,-1] = 0 # No color to last two columns (mean and STD)
vals_color[:,-2] = 0
# scietnific notation for mean and std
for i in np.arange(0,5):
    vals[i,-1] = np.around(my_df.values[i,-1],3)
    vals[i,-2] = np.around(my_df.values[i,-2],3)
  

nrows, ncols = my_df.index.size, my_df.columns.size
fig=plt.figure(figsize=(ncols*wcell+wpad, nrows*hcell+hpad))
ax = fig.add_subplot(111)
ax.axis('off')
#do the table
header = ax.table(cellText=[['']],
                      colLabels=['-- Stratification --'],
                      loc='center'
                      )
header.set_fontsize(13)
the_table=ax.table(cellText=vals, rowLabels=my_df.index, colLabels=None,
                    loc='center', cellColours=cmap(norm(vals_color)), cellLoc='center',
                    bbox=[0, 0, 1, 0.5]
                    )
# change font color to white where needed:
the_table.auto_set_font_size(False)
the_table.set_fontsize(13)
table_props = the_table.properties()
table_cells = table_props['child_artists']
last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
for key, cell in the_table.get_celld().items():
    cell_text = cell.get_text().get_text() 
    if is_number(cell_text) == False:
        pass
    ## elif key[0] == 0: #year's row = no color
    ##     pass
    elif key[1] in last_columns:
         cell._text.set_color('darkslategray')
    elif (np.float(cell_text) <= -2) | (np.float(cell_text) >= 2) :
        cell._text.set_color('white')
    elif (cell_text=='nan'):
        cell._set_facecolor('darkgray')
        cell._text.set_color('darkgray')
plt.savefig("scorecards_s27_strat.png", dpi=300)
os.system('convert -trim scorecards_s27_strat.png scorecards_s27_strat.png')

# French table
my_df.index = [u'strat. hiver', u'strat. printemps', u'strat. été', u'strat. automne', u'strat. annuel']

header = ax.table(cellText=[['']],
                      colLabels=[u'-- Stratification --'],
                      loc='center'
                      )
header.set_fontsize(13)
the_table=ax.table(cellText=vals, rowLabels=my_df.index, colLabels=None,
                    loc='center', cellColours=cmap(norm(vals_color)), cellLoc='center',
                    bbox=[0, 0, 1, 0.5]
                    )
# change font color to white where needed:
the_table.auto_set_font_size(False)
the_table.set_fontsize(13)
table_props = the_table.properties()
table_cells = table_props['child_artists']
last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
for key, cell in the_table.get_celld().items():
    cell_text = cell.get_text().get_text()
    if is_number(cell_text) == False:
        pass
    ## elif key[0] == 0: #year's row = no color
    ##     pass
    elif key[1] in last_columns:
         cell._text.set_color('darkslategray')
    elif (np.float(cell_text) <= -2) | (np.float(cell_text) >= 2) :
        cell._text.set_color('white')
    elif (cell_text=='nan'):
        cell._set_facecolor('darkgray')
        cell._text.set_color('darkgray')
plt.savefig("scorecards_s27_strat_FR.png", dpi=300)
os.system('convert -trim scorecards_s27_strat_FR.png scorecards_s27_strat_FR.png')


#6. Merge all together
# English
#os.system('montage  -gravity west scorecards_s27_T.png scorecards_s27_S.png scorecards_s27_CIL.png scorecards_s27_MLD.png scorecards_s27_strat.png -tile 1x5 -geometry +1+1  -background white scorecards_s27.png') 
os.system('convert scorecards_s27_T.png scorecards_s27_S.png scorecards_s27_CIL.png scorecards_s27_MLD.png scorecards_s27_strat.png  -gravity East -append -geometry +1+1 scorecards_s27.png')

# French
#os.system('montage   -gravity east scorecards_s27_T_FR.png scorecards_s27_S_FR.png scorecards_s27_CIL_FR.png scorecards_s27_MLD_FR.png scorecards_s27_strat_FR.png -tile 1x5 -geometry +1+1  -background white scorecards_s27_FR.png') 
os.system('convert scorecards_s27_T_FR.png scorecards_s27_S_FR.png scorecards_s27_CIL_FR.png scorecards_s27_MLD_FR.png scorecards_s27_strat_FR.png  -gravity East -append -geometry +1+1 scorecards_s27_FR.png')



