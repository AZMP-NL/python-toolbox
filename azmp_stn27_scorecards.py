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

year_clim = [1991, 2020]
years = [1981, 2021]

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

#Set the year after the beginning of 1950
df_temp = df_temp[df_temp.index.year>=1950]

## - Temperature - ##
# stack year, month
ts_stack = df_temp.groupby([(df_temp.index.year),(df_temp.index.month)]).mean()
# rename index
ts_stack.index = ts_stack.index.set_names(['year', 'month'])
# Isolate period for climatology
T_clim_period = df_temp[(df_temp.index.year>=year_clim[0]) & (df_temp.index.year<=year_clim[1])]
# Calculate Monthly clim
monthly_clim_mean = T_clim_period.groupby(T_clim_period.index.month).mean()
monthly_clim_stdv = T_clim_period.groupby(T_clim_period.index.month).std()
ts_monthly_clim = monthly_clim_mean.mean(axis=1)
ts_monthly_std = monthly_clim_stdv.mean(axis=1)
# Tile the climatology for however many years are present
ts_years = len(df_temp.index.year.unique())
monthly_clim_mean = pd.concat([monthly_clim_mean] * ts_years) #<-- good job, I like this!
monthly_clim_stdv = pd.concat([monthly_clim_stdv] * ts_years)
monthly_clim_mean = monthly_clim_mean.iloc[:ts_stack.index.size,:]
monthly_clim_stdv = monthly_clim_stdv.iloc[:ts_stack.index.size,:]
# Set multi-index to clim (using ts_stack index)
monthly_clim_mean.set_index(ts_stack.index, inplace=True)
monthly_clim_stdv.set_index(ts_stack.index, inplace=True)
# Calculate anomalies (mi = miltiindex)
anom_mi = ts_stack-monthly_clim_mean
std_anom_mi = (ts_stack-monthly_clim_mean) / monthly_clim_stdv

# vertically-averaged anomalies to three different groups
T_0_btm_anom_monthly = anom_mi.mean(axis=1)
T_0_btm_std_anom_monthly = std_anom_mi.mean(axis=1)
T_0_50_anom_monthly = anom_mi[anom_mi.columns[(anom_mi.columns<=50)]]
T_0_50_std_anom_monthly = std_anom_mi[std_anom_mi.columns[(std_anom_mi.columns<=50)]]
T_150_btm_anom_monthly = anom_mi[anom_mi.columns[(anom_mi.columns>=150)]]
T_150_btm_std_anom_monthly = std_anom_mi[std_anom_mi.columns[(std_anom_mi.columns>=150)]]

#Convert to yearly
Tave_anom = T_0_btm_anom_monthly.unstack().mean(axis=1)
Tave_anom_std = T_0_btm_std_anom_monthly.unstack().mean(axis=1)
Tsurf_anom = T_0_50_anom_monthly.unstack().mean(axis=1)
Tsurf_anom_std = T_0_50_std_anom_monthly.unstack().mean(axis=1)
Tbot_anom = T_150_btm_anom_monthly.unstack().mean(axis=1)
Tbot_anom_std = T_150_btm_std_anom_monthly.unstack().mean(axis=1)
Tave_anom.index = pd.to_datetime(Tave_anom.index, format='%Y')
Tsurf_anom.index = pd.to_datetime(Tsurf_anom.index, format='%Y')
Tbot_anom.index = pd.to_datetime(Tbot_anom.index, format='%Y')
Tave_anom_std.index = pd.to_datetime(Tave_anom_std.index, format='%Y')
Tsurf_anom_std.index = pd.to_datetime(Tsurf_anom_std.index, format='%Y')
Tbot_anom_std.index = pd.to_datetime(Tbot_anom_std.index, format='%Y')

#Set the year after the beginning of 1950
df_sal = df_sal[df_sal.index.year>=1950]

## - Salinity - ##
# stack year, month
ts_stack = df_sal.groupby([(df_sal.index.year),(df_sal.index.month)]).mean()
# rename index
ts_stack.index = ts_stack.index.set_names(['year', 'month'])
# Isolate period for climatology
S_clim_period = df_sal[(df_sal.index.year>=year_clim[0]) & (df_sal.index.year<=year_clim[1])]
# Calculate Monthly clim
monthly_clim_mean = S_clim_period.groupby(S_clim_period.index.month).mean()
monthly_clim_stdv = S_clim_period.groupby(S_clim_period.index.month).std()
ts_monthly_clim = monthly_clim_mean.mean(axis=1)
ts_monthly_std = monthly_clim_stdv.mean(axis=1)
# Tile the climatology for however many years are present
ts_years = len(df_sal.index.year.unique())
monthly_clim_mean = pd.concat([monthly_clim_mean] * ts_years) #<-- good job, I like this!
monthly_clim_stdv = pd.concat([monthly_clim_stdv] * ts_years)
monthly_clim_mean = monthly_clim_mean.iloc[:ts_stack.index.size,:]
monthly_clim_stdv = monthly_clim_stdv.iloc[:ts_stack.index.size,:]
# Set multi-index to clim (using ts_stack index)
monthly_clim_mean.set_index(ts_stack.index, inplace=True)
monthly_clim_stdv.set_index(ts_stack.index, inplace=True)
# Calculate anomalies (mi = miltiindex)
anom_mi = ts_stack-monthly_clim_mean
std_anom_mi = (ts_stack-monthly_clim_mean) / monthly_clim_stdv

# vertically-averaged anomalies to three different groups
S_0_btm_anom_monthly = anom_mi.mean(axis=1)
S_0_btm_std_anom_monthly = std_anom_mi.mean(axis=1)
S_0_50_anom_monthly = anom_mi[anom_mi.columns[(anom_mi.columns<=50)]]
S_0_50_std_anom_monthly = std_anom_mi[std_anom_mi.columns[(std_anom_mi.columns<=50)]]
S_150_btm_anom_monthly = anom_mi[anom_mi.columns[(anom_mi.columns>=150)]]
S_150_btm_std_anom_monthly = std_anom_mi[std_anom_mi.columns[(std_anom_mi.columns>=150)]]

#Convert to yearly
Save_anom = S_0_btm_anom_monthly.unstack().mean(axis=1)
Save_anom_std = S_0_btm_std_anom_monthly.unstack().mean(axis=1)
Ssurf_anom = S_0_50_anom_monthly.unstack().mean(axis=1)
Ssurf_anom_std = S_0_50_std_anom_monthly.unstack().mean(axis=1)
Sbot_anom = S_150_btm_anom_monthly.unstack().mean(axis=1)
Sbot_anom_std = S_150_btm_std_anom_monthly.unstack().mean(axis=1)
Save_anom.index = pd.to_datetime(Save_anom.index, format='%Y')
Ssurf_anom.index = pd.to_datetime(Ssurf_anom.index, format='%Y')
Sbot_anom.index = pd.to_datetime(Sbot_anom.index, format='%Y')
Save_anom_std.index = pd.to_datetime(Save_anom_std.index, format='%Y')
Ssurf_anom_std.index = pd.to_datetime(Ssurf_anom_std.index, format='%Y')
Sbot_anom_std.index = pd.to_datetime(Sbot_anom_std.index, format='%Y')



# merge anomalies
T_anom = pd.concat([Tave_anom,Tsurf_anom,Tbot_anom], axis=1, keys=['Temp 0-176m','Temp 0-50m','Temp 150-176m'])
T_anom_std = pd.concat([Tave_anom_std,Tsurf_anom_std,Tbot_anom_std], axis=1, keys=['Temp 0-176m','Temp 0-50m','Temp 150-176m'])
T_anom = T_anom[T_anom.index.year>=1947]
T_anom_std = T_anom_std[T_anom_std.index.year>=1947]
S_anom = pd.concat([Save_anom,Ssurf_anom,Sbot_anom], axis=1, keys=['Sal 0-176m','Sal 0-50m','Sal 150-176m ${~}$'])
S_anom_std = pd.concat([Save_anom_std,Ssurf_anom_std,Sbot_anom_std], axis=1, keys=['Sal 0-176m','Sal 0-50m','Sal 150-176m ${~}$'])
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

#### ------------- CIL ---------------- ####
# Load pickled data
df_CIL = pd.read_pickle('S27_CIL_summer_stats.pkl')
df_CIL = df_CIL[df_CIL.index.year>=years[0]]
#df_CIL.dropna(inplace=True)
# Compute anomalies 
CIL_clim_period = df_CIL[(df_CIL.index.year>=year_clim[0]) & (df_CIL.index.year<=year_clim[1])]
# Anomaly
CIL_anom = (df_CIL - CIL_clim_period.mean())
CIL_anom_std = CIL_anom / CIL_clim_period.std()

#### ------------- MLD ---------------- ####
# Load pickled data
mld_monthly = pd.read_pickle('S27_MLD_monthly.pkl')
mld = mld_monthly
# flag some data:
mld[mld.index=='2019-03-15']=np.nan
mld[mld.index=='1980-03-15']=np.nan
mld = mld[mld.index.year>=years[0]]
# limit to after 1990
#mld[mld.index<'1990-01-01']=np.nan
# Annual and seasonal anomalies (average of monthly anomalies rather than anomalies of monthly values))
#stack months
mld_stack = mld.groupby([(mld.index.year),(mld.index.month)]).mean()
mld_unstack = mld_stack.unstack()
# compute clim
mld_clim_period = mld[(mld.index.year>=year_clim[0]) & (mld.index.year<=year_clim[1])]
mld_monthly_stack = mld_clim_period.groupby([(mld_clim_period.index.year),(mld_clim_period.index.month)]).mean()
mld_monthly_clim = mld_monthly_stack.groupby(level=1).mean()
mld_monthly_std = mld_monthly_stack.groupby(level=1).std()
monthly_anom = mld_unstack - mld_monthly_clim 
monthly_stdanom = (mld_unstack - mld_monthly_clim) /  mld_monthly_std
# Seasonal and annual anomalies
mld_winter_anom = monthly_stdanom[[1,2,3]].mean(axis=1)
mld_spring_anom = monthly_stdanom[[4,5,6]].mean(axis=1)
mld_summer_anom = monthly_stdanom[[7,8,9]].mean(axis=1)
mld_fall_anom = monthly_stdanom[[10,11,12]].mean(axis=1)
mld_annual_anom = monthly_stdanom.mean(axis=1)
# Mean values for climatology period (for last columns)
mld_winter_clim = mld_monthly_clim[[1,2,3]]
mld_spring_clim = mld_monthly_clim[[4,5,6]]
mld_summer_clim = mld_monthly_clim[[7,8,9]]
mld_fall_clim = mld_monthly_clim[[10,11,12]]
mld_annual_clim = mld_monthly_clim
# concat clim and anomaly
mld_clim = pd.concat([mld_winter_clim, mld_spring_clim, mld_summer_clim, mld_fall_clim, mld_annual_clim], axis=1, keys=['MLD winter', 'MLD spring', 'MLD summer', 'MLD fall', 'MLD annual'])
mld_anom_std = pd.concat([mld_winter_anom, mld_spring_anom, mld_summer_anom, mld_fall_anom, mld_annual_anom], axis=1, keys=['MLD winter', 'MLD spring', 'MLD summer', 'MLD fall', 'MLD annual'])
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


#### ------------- Stratification (0-50m)---------------- ####
# Load pickled data
strat_monthly = pd.read_pickle('S27_stratif_0-50_monthly.pkl')
strat = strat = strat_monthly
# flag some data:
strat[strat.index=='2019-03-15']=np.nan
strat[strat.index=='1980-03-15']=np.nan
strat = strat[strat.index.year>=years[0]]
# limit to after 1990
#strat[strat.index<'1990-01-01']=np.nan
# Annual and seasonal anomalies (average of monthly anomalies rather than anomalies of monthly values))
#stack months
strat_stack = strat.groupby([(strat.index.year),(strat.index.month)]).mean()
strat_unstack = strat_stack.unstack()
# compute clim
strat_clim_period = strat[(strat.index.year>=year_clim[0]) & (strat.index.year<=year_clim[1])]
strat_monthly_stack = strat_clim_period.groupby([(strat_clim_period.index.year),(strat_clim_period.index.month)]).mean()
strat_monthly_clim = strat_monthly_stack.groupby(level=1).mean()
strat_monthly_std = strat_monthly_stack.groupby(level=1).std()
monthly_anom = strat_unstack - strat_monthly_clim 
monthly_stdanom = (strat_unstack - strat_monthly_clim) /  strat_monthly_std

# Seasonal and annual anomalies
strat_winter_anom = monthly_stdanom[[1,2,3]].mean(axis=1)
strat_spring_anom = monthly_stdanom[[4,5,6]].mean(axis=1)
strat_summer_anom = monthly_stdanom[[7,8,9]].mean(axis=1)
strat_fall_anom = monthly_stdanom[[10,11,12]].mean(axis=1)
strat_annual_anom = monthly_stdanom.mean(axis=1)
# Mean values for climatology period (for last columns)
strat_winter_clim = strat_monthly_clim[[1,2,3]]
strat_spring_clim = strat_monthly_clim[[4,5,6]]
strat_summer_clim = strat_monthly_clim[[7,8,9]]
strat_fall_clim = strat_monthly_clim[[10,11,12]]
strat_annual_clim = strat_monthly_clim
# concat clim and anomaly
strat_clim_0to50 = pd.concat([strat_winter_clim, strat_spring_clim, strat_summer_clim, strat_fall_clim, strat_annual_clim], axis=1, keys=['strat winter', 'strat spring', 'strat summer', 'strat fall', 'strat annual'])
strat_anom_std_0to50 = pd.concat([strat_winter_anom, strat_spring_anom, strat_summer_anom, strat_fall_anom, strat_annual_anom], axis=1, keys=['strat winter', 'strat spring', 'strat summer', 'strat fall', 'strat annual'])



#### ------------- Stratification (10-150m)---------------- ####
# Load pickled data
strat_monthly = pd.read_pickle('S27_stratif_10-150_monthly.pkl')
strat = strat = strat_monthly
# flag some data:
strat[strat.index=='2019-03-15']=np.nan
strat[strat.index=='1980-03-15']=np.nan
strat = strat[strat.index.year>=years[0]]
# limit to after 1990
#strat[strat.index<'1990-01-01']=np.nan
# Annual and seasonal anomalies (average of monthly anomalies rather than anomalies of monthly values))
#stack months
strat_stack = strat.groupby([(strat.index.year),(strat.index.month)]).mean()
strat_unstack = strat_stack.unstack()
# compute clim
strat_clim_period = strat[(strat.index.year>=year_clim[0]) & (strat.index.year<=year_clim[1])]
strat_monthly_stack = strat_clim_period.groupby([(strat_clim_period.index.year),(strat_clim_period.index.month)]).mean()
strat_monthly_clim = strat_monthly_stack.groupby(level=1).mean()
strat_monthly_std = strat_monthly_stack.groupby(level=1).std()
monthly_anom = strat_unstack - strat_monthly_clim 
monthly_stdanom = (strat_unstack - strat_monthly_clim) /  strat_monthly_std

# Seasonal and annual anomalies
strat_winter_anom = monthly_stdanom[[1,2,3]].mean(axis=1)
strat_spring_anom = monthly_stdanom[[4,5,6]].mean(axis=1)
strat_summer_anom = monthly_stdanom[[7,8,9]].mean(axis=1)
strat_fall_anom = monthly_stdanom[[10,11,12]].mean(axis=1)
strat_annual_anom = monthly_stdanom.mean(axis=1)
# Mean values for climatology period (for last columns)
strat_winter_clim = strat_monthly_clim[[1,2,3]]
strat_spring_clim = strat_monthly_clim[[4,5,6]]
strat_summer_clim = strat_monthly_clim[[7,8,9]]
strat_fall_clim = strat_monthly_clim[[10,11,12]]
strat_annual_clim = strat_monthly_clim
# concat clim and anomaly
strat_clim_10to150 = pd.concat([strat_winter_clim, strat_spring_clim, strat_summer_clim, strat_fall_clim, strat_annual_clim], axis=1, keys=['strat winter', 'strat spring', 'strat summer', 'strat fall', 'strat annual'])
strat_anom_std_10to150 = pd.concat([strat_winter_anom, strat_spring_anom, strat_summer_anom, strat_fall_anom, strat_annual_anom], axis=1, keys=['strat winter', 'strat spring', 'strat summer', 'strat fall', 'strat annual'])



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
#table_cells = table_props['child_artists']
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
    elif (float(cell_text) <= -2) | (float(cell_text) >= 2) :
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
#table_cells = table_props['child_artists']
last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
for key, cell in the_table.get_celld().items():
    cell_text = cell.get_text().get_text()
    if is_number(cell_text) == False:
        pass
    elif key[0] == 0: #year's row = no color
        pass
    elif key[1] in last_columns:
         cell._text.set_color('darkslategray')
    elif (float(cell_text) <= -2) | (float(cell_text) >= 2) :
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
#table_cells = table_props['child_artists']
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
    elif (float(cell_text) <= -2) | (float(cell_text) >= 2) :
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
#table_cells = table_props['child_artists']
last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
for key, cell in the_table.get_celld().items():
    cell_text = cell.get_text().get_text()
    if is_number(cell_text) == False:
        pass
    ## elif key[0] == 0: #year's row = no color
    ##     pass
    elif key[1] in last_columns:
         cell._text.set_color('darkslategray')
    elif (float(cell_text) <= -2) | (float(cell_text) >= 2) :
        cell._text.set_color('white')
    elif (cell_text=='nan'):
        cell._set_facecolor('darkgray')
        cell._text.set_color('darkgray')
plt.savefig("scorecards_s27_S_FR.png", dpi=300)
os.system('convert -trim scorecards_s27_S_FR.png scorecards_s27_S_FR.png')


### 3. CIL

#Remove CIL core depth for now
my_df = CIL_anom_std.drop('CIL core depth', axis=1).T
#my_df = CIL_anom_std.T
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
#Add back if you're doing depth
#vals_color[3,:] = vals_color[3,:]*-1

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
#table_cells = table_props['child_artists']
last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
for key, cell in the_table.get_celld().items():
    cell_text = cell.get_text().get_text() 
    if is_number(cell_text) == False:
        pass
    ## elif key[0] == 0: #year's row = no color
    ##     pass
    elif key[1] in last_columns:
         cell._text.set_color('darkslategray')
    elif (float(cell_text) <= -2) | (float(cell_text) >= 2) :
        cell._text.set_color('white')
    elif (cell_text=='nan'):
        cell._set_facecolor('darkgray')
        cell._text.set_color('darkgray')
plt.savefig("scorecards_s27_CIL.png", dpi=300)
os.system('convert -trim scorecards_s27_CIL.png scorecards_s27_CIL.png')

# French table
my_df.index = [u'CIF temp', u'CIF T coeur', u'CIF épaisseur']
#Add back if you're doing depth
#my_df.index = [u'CIF temp', u'CIF T coeur', u'CIF prof coeur', u'CIF épaisseur']

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
#table_cells = table_props['child_artists']
last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
for key, cell in the_table.get_celld().items():
    cell_text = cell.get_text().get_text()
    if is_number(cell_text) == False:
        pass
    ## elif key[0] == 0: #year's row = no color
    ##     pass
    elif key[1] in last_columns:
         cell._text.set_color('darkslategray')
    elif (float(cell_text) <= -2) | (float(cell_text) >= 2) :
        cell._text.set_color('white')
    elif (cell_text=='nan'):
        cell._set_facecolor('darkgray')
        cell._text.set_color('darkgray')
plt.savefig("scorecards_s27_CIL_FR.png", dpi=300)
os.system('convert -trim scorecards_s27_CIL_FR.png scorecards_s27_CIL_FR.png')


### 4. MLD
my_df = mld_anom_std.T
my_df['MEAN'] = mld_clim.mean()
my_df['SD'] =  mld_clim.std()

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
#table_cells = table_props['child_artists']
last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
for key, cell in the_table.get_celld().items():
    cell_text = cell.get_text().get_text() 
    if is_number(cell_text) == False:
        pass
    ## elif key[0] == 0: #year's row = no color
    ##     pass
    elif key[1] in last_columns:
         cell._text.set_color('darkslategray')
    elif (float(cell_text) <= -2) | (float(cell_text) >= 2) :
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
#table_cells = table_props['child_artists']
last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
for key, cell in the_table.get_celld().items():
    cell_text = cell.get_text().get_text()
    if is_number(cell_text) == False:
        pass
    ## elif key[0] == 0: #year's row = no color
    ##     pass
    elif key[1] in last_columns:
         cell._text.set_color('darkslategray')
    elif (float(cell_text) <= -2) | (float(cell_text) >= 2) :
        cell._text.set_color('white')
    elif (cell_text=='nan'):
        cell._set_facecolor('darkgray')
        cell._text.set_color('darkgray')
plt.savefig("scorecards_s27_MLD_FR.png", dpi=300)
os.system('convert -trim scorecards_s27_MLD_FR.png scorecards_s27_MLD_FR.png')


### 5. strat (0-50m)
my_df = strat_anom_std_0to50.T
my_df['MEAN'] = strat_clim_0to50.mean()
my_df['SD'] =  strat_clim_0to50.std()

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
                      colLabels=['-- Stratification (0-50m) --'],
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
#table_cells = table_props['child_artists']
last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
for key, cell in the_table.get_celld().items():
    cell_text = cell.get_text().get_text() 
    if is_number(cell_text) == False:
        pass
    ## elif key[0] == 0: #year's row = no color
    ##     pass
    elif key[1] in last_columns:
         cell._text.set_color('darkslategray')
    elif (float(cell_text) <= -2) | (float(cell_text) >= 2) :
        cell._text.set_color('white')
    elif (cell_text=='nan'):
        cell._set_facecolor('darkgray')
        cell._text.set_color('darkgray')
plt.savefig("scorecards_s27_strat_0-50m.png", dpi=300)
os.system('convert -trim scorecards_s27_strat_0-50m.png scorecards_s27_strat_0-50m.png')

# French table
my_df.index = [u'strat. hiver', u'strat. printemps', u'strat. été', u'strat. automne', u'strat. annuel']

header = ax.table(cellText=[['']],
                      colLabels=[u'-- Stratification (0-50m) --'],
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
#table_cells = table_props['child_artists']
last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
for key, cell in the_table.get_celld().items():
    cell_text = cell.get_text().get_text()
    if is_number(cell_text) == False:
        pass
    ## elif key[0] == 0: #year's row = no color
    ##     pass
    elif key[1] in last_columns:
         cell._text.set_color('darkslategray')
    elif (float(cell_text) <= -2) | (float(cell_text) >= 2) :
        cell._text.set_color('white')
    elif (cell_text=='nan'):
        cell._set_facecolor('darkgray')
        cell._text.set_color('darkgray')
plt.savefig("scorecards_s27_strat_0-50m_FR.png", dpi=300)
os.system('convert -trim scorecards_s27_strat_0-50m_FR.png scorecards_s27_strat_0-50m_FR.png')


### 6. strat (10-150m)
my_df = strat_anom_std_10to150.T
my_df['MEAN'] = strat_clim_10to150.mean()
my_df['SD'] =  strat_clim_10to150.std()

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
                      colLabels=['-- Stratification (10-150m) --'],
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
#table_cells = table_props['child_artists']
last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
for key, cell in the_table.get_celld().items():
    cell_text = cell.get_text().get_text() 
    if is_number(cell_text) == False:
        pass
    ## elif key[0] == 0: #year's row = no color
    ##     pass
    elif key[1] in last_columns:
         cell._text.set_color('darkslategray')
    elif (float(cell_text) <= -2) | (float(cell_text) >= 2) :
        cell._text.set_color('white')
    elif (cell_text=='nan'):
        cell._set_facecolor('darkgray')
        cell._text.set_color('darkgray')
plt.savefig("scorecards_s27_strat_10-150m.png", dpi=300)
os.system('convert -trim scorecards_s27_strat_10-150m.png scorecards_s27_strat_10-150m.png')

# French table
my_df.index = [u'strat. hiver', u'strat. printemps', u'strat. été', u'strat. automne', u'strat. annuel']

header = ax.table(cellText=[['']],
                      colLabels=[u'-- Stratification (10-150m) --'],
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
#table_cells = table_props['child_artists']
last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
for key, cell in the_table.get_celld().items():
    cell_text = cell.get_text().get_text()
    if is_number(cell_text) == False:
        pass
    ## elif key[0] == 0: #year's row = no color
    ##     pass
    elif key[1] in last_columns:
         cell._text.set_color('darkslategray')
    elif (float(cell_text) <= -2) | (float(cell_text) >= 2) :
        cell._text.set_color('white')
    elif (cell_text=='nan'):
        cell._set_facecolor('darkgray')
        cell._text.set_color('darkgray')
plt.savefig("scorecards_s27_strat_10-150m_FR.png", dpi=300)
os.system('convert -trim scorecards_s27_strat_10-150m_FR.png scorecards_s27_strat_10-150m_FR.png')



#6. Merge all together
# English
#os.system('montage  -gravity west scorecards_s27_T.png scorecards_s27_S.png scorecards_s27_CIL.png scorecards_s27_MLD.png scorecards_s27_strat.png -tile 1x5 -geometry +1+1  -background white scorecards_s27.png') 
os.system('convert scorecards_s27_T.png scorecards_s27_S.png scorecards_s27_CIL.png scorecards_s27_MLD.png scorecards_s27_strat_0-50m.png scorecards_s27_strat_10-150m.png  -gravity East -append -geometry +1+1 scorecards_s27.png')

# French
#os.system('montage   -gravity east scorecards_s27_T_FR.png scorecards_s27_S_FR.png scorecards_s27_CIL_FR.png scorecards_s27_MLD_FR.png scorecards_s27_strat_FR.png -tile 1x5 -geometry +1+1  -background white scorecards_s27_FR.png') 
os.system('convert scorecards_s27_T_FR.png scorecards_s27_S_FR.png scorecards_s27_CIL_FR.png scorecards_s27_MLD_FR.png scorecards_s27_strat_0-50m_FR.png scorecards_s27_strat_10-150m_FR.png  -gravity East -append -geometry +1+1 scorecards_s27_FR.png')

# remove some images
os.system('rm scorecards_s27_T* scorecards_s27_S* scorecards_s27_CIL* scorecards_s27_MLD* scorecards_s27_strat*')

