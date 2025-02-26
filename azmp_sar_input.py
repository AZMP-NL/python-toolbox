'''
To generate AZMP score cards for bottom temperature

Uses pickled object generated by azmp_bottom_stats.py

Check /home/cyrf0006/AZMP/SAR_files

Files generated:
BT_2J_fall.dat
BT_3K_fall.dat
BT_3LNO_fall.dat
BT_3LNO_spring.dat
BT_3Ps_spring.dat
CIL_Bonavista_0C_Area.dat
CIL_FlemishCap_0C_Area.dat
CIL_SealIsland_0C_Area.dat
CIL_WhiteBay_0C_Area.dat
S27_Integrated.dat
'''

import numpy as  np
import matplotlib.pyplot as plt
import pandas as pd
import os
#import unicodedata

# Parameters 
path = 'operation_files/'
clim_year = [1991, 2020]
year_min = 1948
stn27_months = [5, 11]

#Years to flag
flag_SI_itp = np.array([1932,1936,1937,1939,1940,1941,1950,1966,1967,1968,1972,1981,1989,2022])
flag_SI_stn = np.arange(1948,1994+1)
flag_SI_stn = np.append(flag_SI_stn, [1998,1999,2000,2001,2003,2004,2006,2007,2022])
flag_BB_itp = np.array([1966,1967,1968,1993,2022])
flag_BB_stn = np.arange(1948,1994+1)
flag_BB_stn = np.append(flag_BB_stn, [1995,1997,1998,1999,2000,2001,2002,2022])
flag_FC_itp = np.array([2022])
flag_FC_stn = np.arange(1948,1994+1)
flag_FC_stn = np.append(flag_FC_stn, [1995,1997,1998,1999,2000,2001,2002,2003,2004,2006,2007,2008,2022])
flag_WB_itp = np.arange(1948,1959+1)
flag_WB_itp = np.append(flag_WB_itp, np.arange(1961,1972+1))
flag_WB_itp = np.append(flag_WB_itp, [2011,2019,2022])
flag_WB_stn = flag_WB_itp

#### -------------1. bottom temperature ---------------- ####
## 2J fall
infile = '~/data/CABOTS/csv_averages/fall_2J_regional_averages.csv'
df = pd.read_csv(infile,index_col=0)
#df.index = pd.to_datetime(df.index) # update index to datetime
# Flag bad years (no or weak sampling):
bad_years = np.array([1995])
for i in bad_years:
    df[df.index==i]=np.nan
# compute std anom
df_clim = df[(df.index>=clim_year[0]) & (df.index<=clim_year[1])]
std_anom = (df-df_clim.mean(axis=0))/df_clim.std(axis=0)
# std anom for temperature
df['std_anom'] = std_anom['Tmean']
# keep only 2 columns
df = df[['Tmean', 'std_anom']]
#df.index = df.index.year
df.to_csv('BT_2J_fall.dat', header=False, sep = ' ', float_format='%.2f')

## 3K fall
infile = '~/data/CABOTS/csv_averages/fall_3K_regional_averages.csv'
df = pd.read_csv(infile,index_col=0)
#df.index = pd.to_datetime(df.index) # update index to datetime
# Flag bad years (no or weak sampling):
bad_years = np.array([])
for i in bad_years:
    df[df.index==i]=np.nan
# compute std anom
df_clim = df[(df.index>=clim_year[0]) & (df.index<=clim_year[1])]
std_anom = (df-df_clim.mean(axis=0))/df_clim.std(axis=0)
# std anom for temperature
df['std_anom'] = std_anom['Tmean']
# keep only 2 columns
df = df[['Tmean', 'std_anom']]
#df.index = df.index.year
df.to_csv('BT_3K_fall.dat', header=False, sep = ' ', float_format='%.2f')

## 3LNO fall
infile = '~/data/CABOTS/csv_averages/fall_3LNO_grandbanks_regional_averages.csv'
df = pd.read_csv(infile,index_col=0)
#df.index = pd.to_datetime(df.index) # update index to datetime
# Flag bad years (no or weak sampling):
bad_years = np.array([2021])
for i in bad_years:
    df[df.index==i]=np.nan
# compute std anom
df_clim = df[(df.index>=clim_year[0]) & (df.index<=clim_year[1])]
std_anom = (df-df_clim.mean(axis=0))/df_clim.std(axis=0)
# std anom for temperature
df['std_anom'] = std_anom['Tmean']
# keep only 2 columns
df = df[['Tmean', 'std_anom']]
#df.index = df.index.year
df.to_csv('BT_3LNO_fall.dat', header=False, sep = ' ', float_format='%.2f')

## 3LNO spring
infile = '~/data/CABOTS/csv_averages/spring_3LNO_grandbanks_regional_averages.csv'
df = pd.read_csv(infile,index_col=0)
#df.index = pd.to_datetime(df.index) # update index to datetime
# Flag bad years (no or weak sampling):
bad_years = np.array([2020, 2021])
for i in bad_years:
    df[df.index==i]=np.nan
# compute std anom
df_clim = df[(df.index>=clim_year[0]) & (df.index<=clim_year[1])]
std_anom = (df-df_clim.mean(axis=0))/df_clim.std(axis=0)
# std anom for temperature
df['std_anom'] = std_anom['Tmean']
# keep only 2 columns
df = df[['Tmean', 'std_anom']]
#df.index = df.index.year
df.to_csv('BT_3LNO_spring.dat', header=False, sep = ' ', float_format='%.2f')

## 3Ps spring 
infile = '~/data/CABOTS/csv_averages/spring_3Ps_regional_averages.csv'
df = pd.read_csv(infile,index_col=0)
#df.index = pd.to_datetime(df.index) # update index to datetime
# Flag bad years (no or weak sampling):
bad_years = np.array([1980, 1981, 1985, 1986, 1987, 1988, 1989, 1990, 1991, 1992, 2006, 2020])
for i in bad_years:
    df[df.index==i]=np.nan
# compute std anom
df_clim = df[(df.index>=clim_year[0]) & (df.index<=clim_year[1])]
std_anom = (df-df_clim.mean(axis=0))/df_clim.std(axis=0)
# std anom for temperature
df['std_anom'] = std_anom['Tmean']
# keep only 2 columns
df = df[['Tmean', 'std_anom']]
#df.index = df.index.year
df.to_csv('BT_3Ps_spring.dat', header=False, sep = ' ', float_format='%.2f')

#### ------------- 2. winter NAO ---------------- ####
## nao_file = '/home/cyrf0006/data/AZMP/indices/data.csv'
## df = pd.read_csv(nao_file, header=1)
## # Set index
## df = df.set_index('Date')
## df.index = pd.to_datetime(df.index, format='%Y%m')
## # Select only DJF
## df_winter = df[(df.index.month==12) | (df.index.month==1) | (df.index.month==2) |  (df.index.month==3)]
## # Start Dec-1950
## df_winter = df_winter[df_winter.index>pd.to_datetime('1950-10-01')]
## # Average 3 consecutive values (DJF average); We loose index.
## df_winter = df_winter.groupby(np.arange(len(df_winter))//4).mean()
## # Reset index using years only
## year_unique = pd.unique(df.index.year)[1:,]
## df_winter = df_winter.iloc[np.arange(0, year_unique.size)] # reduce if last month is december (belongs to following year)
## df_winter.index = year_unique
## df_winter.to_csv('NAO_DJFM.dat', header=False, sep = ' ', float_format='%.2f')

#df_winter = pd.read_pickle('operation_files/NAO_winter.pkl')


#### ------------- 3. CIL ---------------- ####
# see /home/cyrf0006/AZMP/state_reports/ColbourneStuff/CIL_AZMP_SPRING_SUMMER_FALL.xlsx
# These timeseries are already calculated for the IROC (except for WB). 
# Check iroc_CIL_area.py and files in /home/cyrf0006/research/WGOH/IROC <-- Should be removed

# UPDATE 2022 check: 
# azmp_CIL_stats.py
# azmp_CIL_stats_update.pu

df_SI = pd.read_csv('operation_files/CIL_area_SI.csv')
df_WB = pd.read_csv('operation_files/CIL_area_WB.csv')
df_BB = pd.read_csv('operation_files/CIL_area_BB.csv')
df_FC = pd.read_csv('operation_files/CIL_area_FC.csv')
df_SI.rename(columns={df_SI.columns[0]:'year'}, inplace=True)
df_WB.rename(columns={df_WB.columns[0]:'year'}, inplace=True)
df_BB.rename(columns={df_BB.columns[0]:'year'}, inplace=True)
df_FC.rename(columns={df_FC.columns[0]:'year'}, inplace=True)
df_SI.set_index('year', inplace=True)
df_WB.set_index('year', inplace=True)
df_BB.set_index('year', inplace=True)
df_FC.set_index('year', inplace=True) 

# remove problem years
df_SI['station-ID'].loc[flag_SI_stn] = np.nan
df_SI['interp_field'].loc[flag_SI_itp] = np.nan
df_WB['station-ID'].loc[flag_WB_stn] = np.nan
df_WB['interp_field'].loc[flag_WB_itp] = np.nan
df_BB['station-ID'].loc[flag_BB_stn] = np.nan
df_BB['interp_field'].loc[flag_BB_itp] = np.nan
df_FC['station-ID'].loc[flag_FC_stn] = np.nan
df_FC['interp_field'].loc[flag_FC_itp] = np.nan

# cut timeseries
df_SI = df_SI[df_SI.index>=year_min]
df_WB = df_WB[df_WB.index>=year_min]
df_BB = df_BB[df_BB.index>=year_min]
df_FC = df_FC[df_FC.index>=year_min]

# Save timeseries
df_SI['interp_field'].to_csv('CIL_SealIsland_0C_Area.dat', header=False, sep = ' ', float_format='%.2f')
df_WB['interp_field'].to_csv('CIL_WhiteBay_0C_Area.dat', header=False, sep = ' ', float_format='%.2f')
df_BB['interp_field'].to_csv('CIL_Bonavista_0C_Area.dat', header=False, sep = ' ', float_format='%.2f')
df_FC['interp_field'].to_csv('CIL_FlemishCap_0C_Area.dat', header=False, sep = ' ', float_format='%.2f')

df_SI['station-ID'].to_csv('CIL_SealIsland_0C_Area_stationBased.dat', header=False, sep = ' ', float_format='%.2f')
df_WB['station-ID'].to_csv('CIL_WhiteBay_0C_Area_stationBased.dat', header=False, sep = ' ', float_format='%.2f')
df_BB['station-ID'].to_csv('CIL_Bonavista_0C_Area_stationBased.dat', header=False, sep = ' ', float_format='%.2f')
df_FC['station-ID'].to_csv('CIL_FlemishCap_0C_Area_stationBased.dat', header=False, sep = ' ', float_format='%.2f')



#### ------------- 4. Stn 27 ---------------- ####
# see /home/cyrf0006/AZMP/S27/station_27_stratification.xlsx
# Load pickled data
df_temp = pd.read_pickle('operation_files/S27_temperature_monthly.pkl')
df_sal = pd.read_pickle('operation_files/S27_salinity_monthly.pkl')
df_strat_shallow = pd.read_pickle('operation_files/S27_stratif_0-50_monthly.pkl')
df_strat_deep = pd.read_pickle('operation_files/S27_stratif_10-150_monthly.pkl')


# Reduce to summer months and annual mean
df_temp = df_temp[(df_temp.index.month>=stn27_months[0]) & (df_temp.index.month<=stn27_months[1])]
df_sal = df_sal[(df_sal.index.month>=stn27_months[0]) & (df_sal.index.month<=stn27_months[1])]
df_strat_shallow = df_strat_shallow[(df_strat_shallow.index.month>=stn27_months[0]) & (df_strat_shallow.index.month<=stn27_months[1])]
df_strat_deep = df_strat_deep[(df_strat_deep.index.month>=stn27_months[0]) & (df_strat_deep.index.month<=stn27_months[1])]


#Cycle through a shallow and deep output
depth_range = {'shallow': [0,50], 'deep': [10,150]}
for depth in depth_range:

    #Beginning temp and saln
    var_letter = ['T','S']
    d_name = str(depth_range[depth][0])+'_'+str(depth_range[depth][1])

    #Cycle through temperature and salinity
    for x,df in enumerate([df_temp, df_sal]):

        #Flag the bad years
        for y in [1950,1980]:
            df[df.index.year == y] = np.nan
        #Determine the vertically averaged temperature
        ts_stack = df.groupby([(df.index.year),(df.index.month)]).mean()
        ts_stack.index = ts_stack.index.set_names(['year', 'month'])
        #Isolate, calculate the period of climatology
        df_clim_period = df[(df.index.year>=clim_year[0]) & (df.index.year<=clim_year[1])]
        monthly_clim_mean = df_clim_period.groupby(df_clim_period.index.month).mean()
        monthly_clim_stdv = df_clim_period.groupby(df_clim_period.index.month).std()
        #ts_monthly_clim = monthly_clim_mean.mean(axis=1)
        #ts_monthly_std = monthly_clim_stdv.mean(axis=1)
        #Tile the climatology for however many years are present
        years = len(df.index.year.unique())
        month_start = ts_stack.index.values[0][1] - stn27_months[0]
        monthly_clim_mean = pd.concat([monthly_clim_mean] * years).iloc[month_start:]
        monthly_clim_stdv = pd.concat([monthly_clim_stdv] * years).iloc[month_start:]
        monthly_clim_mean.set_index(ts_stack.index, inplace=True)
        monthly_clim_stdv.set_index(ts_stack.index, inplace=True)
        #Calculate the anomalies
        anom_mi = ts_stack-monthly_clim_mean
        std_anom_mi = (ts_stack-monthly_clim_mean) / monthly_clim_stdv

        #Isolate for the depths of interest
        for i in ['0_btm','depth','170_btm']:
            if i == '0_btm':
                name = var_letter[x]+'_'+i
                ts_monthly_clim = monthly_clim_mean[monthly_clim_mean.columns[(monthly_clim_mean.columns>=depth_range[depth][0])*(monthly_clim_mean.columns<=depth_range[depth][1])]].mean(axis=1)
                anom_monthly = anom_mi.mean(axis=1)
                std_anom_monthly = std_anom_mi.mean(axis=1)
            if i == 'depth':
                name = var_letter[x]+'_'+d_name
                ts_monthly_clim = monthly_clim_mean[monthly_clim_mean.columns[(monthly_clim_mean.columns>=depth_range[depth][0])*(monthly_clim_mean.columns<=depth_range[depth][1])]].mean(axis=1)
                anom_monthly = anom_mi[anom_mi.columns[(anom_mi.columns>=depth_range[depth][0])*(anom_mi.columns<=depth_range[depth][1])]].mean(axis=1)
                std_anom_monthly = std_anom_mi[std_anom_mi.columns[(std_anom_mi.columns>=depth_range[depth][0])*(std_anom_mi.columns<=depth_range[depth][1])]].mean(axis=1)
            if i == '170_btm':
                name = var_letter[x]+'_'+i
                ts_monthly_clim = monthly_clim_mean[monthly_clim_mean.columns[(monthly_clim_mean.columns>=170)]].mean(axis=1)
                anom_monthly = anom_mi[anom_mi.columns[(anom_mi.columns>=170)]].mean(axis=1)
                std_anom_monthly = std_anom_mi[std_anom_mi.columns[(std_anom_mi.columns>=170)]].mean(axis=1)
            monthly_stdanom = std_anom_monthly.unstack()
            monthly_anom = anom_monthly.unstack()
            #Calculate the annual anomalies
            anom_std = monthly_stdanom.mean(axis=1)
            anom_std.index = pd.to_datetime(anom_std.index, format='%Y')
            anom = monthly_anom.mean(axis=1) 
            anom.index = pd.to_datetime(anom.index, format='%Y')
            #Annual mean is given by annual anomaly + monthly clim
            annual_mean = anom + ts_monthly_clim.mean()
            exec(name + ' = annual_mean')

    # Statification
    if depth == 'shallow':
        strat_depth = df_strat_shallow*42 # the SAR presents density difference, not stratification
    elif depth == 'deep':
        strat_depth = df_strat_deep*140 # the SAR presents density difference, not stratification

    my_ts = strat_depth
    ts_stack = my_ts.groupby([(my_ts.index.year),(my_ts.index.month)]).mean()
    ts_unstack = ts_stack.unstack()
    # Monthly clim (ts_monthly_clim)
    ts_clim_period = my_ts[(my_ts.index.year>=clim_year[0]) & (my_ts.index.year<=clim_year[1])]
    ts_monthly_stack = ts_clim_period.groupby([(ts_clim_period.index.year),(ts_clim_period.index.month)]).mean()
    ts_monthly_clim = ts_monthly_stack.groupby(level=1).mean()
    ts_monthly_std = ts_monthly_stack.groupby(level=1).std()
    # monthly anom and normalized anom
    monthly_anom = ts_unstack - ts_monthly_clim 
    monthly_stdanom = (ts_unstack - ts_monthly_clim) /  ts_monthly_std

    #See if each year has enough data present
    nom = np.sum(~np.isnan(monthly_anom.values),axis=1)
    monthly_anom.iloc[nom<3] = np.nan
    monthly_stdanom.iloc[nom<3] = np.nan

    # annual normalized anomaly
    anom_std = monthly_stdanom.mean(axis=1)
    anom_std.index = pd.to_datetime(anom_std.index, format='%Y')
    # annual anomaly
    anom = monthly_anom.mean(axis=1) 
    anom.index = pd.to_datetime(anom.index, format='%Y')
    # Re-create annual mean by adding annual anomaly to monthly clim
    annual_mean = anom + ts_monthly_clim.mean()
    exec('strat_'+d_name+' = annual_mean')

    # Merge 
    if depth == 'shallow':
        df_stn27 = pd.concat([T_0_50, T_170_btm, S_0_50, strat_0_50], axis=1, keys=['Temp 0-50m', 'Temp 170-176m', 'Sal 0-50m', 'Strat 5-50m'])
        df_stn27.index = df_stn27.index.year
        df_stn27 = df_stn27[df_stn27.index>=year_min]
        df_stn27.to_csv('S27_shallow_Integrated.dat', header=True, sep = ' ', na_rep='-99', float_format='%.3f')
    elif depth == 'deep':
        df_stn27 = pd.concat([T_10_150, T_170_btm, S_10_150, strat_10_150], axis=1, keys=['Temp 10-150m', 'Temp 170-176m', 'Sal 10-150m', 'Strat 10-150m'])
        df_stn27.index = df_stn27.index.year
        df_stn27 = df_stn27[df_stn27.index>=year_min]
        df_stn27.to_csv('S27_deep_Integrated.dat', header=True, sep = ' ', na_rep='-99', float_format='%.3f')

#Zip the files
os.system('zip SAR_azmp-nl_2024.zip *.dat')
os.system('mv *.dat *.zip SAR_files')
