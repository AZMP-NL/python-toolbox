'''

!!! This is the new script with new working tree.

For ResDoc and presentation on 2023 conditions.

This script can be run from anywhere, but a suggestion:
~/AZMP/state_reports/

The following subfolders will be created:
air_temperature/
bergs/
bottom_temp/
CIL/
climate_index/
LC_transport/

and figures further moved here:
./{year-in-review}/

Jan/Feb 2023
Frederic.Cyr@dfo-mpo.gc.ca &
Jonathan.Coyne@dfo-mpo.gc.ca


'''

import os
import matplotlib.pyplot as plt
import numpy as np
import datetime
import pandas as pd
#Provide the path to where custom packages are saved
import sys
sys.path.append(os.path.expanduser('~/github/AZMP-NL/python-toolbox/azmp_modules'))
import azmp_sections_tools as azst
import azmp_sections_climtools as azsct
import azmp_report_tools as azrt
import azmp_genreport as azgen
import azmp_utils as azu
import cc_tools as cc
import azmp_stn27_newtest as azS27

#Choose a year of interest
yoi = '2023'
#Choose a working directory name
work_name = 'AZMP_v1'

## Preamble (create folders to dump figures and data)
if os.path.isdir('operation_files') != True: os.system('mkdir operation_files')
if os.path.isdir('2023') != True: os.system('mkdir '+yoi)
if os.path.isdir('climate_indices') != True: os.system('mkdir climate_indices')
if os.path.isdir('air_temperature') != True: os.system('mkdir air_temperature')
if os.path.isdir('stn27') != True: os.system('mkdir stn27')
if os.path.isdir('bergs') != True: os.system('mkdir bergs')
if os.path.isdir('bottom_temp') != True: os.system('mkdir bottom_temp')
if os.path.isdir('bottom_saln') != True: os.system('mkdir bottom_saln')
if os.path.isdir('bottom_temp_stats') != True: os.system('mkdir bottom_temp_stats')
if os.path.isdir('AZMP_lines') != True: os.system('mkdir AZMP_lines')
if os.path.isdir('climate_indices') != True: os.system('mkdir climate_indices')


## ---- 2023 update ---- ##
# 1.  NAO, AO, AMO
azgen.nao(int(yoi),'~/data/'+work_name+'/climate_indices/nao_data.csv')
azgen.ao(int(yoi),'~/data/'+work_name+'/climate_indices/ao_data.csv')
azgen.amo(int(yoi),'~/data/'+work_name+'/climate_indices/amo_data.csv')
os.system('cp  NAO_winter_1950-'+yoi+'.png NAO_winter_1950-'+yoi+'_FR.png '+yoi)
os.system('mv *.pkl *.csv operation_files')
os.system('mv *.png climate_indices')

# 2. Air temperature (need to dowload AHCCD and download NUUK update) [Done 2023!]
%my_run azmp_dmi_nuukAirT.py
%my_run azmp_airTemp.py # use this one since 2020 conditions report
os.system('cp air_temp_2023.png air_temp_2023_FR.png air_temp_anom.png air_temp_anom_FR.png air_temp_climate_index.png air_temp_climate_index_FR.png ./2023/')
os.system('mv air_temp_2023.png air_temp_2023_FR.png air_temp_anom.png air_temp_anom_FR.png air_temp_climate_index.png air_temp_climate_index_FR.png air_temperature')
os.system('mv *.pkl operation_files')

%my_run azmp_air_scorecards.py
os.system('cp scorecards_air.png scorecards_air_FR.png ./2023/')
# delete tmp files
os.system('rm scorecards_*Air*.png scorecards_*nao*.png ')
os.system('mv *.png ./air_temperature')
os.system('mv *.csv ./operation_files')

# 3. SSTs  [Received from PSG]
# -> [Received from PSG]

# 4. Station 27 [JON to UPDATE]
#%my_run viking2022.py  # NOT IN 2022
#os.system('cp Viking2022.png Viking2022_FR.png ../2022')
'''
# need to delete /home/cyrf0006/AZMP/state_reports/stn27/stn27_all_casts.nc
os.system('mkdir stn27')
#Remember, azmp_stn27.py needs to be run twice (temperature and salinity)
os.system('python /home/jcoyne/Documents/AZMP-NL_python-toolbox/python-toolbox/azmp_stn27.py')
os.system('python /home/jcoyne/Documents/AZMP-NL_python-toolbox/python-toolbox/azmp_stn27_density.py')
os.system('python /home/jcoyne/Documents/AZMP-NL_python-toolbox/python-toolbox/azmp_stn27_analysis.py')
os.system('python /home/jcoyne/Documents/AZMP-NL_python-toolbox/python-toolbox/azmp_stn27_scorecards.py')
os.system('cp scorecards_s27.png scorecards_s27_FR.png s27_CIL_subplots.png s27_CIL_subplots_FR.png s27_TS_subplots.png s27_TS_subplotsFR.png s27_mld_monthly.png  s27_stratif_monthly_shallow.png s27_stratif_monthly_shallow_FR.png   2021/')
os.system('cp s27_salinity_subplot_2021.png s27_salinity_subplot_2021_FR.png s27_temperature_subplot_2021.png s27_temperature_subplot_2021_FR.png 2021/')
os.system("find ./ -maxdepth 1 -type f | xargs mv -t ./stn27") #move all other files to stn27/
'''

#New version using only functions
#Isolate the stn27 data
file_location = '~/data/AZMP/operation_files/stn27_all_casts.nc'
CASTS_path = '~/data/CASTS/*.nc'
s27_loc = [47.54667,-52.58667]
dc = .025
problem_casts = ['1988_10073164.nc','1990_02183115.nc','1991_11229001.nc']
azS27.stn27_dataisolate(file_location,CASTS_path,s27_loc,dc,problem_casts)
#Determine the temperature and salinity climatology
year_clim = [1991, 2020]
current_year = 2023
df,df_monthly,weekly_clim,df_weekly,df_year,anom = azS27.stn27_climatology(file_location,year_clim,current_year)
#Save the current year monthly average
for variable in df_monthly:
    csv_file = 'monthly_' + variable +  '_' + str(current_year) + '.csv'
    df_monthly[variable].T.to_csv(csv_file, float_format='%.4f')
#Create a temperature and salinity climatology figure
XLIM = [datetime.date(current_year, 1, 1), datetime.date(current_year, 12, 31)]
azS27.stn27_climatology_plot(weekly_clim,XLIM)
#Create a temperature and salinity current year plot
azS27.stn27_currentyear_plot(df_weekly,df_year,current_year,XLIM)
#Create a temperature and salinity anomaly year plot
azS27.stn27_anomaly_plot(anom,current_year,XLIM)
#Convert to subplots and remove individual plots
for variable in anom:
    os.system('montage '+\
        's27_'+variable+'_2023.png '+\
        's27_'+variable+'_clim.png '+\
        's27_'+variable+'_anom_2023.png '+\
        '-tile 1x3 -geometry +10+10  -background white  s27_'+variable+'_subplot_'+str(current_year)+'.png')
    os.system('montage '+\
        's27_'+variable+'_2023_FR.png '+\
        's27_'+variable+'_clim_FR.png '+\
        's27_'+variable+'_anom_2023_FR.png '+\
        '-tile 1x3 -geometry +10+10  -background white  s27_'+variable+'_subplot_'+str(current_year)+'_FR.png')
    os.system('rm '+\
        's27_'+variable+'_2023.png '+\
        's27_'+variable+'_clim.png '+\
        's27_'+variable+'_anom_2023.png')
    os.system('rm '+\
        's27_'+variable+'_2023_FR.png '+\
        's27_'+variable+'_clim_FR.png '+\
        's27_'+variable+'_anom_2023_FR.png')
#Create a station occupation figure
azS27.station_occupations(file_location,current_year,year_clim)

#Isolate the density data
df_rho,df_sig,df_SA,df_CT,ds,Z = azS27.density_calculator(df,file_location)
#Save the MLD
azS27.MLD_calculator(df_SA,df_CT,df_rho,Z)
#Save the stratification
azS27.stratification_calculator(file_location)

#Create anomaly output for temperature
years_flag = [1950,1980]
df = pd.read_pickle('S27_temperature_monthly.pkl')
df,anom_std,anom,annual_mean,ts_monthly_clim = azS27.anomaly_calculator(df,years_flag,current_year,year_clim)
#Plot the temperature anomaly
azS27.anomaly_plotter(anom_std,'temperature')
#Plot the temperature climatology
azS27.climatology_plotter(ts_monthly_clim,annual_mean,'temperature')
#Determine the CIL metrics
cil_temp,cil_core,cil_coredepth,cil_thickness = azS27.CIL_calculator(df)
#Create the CIL figures
azS27.CIL_plotter(cil_temp,'CIL mean temperature','Température moyenne de la CIF','s27_CILtemp_anomaly')
azS27.CIL_plotter(cil_core,'CIL core temperature','Température du coeur de la CIF','s27_CILcore_anomaly')
azS27.CIL_plotter(cil_coredepth,'CIL core depth','Profondeur du coeur de la CIF','s27_CILcoredepth_anomaly')
azS27.CIL_plotter(cil_thickness,'CIL thickness','Épaisseur de la CIF','s27_CILthickness_anomaly')

#Create anomaly output for salinity
df = pd.read_pickle('S27_salinity_monthly.pkl')
df,anom_std,anom,annual_mean,ts_monthly_clim = azS27.anomaly_calculator(df,years_flag,current_year,year_clim)
#Plot the salinity anomaly
azS27.anomaly_plotter(anom_std,'salinity')

#Get the stratification ready for plotting
strat_shallow_path = 'S27_stratif_0-50_monthly.pkl'
strat_deep_path = 'S27_stratif_10-150_monthly.pkl'
strat_monthly_shallow,strat_monthly_deep,anom,anom_std,strat_monthly_clim = azS27.stratification_plotter(
    strat_shallow_path,
    strat_deep_path,
    years_flag,current_year,year_clim)
#Create a stratification barplot
azS27.stratification_barplot(anom_std)
#Create a stratification time series
azS27.stratification_timeseries(anom)
#Create the mean stratification time series
azS27.stratification_timeseries_mean(anom,strat_monthly_clim)
#Create a current year stratification bar plot
azS27.stratification_currentyear_barplot(strat_monthly_shallow,strat_monthly_deep,current_year,year_clim)

#Get the MLD ready for plotting
MLD_path = 'S27_MLD_monthly.pkl'
mld,anom,anom_std = azS27.MLD_processor(MLD_path,years_flag,year_clim,current_year)
#Bar plot the MLD
azS27.MLD_barplot(anom_std)
#Time series plot the MLD 
azS27.MLD_timeseries(anom)
#Bar plot the current year MLD
azS27.MLD_currentyear_barplot(mld,current_year,year_clim=[1991,2020])

#Get the scorecard plot ready
year_plot = [1981,2023]
T_anom,T_anom_std,T_clim_period_annual = azS27.scorecard_TS_processor('S27_temperature_monthly.pkl',year_clim,year_plot,'Temp')
S_anom,S_anom_std,S_clim_period_annual = azS27.scorecard_TS_processor('S27_salinity_monthly.pkl',year_clim,year_plot,'Sal')
CIL_anom,CIL_anom_std,CIL_clim_period = azS27.scorecard_CIL_processor('S27_CIL_summer_stats.pkl',year_clim,year_plot)
mld_clim,mld_anom_std = azS27.scorecard_MLD_processor('S27_MLD_monthly.pkl',year_clim,year_plot)
strat_clim_0to50,strat_anom_std_0to50 = azS27.scorecard_strat_processor('S27_stratif_0-50_monthly.pkl',year_clim,year_plot)
strat_clim_10to150,strat_anom_std_10to150 = azS27.scorecard_strat_processor('S27_stratif_10-150_monthly.pkl',year_clim,year_plot)
#Plot the scorecard
azS27.scorecard_plotter(
    T_anom_std,
    T_clim_period_annual,
    '-- Vertically averaged temperature --',
    u'-- Moyenne verticale de température --',
    'T',years_present=True)
azS27.scorecard_plotter(
    S_anom_std,
    S_clim_period_annual,
    '-- Vertically averaged salinity --',
    u'-- Moyenne verticale de salinité --',
    'S',years_present=False)
azS27.scorecard_plotter(
    CIL_anom_std,
    CIL_clim_period,
    '-- Cold intermediate layer (CIL) properties --',
    u'-- Propriétés de la couche intermédiaire froide (CIF) --',
    'CIL',years_present=False)
azS27.scorecard_plotter(
    mld_anom_std,
    mld_clim,
    '-- Mixed layer depth (MLD) --',
    u'-- Profondeur de la couche de mélange (PCM) --',
    'MLD',years_present=False)
azS27.scorecard_plotter(
    strat_anom_std_0to50,
    strat_clim_0to50,
    '-- Stratification (0-50m) --',
    u'-- Stratification (0-50m) --',
    'strat_0-50m',years_present=False)
azS27.scorecard_plotter(
    strat_anom_std_10to150,
    strat_clim_10to150,
    '-- Stratification (10-150m) --',
    u'-- Stratification (10-150m) --',
    'strat_10-150m',years_present=False)
#Montage all the figures together
os.system('convert scorecards_s27_T.png scorecards_s27_S.png scorecards_s27_CIL.png scorecards_s27_MLD.png scorecards_s27_strat_0-50m.png scorecards_s27_strat_10-150m.png  -gravity East -append -geometry +1+1 scorecards_s27.png')
os.system('convert scorecards_s27_T_FR.png scorecards_s27_S_FR.png scorecards_s27_CIL_FR.png scorecards_s27_MLD_FR.png scorecards_s27_strat_0-50m_FR.png scorecards_s27_strat_10-150m_FR.png  -gravity East -append -geometry +1+1 scorecards_s27_FR.png')
os.system('cp scorecards_s27.png cp scorecards_s27_FR.png ./2023/')
os.system('rm scorecards_s27*.png')


# Clean the files:
os.system('cp s27_mld_monthly.png s27_stratif_monthly_deep.png  s27_stratif_monthly_shallow.png ./2023/')
os.system('cp s27_salinity_subplot_2023.png s27_temperature_subplot_2023.png ./2023/')
os.system('mv s27*.png ./stn27')
os.system('mv *.csv *.pkl ./operation_files')


# 5. Sea Ice (/home/cyrf0006/AZMP/state_reports/ice)
# -> [Received from PSG]


# 6. Icebergs (/home/cyrf0006/AZMP/state_reports/bergs) [FRED to UPDATE]
%my_run azmp_bergs.py
os.system('cp icebergs_climate*.png ../2023')
#os.system('cp bergs_annual_FR.png bergs_annual.png bergs_monthly_FR.png bergs_monthly.png ../2023')
os.system('mv *.png ./bergs')
os.system('mv *.pkl ./operation_files')


# 7. bottom temperature maps (FINISHED/WORKING - 2023)
#Calculate the bottom temperature and salinity climatology
for season in ['spring','summer','fall']:
    azu.get_bottomT_climato(
        INFILES='~/data/CABOTS/CABOTS_bottomstats_'+season+'.nc',
        lonLims=[-63, -45],
        latLims=[42, 58],
        bath_file='~/data/GEBCO/GEBCO_2023_sub_ice_topo.nc',
        year_lims=[1991, 2020],
        season=season,
        h5_outputfile='operation_files/Tbot_climato_'+season+'_0.10.h5'
        )
    azu.get_bottomS_climato(
        INFILES='~/data/CABOTS/CABOTS_bottomstats_'+season+'.nc',
        lonLims=[-63, -45],
        latLims=[42, 58],
        bath_file='~/data/GEBCO/GEBCO_2023_sub_ice_topo.nc',
        year_lims=[1991, 2020],
        season=season,
        h5_outputfile='operation_files/Sbot_climato_'+season+'_0.10.h5'
        )
    print('    -> '+season+' done!')

#Calculate the bottom temperature and salinity for specific years
for season in ['spring','fall']:
    azrt.bottom_temperature(
        season=season,
        year='2023',
        lonLims=[-63, -45],
        latLims=[42, 58],
        climato_file='operation_files/Tbot_climato_'+season+'_0.10.h5',
        netcdf_path='~/data/CABOTS/CABOTS_bottomstats_'+season+'.nc',
        CASTS_path='~/data/CASTS/'
        )
    azrt.bottom_salinity(
        season=season,
        year='2023',
        lonLims=[-63, -45],
        latLims=[42, 58],
        climato_file='operation_files/Sbot_climato_'+season+'_0.10.h5',
        netcdf_path='~/data/CABOTS/CABOTS_bottomstats_'+season+'.nc',
        CASTS_path='~/data/CASTS/'
        )
    print('    -> '+season+' done!')
os.system('cp bottomT_spring2023.png bottomT_spring2023_FR.png bottomT_fall2023.png bottomT_fall2023_FR.png 2023')
os.system('cp bottomS_spring2023.png bottomS_spring2023_FR.png bottomS_fall2023.png bottomS_fall2023_FR.png 2023')
os.system('mv bottomT_*.png bottom_temp_*.png bottom_temp/')
os.system('mv bottomS_*.png bottom_sal_*.png bottom_saln/')


## WE'RE HERE I THINK ##

# For NAFO STACFEN and STACFIS input: [NEED TO DO]
azrt.bottom_temperature(season='summer', year='2022', climato_file='Tbot_climato_SA4_summer_0.10.h5')


# bottom stats and scorecards (arange year+1)
#[need to flag years if coverage insufficient] (FINISHED/WORKING - 2023)
os.system('mkdir bottom_temp_stats')
for season in ['spring','fall']:
    azrt.bottom_stats(
        years=np.arange(1980, 2023+1),
        season=season,
        netcdf_path='~/data/CABOTS/CABOTS_bottomstats_'+season+'.nc'
        )
os.system('mv *.pkl operation_files/')
azrt.bottom_scorecards(years=[1980, 2023], clim_year=[1991, 2020])
os.system('cp scorecards_botT_spring.png scorecards_botT_spring_FR.png scorecards_botT_fall_FR.png scorecards_botT_fall.png 2023')
os.system('mv *.png *.csv bottom_temp_stats/')



# For NAFO STACFEN and STACFIS input (for azmp_composite_index.py): [NEED TO DO]
azrt.bottom_stats(years=np.arange(1980, 2022), season='summer', climato_file='Tbot_climato_SA4_summer_0.10.h5')

# bottom temperature bar plots [need to flag years if coverage insufficient] (FINISHED/WORKING - 2023)
os.system('python /home/jcoyne/Documents/AZMP-NL_python-toolbox/python-toolbox/azmp_bottomT_mean_anomaly.py')
os.system('cp bottomT_anomalies_climateindex.png bottomT_anomalies_climateindex_FR.png 2023')
os.system('mv *.pkl bottom_temp_stats/')
os.system('mv *.png bottom_temp_stats/')


## ---------------  Sections plots ------------- ## (FINISHED/WORKING - 2023)
year = 2023
sections = ['SI', 'BB', 'FC']
seasons = ['summer']
#Ensure that climatologies have been set up
for section in sections:
    for season in seasons:
        azsct.section_clim(
            SECTION=section,
            SEASON=season,
            CLIM_YEAR=[1950,2023],
            dlat=2,
            dlon=2,
            z1=2,
            dz=5,
            dc=0.2,
            CASTS_path='~/data/CASTS/',
            bath_path='~/data/GEBCO/GEBCO_2023_sub_ice_topo.nc'
            )
os.system('mv *.pkl operation_files/')
variables = ['temperature', 'salinity']
for section in sections:
    for season in seasons:
        for var in variables:
            #Ensure that STANDARD_SECTIONS.xlsx has been put in operation_files/
            #Esnure that bathymetry .txt files are in operation_files/
            azst.seasonal_section_plot(
                VAR=var,
                SECTION=section,
                SEASON=season,
                bath_path='~/data/GEBCO/GEBCO_2023_sub_ice_topo.nc',
                CASTS_path='~/data/CASTS/',
                YEAR=year,
                ZMAX=500,
                STATION_BASED=True
                )
            plt.close('all')

        command = 'montage temperature_' + section + '_' + season + '_' + str(year) + '.png salinity_' + section + '_' + season + '_' + str(year) + '.png  -tile 2x1 -geometry +10+10  -background white ' + section + '_stn_' + season + '_' + str(year) + '.png'
        os.system(command)
        command = 'montage temperature_' + section + '_' + season + '_' + str(year) + '_FR.png salinity_' + section + '_' + season + '_' + str(year) + '_FR.png  -tile 2x1 -geometry +10+10  -background white ' + section + '_stn_' + season + '_' + str(year) + '_FR.png'
        os.system(command)

        os.system('cp ' + section + '_stn_' + season + '_' + str(year) + '.png 2023')
        os.system('cp ' + section + '_stn_' + season + '_' + str(year) + '_FR.png 2023')
os.system('mv *.png *.csv *.pkl AZMP_lines/')


# Section CIL [CANNOT UPDATE 2022]
%my_run azmp_CIL_scorecards.py  # update year in script!
os.system('rm scorecards_CIL_SI* scorecards_CIL_BB* scorecards_CIL_FC*')
os.system('cp scorecards_CIL.png scorecards_CIL_FR.png ../2022')

#(FINISHED/WORKING - 2023)
os.system('python /home/jcoyne/Documents/AZMP-NL_python-toolbox/python-toolbox/azmp_CIL_mean_anomaly.py')
os.system('cp CIL_volume_climateindex.png CIL_volume_climateindex_FR.png ../2023')
os.system('mv *.png AZMP_lines/')
os.system('mv *.pkl operation_files/')

## ----------- NLCI ---------------- ##
%my_run azmp_climate_index.py



## ----------- AZMP SAR / IROC ---------------- ##
#%my_run azmp_CIL_stats.py 
%my_run azmp_CIL_stats_update.py # <---- preferred if only an update is needed (need to edit sections)
%my_run azmp_sar_input.py
os.system('cp NL_climate_index_ms_scorecards_FR.png NL_climate_index_ms_scorecards.png NL_climate_index_ms_FR.png NL_climate_index_ms.png ../2022')

## ----------- CSAS DATA ---------------- ##
%my_run csas_crab_stats.py
%my_run azmp_bottomT_habitat.py    
#%my_run NSRF_bottomT.py 
#%my_run NSRF_bottomS.py
#%my_run NSRF_bottomT_habitat.py 
#%my_run azmp_bottomT_shrimp_habitat.py 
#azrt.bottom_stats(years=np.arange(1980, 2020), season='summer') # for crab 4R
azrt.bottom_temperature(season='summer', year='2020') # for crab 4R

# SFAs (clim fill is new in 2022) [done 2022]
#azrt.sfa_bottom_stats(years=np.arange(2006, 2023), season='summer', plot=True, climato_file='Tbot_climato_NSRFx_summer_2006-2021.h5')
azrt.sfa_bottom_stats(years=np.arange(2006, 2023), season='summer', climato_file='Tbot_climato_NSRFx_summer_2006-2021.h5', plot=True, clim_fill=True, plot_biomass=True)
os.system('cp sfa_bottomT_2022.png ../../2022')
azrt.sfa_bottom_scorecards(years=np.arange(2006, 2023), clim_year=[2006, 2020])
os.system('cp scorecards_botT_SFA2-4_summer.png scorecards_botT_SFA2-4_summer_FR.png ../../2022')
# salinity (new in 2022)
azrt.sfa_bottom_stats(years=np.arange(2006, 2023), season='summer', climato_file='Sbot_climato_NSRFx_summer_2006-2021.h5', plot=True, clim_fill=True, plot_biomass=True ,var='salinity')
os.system('cp sfa_bottomS_2022.png ../../2022')
azrt.sfa_bottomS_scorecards(years=np.arange(2006, 2023), clim_year=[2006, 2021])
os.system('cp scorecards_botS_SFA2-4_summer.png scorecards_botS_SFA2-4_summer_FR.png ../../2022')
# In ~/research/lobster:
%my_run lfa_sst_extract.py
%my_run lfa_sst_analysis.py


## ----------- OA data for BGC ResDoc ---------------- ##
os.system('cd /home/cyrf0006/AZMP/oa/NL_obs')
years = [2019, 2020] 
seasons = ['summer', 'fall']
variables = ['Omega_Aragonite_(unitless)', 'pH_Total_(total_scale)', 'Oxygen_Saturation_(%)', 'Dissolved_Oxygen_(mL/L)']
depth = 'bottom'

for year in years:
    print('-> ' + str(year))
    for season in seasons:
        print(' --> ' + season)
        for variable in variables:
            print('  ---> ' + variable)
            cc.seasonal_map_NL(variable, year, season, depth)
            plt.close('all')

# for 2020:
#montage AZMP_OA_2020_summer_OmegaA_surface.png AZMP_OA_2020_summer_pH_surface.png AZMP_OA_2020_summer_DO_perc_surface.png -tile 1x3 -geometry +10+10  -background white AZMP_OA_2020_summer_surface.png

# for 2021:
#montage AZMP_OA_2021_summer_OmegaA_surface.png AZMP_OA_2021_summer_pH_surface.png AZMP_OA_2021_summer_DO_perc_surface.png -tile 1x3 -geometry +10+10  -background white AZMP_OA_2021_summer_surface.png

# Montage 2 years in review:
os.system('montage NL_OA_2019_summer_pH_bottom.png NL_OA_2019_summer_OmegaA_bottom.png NL_OA_2020_summer_pH_bottom.png NL_OA_2020_summer_OmegaA_bottom.png -tile 2x2 -geometry +10+10  -background white NL_OA_summer_2019-2020.png')
os.system('montage NL_OA_2019_summer_pH_bottom.png NL_OA_2019_fall_pH_bottom.png NL_OA_2020_summer_pH_bottom.png NL_OA_2020_fall_pH_bottom.png -tile 2x2 -geometry +10+10  -background white NL_pH_summer-fall_2019-2020.png')
os.system('montage NL_OA_2019_summer_DO_perc_bottom.png NL_OA_2019_fall_DO_perc_bottom.png NL_OA_2020_summer_DO_perc_bottom.png NL_OA_2020_fall_DO_perc_bottom.png -tile 2x2 -geometry +10+10  -background white NL_DO_perc_2019-2020.png')
# Montage for same year summer and fall:
os.system('montage NL_OA_2019_summer_pH_bottom.png NL_OA_2019_fall_pH_bottom.png NL_OA_2019_summer_OmegaA_bottom.png NL_OA_2019_fall_OmegaA_bottom.png -tile 2x2 -geometry +10+10  -background white NL_OA_2019.png')
os.system('montage NL_OA_2020_summer_pH_bottom.png NL_OA_2020_fall_pH_bottom.png NL_OA_2020_summer_OmegaA_bottom.png NL_OA_2020_fall_OmegaA_bottom.png -tile 2x2 -geometry +10+10  -background white NL_OA_2020.png')

# Montage in French:
os.system('montage NL_OA_2019_summer_pH_bottom_FR.png NL_OA_2019_summer_OmegaA_bottom_FR.png NL_OA_2020_summer_pH_bottom_FR.png NL_OA_2020_summer_OmegaA_bottom_FR.png -tile 2x2 -geometry +10+10  -background white NL_OA_summer_2019-2020_FR.png')
os.system('montage NL_OA_2019_summer_pH_bottom_FR.png NL_OA_2019_fall_pH_bottom_FR.png NL_OA_2020_summer_pH_bottom_FR.png NL_OA_2020_fall_pH_bottom_FR.png -tile 2x2 -geometry +10+10  -background white NL_pH_summer-fall_2019-2020_FR.png')
os.system('montage NL_OA_2019_summer_DO_perc_bottom_FR.png NL_OA_2019_fall_DO_perc_bottom_FR.png NL_OA_2020_summer_DO_perc_bottom_FR.png NL_OA_2020_fall_DO_perc_bottom_FR.png -tile 2x2 -geometry +10+10  -background white NL_DO_perc_2019-2020_FR.png')
os.system('montage NL_OA_2019_summer_pH_bottom_FR.png NL_OA_2019_fall_pH_bottom_FR.png NL_OA_2019_summer_OmegaA_bottom_FR.png NL_OA_2019_fall_OmegaA_bottom_FR.png -tile 2x2 -geometry +10+10  -background white NL_OA_2019_FR.png')
os.system('montage NL_OA_2020_summer_pH_bottom_FR.png NL_OA_2020_fall_pH_bottom_FR.png NL_OA_2020_summer_OmegaA_bottom_FR.png NL_OA_2020_fall_OmegaA_bottom_FR.png -tile 2x2 -geometry +10+10  -background white NL_OA_2020_FR.png')


