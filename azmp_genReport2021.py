'''
For ResDoc and presentation on 2021 conditions.
Started Oct. 2021
Frederic.Cyr@dfo-mpo.gc.ca

'''

import os
import matplotlib.pyplot as plt
import numpy as np
import azmp_sections_tools as azst
import azmp_report_tools as azrt  

## ---- 2021 update ---- ##
# 1.  NAO (should be a function with year as input) [Done 2021]
# in /home/cyrf0006/AZMP/state_reports/airTemp
%my_run azmp_nao.py # had to re-write 2021
%my_run azmp_ao.py  # had to re-write 2021
%my_run azmp_amo.py
os.system('cp  NAO_winter_1950-2021.png NAO_winter_1950-2021_FR.png ../2021/')


# 2. Air temperature (need to dowload AHCCD and download NUUK update)
%my_run azmp_dmi_nuukAirT.py
%my_run azmp_airTemp.py # use this one since 2020 conditions report
os.system('cp air_temp_2021.png air_temp_2021_FR.png air_temp_anom.png air_temp_anom_FR.png ../2021/')

%my_run azmp_air_scorecards.py
os.system('cp scorecards_air.png scorecards_air_FR.png ../2021/')

# 3. SSTs  [Received from PSG]
# wget -m ftp://ftp.dfo-mpo.gc.ca/bometrics/noaa/stats/boxes/*.stat 
# (in /home/cyrf0006/data/BIO_remote/bometrics/noaa/stats/boxes)
#%my_run azmp_SSTs.py # to update bometrics data
#%my_run azmp_SSTs_fromExcel.py # To merge with Eugene's historical Excel data
#os.system('cp SST_index.png ../2021/')
#%my_run azmp_SSTs_scorecards.py # to generate monthly anom, scorecards, etc.
#os.system('cp scorecards_sst_yearly.png scorecards_sst_monthly.png ../2021/')


# 4. Station 27 [Done 2021]
%my_run viking2021.py  # done, but not enough data
Os.system('cp Viking2021.png Viking2021_FR.png ../2021')

# need to delete /home/cyrf0006/AZMP/state_reports/stn27/stn27_all_casts.nc
%my_run azmp_stn27.py # ** Note that the script below should be run first to update stn27_all_casts.nc
%my_run azmp_stn27_density.py # <-- ** Run this script first with "append" option
%my_run azmp_stn27_analysis.py
%my_run azmp_stn27_scorecards.py
os.system('cp scorecards_s27.png scorecards_s27_FR.png s27_CIL_subplots.png s27_CIL_subplots_FR.png s27_TS_subplots.png s27_TS_subplotsFR.png s27_mld_monthly.png  s27_stratif_monthly.png s27_stratif_monthly_FR.png   ../2021/')
os.system('cp s27_salinity_subplot_2021.png s27_salinity_subplot_2021_FR.png s27_temperature_subplot_2021.png s27_temperature_subplot_2021_FR.png ../2021/')

# To prepare MS
#%my_run azmp_stn27_climateindex_ms.py 

# 5. Sea Ice (/home/cyrf0006/AZMP/state_reports/ice) [Received from PSG]

#%my_run azmp_ice_index.py
#os.system('cp ice_index.png ice_index_FR.png ../2021/')


# 6. Icebergs (/home/cyrf0006/AZMP/state_reports/bergs) [DONE 2021]
%my_run azmp_bergs.py
os.system('cp bergs_annual_FR.png bergs_annual.png bergs_monthly_FR.png bergs_monthly.png ../2021')


# 7. bottom temperature maps [Done 2021!]
azrt.bottom_temperature(season='spring', year='2021', climato_file='Tbot_climato91_spring_0.10.h5') 
azrt.bottom_temperature(season='fall', year='2021', climato_file='Tbot_climato91_fall_0.10.h5')
os.system('cp bottomT_spring2021.png bottomT_spring2021_FR.png bottomT_fall2021.png bottomT_fall2021_FR.png ../2021')
# For NAFO STACFEN and STACFIS input:
azrt.bottom_temperature(season='summer', year='2021', climato_file='Tbot_climato_SA4_summer_0.10.h5')

# bottom salinity maps
azrt.bottom_salinity(season='spring', year='2021', climato_file='Sbot_climato91_spring_0.10.h5') 
azrt.bottom_salinity(season='fall', year='2021', climato_file='Sbot_climato91_fall_0.10.h5')
os.system('cp bottomS_spring2021.png bottomS_spring2021_FR.png bottomS_fall2021.png bottomS_fall2021_FR.png ../2021')

# bottom stats and scorecards (arange year+1) [Need to find a way to speed this by not re-processing everything!!]
azrt.bottom_stats(years=np.arange(1980, 2022), season='spring') #-> CHECK!
azrt.bottom_stats(years=np.arange(1980, 2022), season='fall')
azrt.bottom_stats(years=np.arange(1980, 2022), season='summer')
azrt.bottom_scorecards(years=[1980, 2022], clim_year=[1991, 2020])
os.system('cp scorecards_botT_spring.png scorecards_botT_spring_FR.png scorecards_botT_fall_FR.png scorecards_botT_fall.png ../2021')
# For NAFO STACFEN and STACFIS input (for azmp_composite_index.py):
#azrt.bottom_stats(years=np.arange(1980, 2022), season='summer', climato_file='Tbot_climato_SA4_summer_0.10.h5')

# bottom temperature bar plots
%my_run azmp_bottomT_mean_anomaly.py # same as previous
os.system('cp mean_anomalies_fall.png mean_anomalies_spring.png ../2021')
os.system('cp mean_anomalies_fall_FR.png mean_anomalies_spring_FR.png ../2021')


## ---------------  Sections plots ------------- ## [Done 2021!]
year = 2021 
sections = ['SI', 'BB', 'FC']
seasons = ['summer']
variables = ['temperature', 'salinity']
for section in sections:
    for season in seasons:
        for var in variables:
            azst.seasonal_section_plot(VAR=var, SECTION=section, SEASON=season, YEAR=year, ZMAX=500, STATION_BASED=True) 
            plt.close('all')

        command = 'montage temperature_' + section + '_' + season + '_' + str(year) + '.png salinity_' + section + '_' + season + '_' + str(year) + '.png  -tile 2x1 -geometry +10+10  -background white ' + section + '_stn_' + season + '_' + str(year) + '.png'
        os.system(command)
        command = 'montage temperature_' + section + '_' + season + '_' + str(year) + '_FR.png salinity_' + section + '_' + season + '_' + str(year) + '_FR.png  -tile 2x1 -geometry +10+10  -background white ' + section + '_stn_' + season + '_' + str(year) + '_FR.png'
        os.system(command)

        os.system('rm temperature*.png salinity*.png')
        os.system('cp ' + section + '_stn_' + season + '_' + str(year) + '.png ../2021')
        os.system('cp ' + section + '_stn_' + season + '_' + str(year) + '_FR.png ../2021')

# Section CIL
%my_run azmp_CIL_scorecards.py  # update year in script!
os.system('rm scorecards_CIL_SI* scorecards_CIL_BB* scorecards_CIL_FC*')
os.system('cp scorecards_CIL.png scorecards_CIL_FR.png ../2021')
%my_run azmp_CIL_mean_anomaly.py # update year in script!
os.system('cp section_CIL_anomaly.png section_CIL_anomaly_FR.png ../2021')

## ----------- AZMP SAR / IROC ---------------- ##
#%my_run azmp_CIL_stats.py 
%my_run azmp_CIL_stats_update.py # <---- preferred if only an update is needed (need to edit sections)
%my_run azmp_sar_input.py


## ----------- CSAS DATA ---------------- ##
%my_run csas_crab_stats.py
%my_run azmp_bottomT_habitat.py    
%my_run NSRF_bottomT.py 
%my_run NSRF_bottomS.py
%my_run NSRF_bottomT_habitat.py 
%my_run azmp_bottomT_shrimp_habitat.py 
#azrt.bottom_stats(years=np.arange(1980, 2020), season='summer') # for crab 4R
azrt.bottom_temperature(season='summer', year='2020') # for crab 4R
# SFAs
azrt.sfa_bottom_stats(years=np.arange(2006, 2021), season='summer', plot=True, climato_file='Tbot_climato_NSRFx_summer_2006-2018.h5')
azrt.sfa_bottom_scorecards(years=np.arange(2006, 2021), clim_year=[2006, 2020])

# In ~/research/lobster:
%my_run lfa_sst_extract.py
%my_run lfa_sst_analysis.py

