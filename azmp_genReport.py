# Preliminary attempt to build an "all inclusive" function
# that would generate all figures for AZMP resDoc.

import os
import matplotlib.pyplot as plt
import numpy as np
import azmp_sections_tools as azst
import azmp_report_tools as azrt  
            

# 1.  NAO (should be a function with year as input)
%my_run azmp_nao.py   
%my_run azmp_ao.py   
%my_run azmp_amo.py
os.system('cp NAO_winter_1950-2020.png ../2020/')

# 2. Air temperature (need to update Excel files and download NUUK update)
%my_run azmp_dmi_nuukAirT.py
%my_run azmp_airTemp_fromExcel.py 
os.system('cp air_temp_2020.png air_temp_2020_FR.png air_temp_anom.png air_temp_anom_FR.png ../2020/')

%my_run azmp_air_scorecards.py
os.system('cp scorecards_air.png scorecards_air_FR.png ../2020/')

# 3. SSTs
# wget -m ftp://ftp.dfo-mpo.gc.ca/bometrics/noaa/stats/boxes/*.stat 
# (in /home/cyrf0006/data/BIO_remote/bometrics/noaa/stats/boxes)
%my_run azmp_SSTs.py # to update bometrics data
%my_run azmp_SSTs_fromExcel.py # To merge with Eugene's historical Excel data
os.system('cp SST_index.png ../2020/')
%my_run azmp_SSTs_scorecards.py # to generate monthly anom, scorecards, etc.
os.system('cp scorecards_sst_yearly.png scorecards_sst_monthly.png ../2020/')


# 4. Station 27
%my_run viking2020.py  
os.system('cp Viking2020.png Viking2020_FR.png ../2020')

# 5. Sea Ice (/home/cyrf0006/AZMP/state_reports/ice)
%my_run azmp_ice_index.py
os.system('cp ice_index.png ice_index_FR.png ../2020/')


# 5. bottom temperature maps
azrt.bottom_temperature(season='spring', year='2020') 
azrt.bottom_temperature(season='fall', year='2020')
# For NAFO STACFEN and STACFIS input:
azrt.bottom_temperature(season='summer', year='2020', climato_file='Tbot_climato_SA4_summer_0.10.h5')

# bottom salinity maps
azrt.bottom_salinity(season='spring', year='2020') 
azrt.bottom_salinity(season='fall', year='2020')

# bottom stats and scorecards
azrt.bottom_stats(years=np.arange(1980, 2020), season='spring')
azrt.bottom_stats(years=np.arange(1980, 2020), season='fall')
azrt.bottom_stats(years=np.arange(1980, 2020), season='summer')
azrt.bottom_scorecards(years=[1980, 2020])
os.system('cp scorecards_botT_spring.png scorecards_botT_spring_FR.png scorecards_botT_fall_FR.png scorecards_botT_fall.png ../2020')
# For NAFO STACFEN and STACFIS input (for azmp_composite_index.py):
azrt.bottom_stats(years=np.arange(1980, 2020), season='summer', climato_file='Tbot_climato_SA4_summer_0.10.h5')

# bottom temperature bar plots
%my_run azmp_bottomT_mean_anomaly.py # same as previous
os.system('cp mean_anomalies_fall.png mean_anomalies_spring.png ../2020')
os.system('cp mean_anomalies_fall_FR.png mean_anomalies_spring_FR.png ../2020')

# Sections plots
year = 2020
sections = ['SI', 'BB', 'FC', 'MB']
sections = ['MB']

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
        os.system('cp ' + section + '_stn_' + season + '_' + str(year) + '.png ../2020')
        os.system('cp ' + section + '_stn_' + season + '_' + str(year) + '_FR.png ../2020')



## ----------- CSAS DATA ---------------- ##
%my_run csas_crab_stats.py
%my_run azmp_bottomT_habitat.py    
%my_run NSRF_bottomT.py 
%my_run NSRF_bottomS.py 
%my_run NSRF_bottomT_habitat.py 
%my_run azmp_bottomT_shrimp_habitat.py 
#azrt.bottom_stats(years=np.arange(1980, 2020), season='summer') # for crab 4R
azrt.bottom_temperature(season='summer', year='2020') # for crab 4R

## ----------- VIKING BUOY ---------------- ##
%my_run viking2017.py # <--- should make a function taking spring/fall and year as input (must run twice now)
os.system('cp Viking2017_sal.png Viking2017_temp.png ~/to_windows/')
os.system('cp Viking2017_sal_FR.png Viking2017_temp_FR.png ~/to_windows/')



## ----------- SECTION PLOT ---------------- ##
# (** or look at azmp_section_report_plot.py)
years = np.arange(1980, 2020)
sections = ['SI', 'BB', 'WB', 'FC', 'SEGB', 'SESPB']

seasons = ['spring', 'summer', 'fall']
variables = ['temperature', 'salinity']

# For MB/BI:
years = np.arange(2000, 2020)
sections = ['MB', 'BI']
seasons = ['summer']
variables = ['temperature', 'salinity']

for year in years:
    for section in sections:
        for season in seasons:
            for var in variables:
                azst.seasonal_section_plot(VAR=var, SECTION=section, SEASON=season, YEAR=year, ZMAX=500, STATION_BASED=True) 
                plt.close('all')

            command = 'montage temperature_' + section + '_' + season + '_' + str(year) + '.png salinity_' + section + '_' + season + '_' + str(year) + '.png  -tile 2x1 -geometry +10+10  -background white ' + section + '_stn_' + season + '_' + str(year) + '.png'                 
            os.system(command)
            os.system('rm temperature*.png salinity*.png')

            for var in variables:            
                azst.seasonal_section_plot(VAR=var, SECTION=section, SEASON=season, YEAR=year, ZMAX=500, STATION_BASED=False) 
                plt.close('all')

            command = 'montage temperature_' + section + '_' + season + '_' + str(year) + '.png salinity_' + section + '_' + season + '_' + str(year) + '.png  -tile 2x1 -geometry +10+10  -background white ' + section + '_itp_' + season + '_' + str(year) + '.png'                 
            os.system(command)
            os.system('rm temperature*.png salinity*.png')

            
## ----------- DENSITY SECTION PLOT ---------------- ##
years = np.arange(1999, 2020)
sections = ['BB']
seasons = ['summer']
variables = ['sigma-t']

for year in years:
    for section in sections:
        for season in seasons:
            for var in variables:
                azst.seasonal_section_plot(VAR=var, SECTION=section, SEASON=season, YEAR=year, ZMAX=2000, STATION_BASED=True) 
                plt.close('all')
