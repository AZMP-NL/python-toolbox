'''
For ResDoc and presentation on 2020 conditions.
Started Feb. 2021
Frederic.Cyr@dfo-mpo.gc.ca

'''

import os
import matplotlib.pyplot as plt
import numpy as np
import azmp_sections_tools as azst
import azmp_report_tools as azrt  

## ----  For new climatologies ---- ##
import azmp_utils as azu
dc = .1
lonLims = [-60, -43] # fish_hab region
latLims = [39, 56]
lonLims = [-60, -45] # FC AZMP report region
latLims = [42, 56]
lonLims = [-63, -45] # include 2H in above
latLims = [42, 58]
lon_reg = np.arange(lonLims[0]+dc/2, lonLims[1]-dc/2, dc)
lat_reg = np.arange(latLims[0]+dc/2, latLims[1]-dc/2, dc)
azu.get_bottomT_climato('/home/cyrf0006/data/dev_database/netCDF/*.nc', lon_reg, lat_reg, year_lims=[1991, 2020], season='spring', h5_outputfile='Tbot_climato91_spring_0.10.h5')
azu.get_bottomT_climato('/home/cyrf0006/data/dev_database/netCDF/*.nc', lon_reg, lat_reg, year_lims=[1991, 2020], season='fall', h5_outputfile='Tbot_climato91_fall_0.10.h5')
azu.get_bottomS_climato('/home/cyrf0006/data/dev_database/netCDF/*.nc', lon_reg, lat_reg, year_lims=[1991, 2020], season='spring', h5_outputfile='Sbot_climato91_spring_0.10.h5')
azu.get_bottomS_climato('/home/cyrf0006/data/dev_database/netCDF/*.nc', lon_reg, lat_reg, year_lims=[1991, 2020], season='fall', h5_outputfile='Sbot_climato91_fall_0.10.h5')
       

## ---- 2020 update ---- ##
# 1.  NAO (should be a function with year as input)
%my_run azmp_nao.py   
%my_run azmp_ao.py   
%my_run azmp_amo.py
os.system('cp  NAO_winter_1950-2020.png NAO_winter_1950-2020_FR.png ../2020/')


# 2. Air temperature (need to dowload AHCCD and download NUUK update)
%my_run azmp_dmi_nuukAirT.py
%my_run azmp_airTemp.py # use this one since 2020 conditions report
os.system('cp air_temp_2020.png air_temp_2020_FR.png air_temp_anom.png air_temp_anom_FR.png ../2020/')

%my_run azmp_air_scorecards.py
os.system('cp scorecards_air.png scorecards_air_FR.png ../2020/')

# 3. SSTs  <--- Now produced by PSG
# wget -m ftp://ftp.dfo-mpo.gc.ca/bometrics/noaa/stats/boxes/*.stat 
# (in /home/cyrf0006/data/BIO_remote/bometrics/noaa/stats/boxes)
%my_run azmp_SSTs.py # to update bometrics data
%my_run azmp_SSTs_fromExcel.py # To merge with Eugene's historical Excel data
os.system('cp SST_index.png ../2020/')
%my_run azmp_SSTs_scorecards.py # to generate monthly anom, scorecards, etc.
os.system('cp scorecards_sst_yearly.png scorecards_sst_monthly.png ../2020/')


# 4. Station 27 [!clim updated!]
#%my_run viking2020.py  # not done
#os.system('cp Viking2020.png Viking2020_FR.png ../2020')

# need to delete /home/cyrf0006/AZMP/state_reports/stn27/stn27_all_casts.nc
%my_run azmp_stn27.py 
%my_run azmp_stn27_density.py # <-- no need to update climato
%my_run azmp_stn27_analysis.py
%my_run azmp_stn27_scorecards.py
os.system('cp scorecards_s27.png scorecards_s27_FR.png s27_CIL_subplots.png s27_CIL_subplotsFR.png s27_TS_subplots.png s27_TS_subplotsFR.png s27_mld_monthly.png  s27_stratif_monthly.png s27_stratif_monthly_FR.png   ../2020/')
os.system('cp s27_salinity_subplot_2020.png s27_salinity_subplot_2020_FR.png s27_temperature_subplot_2020.png s27_temperature_subplot_2020_FR.png ../2020/')


%my_run azmp_stn27_climateindex_ms.py 

# 5. Sea Ice (/home/cyrf0006/AZMP/state_reports/ice)
%my_run azmp_ice_index.py
os.system('cp ice_index.png ice_index_FR.png ../2020/')

# 6. Icebergs (/home/cyrf0006/AZMP/state_reports/bergs) [!clim updated!]
%my_run azmp_bergs.py
os.system('cp bergs_annual_FR.png bergs_annual.png bergs_monthly_FR.png bergs_monthly.png ../2020')


# 7. bottom temperature maps
azrt.bottom_temperature(season='spring', year='2020', climato_file='Tbot_climato91_spring_0.10.h5') 
azrt.bottom_temperature(season='fall', year='2020', climato_file='Tbot_climato91_fall_0.10.h5')
os.system('cp bottomT_spring2020.png bottomT_spring2020_FR.png bottomT_fall2020.png bottomT_fall2020_FR.png ../2020')
# For NAFO STACFEN and STACFIS input:
#azrt.bottom_temperature(season='summer', year='2020', climato_file='Tbot_climato_SA4_summer_0.10.h5')

# bottom salinity maps
azrt.bottom_salinity(season='spring', year='2020', climato_file='Sbot_climato91_spring_0.10.h5') 
azrt.bottom_salinity(season='fall', year='2020', climato_file='Sbot_climato91_fall_0.10.h5')
os.system('cp bottomS_spring2020.png bottomS_spring2020_FR.png bottomS_fall2020.png bottomS_fall2020_FR.png ../2020')

# bottom stats and scorecards (arange year+1)
azrt.bottom_stats(years=np.arange(1980, 2021), season='spring')
azrt.bottom_stats(years=np.arange(1980, 2021), season='fall')
azrt.bottom_stats(years=np.arange(1980, 2021), season='summer')
azrt.bottom_scorecards(years=[1980, 2020], clim_year=[1991, 2020])
os.system('cp scorecards_botT_spring.png scorecards_botT_spring_FR.png scorecards_botT_fall_FR.png scorecards_botT_fall.png ../2020')
# For NAFO STACFEN and STACFIS input (for azmp_composite_index.py):
azrt.bottom_stats(years=np.arange(1980, 2021), season='summer', climato_file='Tbot_climato_SA4_summer_0.10.h5')

# bottom temperature bar plots
%my_run azmp_bottomT_mean_anomaly.py # same as previous
os.system('cp mean_anomalies_fall.png mean_anomalies_spring.png ../2020')
os.system('cp mean_anomalies_fall_FR.png mean_anomalies_spring_FR.png ../2020')


## ---------------  Sections plots ------------- ##
# in 2020, need to update clim:
%my_run azmp_section_clim # after updating each section name / seasons
# (see azmp_genReport.py for other section plots)
year = 2020
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
        os.system('cp ' + section + '_stn_' + season + '_' + str(year) + '.png ../2020')
        os.system('cp ' + section + '_stn_' + season + '_' + str(year) + '_FR.png ../2020')

# Section CIL
%my_run azmp_air_scorecards.py  
os.system('rm scorecards_CIL_SI* scorecards_CIL_BB* scorecards_CIL_FC*')
os.system('cp scorecards_CIL.png scorecards_CIL_FR.png ../2020')
%my_run azmp_CIL_mean_anomaly.py
os.system('cp section_CIL_anomaly.png section_CIL_anomaly_FR.png ../2020')


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


## ----------- VIKING BUOY ---------------- ##
%my_run viking2017.py # <--- should make a function taking spring/fall and year as input (must run twice now)
os.system('cp Viking2017_sal.png Viking2017_temp.png ~/to_windows/')
os.system('cp Viking2017_sal_FR.png Viking2017_temp_FR.png ~/to_windows/')

