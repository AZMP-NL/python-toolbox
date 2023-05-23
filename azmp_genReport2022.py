'''
For ResDoc and presentation on 2022 conditions.

This script is usuall run into:
/home/cyrf0006/AZMP/state_reports/

with subfolders:
airTemp/
bergs/
bottomT/
CIL/
climate_index/
LC_transport/
SSTs/

and figures further moved here:
2022/

Started Jan. 2023
Frederic.Cyr@dfo-mpo.gc.ca

'''

import os
import matplotlib.pyplot as plt
import numpy as np
import azmp_sections_tools as azst
import azmp_report_tools as azrt  
import azmp_genreport as azgen
import cc_tools as cc

## ---- 2022 update ---- ## [DONE 2022x]
# 1.  NAO (should be a function with year as input) [Done 2021]
# in /home/cyrf0006/AZMP/state_reports/airTemp
# previously:
#%my_run azmp_nao.py # had to re-write in 2021
#%my_run azmp_ao.py  # had to re-write in 2021
#%my_run azmp_amo.py
# New from 2022:
azgen.nao(2022)
azgen.ao(2022)
azgen.amo(2022)
os.system('cp  NAO_winter_1950-2022.png NAO_winter_1950-2022_FR.png ../2022/')


# 2. Air temperature (need to dowload AHCCD and download NUUK update) [DONE 2022x]
%my_run azmp_dmi_nuukAirT.py
%my_run azmp_airTemp.py # use this one since 2020 conditions report
os.system('cp air_temp_2022.png air_temp_2022_FR.png air_temp_anom.png air_temp_anom_FR.png air_temp_climate_index.png air_temp_climate_index_FR.png ../2022/')

%my_run azmp_air_scorecards.py
os.system('cp scorecards_air.png scorecards_air_FR.png ../2022/')

# 3. SSTs  [Received from PSG]
# wget -m ftp://ftp.dfo-mpo.gc.ca/bometrics/noaa/stats/boxes/*.stat 
# (in /home/cyrf0006/data/BIO_remote/bometrics/noaa/stats/boxes)
#%my_run azmp_SSTs.py # to update bometrics data
#%my_run azmp_SSTs_fromExcel.py # To merge with Eugene's historical Excel data
#os.system('cp SST_index.png ../2022/')
#%my_run azmp_SSTs_scorecards.py # to generate monthly anom, scorecards, etc.
#os.system('cp scorecards_sst_yearly.png scorecards_sst_monthly.png ../2022/')


# 4. Station 27 [Done 2022; but problems with data]
#%my_run viking2022.py  # NOT IN 2022
#os.system('cp Viking2022.png Viking2022_FR.png ../2022')

# need to delete /home/cyrf0006/AZMP/state_reports/stn27/stn27_all_casts.nc
%my_run azmp_stn27.py # ** Note that the script below should be run first to update stn27_all_casts.nc
%my_run azmp_stn27_density.py # <-- ** Run this script first with "append" option
%my_run azmp_stn27_analysis.py
%my_run azmp_stn27_scorecards.py
os.system('cp scorecards_s27.png scorecards_s27_FR.png s27_CIL_subplots.png s27_CIL_subplots_FR.png s27_TS_subplots.png s27_TS_subplotsFR.png s27_mld_monthly.png  s27_stratif_monthly.png s27_stratif_monthly_FR.png   ../2022/')
os.system('cp s27_salinity_subplot_2022.png s27_salinity_subplot_2022_FR.png s27_temperature_subplot_2022.png s27_temperature_subplot_2022_FR.png ../2022/')

# To prepare MS
#%my_run azmp_stn27_climateindex_ms.py 

# 5. Sea Ice (/home/cyrf0006/AZMP/state_reports/ice) [Received from PSG]

#%my_run azmp_ice_index.py
#os.system('cp ice_index.png ice_index_FR.png ../2022/')


# 6. Icebergs (/home/cyrf0006/AZMP/state_reports/bergs) [DONE 2021]
%my_run azmp_bergs.py
os.system('cp bergs_annual_FR.png bergs_annual.png bergs_monthly_FR.png bergs_monthly.png ../2022')


# 7. bottom temperature maps [Done 2022x]
azrt.bottom_temperature(season='spring', year='2022', climato_file='Tbot_climato_spring_0.10.h5') 
azrt.bottom_temperature(season='fall', year='2022', climato_file='Tbot_climato_fall_0.10.h5')
os.system('cp bottomT_spring2022.png bottomT_spring2022_FR.png bottomT_fall2022.png bottomT_fall2022_FR.png ../2022')

# For NAFO STACFEN and STACFIS input: [NEED TO DO]
azrt.bottom_temperature(season='summer', year='2022', climato_file='Tbot_climato_SA4_summer_0.10.h5')

# bottom salinity maps [Done 2022x]
azrt.bottom_salinity(season='spring', year='2022', climato_file='Sbot_climato_spring_0.10.h5') 
azrt.bottom_salinity(season='fall', year='2022', climato_file='Sbot_climato_fall_0.10.h5')
os.system('cp bottomS_spring2022.png bottomS_spring2022_FR.png bottomS_fall2022.png bottomS_fall2022_FR.png ../2022')

# bottom stats and scorecards (arange year+1)
#[need to flag years if coverage insufficient]
azrt.bottom_stats(years=np.arange(1980, 2023), season='spring') 
azrt.bottom_stats(years=np.arange(1980, 2023), season='fall')
#azrt.bottom_stats(years=np.arange(1980, 2023), season='summer')
azrt.bottom_scorecards(years=[1980, 2023], clim_year=[1991, 2020])
os.system('cp scorecards_botT_spring.png scorecards_botT_spring_FR.png scorecards_botT_fall_FR.png scorecards_botT_fall.png ../2022')

# For NAFO STACFEN and STACFIS input (for azmp_composite_index.py): [NEED TO DO]
azrt.bottom_stats(years=np.arange(1980, 2022), season='summer', climato_file='Tbot_climato_SA4_summer_0.10.h5')

# bottom temperature bar plots [need to flag years if coverage insufficient]
%my_run azmp_bottomT_mean_anomaly.py # same as previous
os.system('cp bottomT_anomalies_climateindex.png bottomT_anomalies_climateindex_FR.png ../2022')


## ---------------  Sections plots ------------- ## [Done 2022!]
year = 2022 
sections = ['SI', 'BB', 'FC']
seasons = ['fall']
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
        os.system('cp ' + section + '_stn_' + season + '_' + str(year) + '.png ../2022')
        os.system('cp ' + section + '_stn_' + season + '_' + str(year) + '_FR.png ../2022')

# Section CIL [CANNOT UPDATE 2022]
%my_run azmp_CIL_scorecards.py  # update year in script!
os.system('rm scorecards_CIL_SI* scorecards_CIL_BB* scorecards_CIL_FC*')
os.system('cp scorecards_CIL.png scorecards_CIL_FR.png ../2022')
%my_run azmp_CIL_mean_anomaly.py # update year in script!
os.system('cp section_CIL_anomaly.png section_CIL_anomaly_FR.png ../2022')

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


