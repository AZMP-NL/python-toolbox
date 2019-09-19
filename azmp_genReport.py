# Preliminary attempt to build an "all inclusive" function
# that would generate all figures for AZMP resDoc.

import os
from seabird.cnv import fCNV
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import azmp_sections_tools as azst
from scipy.interpolate import griddata
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

            
## ---------- BOTTOM TEMPERATURE ---------- ##
%my_run azmp_bottomT.py # <--- should make a function taking spring/fall and year as input (must run twice now)
os.system('cp bottomT_fall2017.png bottomT_spring2017.png ~/to_windows/')
os.system('cp bottomT_fall2017_FR.png bottomT_spring2017_FR.png ~/to_windows/')

%my_run azmp_bottomS.py # <--- Now climato is 1981-2010 while Eugene uses 2000-2015
os.system('cp bottomS_fall2017.png bottomS_spring2017.png ~/to_windows/')
os.system('cp bottomS_fall2017_FR.png bottomS_spring2017_FR.png ~/to_windows/')

%my_run azmp_bottom_stats.py # <--- should make a function taking spring/fall as input (must run twice now)
%my_run azmp_bottomT_scorecards.py
os.system('cp scorecards_botT_spring.png scorecards_botT_fall.png ~/to_windows/')
os.system('cp scorecards_botT_spring_FR.png scorecards_botT_fall_FR.png ~/to_windows/')

%my_run azmp_bottomT_mean_anomaly.py
os.system('cp mean_anomalies_fall.png mean_anomalies_spring.png ~/to_windows/')
os.system('cp mean_anomalies_fall_FR.png mean_anomalies_spring_FR.png ~/to_windows/')


## ----------- VIKING BUOY ---------------- ##
%my_run viking2017.py # <--- should make a function taking spring/fall and year as input (must run twice now)
os.system('cp Viking2017_sal.png Viking2017_temp.png ~/to_windows/')
os.system('cp Viking2017_sal_FR.png Viking2017_temp_FR.png ~/to_windows/')



## ----------- SECTION PLOT ---------------- ##
import numpy as np 
import azmp_sections_tools as azst 
import matplotlib.pyplot as plt
import os
years = np.arange(1980, 2019)
sections = ['SI', 'BB', 'WB', 'FC', 'SEGB', 'SESPB']
seasons = ['spring', 'summer', 'fall']
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

            
