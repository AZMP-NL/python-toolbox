# Preliminary attempt to build an "all inclusive" function
# that would generate all figures for AZMP resDoc.

import os
from seabird.cnv import fCNV
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
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




