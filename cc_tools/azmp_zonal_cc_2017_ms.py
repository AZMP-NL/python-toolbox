# -*- coding: utf-8 -*-
"""

To produce Carbonate manuscript paper.
see cc_tools/seasonal_maps
    
"""
## import os
## import pandas as pd
import matplotlib.pyplot as plt
## from matplotlib.colors import Normalize
## from matplotlib.colors import from_levels_and_colors
## import numpy as np
## import h5py
## import cmocean
## import cmocean.cm as cmo
## import cartopy. crs as ccrs
## import cartopy.feature as cpf
## import matplotlib.ticker as mticker
## from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cc_tools as cc
import cc_variable_list2 as vl


YEAR = 2017
SEASONS = ['spring', 'summer', 'fall']
SEASONS = ['fall']
VARIABLES = ['Temperature_(degC)',
             'Salinity_(psu)',
             'Oxygen_Saturation_(%)',
             'Total_Alkalinity_(umol/kg)',
             'Inorganic_Carbon_(umol/kg)',
             'pCO2_(uatm)',
             'pH_Total_(total_scale)',
             'Omega_Aragonite_(unitless)',
             'Omega_Calcite_(unitless)']
DEPTHS=['surface', 'bottom']



for VARIABLE in VARIABLES:
    for SEASON in SEASONS:
        for DEPTH in DEPTHS:
            print(VARIABLE + ' - ' + str(YEAR) + ' ' + SEASON + ' ' + DEPTH)
            cc.seasonal_map(VARIABLE, YEAR, SEASON, DEPTH)
            plt.close('all')

            
# montage AZMP_OA_2017_spring_Temperature_\(degC\)_surface.png AZMP_OA_2017_spring_Temperature_\(degC\)_bottom.png AZMP_OA_2017_summer_Temperature_\(degC\)_surface.png AZMP_OA_2017_summer_Temperature_\(degC\)_bottom.png AZMP_OA_2017_fall_Temperature_\(degC\)_surface.png AZMP_OA_2017_fall_Temperature_\(degC\)_bottom.png -tile 2x3 -geometry +10+10  -background white  AZMP_OA_2017_Temperature.png

# montage AZMP_OA_2017_spring_Salinity_\(psu\)_surface.png AZMP_OA_2017_spring_Salinity_\(psu\)_bottom.png AZMP_OA_2017_summer_Salinity_\(psu\)_surface.png AZMP_OA_2017_summer_Salinity_\(psu\)_bottom.png AZMP_OA_2017_fall_Salinity_\(psu\)_surface.png AZMP_OA_2017_fall_Salinity_\(psu\)_bottom.png -tile 2x3 -geometry +10+10  -background white  AZMP_OA_2017_Salinity.png

# montage AZMP_OA_2017_spring_Oxygen_Saturation_\(%\)_surface.png AZMP_OA_2017_spring_Oxygen_Saturation_\(%\)_bottom.png AZMP_OA_2017_summer_Oxygen_Saturation_\(%\)_surface.png AZMP_OA_2017_summer_Oxygen_Saturation_\(%\)_bottom.png AZMP_OA_2017_fall_Oxygen_Saturation_\(%\)_surface.png AZMP_OA_2017_fall_Oxygen_Saturation_\(%\)_bottom.png -tile 2x3 -geometry +10+10  -background white  AZMP_OA_2017_Oxygen.png

# montage AZMP_OA_2017_spring_Total_Alkalinity_\(umol-per-kg\)_surface.png AZMP_OA_2017_spring_Total_Alkalinity_\(umol-per-kg\)_bottom.png AZMP_OA_2017_summer_Total_Alkalinity_\(umol-per-kg\)_surface.png AZMP_OA_2017_summer_Total_Alkalinity_\(umol-per-kg\)_bottom.png AZMP_OA_2017_fall_Total_Alkalinity_\(umol-per-kg\)_surface.png AZMP_OA_2017_fall_Total_Alkalinity_\(umol-per-kg\)_bottom.png -tile 2x3 -geometry +10+10  -background white  AZMP_OA_2017_TA.png

# montage AZMP_OA_2017_spring_Inorganic_Carbon_\(umol-per-kg\)_surface.png AZMP_OA_2017_spring_Inorganic_Carbon_\(umol-per-kg\)_bottom.png AZMP_OA_2017_summer_Inorganic_Carbon_\(umol-per-kg\)_surface.png AZMP_OA_2017_summer_Inorganic_Carbon_\(umol-per-kg\)_bottom.png AZMP_OA_2017_fall_Inorganic_Carbon_\(umol-per-kg\)_surface.png AZMP_OA_2017_fall_Inorganic_Carbon_\(umol-per-kg\)_bottom.png -tile 2x3 -geometry +10+10  -background white  AZMP_OA_2017_DIC.png

# montage AZMP_OA_2017_spring_pCO2_\(uatm\)_surface.png AZMP_OA_2017_spring_pCO2_\(uatm\)_bottom.png AZMP_OA_2017_summer_pCO2_\(uatm\)_surface.png AZMP_OA_2017_summer_pCO2_\(uatm\)_bottom.png AZMP_OA_2017_fall_pCO2_\(uatm\)_surface.png AZMP_OA_2017_fall_pCO2_\(uatm\)_bottom.png -tile 2x3 -geometry +10+10  -background white  AZMP_OA_2017_pCO2.png

# montage AZMP_OA_2017_spring_pH_Total_\(total_scale\)_surface.png AZMP_OA_2017_spring_pH_Total_\(total_scale\)_bottom.png AZMP_OA_2017_summer_pH_Total_\(total_scale\)_surface.png AZMP_OA_2017_summer_pH_Total_\(total_scale\)_bottom.png AZMP_OA_2017_fall_pH_Total_\(total_scale\)_surface.png AZMP_OA_2017_fall_pH_Total_\(total_scale\)_bottom.png -tile 2x3 -geometry +10+10  -background white  AZMP_OA_2017_pH.png

# montage AZMP_OA_2017_spring_Omega_Aragonite_\(unitless\)_surface.png AZMP_OA_2017_spring_Omega_Aragonite_\(unitless\)_bottom.png AZMP_OA_2017_summer_Omega_Aragonite_\(unitless\)_surface.png AZMP_OA_2017_summer_Omega_Aragonite_\(unitless\)_bottom.png AZMP_OA_2017_fall_Omega_Aragonite_\(unitless\)_surface.png AZMP_OA_2017_fall_Omega_Aragonite_\(unitless\)_bottom.png -tile 2x3 -geometry +10+10  -background white  AZMP_OA_2017_Omega_A.png

# montage AZMP_OA_2017_spring_Omega_Calcite_\(unitless\)_surface.png AZMP_OA_2017_spring_Omega_Calcite_\(unitless\)_bottom.png AZMP_OA_2017_summer_Omega_Calcite_\(unitless\)_surface.png AZMP_OA_2017_summer_Omega_Calcite_\(unitless\)_bottom.png AZMP_OA_2017_fall_Omega_Calcite_\(unitless\)_surface.png AZMP_OA_2017_fall_Omega_Calcite_\(unitless\)_bottom.png -tile 2x3 -geometry +10+10  -background white  AZMP_OA_2017_Omega_C.png


# Fall 2017 main figures
#montage AZMP_OA_2017_fall_Temperature_surface.png AZMP_OA_2017_fall_Temperature_bottom.png AZMP_OA_2017_fall_Salinity_surface.png AZMP_OA_2017_fall_Salinity_bottom.png AZMP_OA_2017_fall_DO_perc_surface.png AZMP_OA_2017_fall_DO_perc_bottom.png -tile 2x3 -geometry +10+10  -background white Figure06.png

#montage AZMP_OA_2017_fall_TA_surface.png AZMP_OA_2017_fall_TA_bottom.png AZMP_OA_2017_fall_DIC_surface.png AZMP_OA_2017_fall_DIC_bottom.png AZMP_OA_2017_fall_pCO2_surface.png AZMP_OA_2017_fall_pCO2_bottom.png -tile 2x3 -geometry +10+10  -background white Figure07.png

#montage AZMP_OA_2017_fall_pH_surface.png AZMP_OA_2017_fall_pH_bottom.png AZMP_OA_2017_fall_OmegaA_surface.png AZMP_OA_2017_fall_OmegaA_bottom.png AZMP_OA_2017_fall_OmegaC_surface.png AZMP_OA_2017_fall_OmegaC_bottom.png -tile 2x3 -geometry +10+10  -background white Figure08.png

