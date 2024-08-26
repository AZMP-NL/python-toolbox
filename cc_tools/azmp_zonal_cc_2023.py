# -*- coding: utf-8 -*-
"""

To produce AZMP CC SAR figure

For previous years, a script that needed to be edited was used.
Starting in 2023 (2022 reporting), a function from cc_tools is used.

Frederic.Cyr@dfo-mpo.gc.ca
March 2023
"""




import cc_tools as cc
import os

#years = [2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023]
years = [2018, 2023]
#seasons = ['spring','summer','fall']
seasons = ['summer']
#variables = ['Omega_Aragonite_(unitless)', 'pH_Total_(total_scale)', 'Oxygen_Saturation_(%)']
variables = ['Omega_Aragonite_(unitless)', 'Omega_Calcite_(unitless)', 'pH_Total_(total_scale)', 'Oxygen_Saturation_(%)']
#depths = ['surface','bottom']
depths = ['bottom']


for year in years:
    print(str(year))
    for season in seasons:
        print(' ' + season)
        for variable in variables:
            print('  ' + variable )
            for depth in depths:
                print('   ' + depth)
                
                cc.seasonal_map(variable, year, season, depth)
            

# for 2014:
#montage AZMP_OA_2014_fall_OmegaA_surface.png AZMP_OA_2014_fall_pH_surface.png AZMP_OA_2014_fall_DO_perc_surface.png -tile 1x3 -geometry +10+10  -background white AZMP_OA_2014_fall_surface.png

# for 2015:
#montage AZMP_OA_2015_summer_OmegaA_surface.png AZMP_OA_2015_summer_pH_surface.png AZMP_OA_2015_summer_DO_perc_surface.png -tile 1x3 -geometry +10+10  -background white AZMP_OA_2015_summer_surface.png

# for 2017:
#montage AZMP_OA_2017_summer_OmegaA_surface.png AZMP_OA_2017_summer_pH_surface.png AZMP_OA_2017_summer_DO_perc_surface.png -tile 1x3 -geometry +10+10  -background white AZMP_OA_2017_summer_surface.png

# for 2020:
#montage AZMP_OA_2020_summer_OmegaA_surface.png AZMP_OA_2020_summer_pH_surface.png AZMP_OA_2020_summer_DO_perc_surface.png -tile 1x3 -geometry +10+10  -background white AZMP_OA_2020_summer_surface.png

# for 2021:
#montage AZMP_OA_2021_summer_OmegaA_surface.png AZMP_OA_2021_summer_pH_surface.png AZMP_OA_2021_summer_DO_perc_surface.png -tile 1x3 -geometry +10+10  -background white AZMP_OA_2021_summer_surface.png

# for 2018:
#montage AZMP_OA_2018_summer_OmegaA_bottom.png AZMP_OA_2018_summer_pH_bottom.png AZMP_OA_2018_summer_DO_perc_bottom.png -tile 1x3 -geometry +10+10  -background white AZMP_OA_2018_summer_bottom.png

# for 2022:
#montage AZMP_OA_2022_summer_OmegaA_bottom.png AZMP_OA_2022_summer_pH_bottom.png AZMP_OA_2022_summer_DO_perc_bottom.png -tile 1x3 -geometry +10+10  -background white AZMP_OA_2022_summer_bottom.png

#montage AZMP_OA_2023_summer_OmegaA_bottom.png AZMP_OA_2023_summer_pH_bottom.png AZMP_OA_2023_summer_DO_perc_bottom.png -tile 1x3 -geometry +10+10  -background white AZMP_OA_2023_summer_bottom.png
              
# for 2023:
#montage AZMP_OA_2023_summer_OmegaA_bottom.png AZMP_OA_2023_summer_pH_bottom.png AZMP_OA_2023_summer_DO_perc_bottom.png -tile 1x3 -geometry +10+10  -background white AZMP_OA_2023_summer_bottom.png
#montage AZMP_OA_2023_summer_OmegaA_bottom.png AZMP_OA_2023_summer_pH_bottom.png AZMP_OA_2023_summer_DO_perc_bottom.png -tile 1x3 -geometry +10+10  -background white AZMP_OA_2023_summer_bottom.png
#montage AZMP_OA_2023_spring_OmegaA_bottom.png AZMP_OA_2023_spring_pH_bottom.png AZMP_OA_2023_spring_DO_perc_bottom.png -tile 1x3 -geometry +10+10  -background white AZMP_OA_2023_spring_bottom.png                                
#montage AZMP_OA_2021_summer_bottom.png  AZMP_OA_2022_summer_bottom.png -tile 2x1 -geometry +10+10  -background white AZMP_OA_2021-2022_summer_bottom.png

#montage AZMP_OA_2020_summer_bottom.png  AZMP_OA_2022_summer_bottom.png -tile 2x1 -geometry +10+10  -background white AZMP_OA_2020-2022_summer_bottom.png

#montage AZMP_OA_2015_summer_bottom.png  AZMP_OA_2022_summer_bottom.png -tile 2x1 -geometry +10+10  -background white AZMP_OA_2015-2022_summer_bottom.png

#montage AZMP_OA_2014_summer_bottom.png  AZMP_OA_2022_summer_bottom.png -tile 2x1 -geometry +10+10  -background white AZMP_OA_2014-2022_summer_bottom.png

#montage AZMP_OA_2018_summer_bottom.png  AZMP_OA_2023_summer_bottom.png -tile 2x1 -geometry +10+10  -background white AZMP_OA_2018-2023_summer_bottom.png
#montage AZMP_OA_2022_summer_bottom.png  AZMP_OA_2023_summer_bottom.png -tile 2x1 -geometry +10+10  -background white AZMP_OA_2022-2023_summer_bottom.png

#montage AZMP_OA_2023_spring_bottom.png  AZMP_OA_2023_summer_bottom.png  AZMP_OA_2023_summer_bottom.png -tile 3x1 -geometry +10+10  -background white AZMP_OA_2023-spring-summer-summer_bottom.png

os.system('montage AZMP_OA_2018_summer_OmegaA_bottom.png AZMP_OA_2023_summer_OmegaA_bottom.png AZMP_OA_2018_summer_OmegaC_bottom.png AZMP_OA_2023_summer_OmegaC_bottom.png -tile 2x2 -geometry +10+10  -background white SAR_OA_2023_omega.png')

os.system('montage AZMP_OA_2018_summer_pH_bottom.png AZMP_OA_2023_summer_pH_bottom.png AZMP_OA_2018_summer_DO_perc_bottom.png AZMP_OA_2023_summer_DO_perc_bottom.png -tile 2x2 -geometry +10+10  -background white SAR_OA_2023_pH.png')

