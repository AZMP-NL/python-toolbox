# -*- coding: utf-8 -*-
"""

To produce AZMP CC SAR figure

For previous years, a script that needed to be edited was used.
Starting in 2023 (2022 reporting), a function from cc_tools is used.

Frederic.Cyr@dfo-mpo.gc.ca
March 2023
"""




import cc_tools as cc

years = [2017]
#seasons = ['spring','summer','fall']
seasons = ['spring']
variables = ['Omega_Aragonite_(unitless)', 'pH_Total_(total_scale)', 'Oxygen_Saturation_(%)']
depths = ['surface','bottom']


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
#montage AZMP_OA_2014_fall_OmegaA_bottom.png AZMP_OA_2014_fall_pH_bottom.png AZMP_OA_2014_fall_DO_perc_bottom.png -tile 1x3 -geometry +10+10  -background white AZMP_OA_2014_fall_bottom.png

# for 2015:
#montage AZMP_OA_2015_spring_OmegaA_bottom.png AZMP_OA_2015_spring_pH_bottom.png AZMP_OA_2015_spring_DO_perc_bottom.png -tile 1x3 -geometry +10+10  -background white AZMP_OA_2015_spring_bottom.png

# for 2017:
#montage AZMP_OA_2017_spring_OmegaA_bottom.png AZMP_OA_2017_spring_pH_bottom.png AZMP_OA_2017_spring_DO_perc_bottom.png -tile 1x3 -geometry +10+10  -background white AZMP_OA_2017_spring_bottom.png

# for 2020:
#montage AZMP_OA_2020_spring_OmegaA_bottom.png AZMP_OA_2020_spring_pH_bottom.png AZMP_OA_2020_spring_DO_perc_bottom.png -tile 1x3 -geometry +10+10  -background white AZMP_OA_2020_spring_bottom.png

# for 2021:
#montage AZMP_OA_2021_spring_OmegaA_bottom.png AZMP_OA_2021_spring_pH_bottom.png AZMP_OA_2021_spring_DO_perc_bottom.png -tile 1x3 -geometry +10+10  -background white AZMP_OA_2021_spring_bottom.png

# for 2022:
#montage AZMP_OA_2022_spring_OmegaA_bottom.png AZMP_OA_2022_spring_pH_bottom.png AZMP_OA_2022_spring_DO_perc_bottom.png -tile 1x3 -geometry +10+10  -background white AZMP_OA_2022_spring_bottom.png

#montage AZMP_OA_2021_spring_bottom.png  AZMP_OA_2022_spring_bottom.png -tile 2x1 -geometry +10+10  -background white AZMP_OA_2021-2022_spring_bottom.png

#montage AZMP_OA_2020_spring_bottom.png  AZMP_OA_2022_spring_bottom.png -tile 2x1 -geometry +10+10  -background white AZMP_OA_2020-2022_spring_bottom.png

#montage AZMP_OA_2015_spring_bottom.png  AZMP_OA_2022_spring_bottom.png -tile 2x1 -geometry +10+10  -background white AZMP_OA_2015-2022_spring_bottom.png

#montage AZMP_OA_2014_fall_bottom.png  AZMP_OA_2022_fall_bottom.png -tile 2x1 -geometry +10+10  -background white AZMP_OA_2014-2022_fall_bottom.png

#montage AZMP_OA_2017_spring_bottom.png  AZMP_OA_2022_spring_bottom.png -tile 2x1 -geometry +10+10  -background white AZMP_OA_2017-2022_spring_bottom.png
