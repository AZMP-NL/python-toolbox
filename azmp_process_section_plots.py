## ----------- SECTION PLOT ---------------- ##
import numpy as np 
import azmp_sections_tools as azst 
import matplotlib.pyplot as plt
import os

## # Original figure generation:
## years = np.arange(1980, 2019)
## sections = ['SI', 'BB', 'WB', 'FC', 'SEGB']
## seasons = ['spring', 'summer', 'fall']
## variables = ['temperature', 'salinity']


## # Mathilde request:
## years = np.arange(1999, 2019)
## sections = ['SI', 'BB', 'WB']
## seasons = ['summer']
## variables = ['temperature', 'salinity']


## for year in years:
##     for section in sections:
##         for season in seasons:
##             for var in variables:
##                 print(str(year), section, season, var)
##                 azst.seasonal_section_plot(VAR=var, SECTION=section, SEASON=season, YEAR=year, ZMAX=500, STATION_BASED=True) 
##                 plt.close('all')

##             command = 'montage temperature_' + section + '_' + season + '_' + str(year) + '.png salinity_' + section + '_' + season + '_' + str(year) + '.png  -tile 2x1 -geometry +10+10  -background white ' + section + '_stn_' + season + '_' + str(year) + '.png'                 
##             os.system(command)
##             os.system('rm temperature*.png salinity*.png')

            ## for var in variables:     
                
            ##     azst.seasonal_section_plot(VAR=var, SECTION=section, SEASON=season, YEAR=year, ZMAX=500, STATION_BASED=False) 
            ##     plt.close('all')

            ## command = 'montage temperature_' + section + '_' + season + '_' + str(year) + '.png salinity_' + section + '_' + season + '_' + str(year) + '.png  -tile 2x1 -geometry +10+10  -background white ' + section + '_itp_' + season + '_' + str(year) + '.png'                 
            ## os.system(command)
            ## os.system('rm temperature*.png salinity*.png')

            
# SI since 1950
years = np.arange(1999,2020)
sections = ['FC']
seasons = ['spring', 'fall']
variables = ['temperature', 'salinity']

for year in years:
    for section in sections:
        for season in seasons:
            for var in variables:     
                print(str(year), section, season, var)
                azst.seasonal_section_plot(VAR=var, SECTION=section, SEASON=season, YEAR=year, ZMAX=1500, STATION_BASED=True) 
                plt.close('all')

            command = 'montage temperature_' + section + '_' + season + '_' + str(year) + '.png salinity_' + section + '_' + season + '_' + str(year) + '.png  -tile 2x1 -geometry +10+10  -background white ' + section + '_station_' + season + '_' + str(year) + '.png'                 
            os.system(command)
            os.system('rm temperature*.png salinity*.png')



 # Bottles sections
years = np.arange(1999,2019)
sections = ['SEGB']
seasons = ['fall']
variables = ['temperature', 'salinity', 'oxygen', 'NO3', 'PO4', 'SIO']

for year in years:
    for section in sections:
        for season in seasons:
            for var in variables:     
                print(str(year), section, season, var)
                azst.btl_section_plot(VAR=var, SECTION=section, SEASON=season, YEAR=year, ZMAX=500) 
                plt.close('all')

            ## command = 'montage temperature_' + section + '_' + season + '_' + str(year) + '.png salinity_' + section + '_' + season + '_' + str(year) + '.png  oxygen_' + section + '_' + season + '_' + str(year) + '.png -tile 3x1 -geometry +10+10  -background white ' + section + '_btl_TSO_' + season + '_' + str(year) + '.png'                 
            ## os.system(command)
            ## os.system('rm temperature*.png salinity*.png')
           
            ## command = 'montage NO3_' + section + '_' + season + '_' + str(year) + '.png PO4_' + section + '_' + season + '_' + str(year) + '.png SIO_' + section + '_' + season + '_' + str(year) + '.png -tile 3x1 -geometry +10+10  -background white ' + section + '_btl_NPS_' + season + '_' + str(year) + '.png'                 
            ## os.system(command)
            ## os.system('rm temperature*.png salinity*.png')
