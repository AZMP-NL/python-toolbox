'''
General script to generate btl section plots
'''

## ----------- SECTION PLOT ---------------- ##
import numpy as np 
import azmp_sections_tools as azst 
import matplotlib.pyplot as plt
import os

# summer
years = np.arange(1999, 2019)
sections = ['SI', 'BB', 'WB', 'FC']
seasons = ['summer']
variables = ['satO2_perc', 'PO4', 'NO3', 'SIO', 'oxygen']

# spring and fall
years = np.arange(1999, 2019)
sections = ['BB', 'FC', 'SEGB', 'SESPB']
seasons = ['spring', 'fall']
variables = ['satO2_perc', 'PO4', 'NO3', 'SIO', 'oxygen']

# summer
years = np.arange(1999, 2019)
sections = ['SI', 'BB', 'WB', 'FC']
seasons = ['summer']
variables = ['satO2_perc', 'PO4', 'NO3', 'SIO', 'oxygen']

for year in years:
    for section in sections:
        for season in seasons:
            for var in variables:
                print(str(year), section, season, var)
                azst.btl_section_plot(VAR=var, SECTION=section, SEASON=season, YEAR=year, ZMAX=500) 
                plt.close('all')

