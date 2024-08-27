'''
Generation of annual maps and 3D cube for temperature from interpolated CTD casts.

Originally developed for H. Andres

Frederic.Cyr@dfo-mpo.gc.ca
June 2022

Run in /home/cyrf0006/AZMP/state_reports/bottomT/T-SCubes

'''


import azmp_utils as azu
import numpy as np


climato_springT = '/home/cyrf0006/AZMP/state_reports/bottomT/Tbot_climato_spring_0.10.h5'
climato_fallT = '/home/cyrf0006/AZMP/state_reports/bottomT/Tbot_climato_fall_0.10.h5'

climato_springS = '/home/cyrf0006/AZMP/state_reports/bottomT/Sbot_climato_spring_0.10.h5'
climato_fallS = '/home/cyrf0006/AZMP/state_reports/bottomT/Sbot_climato_fall_0.10.h5'


for i in np.arange(1993,  2022):
    print(i)
    year_file = '/home/cyrf0006/data/dev_database/netCDF/' + str(i) + '.nc'

    Tbot_dict = azu.get_bottomT(year_file, 'fall', climato_fallT)
    Tbot_dict = azu.get_bottomT(year_file, 'spring', climato_springT)
    
    Sbot_dict = azu.get_bottomS(year_file, 'fall', climato_fallS)
    Sbot_dict = azu.get_bottomS(year_file, 'spring', climato_springS)  



for i in np.arange(1993,  2022):
    # spring
    Tsh5 = 'Tcube_spring' + str(i) + '.h5'
    Ssh5 = 'Scube_spring' + str(i) + '.h5'
    ncfile_spring = 'T-Scube_spring' + str(i) + '.nc'
    # fall
    Tfh5 = 'Tcube_fall' + str(i) + '.h5'
    Sfh5 = 'Scube_fall' + str(i) + '.h5'
    ncfile_fall ='T-Scube_fall' + str(i) + '.nc'
    
    azu.h5cubes_2nc(Tsh5, Ssh5, ncfile_spring)
    azu.h5cubes_2nc(Tfh5, Sfh5, ncfile_fall)
