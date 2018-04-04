"""A script to loop on pfile lists and generate multiple (yearly) netCDF files

----------
Frederic.Cyr@dfo-mpo.gc.ca, October 2017

"""

# NOTE: azmp*.list were generated with bash command (folder ~/research/AZMP_database):
# $ for i in `seq 1950 2016`; do ls data/*.p$i > azmp$i.list; done


import pfile_tools as p
import glob
import os

lists = glob.glob('*.list')

for yearfile in lists:
    outfile = os.path.splitext(yearfile)[0] + '.nc'
    p.pfiles_to_netcdf(yearfile, outfile, zbin=5, zmax=6000)
    print ' -> ' + outfile + ' done!'
    expr = 'mv ' + yearfile + ' ./list_done'
    os.system(expr)

## To generate the lists:
## import numpy as np
## for i in np.arange(2017, 2018):
##     #expr = 'ls ./data2/*.p' + np.str(i) + ' > ' + np.str(i) + '.list'
##     # next command takes more time but is better for large lists
##     expr = 'find data2 -type f -name *.p' + np.str(i) + ' > ' + np.str(i) + '.list'
##     print expr
##     os.system(expr)
    
## ## To correct 24:00 time
## import numpy as np
## for i in np.arange(1910, 2017):
##     expr = './data2/*.p' + np.str(i)
##     my_list = glob.glob(expr)
##     for file in my_list:
##         grep_call = 'grep -l 24:00 ' + file
##         out = os.system(grep_call)
##         if out == 0:
##             sed_call = 'sed s/24:00/00:00/ ' + file + ' > ./tmp.tmp'
##             os.system(sed_call)
##             mv_call = 'mv ./tmp.tmp ' + file
##             os.system(mv_call)
##             print file

## To generate the lists from DevClone:
## import numpy as np
## for i in np.arange(1910, 2018):
##     expr = 'find /media/cyrf0006/DevClone/data/seabird/ -type f -name *.p' + np.str(i) + ' > ' + np.str(i) + '.list'
##     print expr
##     os.system(expr)
    
## ## To correct 24:00 time over all files
import numpy as np
lists = glob.glob('*.list')
for yearfile in lists:
    with open(yearfile) as f:
        for file in f:
            grep_call = 'grep -l 24:00 ' + file
            print grep_call
            out = os.system(grep_call)
            if out == 0:                                
                keyboard

                
                sed_call = 'sed s/24:00/00:00/ ' + file + ' > ./tmp.tmp'
                os.system(sed_call)
                mv_call = 'mv ./tmp.tmp ' + file
                os.system(mv_call)
                print file

                

            
