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

os.system('echo " --------- New run ------ " >> .netcdfgen_log.txt')

for yearfile in lists:
    outfile = os.path.splitext(yearfile)[0] + '.nc'
    p.pfiles_to_netcdf(yearfile, outfile, zbin=5, zmax=2000)
    print ' -> ' + outfile + ' done!'
    expr = 'mv ' + yearfile + ' ./list_done'
    os.system(expr)

## To generate the lists:
## import numpy as np
## for i in np.arange(1912, 2018):
##     #expr = 'ls ./data2/*.p' + np.str(i) + ' > ' + np.str(i) + '.list'
##     # next command takes more time but is better for large lists
##     expr = 'find data2 -type f -name *.p' + np.str(i) + ' > ' + np.str(i) + '.list'
##     print expr
##     os.system(expr)

# too long...
#rsync -avL --ignore-existing cyrf@nl-bpo-dev.ent.dfo-mpo.ca:/data/seabird/symlinks_pfiles/* /home/cyrf0006/data/symlinks_pfiles/. 
