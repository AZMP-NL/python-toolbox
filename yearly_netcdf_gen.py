"""A script to loop on pfile lists and generate multiple (yearly) netCDF files

----------
Frederic.Cyr@dfo-mpo.gc.ca, October 2017

"""

# NOTE: azmp*.list were generated with bash command (folder ~/research/AZMP_database):
# $ for i in `seq 1950 2016`; do ls data/*.p$i > azmp$i.list; done


import pfile_tools as p
import glob
from os import path


lists = glob.glob('*.list')

for yearfile in lists:
    outfile = path.splitext(yearfile)[0] + '.nc'
    p.pfiles_to_netcdf(yearfile, outfile)
    print ' -> ' + outfile + ' done!'
    
    
