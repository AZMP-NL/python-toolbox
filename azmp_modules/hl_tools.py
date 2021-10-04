"""Some tools to process Headlands project data


References

ADD reference to Headland project


----------

Atlantic Zone Monitoring Program @NAFC:
https://azmp-nl.github.io/

"""

__author__ = 'Frederic.Cyr@dfo-mpo.gc.ca'
__version__ = '0.1'

## to do:
# - option for dataFrame with z as columns
# - export multiple dataframe into a netCDF or H5
# - From header, xtract lat, lon, station, line, trip, ship and other attributes (meteo?)

import re
import pandas as pd
import pfile_tools
import numpy as np
import time as tt
import netCDF4 as nc
import os
from sys import version_info

def rpf_eoh():
    """End-of-header key for .RPF files

    """
    # end-of-header key
    eoh = '-- DATA --'

    return eoh


def rpf_header(filename):
    """Returns the header of a given RPF file

    """

    eoh = rpf_eoh()

    # Read header
    header = []
    with open(filename, 'r') as td:
        for line in td:
        
            if re.match(eoh, line): # end-of-header            
                break
            else:
                header.append(line)

    return header

def rpf_header_to_pandas(filename):
    """Returns RPF header in a pandas series

    """

    header = rpf_header(filename)    

    for i in header:
        my_str = i.strip()
        my_str = re.sub(',','', my_str)
        my_str = columns.split(' ')
        
    # clean columns
    columns = re.sub('\n','', columns)
    columns = columns.strip()
    columns = re.sub(' +',' ', columns)
    columns = columns.split(' ')
    # HERE!!!
    
    
    return df

def rpf_to_dataframe(filename):
    """Reads a thermistor .RPF given in 'filename' as returns a Pandas.DataFrame with
       data and depth columns

    """

    header = rpf_header(filename)    
    df = pd.read_csv(filename, sep='\s+',  parse_dates={'datetime': [0, 1]}, header=np.shape(header)[0]+1)
    df = df.set_index('datetime')
    df.columns = ['temperature']
    df = df.replace(9999.99, np.NaN)
        
    return df

