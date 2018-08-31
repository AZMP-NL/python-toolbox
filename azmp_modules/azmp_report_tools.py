"""Some tools to generate AZMP ResDocs

Contains:
- build_bottomT_climato(infile, year_lims)

----------

Atlantic Zone Monitoring Program @NAFC:
https://azmp-nl.github.io/

"""

__author__ = 'Frederic.Cyr@dfo-mpo.gc.ca'
__version__ = '0.1'

import numpy as np
#import time as tt
import os
import sys
#from sys import version_info
#from shapely.geometry import Point
#from shapely.geometry.polygon import Polygon
import netCDF4
import xarray as xr
from mpl_toolkits.basemap import Basemap
#import openpyxl, pprint
import pandas as pd
from scipy.interpolate import griddata
from scipy.interpolate import interp1d  # to remove NaNs in profiles
import h5py


def add_path(PATH):
    """ Since this module uses (e.g.) bathymetry data, the path to the file must be specified if not already permanent.

    Usage ex: (for bathymetry data)
    import azmp_report_tools as az_r
    az_r.add_path('/home/cyrf0006/data/GEBCO/')
    az_r.add_path('/home/cyrf0006/data/dev_database/')

    ** Turns out to be a useless function...
    """

    sys.path.append(PATH)  # or .insert(0, YOUR_PATH) may give higher priority
    
