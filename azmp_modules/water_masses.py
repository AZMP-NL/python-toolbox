"""Some tools for NW Atlantic water masses

References:
Initial script based on  T-S properties according to a script from Brian Petrie.
References for these water masses are however unclear and should be referenced.
-> see /home/cyrf0006/research/Pepin_stuff/WaterMass_Summary\ \(Brian\ Petrie\).xls
----------

Atlantic Zone Monitoring Program @NAFC:
https://azmp-nl.github.io/

"""

__author__ = 'Frederic.Cyr@dfo-mpo.gc.ca'
__version__ = '0.1'

import numpy as np
import time as tt
import os
from sys import version_info
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon


def water_masses_def():
    """ Will load water masses definition in memory (lists)

    """
    # NorthAtlantic Water
    NAtlW = [
        (36.18, 20),
        (34.85, 6),
        (35.15, 6),
        (36.48, 20),
        (36.18, 20)
        ]

    # ? Western Slope Water ?
    WSW	 = [
        (34, 6),
        (34.85, 6),
        (36.18, 20),
        (35.33, 20),
        (34,	6),
        ]

    # Newfoundland slope water        	
    NLSW = [
        (35.15, 6),
        (34, 6),
        (32.7, 2),
        (32.7, 0),
        (33, 0),
        (35.15, 4),
        (35.15, 6),
        ]
		
    # Gulf of St. Lawrence Water
    GSLW = [
        (32.7, 2),
        (32.7, 0),
        (32, -0.8),
        (32, -2),
        (30.75, -2),
        (30.75, 0),
        (32.7, 2),
        ]

    # Inner Labrador Current        
    InLC = [
        (32, -2),
        (32, -0.8),
        (32.7, 0),
        (33.5, 0),
        (33.5, -2),
        (32,-2),
        ]	
		
    # Surface Shelf
    SurfShlf = [
        (28, 20),
        (28, -2),
        (30.75, -2),
        (30.75, 0),
        (32.7, 2),
        (33, 2.92),
        (33, 20),
        (28, 20),
        ]

    # Surface slope water
    SurfSlp = [
        (33, 2.92),
        (33, 20),
        (35.33, 20),
        (34, 6),
        (33, 2.92),
        ]
    
    dict = {}
    dict['NAtlW'] = NAtlW
    dict['WSW'] = WSW
    dict['NLSW'] = NLSW
    dict['GSLW'] = GSLW
    dict['InLC'] = InLC
    dict['SurfShlf'] = SurfShlf
    dict['SurfSlp'] = SurfSlp

    return dict


def water_mass_id(T,S):
    """Returns the water mass associated to the given T-S properties.

    """

    wm =  water_masses_def()

    point = Point(S,T)
    wm_id = ''
    for name in wm.keys():
        polygon = Polygon(wm.get(name))
        if polygon.contains(point):
            wm_id = name
    
    return wm_id
        

def get_water_mass(T,S, wm_name):
    """Given a list of T-S properties and the name of a water mass, the function returns the indices that fits in the definition.

    """
    
    indices = []
    for idx, tt in enumerate(T):
        if water_mass_id(T[idx],S[idx]) == wm_name:
            indices.append(idx)
    
    return indices
        
