# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 09:34:09 2019
add some comments
@author: gibbo
"""
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.colors import from_levels_and_colors
import matplotlib as mpl
import numpy as np
import seaborn as sns
import datetime
from scipy import stats
import water_masses as wm
import seawater as swx
import cmocean
import cmocean.cm as cmo
import cartopy. crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.feature as cpf
from cartopy.mpl.geoaxes import GeoAxes
import matplotlib.ticker as mticker
import os
import netCDF4
import io
import h5py
import sys

########import variable_parameters as var_params

#my_variable='pH'

def variable_parameters(my_variable):

    print(my_variable)
    if my_variable == 'pH':
        return [17, 7.5, 8.3, 7.8, cmo.dense, [7.5, 7.7, 7.9, 8.1, 8.3], r'$\rm pH_{T,is}$', 'both']

#    elif my_variable == 'temperature':
#        return [25, -2, 22, 12, plt.cm.plasma, [-2, 2, 6, 10, 14, 18, 22], r'$\theta\ (^\circ$C)', 'both']

    elif my_variable == 'temperature':
        return [25, -2, 26, 12, plt.cm.plasma, [-2, 2, 6, 10, 14, 18, 22, 26], r'$\theta\ (^\circ$C)', 'both']

    elif my_variable == 'Omega_A':
        return [21, 0, 4, 1, plt.cm.seismic_r, [0, 1, 2, 3, 4], r'$\Omega_{\rm{arg}}$', 'both']

    elif my_variable == 'salinity':
        return [21, 21, 37, 29, cmo.haline, [21, 25, 29, 33, 37], r'Salinity', 'both']

#    elif my_variable == 'salinity':
#        return [21, 27, 37, 32, cmo.haline, [27, 29, 31, 33, 35, 37], r'Salinity', 'both']
    
    elif my_variable == 'satO2_perc':
        return [13, 0, 120, 30, plt.cm.seismic_r, [0,20,40,60,80,100, 120], r'$\rm O_{2,sat}$ (%)', 'both']

#    elif my_variable == 'satO2_perc':
#        return [11, 50, 100, 75, plt.cm.plasma, [50,60,70,80,90,100], r'$\rm O_{2,sat}$ (%)', 'both']

    elif my_variable == 'TA':
        return [20, 1750, 2450, 2100, plt.cm.plasma, [1800, 1900, 2000,2100,2200,2300,2400], r'TA $(\rm \mu $mol/kg)', 'both']

    elif my_variable == 'TIC':
        return [20, 1750, 2350, 2050, plt.cm.plasma, [1800,1900,2000,2100,2200,2300], r'DIC $(\rm \mu $mol/kg)', 'both']
    
    elif my_variable == 'pCO2':
        return [20, 250, 1350, 800, plt.cm.plasma, [250,500,750,1000,1250], r'$\rm pCO_{2} (\rm \mu $atm)', 'both']
    
    elif my_variable == 'chla':
        return [11, 0, 2.5, 1.25, cmo.algae, [0,0.5,1,1.5,2,2.5], r'Chlorophyll a (mg/m$^{3}$))', 'both']
    
    elif my_variable == 'Omega_C':
        return [21, 0, 5, 1, plt.cm.seismic_r, [0, 1, 2, 3, 4, 5], r'$\Omega_{\rm{cal}}$', 'both']
    
#    elif my_variable == 'depth':
#        return [26, 0, 250, 125, cmo.deep, [0,50,100,150,200,250], r'Depth (m)', 'max']
    
    elif my_variable == 'depth':
        return [26, 0, 250, 125, cmo.deep, [0,50,100,150,200,250], r'Depth (m)', 'max']
            
    elif my_variable == 'region_number':
        return [26, 0, 300, 250, plt.cm.plasma, [0,100,200,300], r'Depth (m)', 'max']
    
        
    elif my_variable == 'kml':
        return [100, 0, 2, 1, plt.cm.plasma, [0,1,2], r'Depth (m)', 'max']
    
    else:
        print ('nothing')
#var_y=variable_parameters('pH')

#for i, v in zip(plt_var, axis):
#    
#    if i == 'pH':
#        num_levels = 17
#        vmin = 7.5
#        vmax = 8.3
#        midpoint = 7.8
#        colors = cmo.dense
#        ticks = [7.5, 7.6, 7.7, 7.8, 7.9, 8.0, 8.1, 8.2, 8.3]
#        axis_label = r'$\rm pH_{T,is}$'
#        extent = 'both'
#    
#    if i == 'Omega_A':
#        num_levels = 21
#        vmin = 0
#        vmax = 4
#        midpoint = 1
#        colors = plt.cm.seismic_r
#        ticks = [0, 1, 2, 3, 4]
#        axis_label = r'$\Omega_{\rm{arg}}$'
#        extent = 'both'
#    
#    if i == 'temperature':
#        num_levels = 25
#        vmin = -2
#        vmax = 22
#        midpoint = 12
#        colors = plt.cm.plasma
#        ticks = [-2, 2, 6, 10, 14, 18, 22]
#        axis_label = r'Temperature ($^\circ$C)'
#        extent = 'both'
#    
#    if i == 'salinity':
#        num_levels = 21
#        vmin = 22#18
#        vmax = 38#38
#        midpoint = 32#28
#        colors = cmo.haline
#        ticks = [27, 29, 31, 33, 35, 37]
#        ticks = [18,22,26,30,34,38]
#        axis_label = r'Salinity'
#        extent = 'both'
#    
#    if i == 'satO2_perc':
#        num_levels = 11
#        vmin = 0
#        vmax = 100
#        midpoint = 30
#        colors = plt.cm.seismic_r
#        ticks = [0,20,40,60,80,100]
#        axis_label = r'$\rm O_{2,sat}$ (%)'
#        extent = 'max'
#        
#    if i == 'TA':
#        num_levels = 20
#        vmin = 1800
#        vmax = 2500
#        midpoint = 2125
#        colors = plt.cm.plasma
#        ticks = [1800, 1900, 2000,2100,2200,2300,2400]
#        axis_label = r'Total Alkalinity $(\rm \mu $mol/kg)'
#        extent = 'both'
#        
#    if i == 'TIC':
#        num_levels = 20
#        vmin = 1800
#        vmax = 2400
#        midpoint = 2100
#        colors = plt.cm.plasma
#        ticks = [1800,1900,2000,2100,2200,2300, 2400]
#        axis_label = r'DIC $(\rm \mu $mol/kg)'
#        extent = 'both'
#
#    if i == 'nTA':
#        num_levels = 20
#        vmin = 2100
#        vmax = 2400
#        midpoint = 2250
#        colors = plt.cm.plasma
#        ticks = [2100,2200,2300,2400]
#        axis_label = r'nTA $(\rm \mu $mol/kg)'
#        extent = 'both'
#        
#    if i == 'nTIC':
#        num_levels = 20
#        vmin = 1900
#        vmax = 2300
#        midpoint = 2100
#        colors = plt.cm.plasma
#        ticks = [1900,2000,2100,2200,2300]
#        axis_label = r'nDIC $(\rm \mu $mol/kg)'
#        extent = 'both'        
#
#    if i == 'pCO2':
#        num_levels = 20
#        vmin = 250
#        vmax = 1350
#        midpoint = 800
#        colors = plt.cm.plasma
#        ticks = [250,500,750,1000,1250]
#        axis_label = r'$\rm pCO_{2} (\rm \mu $atm)'
#        extent = 'both'
#
#    if i == 'pCO2T':
#        num_levels = 20
#        vmin = 100
#        vmax = 600
#        midpoint = 610
#        colors = plt.cm.plasma
#        ticks = []
#        axis_label = r'$\rm pCO_{2} (T_{mean})(\rm \mu $atm)'
#        extent = 'both'
#
#    if i == 'pCO2Tdiff':
#        num_levels = 20
#        vmin = -300
#        vmax = 300
#        midpoint = 0
#        colors = plt.cm.seismic_r
#        ticks = []
#        axis_label = r'$\rm pCO_{2} (T-T_{mean})(\rm \mu $atm)'
#        extent = 'both'
#        
#    if i == 'depth':
#        num_levels = 26
#        vmin = 0
#        vmax = 300
#        midpoint = 250
#        colors = cmo.deep
#        ticks = [0,100,200,300]
#        axis_label = r'Depth (m)'    
#        extent = 'max'
#        
#    if i == 'TADIC':
#        num_levels = 11
#        vmin = 0.95
#        vmax = 1.15
#        midpoint = 1.1
#        colors = plt.cm.seismic_r
#        ticks = [0.8,0.9,1,1.1,1.2]
#        axis_label = r'TA:DIC'
#        extent = 'both'
#
#    if i == 'oxygen':  #/44.661 ##change units of GSL data from umol/L to ml/L (/44.661)
#        num_levels = 11
#        vmin = 3.3
#        vmax = 8.3
#        midpoint = 5.75
#        colors = plt.cm.seismic_r
#        ticks = [3.3, 4.3, 5.3, 6.3, 7.3, 8.3]
#        axis_label = r'$\rm O_{2}(\rm \mu $mol/kg)' 
#        extent = 'both'
#        
#    if i == 'chla':
#        df = df.dropna()
#        num_levels = 11
#        vmin = 0
#        vmax = 2.5
#        midpoint = 1.25
#        colors = cmo.algae
#        ticks = [0,0.5,1,1.5,2,2.5]
#        axis_label = r'Chlorophyll a (mg/m$^{3}$))'
#    
#    if i == 'NO3':
#        df = df.dropna()
#        num_levels = 13
#        vmin = 0
#        vmax = 25
#        midpoint = 12.5
#        colors = cmo.algae
#        ticks = [0,5,10,15,20,25]
#        axis_label = r'NO3 (mg/m$^{3}$))'
#        
#    locals()['vmin_'+(v)] = vmin
#    locals()['vmax_'+(v)] = vmax
#    locals()['axis_label_'+(v)] = axis_label
#    locals()['ticks_'+(v)] = ticks

