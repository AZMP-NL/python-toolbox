# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 09:34:09 2019
add some comments
@author: gibbo
"""
import matplotlib.pyplot as plt
import cmocean
import cmocean.cm as cmo


def variable_parameters(my_variable):

    print(my_variable)
    if my_variable == 'pH_Total_(total_scale)':
        return [7, 7.6, 8.2, 7.9, cmo.thermal_r, [7.6, 7.7, 7.8, 7.9, 8.0, 8.1, 8.2], r'$\rm pH_{T,is}$', 'both', -.5, .5, 11]
        #return [17, 7.4, 8.2, 7.7, cmo.halin_r, [7.4, 7.6, 7.8, 8.0, 8.2], r'$\rm pH_{T,is}$', 'both']
#        return [17, 7.5, 8.3, 7.8, cmo.thermal, [7.5, 7.7, 7.9, 8.1, 8.3], r'$\rm pH_{T,is}$', 'both']

#    elif my_variable == 'temperature':
#        return [25, -2, 22, 12, plt.cm.plasma, [-2, 2, 6, 10, 14, 18, 22], r'$\theta\ (^\circ$C)', 'both']

    elif my_variable == 'Temperature_(degC)':
        return [25, -2, 26, 12, cmo.thermal, [-2, 2, 6, 10, 14, 18, 22, 26], r'$\theta\ (^\circ$C)', 'both']

    elif my_variable == 'Omega_Aragonite_(unitless)':
        return [21, 0, 4, 1, plt.cm.seismic_r, [0, 1, 2, 3, 4], r'$\Omega_{\rm{arg}}$', 'both', -.25, .25, 11]

    elif my_variable == 'Salinity_(psu)':
        return [21, 21, 37, 29, cmo.haline, [21, 25, 29, 33, 37], r'S', 'both', -1, 1, 11]

    ## elif my_variable == 'Oxygen_Saturation_(%)':
    ##     return [11, 0, 100, 30, plt.cm.seismic_r, [0,20,40,60,80,100], r'$\rm O_{2,sat}$ (%)', 'both', -10, 10, 11]

    elif my_variable == 'Oxygen_Saturation_(%)':
        return [9, 60, 100, 80, cmo.ice, [60,70,80,90,100], r'$\rm O_{2,sat}$ (%)', 'both', -10, 10, 11]
    elif my_variable == 'Dissolved_Oxygen_(mL/L)':
        return [16, 6, 9, 7.5, cmo.thermal, [6, 7, 8, 9], r'$\rm O_{2} (mL\,L^{-1})$', 'both', -1, 1, 11]
    elif my_variable == 'AOU':
        return [11, 0, 2, 1, cmo.ice, [0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75], r'$\rm AOU (mL\,L^{-1})$', 'both', -1, 1, 11]

    elif my_variable == 'Total_Alkalinity_(umol/kg)':
        return [20, 1750, 2450, 2100, cmo.thermal_r, [1800, 1900, 2000,2100,2200,2300,2400], r'TA $(\rm \mu $mol/kg)', 'both', -20, 20, 11]

    elif my_variable == 'Inorganic_Carbon_(umol/kg)':
        return [20, 1750, 2350, 2050, cmo.thermal_r, [1800,1900,2000,2100,2200,2300], r'DIC $(\rm \mu $mol/kg)', 'both', -50, 50, 11]
    
    elif my_variable == 'pCO2_(uatm)':
        return [20, 250, 1350, 800, cmo.thermal_r, [250,500,750,1000,1250], r'$\rm pCO_{2} (\rm \mu $atm)', 'both']
    
    elif my_variable == 'chla':
        return [11, 0, 2.5, 1.25, cmo.algae, [0,0.5,1,1.5,2,2.5], r'Chlorophyll a (mg/m$^{3}$))', 'both']    
    elif my_variable == 'Silicate_Concentration_(mmol/m3)':
        return [11, 0, 10, 5, plt.cm.seismic_r, [0, 2, 4, 6, 8], r'$\rm [SiO] (mmol m^{-3})$', 'both', -20, 20, 11]
    
    elif my_variable == 'Nitrate_Concentration_(mmol/m3)':
        return [13, 0, 12, 6, plt.cm.seismic_r, [0, 2, 4, 6, 8, 10], r'$\rm [SiO] (mmol m^{-3})$', 'both', -20, 20, 11]

    elif my_variable == 'Phosphate_Concentration_(mmol/m3)':
        return [71, 0, 1.5, .75, plt.cm.seismic_r, [0, .25, 0.5, .75, 1], r'$\rm [SiO] (mmol m^{-3})$', 'both', -20, 20, 11]
    
    elif my_variable == 'Depth_(dbar)':
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

