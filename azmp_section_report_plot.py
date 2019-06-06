'''
Pickled climatologies are generated by azmp_section_clim.py



'''
import os
import netCDF4
import xarray as xr
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import pandas as pd
import numpy as  np
from scipy.interpolate import interp1d  # to remove NaNs in profiles
from scipy.interpolate import griddata
import azmp_sections_tools as azst
import cmocean


## ---- Region parameters ---- ## <-------------------------------Would be nice to pass this in a config file '2017.report'
VAR = 'temperature'
SECTION = 'FC'
SEASON = 'summer'
YEAR = 2018

# derived parameters
if VAR == 'temperature':
    v = np.arange(-2,11,1)
    v_anom = np.linspace(-3.5, 3.5, 8)
    CMAP = cmocean.cm.thermal
elif VAR == 'salinity':
    v = np.arange(29,36,.5)
    v_anom = np.linspace(-1.5, 1.5, 16)
    CMAP = cmocean.cm.haline
else:
    v = 10
    v_anom = 10
    
SECTION_BATHY = SECTION

    
    
# CIL surface (Note that there is a bias because )
def area(vs):
    a = 0
    x0,y0 = vs[0]
    for [x1,y1] in vs[1:]:
        dx = x1-x0
        dy = y1-y0
        a += 0.5*(y0*dx - x0*dy)
        x0 = x1
        y0 = y1
    return a


## ---- Get this year's section ---- ## 
df_section_stn, df_section_itp = azst.get_section(SECTION, YEAR, SEASON, VAR)

## ---- Get climatology ---- ## 
clim_name = 'df_' + VAR + '_' + SECTION + '_' + SEASON + '_clim.pkl' 
df_clim = pd.read_pickle(clim_name)

## ---- Retrieve bathymetry using function ---- ##
bathymetry = azst.section_bathymetry(SECTION_BATHY)

## ---  ---- ## 
df_anom =df_section_itp.reset_index(level=1, drop=True) - df_clim
#df_section_itp = df_section_itp.reset_index(level=0, drop=True)

## ---- plot Figure ---- ##
XLIM = df_section_itp.index[-1][1]
fig = plt.figure()
# ax1
ax = plt.subplot2grid((3, 1), (0, 0))
c = plt.contourf(df_section_itp.index.droplevel(0), df_section_itp.columns, df_section_itp.T, v, cmap=CMAP, extend='max')
if VAR == 'temperature':
    c_cil_itp = plt.contour(df_section_itp.index.droplevel(0), df_section_itp.columns, df_section_itp.T, [0,], colors='k', linewidths=2)
ax.set_ylim([0, 400])
ax.set_xlim([0,  XLIM])
ax.set_ylabel('Depth (m)', fontWeight = 'bold')
ax.invert_yaxis()
Bgon = plt.Polygon(bathymetry,color=np.multiply([1,.9333,.6667],.4), alpha=0.8)
ax.add_patch(Bgon)
plt.colorbar(c)
ax.xaxis.label.set_visible(False)
ax.tick_params(labelbottom='off')
ax.set_title(VAR + 'for section ' + SECTION + ' - ' + SEASON + ' ' + str(YEAR))

# ax2
ax2 = plt.subplot2grid((3, 1), (1, 0))
c = plt.contourf(df_section_itp.index.droplevel(0), df_clim.columns, df_clim.T, v, cmap=CMAP, extend='max')
if VAR == 'temperature':
    c_cil_itp = plt.contour(df_section_itp.index.droplevel(0), df_section_itp.columns, df_clim.T, [0,], colors='k', linewidths=2)
ax2.set_ylim([0, 400])
ax2.set_xlim([0,  XLIM])
ax2.set_ylabel('Depth (m)', fontWeight = 'bold')
ax2.invert_yaxis()
Bgon = plt.Polygon(bathymetry,color=np.multiply([1,.9333,.6667],.4), alpha=0.8)
ax2.add_patch(Bgon)
plt.colorbar(c)
ax2.xaxis.label.set_visible(False)
ax2.tick_params(labelbottom='off')
ax2.set_title('1981-2010 climatology')



# ax3
v=20
ax3 = plt.subplot2grid((3, 1), (2, 0))
c = plt.contourf(df_section_itp.index.droplevel(0), df_section_itp.columns, df_anom.T, v_anom, cmap=cmocean.cm.balance, extend='both')
ax3.set_ylim([0, 400])
ax3.set_xlim([0,  XLIM])
ax3.set_ylabel('Depth (m)', fontWeight = 'bold')
ax3.set_xlabel('Distance (km)', fontWeight = 'bold')
ax3.invert_yaxis()
Bgon = plt.Polygon(bathymetry,color=np.multiply([1,.9333,.6667],.4), alpha=0.8)
ax3.add_patch(Bgon)
plt.colorbar(c)
ax3.set_title(r'Anomaly')

fig.set_size_inches(w=8,h=12)
fig_name = VAR + '_' + SECTION + '_' + SEASON + '_' + str(YEAR) + '.png' 
fig.savefig(fig_name, dpi=200)
os.system('convert -trim ' + fig_name + ' ' + fig_name)


#montage temperature_BB_summer_2018.png salinity_BB_summer_2018.png  -tile 2x1 -geometry +10+10  -background white BB_summer_2018.png 
#montage temperature_SI_summer_2018.png salinity_SI_summer_2018.png  -tile 2x1 -geometry +10+10  -background white SI_summer_2018.png 
#montage temperature_FC_summer_2018.png salinity_FC_summer_2018.png  -tile 2x1 -geometry +10+10  -background white FC_summer_2018.png 
