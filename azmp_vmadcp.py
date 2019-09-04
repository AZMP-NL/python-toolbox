#wget https://www.ncdc.noaa.gov/teleconnections/nao/data.csv
# check in : /home/cyrf0006/AZMP/annual_meetings/2019

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import datetime
import matplotlib.dates as mdates
import pandas as pd
import os
from sys import version_info
from scipy.interpolate import griddata
import azmp_sections_tools as azst
import cmocean

# Adjust fontsize/weight
font = {'family' : 'DejaVu Sans',
        'weight' : 'bold',
        'size'   : 14}
plt.rc('font', **font)


vmadcp_file = '/home/cyrf0006/data/AZMP/VMADCP/HUDSON/GRID/hud112_fc.grid'


# Load using pandas
df = pd.read_csv(vmadcp_file, sep='\t', parse_dates={'datetime': ['year-mon-day', 'hr:min:sec']}) 
df = df.set_index('datetime')
df.index = pd.to_datetime(df.index)



# Grid data
U = df['U_Vel(cm/s)']
V = df['V_Vel(cm/s)']
Z = df['depth(m)']*-1
X = df['distance(5Km_bin)']

x = X.unique()
z = Z.unique()
x_grid, z_grid = np.meshgrid(x, z)
Ugrid = griddata((X, Z), U, (x_grid, z_grid), method='linear')
Vgrid = griddata((X, Z), V, (x_grid, z_grid), method='linear')

## ---- Retrieve bathymetry using function ---- ##
bathymetry = azst.section_bathymetry('FC')


## ---- Figure ---- ##
VV = np.linspace(-80, 80, 17)
## ---- plot Figure ---- ##
fig = plt.figure()
# ax1
ax = plt.subplot2grid((2, 1), (0, 0))
c = plt.contourf(x, z, Ugrid, VV, cmap=cmocean.cm.balance, extend='both')
#plt.contourf(df_section_itp.index.droplevel(0), df_section_itp.columns, df_section_itp.T, v, cmap=CMAP, extend='max')
#ax.set_ylim([0, 400])
#ax.set_xlim([0,  XLIM])
ax.set_ylabel('Depth (m)', fontWeight = 'bold')
ax.invert_yaxis()
Bgon = plt.Polygon(bathymetry,color=np.multiply([1,.9333,.6667],.4), alpha=1)
ax.add_patch(Bgon)
cb = plt.colorbar(c)
cb.set_label('U (cm/s)', fontsize=13, fontweight='bold')
ax.xaxis.label.set_visible(False)
ax.tick_params(labelbottom='off')
#ax.set_title(VAR + ' for section ' + SECTION + ' - ' + SEASON + ' ' + str(YEAR))
plt.title(df.index[0].strftime('%B %Y'))

# ax2
ax2 = plt.subplot2grid((2, 1), (1, 0))
c = plt.contourf(x, z, Vgrid, VV, cmap=cmocean.cm.balance, extend='both')
#ax2.set_ylim([0, 400])
#ax2.set_xlim([0,  XLIM])
ax2.set_ylabel('Depth (m)', fontWeight = 'bold')
ax2.invert_yaxis()
Bgon = plt.Polygon(bathymetry,color=np.multiply([1,.9333,.6667],.4), alpha=1)
ax2.add_patch(Bgon)
cb = plt.colorbar(c)
cb.set_label('V (cm/s)', fontsize=13, fontweight='bold')

plt.show()




