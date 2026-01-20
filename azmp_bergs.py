'''
AZMP reporting - Iceberg counts from :
Iceberg data from USCG International Ice Patrol
Contact  Michael Hicks USCG
Hicks, Michael R CIV <Michael.R.Hicks@uscg.mil>
The Iceberg Season runs fron Oct 1 to Sept 30

data in /home/cyrf0006/data/AZMP/ColbourneStuff/NL_ICE_BERG_NUMBERS_DATA_1900_2017.xlsx

(script ran in /home/cyrf0006/AZMP/state_report/bergs)


Frederic.Cyr@dfo-mpo.gc.ca - June 2019

Modifications:

- Jan 2021 to now take data from NSIDC 

'''

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import seaborn as sns
 
# Adjust fontsize/weight
font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 10}
plt.rc('font', **font)

clim_year = [1991, 2020]
current_year = 2025

## ----  Prepare the data ---- ##
# Legacy - Load from Excel sheets
#df = pd.read_excel('/home/cyrf0006/data/AZMP/ColbourneStuff/NL_ICE_BERG_NUMBERS_DATA_1900_2019.xlsx', header=5, index_col='YEAR') 
#df = df.drop(columns = 'TOT SEASON')

# Since 2021, load from NSIDC (https://nsidc.org/data/G10028/versions/1)
# This file is updated manually using the IIP 202* North Atlantic Iceberg Summary (Contact Alexis Denton)
df = pd.read_csv('~/github/AZMP-NL/external_data/NSIDC/Icebergs/G10028_Icebergs_South_of_48N.csv', index_col='YEAR') 
df = df.drop(columns = 'TOTAL')

# Rename columns
df.rename(columns={'OCT':'Oct', 'NOV':'Nov', 'DEC':'Dec', 'JAN':'Jan', 'FEB':'Feb', 'MAR':'Mar', 'APR':'Apr', 'MAY':'May', 'JUN':'Jun', 'JUL':'Jul', 'AUG':'Aug', 'SEP':'Sep'}, inplace=True)
# Stack months under Years (pretty cool!)
#df = df.stack() 

# Transform to a series with values based the 15th of each month (had to convert years to string)
#df_BB.index = pd.to_datetime('15-' + df_BB.index.get_level_values(1) + '-' + df_BB.index.get_level_values(0).values.astype(np.str))

# Annual mean
# We now get the annual number from Section 2.2 of the Ice and Environmental Conditions in Ice Year 202* report
df_annual = df.sum(axis=1)
# Temporary:
df_annual.loc[2022] = 58
df_annual.loc[2023] = 385
df_annual.loc[2024] = 22
df_annual.loc[2025] = 250
df_annual.to_pickle('bergs_annual.pkl')
df_annual_clim = df_annual[(df_annual.index>=clim_year[0]) & (df_annual.index<=clim_year[1])]
df_annual_anom = df_annual - df_annual_clim.mean()
df_annual_std_anom = df_annual_anom/df_annual_clim.std()
df_annual_std_anom.to_pickle('bergs_std_anom.pkl')


# Monthly mean
'''
df_monthly = df[df.index==current_year]
df_monthly_clim = df[(df.index>=clim_year[0]) & (df.index<=clim_year[1])]
df_monthly_std = df_monthly_clim.std(axis=0)
df_monthly_clim = df_monthly_clim.mean(axis=0)


## ---- plot monthly ---- ##
ind = np.arange(len(df_monthly.keys()))  # the x locations for the groups
width = 0.35  # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(ind - width/2, df_monthly_clim.values, width, yerr=df_monthly_std.values*.5,
                label='1991-2020')
rects2 = ax.bar(ind + width/2, np.squeeze(df_monthly.values), width, yerr=None,
                label=str(current_year))

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Counts')
#ax.set_title('Number of icebergs')
ax.set_xticks(ind)
ax.set_xticklabels(df_monthly_clim.index)
ax.legend()
ax.yaxis.grid() # horizontal lines
plt.ylim([0, 400])

fig.set_size_inches(w=6,h=3)
fig_name = 'bergs_monthly.png'
#plt.annotate('data source: www.ncdc.noaa.gov/teleconnections/', xy=(.58, .01), xycoords='figure fraction', annotation_clip=False, FontSize=12)
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)

# Save French Figure
french_months = ['oct', 'nov', 'déc', 'jan', 'fev', 'mar', 'avr', 'mai', 'juin', 'juil', 'aou', 'sep']
ax.set_ylabel('Nombre d\'icebergs')
ax.set_xticklabels(french_months, rotation='horizontal')
fig_name = 'bergs_monthly_FR.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)
'''

## ---- plot annual ---- ##
width = 0.75  # the width of the bars

fig, ax = plt.subplots()
ax.bar(df_annual.index, df_annual.values, width)

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Counts')
#ax.set_title('Number of icebergs')
#ax.set_xticks(ind)
#ax.set_xticklabels(df_annual.index)
#ax.legend()
#ax.yaxis.grid() # horizontal lines
#plt.ylim([0, 330])
plt.grid()
ax.axhspan(df_annual_clim.mean()-df_annual_clim.std()/2, df_annual_clim.mean()+df_annual_clim.std()/2, alpha=0.25, color='gray')

fig.set_size_inches(w=6,h=3)
fig_name = 'bergs_annual.png'
#plt.annotate('data source: www.ncdc.noaa.gov/teleconnections/', xy=(.58, .01), xycoords='figure fraction', annotation_clip=False, FontSize=12)
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)

# Save French Figure
french_months = ['oct', 'nov', 'déc', 'jan', 'fev', 'mar', 'avr', 'mai', 'juin', 'juil', 'aou', 'sep']
ax.set_ylabel('Nombre d\'icebergs')
fig_name = 'bergs_annual_FR.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)



## ---- plot annual normalized (with scorecards) ---- ##
# preamble
from matplotlib.colors import from_levels_and_colors
# Build the colormap
vmin = -3.49
vmax = 3.49
midpoint = 0
levels = np.linspace(vmin, vmax, 15)
midp = np.mean(np.c_[levels[:-1], levels[1:]], axis=1)
colvals = np.interp(midp, [vmin, midpoint, vmax], [-1, 0., 1])
normal = plt.Normalize(-3.49, 3.49)
reds = plt.cm.Reds(np.linspace(0,1, num=7))
blues = plt.cm.Blues_r(np.linspace(0,1, num=7))
whites = [(1,1,1,1)]*2
colors = np.vstack((blues[0:-1,:], whites, reds[1:,:]))
colors = np.concatenate([[colors[0,:]], colors, [colors[-1,:]]], 0)
cmap, norm = from_levels_and_colors(levels, colors, extend='both')
cmap_r, norm_r = from_levels_and_colors(levels, np.flipud(colors), extend='both')
# Common parameters
hcell, wcell = 0.5, 0.6
hpad, wpad = 0, 0

     
fig, ax = plt.subplots() 
ax.bar(df_annual.index, df_annual.values, width)
ax.set_ylabel('Counts')
plt.xlim([1899.5, current_year+0.5])
plt.grid()
ax.axhspan(df_annual_clim.mean()-df_annual_clim.std()/2, df_annual_clim.mean()+df_annual_clim.std()/2, alpha=0.25, color='gray')

std_anom = (df_annual - df_annual_clim.mean()) / df_annual_clim.std()

ax.set_title('Annual icebergs count')
colors = cmap(normal(std_anom.values*-1))
cell_text = [std_anom.values.round(1)]
the_table = ax.table(cellText=cell_text,
        rowLabels=['Icebergs subindex'],
        colLabels=None,
        cellColours = [colors],
        cellLoc = 'center', rowLoc = 'center',
        loc='bottom', bbox=[0, -0.1, 1, 0.05])
the_table.auto_set_font_size (False)
the_table.set_fontsize(6)

for key, cell in the_table.get_celld().items():
    if key[1] == -1:
        cell.set_linewidth(0)
        cell.set_fontsize(10)
    else:
        cell._text.set_rotation(90)

fig.set_size_inches(w=13,h=8)
fig_name = 'icebergs_climate_index.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim -bordercolor White -border 10x10 ' + fig_name + ' ' + fig_name)

## In French ##
fig, ax = plt.subplots() 
ax.bar(df_annual.index, df_annual.values, width)
ax.set_ylabel("Nombre d'icebergs")
plt.xlim([1899.5, 2021.5])
plt.grid()
ax.axhspan(df_annual_clim.mean()-df_annual_clim.std()/2, df_annual_clim.mean()+df_annual_clim.std()/2, alpha=0.25, color='gray')

std_anom = (df_annual - df_annual_clim.mean()) / df_annual_clim.std()

ax.set_title("Décompte annuel d'icebergs")
colors = cmap(normal(std_anom.values*-1))
cell_text = [std_anom.values.round(1)]
the_table = ax.table(cellText=cell_text,
        rowLabels=['Sous-indice Icebergs'],
        colLabels=None,
        cellColours = [colors],
        cellLoc = 'center', rowLoc = 'center',
        loc='bottom', bbox=[0, -0.1, 1, 0.05])
the_table.auto_set_font_size (False)
the_table.set_fontsize(6)

for key, cell in the_table.get_celld().items():
    if key[1] == -1:
        cell.set_linewidth(0)
        cell.set_fontsize(10)
    else:
        cell._text.set_rotation(90)

fig.set_size_inches(w=13,h=8)
fig_name = 'icebergs_climate_index_FR.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim -bordercolor White -border 10x10 ' + fig_name + ' ' + fig_name)
