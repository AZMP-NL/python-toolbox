# -*- coding: utf-8 -*-
'''
To generate AZMP air variables scorecards

Uses these pickled DataFrames:

1. NAO (from azmp_nao.py):
/home/cyrf0006/AZMP/state_reports/NAO/NAO_annual.pkl
/home/cyrf0006/AZMP/state_reports/NAO/NAO_winter.pkl

2. Air temperature (from azmp_airTemp_fromExcel.py):
/home/cyrf0006/AZMP/state_reports/airTemp/airT_monthly.pkl



'''

import numpy as  np
import matplotlib.pyplot as plt
import pandas as pd
import os
import unicodedata
from matplotlib.colors import from_levels_and_colors

clim_year = [1981, 2010]
years = [1980, 2019]

def is_number(s):
    #https://www.pythoncentral.io/how-to-check-if-a-string-is-a-number-in-python-including-unicode/
    try:
        float(s)
        return True
    except ValueError:
        pass 
    try:
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass 
    return False


#### ------------- NAO & AO ---------------- ####
nao_winter = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/NAO/NAO_winter.pkl')
ao = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/NAO/AO_annual.pkl')
amo = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/NAO/AMO_annual.pkl')
ao.to_csv('ao_annual.csv', float_format='%.2f')
amo.to_csv('amo_annual.csv', float_format='%.2f')

# restrict years
nao_winter = nao_winter[(nao_winter.index>=years[0]) & (nao_winter.index<=years[1])]
ao = ao[(ao.index>=years[0]) & (ao.index<=years[1])]
amo = amo[(amo.index>=years[0]) & (amo.index<=years[1])]
# Rename columns
nao_winter = nao_winter.rename(columns={'Value': r'  $\rm NAO_{winter }$'})
ao = ao.rename(columns={'Value': '    AO'})
amo.name='   AMO' 
df_indices = pd.concat([nao_winter, ao, amo], axis=1)
df_ind_clim = df_indices[(df_indices.index>=clim_year[0]-1) & (df_indices.index<=clim_year[1])]
#df_ind_anom = (df_indices - df_ind_clim.mean()) / df_ind_clim.std()
#df_indices = df_ind_anom


#### ------------- Air Temperature ---------------- ####
df = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/airTemp/airT_monthly.pkl')
# restrict years (years[0]-1 for winter...)
df = df[(df.index.year>=years[0]-1) & (df.index.year<=years[1])]

# 1. Annual mean !!! Verify if better to do average of monthly anomalies for next year (see azmp_stn27_scorecards for MLD and Stratif)
df_annual = df.resample('As').mean()
df_annual = df_annual[df_annual.index.year>=years[0]]
clim_annual = df_annual[(df_annual.index.year>=clim_year[0]) & (df_annual.index.year<=clim_year[1])].mean()
std_annual = df_annual[(df_annual.index.year>=clim_year[0]) & (df_annual.index.year<=clim_year[1])].std()
std_anom_annual = (df_annual - clim_annual)/std_annual
std_anom_annual.index = std_anom_annual.index.year



# 2. Winter only (see azmp_nao.py for explanations)
df_winter = df[(df.index.month==12) | (df.index.month==1) | (df.index.month==2)]
df_winter = df_winter[df_winter.index>pd.to_datetime(str(years[0]-1) + '-10-01')] # start in December
df_winter = df_winter.groupby(np.arange(len(df_winter))//3).mean()
year_unique = pd.unique(df.index.year)[1:,]
df_winter = df_winter.iloc[np.arange(0, year_unique.size)] # reduce if last month is december (belongs to following year)
df_winter.index = df_annual.index
clim_winter = df_winter[(df_winter.index.year>=clim_year[0]) & (df_winter.index.year<=clim_year[1])].mean()
std_winter = df_winter[(df_winter.index.year>=clim_year[0]) & (df_winter.index.year<=clim_year[1])].std()
std_anom_winter = (df_winter - clim_winter)/std_winter
std_anom_winter.index = std_anom_winter.index.year



#### ------------- Build the scorecard ---------------- ####
# preamble
year_list = df_annual.index.year.astype('str')
year_list = [i[2:4] for i in year_list] # 2-digit year
year_list.append(r'$\rm \overline{x}$') # add 2 extra columns
year_list.append(r'sd')
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


# 1. NAO
#my_df = nao_winter.T
my_df = df_indices.T
my_df['MEAN'] = df_ind_clim.mean()
my_df['SD'] =   df_ind_clim.std()
#my_df = my_df.rename({'Value': r'  $\rm NAO_{winter }$'})
my_df.rename(columns={'MEAN': r'$\rm \overline{x}$', 'SD': r'sd'}, inplace=True)

# Get text values +  cell color
vals = np.around(my_df.values,1)
vals[vals==-0.] = 0.
vals_color = vals.copy()
vals_color[:,-1] = 0 # No color to last two columns (mean and STD)
vals_color[:,-2] = 0
# Reverse AMO values (colormap already reverse for AO and NAO)
vals_color[2,:] = vals_color[2,:]*-1

# Change for text
vals[:,-1] = [np.nan, np.nan, np.nan]
vals[:,-2] =   [np.nan, np.nan, np.nan]

nrows, ncols = my_df.index.size+1, my_df.columns.size
fig=plt.figure(figsize=(ncols*wcell+wpad, nrows*hcell+hpad))
ax = fig.add_subplot(111)
ax.axis('off')
#do the table
header = ax.table(cellText=[['']],
                      colLabels=['-- Climate indices --'],
                      loc='center'
                      )
header.set_fontsize(13)
the_table=ax.table(cellText=vals, rowLabels=my_df.index, colLabels=year_list,
                    loc='center', cellColours=cmap_r(norm_r(vals_color)), cellLoc='center',
                    bbox=[0, 0, 1, 0.5]
                    )
# change font color to white where needed:
the_table.auto_set_font_size(False)
the_table.set_fontsize(13)
table_props = the_table.properties()
table_cells = table_props['child_artists']
last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
for key, cell in the_table.get_celld().items():
    cell_text = cell.get_text().get_text() 
    if is_number(cell_text) == False:
        pass
    elif key[0] == 0: #year's row = no color
        pass
    elif key[1] in last_columns:
         cell._text.set_color('darkslategray')
    elif (np.float(cell_text) <= -2) | (np.float(cell_text) >= 2) :
        cell._text.set_color('white')

    if str(cell_text)=='nan':
        cell._set_facecolor('darkgray')
        cell._text.set_color('darkgray')
plt.savefig("scorecards_nao.png", dpi=300)
os.system('convert -trim scorecards_nao.png scorecards_nao.png')

# French table
my_df = my_df.rename({r'  $\rm NAO_{winter }$' : r'   $\rm ONA_{hiver}$'})
my_df = my_df.rename({r'    AO' : r'    OA'})
my_df = my_df.rename({r'   AMO' : r'   OMA'})
year_list[-1] = u'ET'

header = ax.table(cellText=[['']],
                      colLabels=['-- Indices climatiques --'],
                      loc='center'
                      )
header.set_fontsize(13)
the_table=ax.table(cellText=vals, rowLabels=my_df.index, colLabels=year_list,
                    loc='center', cellColours=cmap(norm(vals_color)), cellLoc='center',
                    bbox=[0, 0, 1, 0.5]
                    )
# change font color to white where needed:
the_table.auto_set_font_size(False)
the_table.set_fontsize(13)
table_props = the_table.properties()
table_cells = table_props['child_artists']
last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
for key, cell in the_table.get_celld().items():
    cell_text = cell.get_text().get_text()
    if is_number(cell_text) == False:
        pass
    elif key[0] == 0: #year's row = no color
        pass
    elif key[1] in last_columns:
         cell._text.set_color('darkslategray')
    elif (np.float(cell_text) <= -2) | (np.float(cell_text) >= 2) :
        cell._text.set_color('white')

    if str(cell_text)=='nan':
        cell._set_facecolor('darkgray')
        cell._text.set_color('darkgray')

plt.savefig("scorecards_nao_FR.png", dpi=300)
os.system('convert -trim scorecards_nao_FR.png scorecards_nao_FR.png')



# 2. Winter AirT
my_df = std_anom_winter.T
my_df['MEAN'] = clim_winter
my_df['SD'] = std_winter
my_df.rename(columns={'MEAN': r'$\rm \overline{x}$', 'SD': r'sd'}, inplace=True)

# Get text values +  cell color
vals = np.around(my_df.values,1)
vals[vals==-0.] = 0.
vals_color = vals.copy()
vals_color[:,-1] = 0 # No color to last two columns (mean and STD)
vals_color[:,-2] = 0

nrows, ncols = my_df.index.size+1, my_df.columns.size
fig=plt.figure(figsize=(ncols*wcell+wpad, nrows*hcell+hpad))
ax = fig.add_subplot(111)
ax.axis('off')
#do the table
header = ax.table(cellText=[['']],
                      colLabels=['-- Winter Air Temperature --'],
                      loc='center'
                      )
header.set_fontsize(13)
the_table=ax.table(cellText=vals, rowLabels=my_df.index, colLabels=None,
                    loc='center', cellColours=cmap(norm(vals_color)), cellLoc='center',
                    bbox=[0, 0, 1, 0.5]
                    )
# change font color to white where needed:
the_table.auto_set_font_size(False)
the_table.set_fontsize(13)
table_props = the_table.properties()
table_cells = table_props['child_artists']
last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
for key, cell in the_table.get_celld().items():
    cell_text = cell.get_text().get_text() 
    if is_number(cell_text) == False:
        pass
    ## elif key[0] == 0: # <-- No year, remove
    ##     pass
    elif key[1] in last_columns:
         cell._text.set_color('darkslategray')
    elif (np.float(cell_text) <= -2) | (np.float(cell_text) >= 2) :
        cell._text.set_color('white')
plt.savefig("scorecards_winterAirT.png", dpi=300)
os.system('convert -trim scorecards_winterAirT.png scorecards_winterAirT.png')


# French table
year_list[-1] = u'ET'

header = ax.table(cellText=[['']],
                      colLabels=['-- Temperature hivernale de l\'air --'],
                      loc='center'
                      )
header.set_fontsize(13)
the_table=ax.table(cellText=vals, rowLabels=my_df.index, colLabels=None,
                    loc='center', cellColours=cmap(norm(vals_color)), cellLoc='center',
                    bbox=[0, 0, 1, 0.5]
                    )
# change font color to white where needed:
the_table.auto_set_font_size(False)
the_table.set_fontsize(13)
table_props = the_table.properties()
table_cells = table_props['child_artists']
last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
for key, cell in the_table.get_celld().items():
    cell_text = cell.get_text().get_text()
    if is_number(cell_text) == False:
        pass
    ## elif key[0] == 0: #year's row = no color
    ##     pass
    elif key[1] in last_columns:
         cell._text.set_color('darkslategray')
    elif (np.float(cell_text) <= -2) | (np.float(cell_text) >= 2) :
        cell._text.set_color('white')
plt.savefig("scorecards_winterAirT_FR.png", dpi=300)
os.system('convert -trim scorecards_winterAirT_FR.png scorecards_winterAirT_FR.png')






# 3. Annual AirT
my_df = std_anom_annual.T
my_df['MEAN'] = clim_annual
my_df['SD'] = std_annual
my_df.rename(columns={'MEAN': r'$\rm \overline{x}$', 'SD': r'sd'}, inplace=True)

# Get text values +  cell color
vals = np.around(my_df.values,1)
vals[vals==-0.] = 0.
vals_color = vals.copy()
vals_color[:,-1] = 0 # No color to last two columns (mean and STD)
vals_color[:,-2] = 0

nrows, ncols = my_df.index.size+1, my_df.columns.size
fig=plt.figure(figsize=(ncols*wcell+wpad, nrows*hcell+hpad))
ax = fig.add_subplot(111)
ax.axis('off')
#do the table
header = ax.table(cellText=[['']],
                      colLabels=['-- Annual Air Temperature --'],
                      loc='center'
                      )
header.set_fontsize(13)
the_table=ax.table(cellText=vals, rowLabels=my_df.index, colLabels=None,
                    loc='center', cellColours=cmap(norm(vals_color)), cellLoc='center',
                    bbox=[0, 0, 1, 0.5]
                    )
# change font color to white where needed:
the_table.auto_set_font_size(False)
the_table.set_fontsize(13)
table_props = the_table.properties()
table_cells = table_props['child_artists']
last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
for key, cell in the_table.get_celld().items():
    cell_text = cell.get_text().get_text() 
    if is_number(cell_text) == False:
        pass
    ## elif key[0] == 0: # <-- No year, remove
    ##     pass
    elif key[1] in last_columns:
         cell._text.set_color('darkslategray')
    elif (np.float(cell_text) <= -2) | (np.float(cell_text) >= 2) :
        cell._text.set_color('white')
plt.savefig("scorecards_annualAirT.png", dpi=300)
os.system('convert -trim scorecards_annualAirT.png scorecards_annualAirT.png')


# French table
year_list[-1] = u'ET'

header = ax.table(cellText=[['']],
                      colLabels=['-- Temperature annuelle de l\'air --'],
                      loc='center'
                      )
header.set_fontsize(13)
the_table=ax.table(cellText=vals, rowLabels=my_df.index, colLabels=None,
                    loc='center', cellColours=cmap(norm(vals_color)), cellLoc='center',
                    bbox=[0, 0, 1, 0.5]
                    )
# change font color to white where needed:
the_table.auto_set_font_size(False)
the_table.set_fontsize(13)
table_props = the_table.properties()
table_cells = table_props['child_artists']
last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
for key, cell in the_table.get_celld().items():
    cell_text = cell.get_text().get_text()
    if is_number(cell_text) == False:
        pass
    ## elif key[0] == 0: #year's row = no color
    ##     pass
    elif key[1] in last_columns:
         cell._text.set_color('darkslategray')
    elif (np.float(cell_text) <= -2) | (np.float(cell_text) >= 2) :
        cell._text.set_color('white')
plt.savefig("scorecards_annualAirT_FR.png", dpi=300)
os.system('convert -trim scorecards_annualAirT_FR.png scorecards_annualAirT_FR.png')



#4. Merge all together
# English
os.system('montage  scorecards_nao.png scorecards_winterAirT.png scorecards_annualAirT.png -tile 1x3 -geometry +1+1  -background white  scorecards_air.png') 
# French
os.system('montage  scorecards_nao_FR.png scorecards_winterAirT_FR.png scorecards_annualAirT_FR.png -tile 1x3 -geometry +1+1  -background white  scorecards_air_FR.png') 


