# -*- coding: utf-8 -*-
'''




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

#### ---- Load the data and compute anomalies ---- ####
df_SI = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/sections_plots/CIL/df_CIL_SI_summer.pkl')
df_BB = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/sections_plots/CIL/df_CIL_BB_summer.pkl')
df_FC = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/sections_plots/CIL/df_CIL_FC_summer.pkl')

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
#cmap_r, norm_r = from_levels_and_colors(levels, np.flipud(colors), extend='both')
# Common parameters
hcell, wcell = 0.6, 0.6
hpad, wpad = 1, 10


#### ------------- Scorecard 1 - SI ---------------- ####
df_years = df_SI[(df_SI.index>=years[0]) & (df_SI.index<=years[1])]
df_clim = df_SI[(df_SI.index>=clim_year[0]) & (df_SI.index<=clim_year[1])]
std_anom = (df_years - df_clim.mean()) / df_clim.std()
std_anom_FR = std_anom.copy()

my_df = std_anom.T
my_df['MEAN'] = df_clim.mean()
my_df['SD'] = df_clim.std()
my_df.rename(index={
    'vol_itp' : r'CIL area ($\rm km^2$)',
    'core_itp':r'CIL core ($\rm ^{\circ}C$)',
    'core_depth_itp':'core depth (m)'
    }, inplace=True)

my_df.rename(columns={
    'MEAN' : r'$\rm \overline{x}$',
    'SD':u'sd'
    }, inplace=True)

year_list = std_anom.index.astype('str')
year_list = [i[2:4] for i in year_list] # 2-digit year
year_list.append(r'$\rm \overline{x}$') # add 2 extra columns
year_list.append(r'sd')

# Cell parameters
hcell, wcell = 0.6, 0.6
hpad, wpad = 1, 0

# Get text values +  cell color
vals = np.around(my_df.values,1)
vals[vals==-0.] = 0.
vals_color = vals.copy()
vals_color[:,-1] = 0 # No color to last two columns (mean and STD)
vals_color[:,-2] = 0
# Reverse colorbar for Area and Depth
vals_color[0,:] = -vals_color[0,:]
vals_color[2,:] = -vals_color[2,:]

nrows, ncols = my_df.index.size+1, my_df.columns.size # <--------- +1 because years
fig=plt.figure(figsize=(ncols*wcell+wpad, nrows*hcell+hpad))
ax = fig.add_subplot(111)
ax.axis('off')
#do the table
header = ax.table(cellText=[['']],
                      colLabels=['-- Seal Island section --'],
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
    elif key[0] == 0: # Years are white
        pass
    elif key[1] in last_columns:
         cell._text.set_color('darkslategray')
    elif (np.float(cell_text) <= -2) | (np.float(cell_text) >= 2) :
        cell._text.set_color('white')
    elif (cell_text=='nan'):
        cell._set_facecolor('darkgray')
        cell._text.set_color('darkgray')
plt.savefig('scorecards_CIL_SI.png', dpi=300)
os.system('convert -trim scorecards_CIL_SI.png scorecards_CIL_SI.png')

###  French table ##
my_df = std_anom_FR.T
my_df['MEAN'] = df_clim.mean()
my_df['SD'] = df_clim.std()
my_df.rename(index={
    'vol_itp' : r'surface CIF ($\rm km^2$)',
    'core_itp':r'coeur CIF ($\rm ^{\circ}C$)',
    'core_depth_itp':'coeur prof. (m)'
    }, inplace=True)

my_df.rename(columns={
    'MEAN' : r'$\rm \overline{x}$',
    'SD':u'e-t'
    }, inplace=True)

year_list = std_anom.index.astype('str')
year_list = [i[2:4] for i in year_list] # 2-digit year
year_list.append(r'$\rm \overline{x}$') # add 2 extra columns
year_list.append(r'sd')

# Cell parameters
hcell, wcell = 0.6, 0.6
hpad, wpad = 1, 0

# Get text values +  cell color
vals = np.around(my_df.values,1)
vals[vals==-0.] = 0.
vals_color = vals.copy()
vals_color[:,-1] = 0 # No color to last two columns (mean and STD)
vals_color[:,-2] = 0
# Reverse colorbar for Area and Depth
vals_color[0,:] = -vals_color[0,:]
vals_color[2,:] = -vals_color[2,:]

nrows, ncols = my_df.index.size+1, my_df.columns.size
fig=plt.figure(figsize=(ncols*wcell+wpad, nrows*hcell+hpad))
ax = fig.add_subplot(111)
ax.axis('off')
#do the table
header = ax.table(cellText=[['']],
                      colLabels=['-- section Seal Island --'],
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
    elif key[0] == 0: # Years are white
        pass
    elif key[1] in last_columns:
         cell._text.set_color('darkslategray')
    elif (np.float(cell_text) <= -2) | (np.float(cell_text) >= 2) :
        cell._text.set_color('white')
    elif (cell_text=='nan'):
        cell._set_facecolor('darkgray')
        cell._text.set_color('darkgray')
plt.savefig('scorecards_CIL_SI_FR.png', dpi=300)
os.system('convert -trim scorecards_CIL_SI_FR.png scorecards_CIL_SI_FR.png')



#### ------------- Scorecard 2 - BB ---------------- ####
df_years = df_BB[(df_BB.index>=years[0]) & (df_BB.index<=years[1])]
df_clim = df_BB[(df_BB.index>=clim_year[0]) & (df_BB.index<=clim_year[1])]
std_anom = (df_years - df_clim.mean()) / df_clim.std()
std_anom_FR = std_anom.copy()


my_df = std_anom.T
my_df['MEAN'] = df_clim.mean()
my_df['SD'] = df_clim.std()
my_df.rename(index={
    'vol_itp' : r'CIL area ($\rm km^2$)',
    'core_itp':r'CIL core ($\rm ^{\circ}C$)',
    'core_depth_itp':'core depth (m)'
    }, inplace=True)

my_df.rename(columns={
    'MEAN' : r'$\rm \overline{x}$',
    'SD':u'sd'
    }, inplace=True)


# Cell parameters
hcell, wcell = 0.6, 0.6
hpad, wpad = 1, 0

# Get text values +  cell color
vals = np.around(my_df.values,1)
vals[vals==-0.] = 0.
vals_color = vals.copy()
vals_color[:,-1] = 0 # No color to last two columns (mean and STD)
vals_color[:,-2] = 0
# Reverse colorbar for Area and Depth
vals_color[0,:] = -vals_color[0,:]
vals_color[2,:] = -vals_color[2,:]

nrows, ncols = my_df.index.size, my_df.columns.size # <--------- remove +1 because no year
fig=plt.figure(figsize=(ncols*wcell+wpad, nrows*hcell+hpad))
ax = fig.add_subplot(111)
ax.axis('off')
#do the table
header = ax.table(cellText=[['']],
                      colLabels=['-- Bonavista section --'],
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
    #elif key[0] == 0: # Years are white
    #    pass
    elif key[1] in last_columns:
         cell._text.set_color('darkslategray')
    elif (np.float(cell_text) <= -2) | (np.float(cell_text) >= 2) :
        cell._text.set_color('white')
    elif (cell_text=='nan'):
        cell._set_facecolor('darkgray')
        cell._text.set_color('darkgray')
plt.savefig('scorecards_CIL_BB.png', dpi=300)
os.system('convert -trim scorecards_CIL_BB.png scorecards_CIL_BB.png')

###  French table ##
my_df = std_anom_FR.T
my_df['MEAN'] = df_clim.mean()
my_df['SD'] = df_clim.std()
my_df.rename(index={
    'vol_itp' : r'surface CIF ($\rm km^2$)',
    'core_itp':r'coeur CIF ($\rm ^{\circ}C$)',
    'core_depth_itp':'coeur prof. (m)'
    }, inplace=True)

my_df.rename(columns={
    'MEAN' : r'$\rm \overline{x}$',
    'SD':u'e-t'
    }, inplace=True)

year_list = std_anom.index.astype('str')
year_list = [i[2:4] for i in year_list] # 2-digit year
year_list.append(r'$\rm \overline{x}$') # add 2 extra columns
year_list.append(r'sd')

# Cell parameters
hcell, wcell = 0.6, 0.6
hpad, wpad = 1, 0

# Get text values +  cell color
vals = np.around(my_df.values,1)
vals[vals==-0.] = 0.
vals_color = vals.copy()
vals_color[:,-1] = 0 # No color to last two columns (mean and STD)
vals_color[:,-2] = 0
# Reverse colorbar for Area and Depth
vals_color[0,:] = -vals_color[0,:]
vals_color[2,:] = -vals_color[2,:]

nrows, ncols = my_df.index.size, my_df.columns.size
fig=plt.figure(figsize=(ncols*wcell+wpad, nrows*hcell+hpad))
ax = fig.add_subplot(111)
ax.axis('off')
#do the table
header = ax.table(cellText=[['']],
                      colLabels=['-- section Bonavista --'],
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
    #elif key[0] == 0: # Years are white
    #    pass
    elif key[1] in last_columns:
         cell._text.set_color('darkslategray')
    elif (np.float(cell_text) <= -2) | (np.float(cell_text) >= 2) :
        cell._text.set_color('white')
    elif (cell_text=='nan'):
        cell._set_facecolor('darkgray')
        cell._text.set_color('darkgray')
plt.savefig('scorecards_CIL_BB_FR.png', dpi=300)
os.system('convert -trim scorecards_CIL_BB_FR.png scorecards_CIL_BB_FR.png')




#### ------------- Scorecard 3 - FC ---------------- ####
df_years = df_FC[(df_FC.index>=years[0]) & (df_FC.index<=years[1])]
df_clim = df_FC[(df_FC.index>=clim_year[0]) & (df_FC.index<=clim_year[1])]
std_anom = (df_years - df_clim.mean()) / df_clim.std()
std_anom_FR = std_anom.copy()

my_df = std_anom.T
my_df['MEAN'] = df_clim.mean()
my_df['SD'] = df_clim.std()
my_df.rename(index={
    'vol_itp' : r'CIL area ($\rm km^2$)',
    'core_itp':r'CIL core ($\rm ^{\circ}C$)',
    'core_depth_itp':'core depth (m)'
    }, inplace=True)

my_df.rename(columns={
    'MEAN' : r'$\rm \overline{x}$',
    'SD':u'sd'
    }, inplace=True)

year_list = std_anom.index.astype('str')
year_list = [i[2:4] for i in year_list] # 2-digit year
year_list.append(r'$\rm \overline{x}$') # add 2 extra columns
year_list.append(r'sd')

# Cell parameters
hcell, wcell = 0.6, 0.6
hpad, wpad = 1, 0

# Get text values +  cell color
vals = np.around(my_df.values,1)
vals[vals==-0.] = 0.
vals_color = vals.copy()
vals_color[:,-1] = 0 # No color to last two columns (mean and STD)
vals_color[:,-2] = 0
# Reverse colorbar for Area and Depth
vals_color[0,:] = -vals_color[0,:]
vals_color[2,:] = -vals_color[2,:]

nrows, ncols = my_df.index.size, my_df.columns.size
fig=plt.figure(figsize=(ncols*wcell+wpad, nrows*hcell+hpad))
ax = fig.add_subplot(111)
ax.axis('off')
#do the table
header = ax.table(cellText=[['']],
                      colLabels=['-- Flemish Cap section --'],
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
    #elif key[0] == 0: # Years are white
    #    pass
    elif key[1] in last_columns:
         cell._text.set_color('darkslategray')
    elif (np.float(cell_text) <= -2) | (np.float(cell_text) >= 2) :
        cell._text.set_color('white')
    elif (cell_text=='nan'):
        cell._set_facecolor('darkgray')
        cell._text.set_color('darkgray')
plt.savefig('scorecards_CIL_FC.png', dpi=300)
os.system('convert -trim scorecards_CIL_FC.png scorecards_CIL_FC.png')

###  French table ##
my_df = std_anom_FR.T
my_df['MEAN'] = df_clim.mean()
my_df['SD'] = df_clim.std()
my_df.rename(index={
    'vol_itp' : r'surface CIF ($\rm km^2$)',
    'core_itp':r'coeur CIF ($\rm ^{\circ}C$)',
    'core_depth_itp':'coeur prof. (m)'
    }, inplace=True)

my_df.rename(columns={
    'MEAN' : r'$\rm \overline{x}$',
    'SD':u'e-t'
    }, inplace=True)

year_list = std_anom.index.astype('str')
year_list = [i[2:4] for i in year_list] # 2-digit year
year_list.append(r'$\rm \overline{x}$') # add 2 extra columns
year_list.append(r'sd')

# Cell parameters
hcell, wcell = 0.6, 0.6
hpad, wpad = 1, 0

# Get text values +  cell color
vals = np.around(my_df.values,1)
vals[vals==-0.] = 0.
vals_color = vals.copy()
vals_color[:,-1] = 0 # No color to last two columns (mean and STD)
vals_color[:,-2] = 0
# Reverse colorbar for Area and Depth
vals_color[0,:] = -vals_color[0,:]
vals_color[2,:] = -vals_color[2,:]

nrows, ncols = my_df.index.size, my_df.columns.size
fig=plt.figure(figsize=(ncols*wcell+wpad, nrows*hcell+hpad))
ax = fig.add_subplot(111)
ax.axis('off')
#do the table
header = ax.table(cellText=[['']],
                      colLabels=['-- section Bonnet Flamand --'],
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
    #elif key[0] == 0: # Years are white
    #    pass
    elif key[1] in last_columns:
         cell._text.set_color('darkslategray')
    elif (np.float(cell_text) <= -2) | (np.float(cell_text) >= 2) :
        cell._text.set_color('white')
    elif (cell_text=='nan'):
        cell._set_facecolor('darkgray')
        cell._text.set_color('darkgray')
plt.savefig('scorecards_CIL_FC_FR.png', dpi=300)
os.system('convert -trim scorecards_CIL_FC_FR.png scorecards_CIL_FC_FR.png')



#4. Merge all together
# English
os.system('montage  scorecards_CIL_SI.png scorecards_CIL_BB.png scorecards_CIL_FC.png -tile 1x3 -geometry +1+1  -background white  scorecards_CIL.png') 
# French
os.system('montage  scorecards_CIL_SI_FR.png scorecards_CIL_BB_FR.png scorecards_CIL_FC_FR.png -tile 1x3 -geometry +1+1  -background white  scorecards_CIL_FR.png') 



