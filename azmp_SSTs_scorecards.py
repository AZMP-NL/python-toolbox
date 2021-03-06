# -*- coding: utf-8 -*-
'''
To generate AZMP SST scorecards

Uses this pickled DataFrame:

/home/cyrf0006/AZMP/state_reports/SSTs/SSTs_merged_monthly.pkl
generated by from azmp_sst_scorecards.py


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
df_month = pd.read_pickle('/home/cyrf0006/AZMP/state_reports/SSTs/SSTs_merged_monthly.pkl')
df_year = df_month.resample('As').mean()
df_clim_month = df_month[(df_month.index.year>=clim_year[0]) & (df_month.index.year<=clim_year[1])]
df_clim_year = df_clim_month.resample('As').mean()

#### ------------- Scorecard 1 - Monthly ---------------- ####
# anomalies
df_clim_month['month'] = df_clim_month.index.month   
monthly_clim = df_clim_month.groupby('month').mean()
monthly_std = df_clim_month.groupby('month').std()
this_year = df_month[(df_month.index.year==years[1])]
this_year.index = this_year.index.month
this_year.index.name = 'month'
monthly_std_anom = (this_year - monthly_clim) / monthly_std
# Reorder columns
cols = [
    'Greenland_Shelf',
    'North_Central_Labrador_Sea',
    'Hudson_Strait',
    'Central_Labrador_Sea',
    'Bravo',
    'Hamilton_Bank',
    'St.Anthony_Basin',
    'Northeast_Nfld_Shelf',
    'Orphan_Knoll',
    'Avalon_Channel',
    'Hybernia',
    'Flemish_Pass', 
    'Flemish_Cap',
    'Green-St._Pierre_Bank'
    ]
monthly_std_anom = monthly_std_anom[cols]

# Rename columns
monthly_std_anom_FR = monthly_std_anom.copy()
monthly_std_anom.rename(columns={
    'Avalon_Channel' : 'Avalon Channel (AC)',
    'Bravo':'Bravo (BRA)',
    'Central_Labrador_Sea':'Cent. Lab. Sea (CLS)',
    'Flemish_Cap' : 'Flemish Cap (FC)',
    'Flemish_Pass' : 'Flemish Pass (FP)',
    'Greenland_Shelf' : 'Greenland Shelf (GS)',
    'Green-St._Pierre_Bank':'St.Pierre Bank (SPB)' ,
    'Hamilton_Bank':'Hamilton Bank (HB)',
    'Hudson_Strait': 'Hudson Strait (HS)',
    'Hybernia' : 'Hibernia (HIB)',
    'North_Central_Labrador_Sea':'NC Lab. Sea (NCLS)',
    'Northeast_Nfld_Shelf':'NE NF shelf (NENL)',
    'Orphan_Knoll':'Orphan Knoll (OK)',
    'St.Anthony_Basin':'St. Anthony B. (SAB)'
    }, inplace=True)

monthly_std_anom_FR.rename(columns={
    'Avalon_Channel' : 'Chenal Avalon (AC)',
    'Bravo':'Bravo (BRA)',
    'Central_Labrador_Sea':'Mer Lab. Centrale (CLS)',
    'Flemish_Cap' : 'Bonnet Flammand (FC)',
    'Flemish_Pass' : 'Passe Flamande (FP)',
    'Greenland_Shelf' : 'Plateau Groenland (GS)',
    'Green-St._Pierre_Bank':'Banc St.Pierre (SPB)' ,
    'Hamilton_Bank':'Banc Hamilton (HB)',
    'Hudson_Strait': 'Detroit Hudson (HS)',
    'Hybernia' : 'Hibernia (HIB)',
    'North_Central_Labrador_Sea':'Mer Lab. Nord Cent. (NCLS)',
    'Northeast_Nfld_Shelf':'Plateau T-N NE (NENL)',
    'Orphan_Knoll':'Orphan Knoll (OK)',
    'St.Anthony_Basin':'Basin St. Anthony (SAB)'
    }, inplace=True)
                                      

# preamble
month_list = ['J','F','M','A','M','J','J','A','S','O','N','D']
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
# Common parameters
hcell, wcell = 0.6, 0.6
hpad, wpad = 1, 10

my_df = monthly_std_anom.T

# Get text values +  cell color
vals = np.around(my_df.values,1)
vals[vals==-0.] = 0.
vals_color = vals.copy()
#vals_colors[np.isnan(vals_colors)] = 0.

nrows, ncols = my_df.index.size+1, my_df.columns.size
fig=plt.figure(figsize=(ncols*wcell+wpad, nrows*hcell+hpad))
ax = fig.add_subplot(111)
ax.axis('off')
#do the table
header = ax.table(cellText=[['']],
                      colLabels=['-- 2019 Monthly Sea Surface Temperature anomalies --'],
                      loc='center'
                      )
header.set_fontsize(13)
the_table=ax.table(cellText=vals, rowLabels=my_df.index, colLabels=month_list,
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
        cell._text.set_font_properties('normal') # Not really bold, but smaller thus fits better.
    elif key[0] == 0: #year's row = no color
        pass
    ## elif key[1] in last_columns:
    ##      cell._text.set_color('darkslategray')
    elif (np.float(cell_text) <= -2.) | (np.float(cell_text) >= 2.) :
        cell._text.set_color('white')
    elif (cell_text=='nan'):
        cell._set_facecolor('darkgray')
        cell._text.set_color('darkgray')
cellDict = the_table.get_celld()
for x in range(1, len(my_df.index)+1): # right-align row label
    cellDict[(x,-1)]._loc = 'right'
    
plt.savefig("scorecards_sst_monthly.png", dpi=300)
os.system('convert -trim scorecards_sst_monthly.png scorecards_sst_monthly.png')

## French table
my_df = monthly_std_anom_FR.T
fig=plt.figure(figsize=(ncols*wcell+wpad, nrows*hcell+hpad))
ax = fig.add_subplot(111)
ax.axis('off')
#do the table
header = ax.table(cellText=[['']],
                      colLabels=['-- Anomalies de température de surface mensuelle en 2019 --'],
                      loc='center'
                      )
header.set_fontsize(13)
the_table=ax.table(cellText=vals, rowLabels=my_df.index, colLabels=month_list,
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
        cell._text.set_font_properties('normal') # Not really bold, but smaller thus fits better.
    elif key[0] == 0: #year's row = no color
        pass
    ## elif key[1] in last_columns:
    ##      cell._text.set_color('darkslategray')
    elif (np.float(cell_text) <= -2.) | (np.float(cell_text) >= 2.) :
        cell._text.set_color('white')
    elif (cell_text=='nan'):
        cell._set_facecolor('darkgray')
        cell._text.set_color('darkgray')
cellDict = the_table.get_celld()
for x in range(1, len(my_df.index)+1): # right-align row label
    cellDict[(x,-1)]._loc = 'right'
    
plt.savefig("scorecards_sst_monthly_FR.png", dpi=300)
os.system('convert -trim scorecards_sst_monthly_FR.png scorecards_sst_monthly_FR.png')



#### ------------- Scorecard 2 - Annual ---------------- ####
## Yearly anomalies
yearly_std_anom = (df_year - df_clim_year.mean()) / df_clim_year.std()
yearly_std_anom.index = yearly_std_anom.index.year
# Reorder columns
yearly_std_anom = yearly_std_anom[cols]
yearly_std_anom_FR = yearly_std_anom.copy()

my_df = yearly_std_anom.T
my_df['MEAN'] = df_clim_year.mean()
my_df['SD'] = df_clim_year.std()
my_df.rename(index={
    'Avalon_Channel' : 'Avalon Channel (AC)',
    'Bravo':'Bravo (BRA)',
    'Central_Labrador_Sea':'Cent. Lab. Sea (CLS)',
    'Flemish_Cap' : 'Flemish Cap (FC)',
    'Flemish_Pass' : 'Flemish Pass (FP)',
    'Greenland_Shelf' : 'Greenland Shelf (GS)',
    'Green-St._Pierre_Bank':'St.Pierre Bank (SPB)' ,
    'Hamilton_Bank':'Hamilton Bank (HB)',
    'Hudson_Strait': 'Hudson Strait (HS)',
    'Hybernia' : 'Hibernia (HIB)',
    'North_Central_Labrador_Sea':'NC Lab. Sea (NCLS)',
    'Northeast_Nfld_Shelf':'NE NF shelf (NENL)',
    'Orphan_Knoll':'Orphan Knoll (OK)',
    'St.Anthony_Basin':'St. Anthony B. (SAB)'
    }, inplace=True)

my_df.rename(columns={
    'MEAN' : r'$\rm \overline{x}$',
    'SD':u'sd'
    }, inplace=True)

year_list = yearly_std_anom.index.astype('str')
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

nrows, ncols = my_df.index.size+1, my_df.columns.size
fig=plt.figure(figsize=(ncols*wcell+wpad, nrows*hcell+hpad))
ax = fig.add_subplot(111)
ax.axis('off')
#do the table
header = ax.table(cellText=[['']],
                      colLabels=['-- Sea Surface temperature anomalies --'],
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
plt.savefig("scorecards_sst_yearly.png", dpi=300)
os.system('convert -trim scorecards_sst_yearly.png scorecards_sst_yearly.png')


# French table
my_df = yearly_std_anom_FR.T
my_df['MEAN'] = df_clim_year.mean()
my_df['SD'] = df_clim_year.std()
my_df.rename(index={
    'Avalon_Channel' : 'Chenal Avalon (AC)',
    'Bravo':'Bravo (BRA)',
    'Central_Labrador_Sea':'Mer Lab. Centrale (CLS)',
    'Flemish_Cap' : 'Bonnet Flammand (FC)',
    'Flemish_Pass' : 'Passe Flamande (FP)',
    'Greenland_Shelf' : 'Plateau Groenland (GS)',
    'Green-St._Pierre_Bank':'Banc St.Pierre (SPB)' ,
    'Hamilton_Bank':'Banc Hamilton (HB)',
    'Hudson_Strait': 'Detroit Hudson (HS)',
    'Hybernia' : 'Hibernia (HIB)',
    'North_Central_Labrador_Sea':'Mer Lab. Nord Cent. (NCLS)',
    'Northeast_Nfld_Shelf':'Plateau T-N NE (NENL)',
    'Orphan_Knoll':'Orphan Knoll (OK)',
    'St.Anthony_Basin':'Basin St. Anthony (SAB)'
    }, inplace=True)

year_list = yearly_std_anom.index.astype('str')
year_list = [i[2:4] for i in year_list] # 2-digit year
year_list.append(r'$\rm \overline{x}$') # add 2 extra columns
year_list.append(r'ET')

header = ax.table(cellText=[['']],
                      colLabels=['-- Anomalies de température de surface --'],
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
plt.savefig("scorecards_sst_yearly_FR.png", dpi=300)
os.system('convert -trim scorecards_sst_yearly_FR.png scorecards_sst_yearly_FR.png')


