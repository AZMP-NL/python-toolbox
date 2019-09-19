'''
AZMP reporting - Iceberg counts from :
Iceberg data from USCG International Ice Patrol
Contact  Michael Hicks USCG
Hicks, Michael R CIV <Michael.R.Hicks@uscg.mil>
The Iceberg Season runs fron Oct 1 to Sept 30

data in /home/cyrf0006/data/AZMP/ColbourneStuff/NL_ICE_BERG_NUMBERS_DATA_1900_2017.xlsx

(script ran in /home/cyrf0006/AZMP/state_report/bergs)


Frederic.Cyr@dfo-mpo.gc.ca - June 2019

'''


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import seaborn as sns
 
# Adjust fontsize/weight
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 10}
plt.rc('font', **font)

clim_year = [1981, 2010]
current_year = 2018

## ----  Prepare the data ---- ##
# load from Excel sheets
df = pd.read_excel('/home/cyrf0006/data/AZMP/ColbourneStuff/NL_ICE_BERG_NUMBERS_DATA_1900_2017.xlsx', header=5, index_col='YEAR') 
df = df.drop(columns = 'TOT SEASON')

# Stack months under Years (pretty cool!)
#df = df.stack() 

# Transform to a series with values based the 15th of each month (had to convert years to string)
#df_BB.index = pd.to_datetime('15-' + df_BB.index.get_level_values(1) + '-' + df_BB.index.get_level_values(0).values.astype(np.str))

# Annual mean
df_annual = df.sum(axis=1)
df_annual.to_pickle('bergs_annual.pkl')
df_annual_clim = df_annual[(df_annual.index>=clim_year[0]) & (df_annual.index<=clim_year[1])]
df_annual_anom = df_annual - df_annual_clim.mean()
df_annual_std_anom = df_annual_anom/df_annual_clim.std()



# Monthly mean
df_monthly = df[df.index==current_year]
df_monthly_clim = df[(df.index>=clim_year[0]) & (df.index<=clim_year[1])]
df_monthly_std = df_monthly_clim.std(axis=0)
df_monthly_clim = df_monthly_clim.mean(axis=0)


## ---- plot monthly ---- ##
ind = np.arange(len(df_monthly.keys()))  # the x locations for the groups
width = 0.35  # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(ind - width/2, df_monthly_clim.values, width, yerr=df_monthly_std.values*.5,
                label='1981-2010')
rects2 = ax.bar(ind + width/2, np.squeeze(df_monthly.values), width, yerr=None,
                label=str(current_year))

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Counts')
#ax.set_title('Number of icebergs')
ax.set_xticks(ind)
ax.set_xticklabels(df_monthly_clim.index)
ax.legend()
ax.yaxis.grid() # horizontal lines
plt.ylim([0, 330])

fig.set_size_inches(w=6,h=3)
fig_name = 'bergs_monthly.png'
#plt.annotate('data source: www.ncdc.noaa.gov/teleconnections/', xy=(.58, .01), xycoords='figure fraction', annotation_clip=False, FontSize=12)
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)


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
