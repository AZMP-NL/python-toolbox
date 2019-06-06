'''
AZMP script to read Headlands thermograph data from D. Senciall's .rpf format

data in : /home/cyrf0006/data/Headlands/
process in : /home/cyrf0006/AZMP/Headlands

Frederic.Cyr@dfo-mpo.gc.ca - February 2019

'''

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import datetime
import pandas as pd
import os
import hl_tools as hlt

# Adjust fontsize/weight
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 18}
plt.rc('font', **font)

# Generate the list
#infiles = 'comfortcove.list'
#os.system('ls ~/data/Headlands/Comfort_Cove/*.rpf > ' + infiles)
infiles = 'oldbonaventure.list'
os.system('ls ~/data/Headlands/Old_Bonaventure/*.rpf > ' + infiles)
infiles = 'arnoldscove.list'
os.system('ls ~/data/Headlands/ArnoldsCove/*.rpf > ' + infiles)
filelist = np.genfromtxt(infiles, dtype=str)
filelist = np.reshape(filelist, filelist.size) 

dfs = []
for fname in filelist:

    header = hlt.rpf_header(fname)      
    df = hlt.rpf_to_dataframe(fname)

    #df['depth'] = 
    
        
    print fname, df.max()
    dfs.append(df)
    
# concatenate all data    
df_all = pd.concat(dfs, axis=0)
df_all = df_all.sort_index()

# A check plot:
# Remove <1980
df_all = df_all[df_all.index.year>1980]
df_all.temperature.plot()
plt.show()


# monthly average
df_monthly = df_all.resample('M').mean()
df_monthly.plot()
plt.show()


#df_monthly.to_csv('comfort_cove_thermograph_1989-2017_monthly.csv')


# June-July only:
df_summer = df_monthly[(df_monthly.index.month>=6) & (df_monthly.index.month<=7)]
df_summer = df_summer.resample('As').mean()
df_summer.index = df_summer.index.year
#df_summer.to_csv('comfort_cove_thermograph_1989-2017_Jnue-July.csv')



## ---- plot summer data in anomalies ---- ##
anom = (df_summer - df_summer.mean()) / df_summer.std()
df1 = anom[anom<0]
df2 = anom[anom>0]
fig = plt.figure(4)
fig.clf()
width = .9
p1 = plt.bar(df1.index, np.squeeze(df1.values), width, alpha=0.8, color='steelblue')
p2 = plt.bar(df2.index, np.squeeze(df2.values), width, bottom=0, alpha=0.8, color='indianred')
plt.ylabel('Standardized Anomaly')
plt.title('Comfort Cove temperature (June-July)')
plt.grid()
fig.set_size_inches(w=15,h=9)
fig_name = 'anomalie_toberenamed.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)




## ---- annual curve ---- ##
df_all['woy'] = df_all.index.weekofyear
weekly_clim = df_all.groupby('woy').mean()
weekly_std = df_all.groupby('woy').std()
df_2018 = df_all[df_all.index.year>=2018]
weekly_2018 = df_2018.groupby('woy').mean()

fig = plt.figure(4)
fig.clf()
plt.plot(weekly_clim.index, weekly_clim.values, linewidth=3)
plt.plot(weekly_2018.index, weekly_2018.values, linewidth=3)
plt.fill_between(weekly_clim.index, np.squeeze(weekly_clim.values+weekly_std.values),np.squeeze(weekly_clim.values-weekly_std.values), facecolor='steelblue', interpolate=True , alpha=.3)
#plt.plot(weekly_clim.index, weekly_clim.values, linewidth=3)
#plt.plot(weekly_2018.index, weekly_2018.values, linewidth=3)
plt.ylabel(r'T ($^{\circ}C$)')
plt.xlabel('Week of the year')
plt.title('Comfort Cove temperature')
plt.xlim([0,53])
plt.ylim([-2,18])
plt.grid()
plt.legend(['1989-2018 average', '2018'])
fig.set_size_inches(w=15,h=9)
fig_name = 'annual_cycle_toberenamed.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)
