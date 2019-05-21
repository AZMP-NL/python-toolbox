'''
A script to read and plot data from Ice EXtent Excel sheet south of 55N on NL shelves
'''
# A first test to read Excel nutrient file and export to Pandas.

# Check in:
#  /home/cyrf0006/research/AZMP_database/biochem

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import datetime
from scipy import stats

# Adjust fontsize/weight
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 18}
plt.rc('font', **font)

## ----  Load ice data ---- ##
df = pd.read_excel('/home/cyrf0006/AZMP/state_reports/ColbourneStuff/NL_SEA_ICE_EXTENT_DATA_1963_2017.xlsx', header=4)
df_time = df[['YEAR', 'MONTH']]
df_time['DAYS'] = pd.Series(np.repeat(15, df.index.size))
df.index = pd.to_datetime(df_time)
df = df['km2']
df = df/1000


df_annual = df.resample('As').mean()

## ---- plot ice extent (for Max) ---- ##
fig = plt.figure(1)
plt.clf()
plt.plot(df.index, df.values, 'k-', linewidth=1)
#plt.plot(df.index, df.rolling(24,center=True).mean(), 'r-', linewidth=3)
plt.plot(df_annual.index, df_annual.values, 'r-', linewidth=3)
plt.ylabel(r'$\rm \times 1000\,km^2$', fontsize=15, fontweight='bold')
plt.xlabel('Year', fontsize=15, fontweight='bold')
plt.title(r'Ice Extent south of 55$^{\circ}$N', fontsize=15, fontweight='bold')
plt.xlim([pd.Timestamp('1970-01-01'), pd.Timestamp('2017-01-01')])
#plt.ylim([-1, 3])
plt.xticks(pd.date_range('1970-01-01', periods=6, freq='10Y'))
plt.legend(['monthly', 'annual'])
plt.grid('on')


fig.set_size_inches(w=9,h=6)
fig_name = 'Ice_extent_NL_1970-2017.png'
fig.set_dpi(300)
fig.savefig(fig_name)

