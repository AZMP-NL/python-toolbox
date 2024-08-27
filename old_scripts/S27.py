'''
To explore relationship between Stn 27 stratification and AOO, etc.


Frederic.Cyr@dfo-mpo.gc.ca - March 2018

'''


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import datetime
import matplotlib.dates as mdates

# Adjust fontsize/weight
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 18}
plt.rc('font', **font)


### ---- S27 data ---- ###
# Load data
S27_file = '/home/cyrf0006/research/AZMP_stateReports/ColbourneStuff/S27_Integrated.dat'
AO_file = '/home/cyrf0006/research/AZMP_stateReports/ColbourneStuff/ARCTIC_OSCILLATION.xlsx'

# Read S27 and drop line 2 of file
df = pd.read_csv(S27_file, delimiter=r"\s+", header=0)
df = df.drop(df.index[0])

# Set index
df = df.set_index('Year')
df.index = pd.to_datetime(df.index, format='%Y')
df_27 = df


### ---- AO data ---- ###
df = pd.read_excel(AO_file, header=5)
df = df.set_index(df.columns[0])
df.index = df.index.set_names('Year')
df = df.drop(df.index[-3:,])
df.index = np.int64(df.index)
df.index = pd.to_datetime(df.index, format='%Y')


### ---- plot compa ---- ###
plt.plot(df['DEC-FEB'])
plt.plot(df_27.index, np.float64(df_27['Strat'])-np.float64(df_27['Strat']['1981':'2010']).mean())
plt.show()
