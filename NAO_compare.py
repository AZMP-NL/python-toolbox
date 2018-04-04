
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import datetime
import matplotlib.dates as mdates

### --- T-S --- ###
df_sal_season = pd.read_pickle('salinity_1948-2017.pkl')
meanS = df_sal_season.iloc[:,df_sal_season.columns<300].mean(axis=1)
df_all = pd.read_pickle('temp_summer_1948-2017.pkl')
CIL = df_all.min(axis=1)

### --- Nutrients --- ###

## BB only
df_BB_nut = pd.read_pickle('BB_nutrients.pkl')
# depth range
idx = ((df_BB_nut['depth']>0) & (df_BB_nut['depth']<=125))
df_BB_nut_surf = df_BB_nut[idx]

idx = ((df_BB_nut['depth']>200) & (df_BB_nut['depth']<=350))
df_BB_nut_btm = df_BB_nut[idx]

df_season_surf = df_BB_nut_surf.resample('4M').mean()
df_season_btm = df_BB_nut_btm.resample('4M').mean()

surf = df_season_surf['NO3']/df_season_surf['PO4']
btm = df_season_btm['NO3']/df_season_btm['PO4']
btm_SIO = df_season_btm['SIO']
btm_PO4 = df_season_btm['PO4']

## All stations
df = pd.read_excel('~/github/AZMP-NL/data/AZMP_Nutrients_1999_2016_Good_Flags_Only.xlsx')
# Set date as index
df = df.set_index('sas_date')
# Drop other time-related columns
df = df.drop(['Day', 'Month', 'Year'], axis=1)

idx = ((df['depth']>125) & (df['depth']<=350))
df = df[idx]

R1 = df['NO3']/df['PO4']
R1 = R1.replace([np.inf, -np.inf], np.nan)
R1 = R1.resample('AS').mean()   
R2 = df['NO3']/df['SIO']
R2 = R2.replace([np.inf, -np.inf], np.nan)
R2 = R2.resample('AS').mean()   
R1.drop(R1.index[R1.index.year==2015], inplace=True)
R2.drop(R2.index[R2.index.year==2015], inplace=True)



### --- NAO index --- ###
# Load data
nao_file = '/home/cyrf0006/research/AZMP_stateReports/data/indices/data.csv'
df = pd.read_csv(nao_file, header=1)

# Set index
df = df.set_index('Date')
df.index = pd.to_datetime(df.index, format='%Y%m')

# Select only DJF
df_winter = df[(df.index.month==12) | (df.index.month==1) | (df.index.month==2)]

# Start Dec-1950
df_winter = df_winter[df_winter.index>pd.to_datetime('1950-10-01')]
df_winter.resample('3M').mean()

# Average 3 consecutive values (DJF average); We loose index.
df_winter = df_winter.groupby(np.arange(len(df_winter))//3).mean()

# Reset index using years only
df_winter.index = pd.unique(df.index.year)[1:,]
df_winter.index = pd.to_datetime(df_winter.index, format='%Y')

### --- AOO index --- ###
# Load data
aoo_file = '/home/cyrf0006/Data/AZMP/Indices/AOO.xlsx'
df_aoo = pd.read_excel(aoo_file)

# Set index
df_aoo = df_aoo.set_index('Year')
df_aoo.index = pd.to_datetime(df_aoo.index)

### --- AMO index --- ###
# Load data
amo_file = '/home/cyrf0006/Data/AZMP/Indices/AMO_index.txt'
df_amo = pd.read_csv(amo_file, header=None, delimiter=r"\s+")
df_amo = df_amo.rename(columns={ df_amo.columns[0]: "Year" })

df_amo = df_amo.set_index('Year')
#df_amo.index = pd.to_datetime(df_amo.index)

# Stack months under Years (pretty cool!)
df_amo = df_amo.stack() 

# Transform to a series with values based the 15th of each month
df_amo.index = pd.to_datetime({'year':df_amo.index.get_level_values(0), 'month':df_amo.index.get_level_values(1), 'day':15})
df_amo = df_amo[df_amo>-10] # remove -99.99 data


### --- plot --- ### (Physics)
fig = plt.figure(1)
fig.clf()
ax = fig.add_subplot(111)

Sanomaly = (meanS-meanS.mean())/meanS.std()
ratio = (surf-surf.mean())/surf.std()
CIL_anom = (CIL-CIL.mean())/CIL.std()
AOO_anom = (df_aoo-df_aoo.mean())/df_aoo.std()
AMO_anom = (df_amo-df_amo.mean())/df_amo.std()
NAO_anom = (df_winter-df_winter.mean())/df_winter.std()

#plt.plot(df_winter.index, df_winter.rolling(window=12, center=True).mean())
plt.plot(NAO_anom.index, NAO_anom.rolling(window=12, center=True).mean(), 'k')
#plt.plot(NAO_anom.index, NAO_anom.rolling(window=12, center=True).mean(), '--')
#plt.plot(AOO_anom.index, AOO_anom.rolling(window=5, center=True).mean())
plt.plot(AMO_anom.index, AMO_anom.rolling(window=60, center=True).mean(), '--k')
#plt.plot(df_amo.index, df_amo.rolling(window=60, center=True).mean())
plt.plot(CIL.index, CIL_anom.rolling(window=5, center=True).mean())
plt.plot(meanS.index, Sanomaly.rolling(window=12, center=True).mean(), 'r')
#plt.plot(df_season_surf.index, ratio.rolling(window=3, center=True).mean())
plt.xlim([pd.to_datetime('1948'), pd.to_datetime('2017')])
plt.ylim([-2,2])
plt.grid()
#plt.legend(['NAO_w','AOO','AMO','CIL','S'])
plt.legend([r'$\rm NAO_{winter}$','AMO','CIL','S'])
plt.xlabel('Year')
plt.ylabel(r'Standardized anom.')
plt.title('multi-index comparison')
fig.set_size_inches(w=10,h=6)
fig_name = 'multi_index_physics.png'
fig.set_dpi(300)
fig.savefig(fig_name)


### --- plot --- ### (Biochem 1)
fig = plt.figure(1)
fig.clf()
ax = fig.add_subplot(111)

Sanomaly = (meanS-meanS.mean())/meanS.std()
#ratio = (surf-surf.mean())/surf.std()
ratio = (btm-btm.mean())/btm.std()
CIL_anom = (CIL-CIL.mean())/CIL.std()
SIO = (btm_SIO-btm_SIO.mean())/btm_SIO.std()
PO4 = (btm_PO4-btm_PO4.mean())/btm_PO4.std()


plt.plot(NAO_anom.index, NAO_anom.rolling(window=3, center=True).mean(), 'k')
#plt.plot(meanS.index, Sanomaly.rolling(window=3, center=True).mean())
plt.plot(AOO_anom.index, AOO_anom.rolling(window=2, center=True).mean(), '--k')
#plt.plot(CIL.index, CIL_anom.rolling(window=5, center=True).mean())
plt.plot(df_season_btm.index, ratio.rolling(window=3, center=True).mean())
plt.plot(df_season_btm.index, SIO.rolling(window=6, center=True).mean())
#plt.plot(df_season_btm.index, PO4.rolling(window=3, center=True).mean())
plt.xlim([pd.to_datetime('1999'), pd.to_datetime('2017')])
plt.ylim([-2,2])
plt.grid()
plt.legend(['NAO','AOO',r'$NO3/PO4_{btm}$',r'$SiO_{btm}$'])
plt.xlabel('Year')
plt.ylabel(r'Standardized anom.')
plt.title('multi-index comparison')
fig.set_size_inches(w=10,h=6)
fig_name = 'multi_index.png'
fig.set_dpi(300)
fig.savefig(fig_name)




### --- plot --- ### (Biochem 2 - Ratios)
fig = plt.figure(1)
fig.clf()
ax = fig.add_subplot(111)

ratio1 = (R1-R1.mean())/R1.std()
ratio2 = (R2-R2.mean())/R2.std()

#plt.bar(ratio1.index, ratio1.rolling(window=1, center=True).mean(),width=w,color='b',align='center')
#plt.bar(ratio2.index, ratio2.rolling(window=1, center=True).mean(),width=w,color='g',align='center')
#plt.plot(ratio1.index, ratio1.rolling(window=1, center=True).mean())
#plt.plot(ratio2.index, ratio2.rolling(window=1, center=True).mean())
plt.plot(ratio1.index, ratio1)
plt.plot(ratio2.index, ratio2)
plt.plot(NAO_anom.index, NAO_anom.rolling(window=3, center=True).mean()*0+100, 'k')
plt.plot(AOO_anom.index, AOO_anom.rolling(window=3, center=True).mean(), '--k')

#plt.plot(R1.index, R1)
#plt.plot(R2.index, R2)
plt.xlim([pd.to_datetime('1999'), pd.to_datetime('2017')])
plt.ylim([-2,3])
plt.grid()
#plt.legend(['NAO','AOO',r'[$PO_4/NO_3]$',r'[$SiO_3/NO_3$]'])
plt.legend([r'$[NO_3]/[PO_4]$',r'$[NO_3]/[SiO_3]$','NAO','AOO'])
plt.xlabel('Year')
plt.ylabel(r'Standardized anom.')
plt.title('multi-index comparison')
fig.set_size_inches(w=10,h=6)
fig_name = 'multi_index_ratiosBTM_AOO.png'
fig.set_dpi(300)
fig.savefig(fig_name)



### --- plot --- ### (Index Only for Ali)
fig = plt.figure(1)
fig.clf()
ax = fig.add_subplot(111)

plt.plot(NAO_anom.index, NAO_anom.rolling(window=3, center=True).mean())
plt.plot(AOO_anom.index, AOO_anom.rolling(window=3, center=True).mean())
plt.plot(AMO_anom.index, AMO_anom.rolling(window=3, center=True).mean())
plt.xlim([pd.to_datetime('1993'), pd.to_datetime('2015')])
plt.ylim([-2,2])
plt.grid()
#plt.legend(['NAO','AOO',r'[$PO_4/NO_3]$',r'[$SiO_3/NO_3$]'])
plt.legend(['NAO','AOO','AMO'])
plt.xlabel('Year')
plt.ylabel(r'Standardized anom.')
plt.title('multi-index comparison')
fig.set_size_inches(w=10,h=6)
fig_name = 'Index_all.png'
fig.set_dpi(300)
fig.savefig(fig_name)

