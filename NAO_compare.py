
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import datetime
import matplotlib.dates as mdates

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 14}

plt.rc('font', **font)

### --- T-S --- ###
df_sal_season = pd.read_pickle('/home/cyrf0006/AZMP/dig_project/salinity_1948-2017.pkl')
meanS = df_sal_season.iloc[:,df_sal_season.columns<300].mean(axis=1)
meanS2 = df_sal_season.iloc[:,(df_sal_season.columns<=200) & (df_sal_season.columns>=100)].mean(axis=1)
meanS1 = df_sal_season.iloc[:,(df_sal_season.columns<=100) & (df_sal_season.columns>=0)].mean(axis=1)
meanS[meanS<20]=np.nan
meanS1[meanS1<20]=np.nan
meanS2[meanS2<20]=np.nan
df_all = pd.read_pickle('/home/cyrf0006/AZMP/dig_project/temp_summer_1948-2017.pkl')
CIL = df_all.min(axis=1)

df_temp_surf = pd.read_pickle('/home/cyrf0006/AZMP/dig_project/historical_SST_25m_1948-2017.pkl')
SST = df_temp_surf.iloc[:,df_temp_surf.columns<10].mean(axis=1)
#SST = SST[SST.index.year>1995]
SST_spring = SST[(SST.index.month>=4) & (SST.index.month<=8)]
SST = SST.resample('AS').mean()
SST_spring = SST_spring.resample('AS').mean()
#HERRE!!!
### --- Nutrients --- ###

## BB only
df_BB_nut = pd.read_pickle('/home/cyrf0006/AZMP/database/BB_nutrients.pkl')
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
df = pd.read_excel('/home/cyrf0006/github/AZMP-NL/data/AZMP_Nutrients_1999_2016_Good_Flags_Only.xlsx')
# Set date as index
df = df.set_index('sas_date')
# Drop other time-related columns
df = df.drop(['Day', 'Month', 'Year'], axis=1)

#idx = ((df['depth']>125) & (df['depth']<=350))

## Nutrients 50-150m ##
idx = ((df['depth']>50) & (df['depth']<=150))
df = df[idx]

R1 = df['NO3']/df['PO4']
R1 = R1.replace([np.inf, -np.inf], np.nan)
R1 = R1.resample('AS').mean()   
R2 = df['NO3']/df['SIO']
R2 = R2.replace([np.inf, -np.inf], np.nan)
R2 = R2.resample('AS').mean()   
R1.drop(R1.index[R1.index.year==2015], inplace=True)
R2.drop(R2.index[R2.index.year==2015], inplace=True)

# nutrient timeseries
NO3_deep = df['NO3'].resample('AS').mean()
SIO_deep = df['SIO'].resample('AS').mean()
PO4_deep = df['PO4'].resample('AS').mean()
O2_deep = df['oxygen'].resample('AS').mean()


### --- S27 stratification --- ###
s27_file = '/home/cyrf0006/AZMP/state_reports/ColbourneStuff/S27_Integrated.dat'
df_s27 = pd.read_csv(s27_file, header=1, delimiter=r"\s+")
df_s27 = df_s27.rename(columns={'----':'Year', '0-50m':'T', '175m':'botT', '0-50m.1':'S', '0-50m.2':'Strat' })
df_s27 = df_s27.set_index('Year')
df_s27.index = pd.to_datetime(df_s27.index, format='%Y')
df_s27_strat = df_s27['Strat']


### --- NAO index --- ###
# Load data
nao_file = '/home/cyrf0006/data/AZMP/indices/data.csv'
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
df_winter = df_winter[(df_winter.index.year>1950) & (df_winter.index.year<2018)]


### --- AOO index --- ###
# Load data
aoo_file = '/home/cyrf0006/data/AZMP/indices/AOO.xlsx'
df_aoo = pd.read_excel(aoo_file)

# Set index
df_aoo = df_aoo.set_index('Year')
df_aoo.index = pd.to_datetime(df_aoo.index)

### --- AMO index --- ###
# Load data
amo_file = '/home/cyrf0006/data/AZMP/indices/amon.us.data'
#amo_file = '/home/cyrf0006/data/AZMP/indices/AMO_index.txt'names=['ID','CODE']
#col_names = ['Year','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
col_names = ['Year','01','02','03','04','05','06','07','08','09','10','11','12']
#df_amo = pd.read_csv(amo_file, header=1, delimiter=r"\s+", index_col=[0])
df_amo = pd.read_csv(amo_file, header=1, delimiter=r"\s+", names=col_names)
df_amo = df_amo.drop(df_amo.tail(4).index) # drop NOAA footer
#df_amo = df_amo.rename(columns={ df_amo.columns[0]: "Year" })

df_amo = df_amo.set_index('Year')
#df_amo.index = pd.to_datetime(df_amo.index)

# Stack months under Years (pretty cool!)
df_amo = df_amo.stack() 

# Transform to a series with values based the 15th of each month
df_amo.index = pd.to_datetime({'year':df_amo.index.get_level_values(0), 'month':df_amo.index.get_level_values(1), 'day':15})
df_amo = df_amo[(df_amo.index.year>1950) & (df_amo.index.year<2018)]

#df_amo = df_amo[df_amo>-10] # remove -99.99 data
df_amo = pd.to_numeric(df_amo) # default dtype is object
df_amo = df_amo.resample('As').mean() 


## ### --- AMO index --- ###
## # Load data
## amo_file = '/home/cyrf0006/data/AZMP/indices/amon.us.data'
## df_amo = pd.read_csv(amo_file, header=None, delimiter=r"\s+")
## df_amo = df_amo.rename(columns={ df_amo.columns[0]: "Year" })

## df_amo = df_amo.set_index('Year')
## #df_amo.index = pd.to_datetime(df_amo.index)

## # Stack months under Years (pretty cool!)
## df_amo = df_amo.stack() 

## # Transform to a series with values based the 15th of each month
## df_amo.index = pd.to_datetime({'year':df_amo.index.get_level_values(0), 'month':df_amo.index.get_level_values(1), 'day':15})
## df_amo = df_amo[df_amo>-10] # remove -99.99 data


### --- plot --- ### (Physics)
Sanomaly = (meanS-meanS.mean())/meanS.std()
Sanomaly1 = (meanS1-meanS1.mean())/meanS1.std()
Sanomaly2 = (meanS2-meanS2.mean())/meanS2.std()
ratio = (surf-surf.mean())/surf.std()
CIL_anom = (CIL-CIL.mean())/CIL.std()
AOO_anom = (df_aoo-df_aoo.mean())/df_aoo.std()
AMO_anom = (df_amo-df_amo.mean())/df_amo.std()
NAO_anom = (df_winter-df_winter.mean())/df_winter.std()
SST_anom = (SST-SST.mean())/SST.std()
SST_anom_spring = (SST_spring-SST_spring.mean())/SST_spring.std()
S27_anom = (df_s27_strat-df_s27_strat.mean())/df_s27_strat.std()



fig = plt.figure(1)
fig.clf()
ax = fig.add_subplot(111)
#plt.plot(df_winter.index, df_winter.rolling(window=12, center=True).mean())
plt.plot(NAO_anom.index, NAO_anom.rolling(window=12, center=True).mean(), 'k')
plt.plot(AMO_anom.index, AMO_anom.rolling(window=60, center=True).mean(), '--k')
plt.plot(AOO_anom.index, AOO_anom.rolling(window=5, center=True).mean(), color='gray')
plt.plot(CIL.index, CIL_anom.rolling(window=5, center=True).mean())
plt.plot(meanS.index, Sanomaly.rolling(window=12, center=True).mean(), 'r')
#plt.plot(df_season_surf.index, ratio.rolling(window=3, center=True).mean())
plt.xlim([pd.to_datetime('1948'), pd.to_datetime('2017')])
plt.xticks(pd.date_range('1950-01-01', periods=7, freq='10Y'))
plt.ylim([-2,2])
plt.grid()
#plt.legend(['NAO_w','AOO','AMO','CIL','S'])
plt.legend([r'$\rm NAO_{winter}$','AMO','AOO','CIL','S'])
plt.xlabel('Year', fontweight='bold')
plt.ylabel(r'Standardized anom.', fontweight='bold')
plt.title('multi-index comparison')
fig.set_size_inches(w=10,h=6)
fig_name = 'multi_index_physics.png'
fig.set_dpi(300)
fig.savefig(fig_name)

# T-S vs NAO
fig = plt.figure(1)
fig.clf()
ax = fig.add_subplot(111)
plt.plot(CIL.index, CIL_anom.rolling(window=5, center=True).mean(), linewidth=3)
plt.plot(meanS.index, Sanomaly.rolling(window=12, center=True).mean(), 'r', linewidth=3)
plt.plot(NAO_anom.index, NAO_anom.rolling(window=12, center=True).mean(), '--k', linewidth=5)
#plt.plot(df_season_surf.index, ratio.rolling(window=3, center=True).mean())
plt.xlim([pd.to_datetime('1948'), pd.to_datetime('2019')])
plt.xticks(pd.date_range('1950-01-01', periods=7, freq='10Y'))
plt.ylim([-2,2])
plt.grid()
plt.legend(['CIL','S', r'$\rm NAO_{winter}$'], loc=2)
plt.xlabel('Year', fontweight='bold')
plt.ylabel(r'Standardized anom.', fontweight='bold')
plt.title('multi-index comparison')

plt.annotate("", xy=(pd.to_datetime('2017-01-01'), 1.75), xytext=(pd.to_datetime('2017-01-01'), .25), arrowprops=dict(arrowstyle="->"))
plt.annotate("", xy=(pd.to_datetime('2017-01-01'), -1.75), xytext=(pd.to_datetime('2017-01-01'), -.25), arrowprops=dict(arrowstyle="->", facecolor='black'))
plt.text(pd.to_datetime('2017-01-01'), 1.75, 'colder winter', horizontalalignment='right', verticalalignment='bottom')
plt.text(pd.to_datetime('2017-01-01'), -1.75, 'milder winter', horizontalalignment='right', verticalalignment='top')

fig.set_size_inches(w=10,h=6)
fig_name = 'TS-vs-NAO.png'
fig.set_dpi(300)
fig.savefig(fig_name)

# T-S vs AOO
fig = plt.figure(1)
fig.clf()
ax = fig.add_subplot(111)
plt.plot(CIL.index, CIL_anom.rolling(window=5, center=True).mean(), linewidth=3)
plt.plot(meanS.index, Sanomaly.rolling(window=12, center=True).mean(), 'r', linewidth=3)
plt.plot(AOO_anom.index, AOO_anom.rolling(window=5, center=True).mean(), '--k', linewidth=5)
#plt.plot(df_season_surf.index, ratio.rolling(window=3, center=True).mean())
plt.xlim([pd.to_datetime('1948'), pd.to_datetime('2019')])
plt.xticks(pd.date_range('1950-01-01', periods=7, freq='10Y'))
plt.ylim([-2,2])
plt.grid()
plt.legend(['CIL','S', 'AOO'], loc=2)
plt.xlabel('Year', fontweight='bold')
plt.ylabel(r'Standardized anom.', fontweight='bold')
plt.title('multi-index comparison')

plt.annotate("", xy=(pd.to_datetime('2017-01-01'), 1.75), xytext=(pd.to_datetime('2017-01-01'), .25), arrowprops=dict(arrowstyle="->"))
plt.annotate("", xy=(pd.to_datetime('2017-01-01'), -1.75), xytext=(pd.to_datetime('2017-01-01'), -.25), arrowprops=dict(arrowstyle="->", facecolor='black'))
plt.text(pd.to_datetime('2017-01-01'), 1.75, 'less Arctic export', horizontalalignment='right', verticalalignment='bottom')
plt.text(pd.to_datetime('2017-01-01'), -1.75, 'more Arctic export', horizontalalignment='right', verticalalignment='top')

fig.set_size_inches(w=10,h=6)
fig_name = 'TS-vs-AOO.png'
fig.set_dpi(300)
fig.savefig(fig_name)

# T-S vs AMO
fig = plt.figure(1)
fig.clf()
ax = fig.add_subplot(111)
plt.plot(CIL.index, CIL_anom.rolling(window=5, center=True).mean(), linewidth=3)
plt.plot(meanS.index, Sanomaly.rolling(window=12, center=True).mean(), 'r', linewidth=3)
plt.plot(AMO_anom.index, AMO_anom.rolling(window=5, center=True).mean(), '--k', linewidth=5)
#plt.plot(S27_anom.index, S27_anom.rolling(window=5, center=True).mean(), 'm', linewidth=3)
#plt.plot(df_season_surf.index, ratio.rolling(window=3, center=True).mean())
plt.xlim([pd.to_datetime('1948'), pd.to_datetime('2019')])
plt.xticks(pd.date_range('1950-01-01', periods=7, freq='10Y'))
plt.ylim([-2,2])
plt.grid()
plt.legend(['CIL','S', 'AMO'], loc=2)
plt.xlabel('Year', fontweight='bold')
plt.ylabel(r'Standardized anom.', fontweight='bold')
plt.title('multi-index comparison')

plt.annotate("", xy=(pd.to_datetime('2017-01-01'), 1.75), xytext=(pd.to_datetime('2017-01-01'), .25), arrowprops=dict(arrowstyle="->"))
plt.annotate("", xy=(pd.to_datetime('2017-01-01'), -1.75), xytext=(pd.to_datetime('2017-01-01'), -.25), arrowprops=dict(arrowstyle="->", facecolor='black'))
plt.text(pd.to_datetime('2017-01-01'), 1.75, 'warmer SST / stronger MOC', horizontalalignment='right', verticalalignment='bottom')
plt.text(pd.to_datetime('2017-01-01'), -1.75, 'cooler SST / weaker MOC', horizontalalignment='right', verticalalignment='top')

fig.set_size_inches(w=10,h=6)
fig_name = 'TS-vs-AMO.png'
fig.set_dpi(300)
fig.savefig(fig_name)

# SST vs AMO
fig = plt.figure(1)
fig.clf()
ax = fig.add_subplot(111)
plt.plot(SST_spring.index, SST_anom_spring.rolling(window=3, center=True).mean(), linewidth=3)
plt.plot(AMO_anom.index, AMO_anom.rolling(window=5, center=True).mean(), '--k', linewidth=5)
#plt.plot(df_season_surf.index, ratio.rolling(window=3, center=True).mean())
plt.xlim([pd.to_datetime('1948'), pd.to_datetime('2019')])
plt.xticks(pd.date_range('1950-01-01', periods=7, freq='10Y'))
plt.ylim([-2,2])
plt.grid()
plt.legend(['SST', 'AMO'], loc=2)
plt.xlabel('Year', fontweight='bold')
plt.ylabel(r'Standardized anom.', fontweight='bold')
plt.title('multi-index comparison')

plt.annotate("", xy=(pd.to_datetime('2017-01-01'), 1.75), xytext=(pd.to_datetime('2017-01-01'), .25), arrowprops=dict(arrowstyle="->"))
plt.annotate("", xy=(pd.to_datetime('2017-01-01'), -1.75), xytext=(pd.to_datetime('2017-01-01'), -.25), arrowprops=dict(arrowstyle="->", facecolor='black'))
plt.text(pd.to_datetime('2017-01-01'), 1.75, 'warmer SST / stronger MOC', horizontalalignment='right', verticalalignment='bottom')
plt.text(pd.to_datetime('2017-01-01'), -1.75, 'cooler SST / weaker MOC', horizontalalignment='right', verticalalignment='top')

fig.set_size_inches(w=10,h=6)
fig_name = 'SST-vs-AMO.png'
fig.set_dpi(300)
fig.savefig(fig_name)

# SST vs AOO
fig = plt.figure(1)
fig.clf()
ax = fig.add_subplot(111)
plt.plot(SST_spring.index, SST_anom_spring.rolling(window=3, center=True).mean(), linewidth=3)
#plt.plot(meanS.index, Sanomaly.rolling(window=12, center=True).mean(), 'r', linewidth=3)
plt.plot(AOO_anom.index, AOO_anom.rolling(window=5, center=True).mean(), '--k', linewidth=5)
#plt.plot(df_season_surf.index, ratio.rolling(window=3, center=True).mean())
plt.xlim([pd.to_datetime('1948'), pd.to_datetime('2019')])
plt.xticks(pd.date_range('1950-01-01', periods=7, freq='10Y'))
plt.ylim([-2,2])
plt.grid()
#plt.legend(['SST', r'$\rm S_{0-200}$', 'AOO'], loc=2)
plt.legend(['SST', 'AOO'], loc=2)
plt.xlabel('Year', fontweight='bold')
plt.ylabel(r'Standardized anom.', fontweight='bold')
plt.title('multi-index comparison')

plt.annotate("", xy=(pd.to_datetime('2017-01-01'), 1.75), xytext=(pd.to_datetime('2017-01-01'), .25), arrowprops=dict(arrowstyle="->"))
plt.annotate("", xy=(pd.to_datetime('2017-01-01'), -1.75), xytext=(pd.to_datetime('2017-01-01'), -.25), arrowprops=dict(arrowstyle="->", facecolor='black'))
plt.text(pd.to_datetime('2017-01-01'), 1.75, 'less Arctic export', horizontalalignment='right', verticalalignment='bottom')
plt.text(pd.to_datetime('2017-01-01'), -1.75, 'more Arctic export', horizontalalignment='right', verticalalignment='top')

fig.set_size_inches(w=10,h=6)
fig_name = 'SST-vs-AOO.png'
fig.set_dpi(300)
fig.savefig(fig_name)

# S, N2 vs AOO
fig = plt.figure(1)
fig.clf()
ax = fig.add_subplot(111)
plt.plot(meanS.index, Sanomaly.rolling(window=12, center=True).mean(), 'r', linewidth=3)
plt.plot(S27_anom.index, S27_anom.rolling(window=5, center=True).mean(), 'm', linewidth=3)
#plt.plot(meanS.index, Sanomaly.rolling(window=12, center=True).mean(), 'r', linewidth=3)
plt.plot(AOO_anom.index, AOO_anom.rolling(window=5, center=True).mean(), '--k', linewidth=5)
#plt.plot(df_season_surf.index, ratio.rolling(window=3, center=True).mean())
plt.xlim([pd.to_datetime('1948'), pd.to_datetime('2019')])
plt.xticks(pd.date_range('1950-01-01', periods=7, freq='10Y'))
plt.ylim([-2,2])
plt.grid()
#plt.legend(['SST', r'$\rm S_{0-200}$', 'AOO'], loc=2)
plt.legend([r'$S_{0-200m}$', r'$Strat_{0-50m}$', 'AOO'], loc=2)
plt.xlabel('Year', fontweight='bold')
plt.ylabel(r'Standardized anom.', fontweight='bold')
plt.title('multi-index comparison')

plt.annotate("", xy=(pd.to_datetime('2017-01-01'), 1.75), xytext=(pd.to_datetime('2017-01-01'), .25), arrowprops=dict(arrowstyle="->"))
plt.annotate("", xy=(pd.to_datetime('2017-01-01'), -1.75), xytext=(pd.to_datetime('2017-01-01'), -.25), arrowprops=dict(arrowstyle="->", facecolor='black'))
plt.text(pd.to_datetime('2017-01-01'), 1.75, 'less Arctic export', horizontalalignment='right', verticalalignment='bottom')
plt.text(pd.to_datetime('2017-01-01'), -1.75, 'more Arctic export', horizontalalignment='right', verticalalignment='top')

fig.set_size_inches(w=10,h=6)
fig_name = 'N2-vs-AOO.png'
fig.set_dpi(300)
fig.savefig(fig_name)



### --- plot --- ### (Biochem 1)
fig = plt.figure(1)
fig.clf()
ax = fig.add_subplot(111)

#Sanomaly = (meanS-meanS.mean())/meanS.std()
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


### ---- Nutrients ---- ###
SIO_deep_anom = (SIO_deep-SIO_deep.mean())/SIO_deep.std()
NO3_deep_anom = (NO3_deep-NO3_deep.mean())/NO3_deep.std()
PO4_deep_anom = (PO4_deep-PO4_deep.mean())/PO4_deep.std()
O2_deep_anom = (O2_deep-O2_deep.mean())/O2_deep.std()
df_aoo_nut = df_aoo[df_aoo.index.year>1997]
df_amo_nut = df_amo[df_amo.index.year>1997]
AOO_anom_nut = (df_aoo_nut-df_aoo_nut.mean())/df_aoo_nut.std()
AMO_anom_nut = (df_amo_nut-df_amo_nut.mean())/df_amo_nut.std()

meanS_nut = meanS[meanS.index.year>1997]
Sanomaly_nut = (meanS_nut-meanS_nut.mean())/meanS_nut.std()

fig = plt.figure(4)
fig.clf()
ax = fig.add_subplot(111)
plt.plot(NO3_deep_anom.index, NO3_deep_anom.rolling(window=3, center=True).mean(), linewidth=3)
plt.plot(PO4_deep_anom.index, PO4_deep_anom.rolling(window=3, center=True).mean(), linewidth=3)
plt.plot(SIO_deep_anom.index, SIO_deep_anom.rolling(window=3, center=True).mean(), linewidth=3)
plt.ylim([-1.75,1.75])
plt.grid()
plt.xlim([pd.to_datetime('1998'), pd.to_datetime('2017')])
plt.legend(['NO3','PO4','SiO'], loc=1)
plt.xlabel('Year', fontsize=16, fontweight='bold')
plt.ylabel(r'Standardized anom.', fontsize=16, fontweight='bold')
fig.set_size_inches(w=10,h=6)
fig_name = 'nutrients_alone.png'
fig.set_dpi(300)
fig.savefig(fig_name)


fig = plt.figure(4)
fig.clf()
ax = fig.add_subplot(111)
plt.plot(Sanomaly_nut.index, Sanomaly_nut.rolling(window=12, center=True).mean(), linewidth=3)
#plt.plot(Sanomaly.index, Sanomaly.rolling(window=12, center=True).mean())
plt.plot(NO3_deep_anom.index, NO3_deep_anom.rolling(window=3, center=True).mean(), linewidth=3)
plt.plot(SIO_deep_anom.index, SIO_deep_anom.rolling(window=3, center=True).mean(), linewidth=3)
plt.plot(AOO_anom_nut.index, AOO_anom_nut.rolling(window=2, center=True).mean(), '--k', linewidth=3)
plt.plot(S27_anom.index, S27_anom.rolling(window=2, center=True).mean(), '--m', linewidth=3)
#plt.plot(O2_deep_anom.index, O2_deep_anom.rolling(window=3, center=True).mean())
plt.ylim([-1.75,1.75])
plt.grid()
plt.xlim([pd.to_datetime('1998'), pd.to_datetime('2017')])
plt.legend(['S','NO3','SiO','AOO', 'Strat'], loc=2)
plt.xlabel('Year', fontsize=16, fontweight='bold')
plt.ylabel(r'Standardized anom.', fontsize=16, fontweight='bold')

plt.annotate("", xy=(pd.to_datetime('2016-01-01'), 1.45), xytext=(pd.to_datetime('2016-01-01'), .25), arrowprops=dict(arrowstyle="->"))
plt.annotate("", xy=(pd.to_datetime('2016-01-01'), -1.45), xytext=(pd.to_datetime('2016-01-01'), -.25), arrowprops=dict(arrowstyle="->", facecolor='black'))
plt.text(pd.to_datetime('2016-01-01'), 1.5, 'less Arctic export', horizontalalignment='right', verticalalignment='bottom')
plt.text(pd.to_datetime('2016-01-01'), -1.5, 'more Arctic export', horizontalalignment='right', verticalalignment='top')
#plt.title('multi-index comparison')
fig.set_size_inches(w=10,h=6)
fig_name = 'nutrients_AOO.png'
fig.set_dpi(300)
fig.savefig(fig_name)


fig = plt.figure(4)
fig.clf()
ax = fig.add_subplot(111)
plt.plot(CIL_anom.index, CIL_anom.rolling(window=3, center=True).mean(), linewidth=3)
plt.plot(Sanomaly_nut.index, Sanomaly_nut.rolling(window=12, center=True).mean(), linewidth=3)
#plt.plot(Sanomaly.index, Sanomaly.rolling(window=12, center=True).mean())
plt.plot(NO3_deep_anom.index, NO3_deep_anom.rolling(window=3, center=True).mean(), linewidth=3)
plt.plot(SIO_deep_anom.index, SIO_deep_anom.rolling(window=3, center=True).mean(), linewidth=3)
plt.plot(AMO_anom_nut.index, AMO_anom_nut.rolling(window=3, center=True).mean(), '--k', linewidth=3)
#plt.plot(O2_deep_anom.index, O2_deep_anom.rolling(window=3, center=True).mean())
plt.ylim([-1.75,1.75])
plt.grid()
plt.xlim([pd.to_datetime('1998'), pd.to_datetime('2019')])
plt.legend(['CIL', 'S','NO3','SiO','AMO'], loc=2)
plt.xlabel('Year', fontsize=16, fontweight='bold')
plt.ylabel(r'Standardized anom.', fontsize=16, fontweight='bold')

plt.annotate("", xy=(pd.to_datetime('2018-01-01'), 1.5), xytext=(pd.to_datetime('2018-01-01'), .25), arrowprops=dict(arrowstyle="->"))
plt.annotate("", xy=(pd.to_datetime('2018-01-01'), -1.5), xytext=(pd.to_datetime('2018-01-01'), -.25), arrowprops=dict(arrowstyle="->", facecolor='black'))
plt.text(pd.to_datetime('2018-01-01'), 1.5, 'warmer SST / stronger MOC', horizontalalignment='right', verticalalignment='bottom')
plt.text(pd.to_datetime('2018-01-01'), -1.5, 'cooler SST / weaker MOC', horizontalalignment='right', verticalalignment='top')
fig.set_size_inches(w=10,h=6)
fig_name = 'nutrients_AMO.png'
fig.set_dpi(300)
fig.savefig(fig_name)


fig = plt.figure(4)
fig.clf()
ax = fig.add_subplot(111)
plt.plot(Sanomaly_nut.index, Sanomaly_nut.rolling(window=12, center=True).mean(), linewidth=3)
#plt.plot(Sanomaly.index, Sanomaly.rolling(window=12, center=True).mean())
plt.plot(NO3_deep_anom.index, NO3_deep_anom.rolling(window=3, center=True).mean(), linewidth=3)
plt.plot(SIO_deep_anom.index, SIO_deep_anom.rolling(window=3, center=True).mean(), linewidth=3)
plt.plot(AOO_anom_nut.index, AOO_anom_nut.rolling(window=2, center=True).mean(), '--k', linewidth=3)
#plt.plot(O2_deep_anom.index, O2_deep_anom.rolling(window=3, center=True).mean())
plt.ylim([-1.75,1.75])
plt.grid()
plt.xlim([pd.to_datetime('1998'), pd.to_datetime('2017')])
plt.legend(['S','NO3','SiO','AOO'], loc=2)
plt.xlabel('Year', fontsize=16, fontweight='bold')
plt.ylabel(r'Standardized anom.', fontsize=16, fontweight='bold')

plt.annotate("", xy=(pd.to_datetime('2016-01-01'), 1.45), xytext=(pd.to_datetime('2016-01-01'), .25), arrowprops=dict(arrowstyle="->"))
plt.annotate("", xy=(pd.to_datetime('2016-01-01'), -1.45), xytext=(pd.to_datetime('2016-01-01'), -.25), arrowprops=dict(arrowstyle="->", facecolor='black'))
plt.text(pd.to_datetime('2016-01-01'), 1.5, 'less Arctic export', horizontalalignment='right', verticalalignment='bottom')
plt.text(pd.to_datetime('2016-01-01'), -1.5, 'more Arctic export', horizontalalignment='right', verticalalignment='top')
#plt.title('multi-index comparison')
fig.set_size_inches(w=10,h=6)
fig_name = 'nutrients_AOO.png'
fig.set_dpi(300)
fig.savefig(fig_name)

### ---- Salinity ---- ###
fig = plt.figure(6)
fig.clf()
ax = fig.add_subplot(111)
plt.plot(NAO_anom.index, NAO_anom.rolling(window=12, center=True).mean(), 'k')
plt.plot(AOO_anom.index, AOO_anom.rolling(window=5, center=True).mean(), color='gray')
plt.plot(meanS.index, Sanomaly.rolling(window=12, center=True).mean(), 'r')
plt.xlim([pd.to_datetime('1948'), pd.to_datetime('2017')])
plt.ylim([-2,2])
plt.grid()
#plt.legend(['NAO_w','AOO','AMO','CIL','S'])
plt.legend([r'$\rm NAO_{winter}$','AOO','S'])
plt.xlabel('Year', fontweight='bold')
plt.ylabel(r'Standardized anom.', fontweight='bold')
plt.title('multi-index comparison')
fig.set_size_inches(w=10,h=6)
fig_name = 'aoo_sal.png'
fig.set_dpi(300)
fig.savefig(fig_name)


### ---- Temp & Salinity - 3 depth range ---- ###
fig = plt.figure(7)
fig.clf()
ax = fig.add_subplot(111)
plt.plot(CIL.index, CIL_anom.rolling(window=5, center=True).mean(), linewidth=3)
plt.plot(meanS.index, Sanomaly.rolling(window=12, center=True).mean(), 'g')
plt.plot(meanS1.index, Sanomaly1.rolling(window=12, center=True).mean(), 'r')
plt.plot(meanS2.index, Sanomaly2.rolling(window=12, center=True).mean(), '--r')
plt.xlim([pd.to_datetime('1948'), pd.to_datetime('2017')])
plt.ylim([-2,2])
plt.grid()
#plt.legend(['NAO_w','AOO','AMO','CIL','S'])
plt.legend([r'CIL',r'$S_{0-300m}$',r'$S_{0-100m}$', r'$S_{100-200m}$'])
plt.xlabel('Year', fontweight='bold')
plt.ylabel(r'Standardized anom.', fontweight='bold')
fig.set_size_inches(w=10,h=6)
fig_name = 'T-S1-S2.png'
fig.set_dpi(300)
fig.savefig(fig_name)
