# -*- coding: utf-8 -*-

###########CSL changed from SS to GSL in here#########
###########for ease of importing, data files have been altered due to formatting#########

"""
Created on Wed Feb  6 09:43:50 2019

@author: GIBBO
"""
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np
import seaborn as sns
import datetime
from scipy import stats
import seawater as swx
import cmocean
import co2sys
from sklearn.impute import SimpleImputer

###Adjust fontsize/weight
font = {'family' : 'DejaVu Sans',
        'weight' : 'bold',
        'size'   : 12}
plt.rc('font', **font)
pd.set_option('display.max_rows', 500)

################### Load BIOChem data into a Pandas DataFrame #####################################
################### the GSL 2014 dataset is also in this dataset #############################
df = pd.read_csv('C:\Users\gibbo\Documents\data\AZMP_OA\Discrete_Data_ 20190305.csv', delimiter=',')
df['timestamp'] = df.year.astype(str) + '-' + df.month.astype(str) + '-' + df.day.astype(str) + ' ' + (df.time.astype(str)).str.zfill(4)
df['timestamp'] = pd.to_datetime(df['timestamp'], format='%Y-%m-%d %H:%M')
#df.index = pd.to_datetime(df.timestamp)
##################### rename columns########################
df = df.rename(columns={'mission_descriptor' : 'TripID'})
df = df.rename(columns={'start_depth' : 'depth'})
df = df.rename(columns={'collector_station_name' : 'StationID'})
df = df.rename(columns={'Alkalinity_umol/kg' : 'TA'})
#df = df.rename(columns={'Chl_a_Holm-Hansen_F' : 'chla'})
df = df.rename(columns={'NO2NO3_Tech_F' : 'NO3'})
df = df.rename(columns={'O2_CTD_mLL' : 'O2CTD'})
df = df.rename(columns={'O2_Winkler_Auto' : 'O2'})
df = df.rename(columns={'PO4_Tech_F' : 'PO4'})
df = df.rename(columns={'SiO4_Tech_F' : 'SiO'})
df = df.rename(columns={'Temp_CTD_1968' : 'temperature'})
df = df.rename(columns={'Salinity_CTD' : 'salinity'})
df = df.rename(columns={'pH_invitro_total_f' : 'pH_25'}) ###this is pH measured at 25 degrees -not from CTD, used to calculate pH insitu
df = df.rename(columns={'pH_invitro_temp' : 'temp_pH'})  ###temperature at which the pH was measured (~25C)

################## combine BIO and IML (GSL 2014) columns (they had different headers) ###################
df.TA.update(df.ALKW_01)
df.TIC.update(df.TICW_01)
#df.chla.update(df.Chl_a_Welschmeyer_sF)
df.NO3.update(df.NO2NO3_Tech_SF)
df.NO3.update(df.NO2NO3_Tech_Fsh)
df.PO4.update(df.PO4_Tech_SF)
df.PO4.update(df.PO4_Tech_Fsh)
df.SiO.update(df.SiO4_Tech_Fsh)
df.SiO.update(df.SiO4_Tech_SF)
df.temperature.update(df.Temp_CTD_1990)
df.pH_25.update(df.pH_invitro_total_p)
df['Region'] = 'SS'
df.loc[(df['collector_event_id'].str.contains('IML')), 'Region'] = 'GSL'

df['O2GSL'] = df.loc[df['Region'].str.contains('GSL')]['O2']/44.661 ##change units of GSL data from umol/L to ml/L (/44.661)
df.O2.update(df.O2GSL)

################upload the BIO flag data ############################
bio_flags = pd.read_excel('C:\Users\gibbo\Documents\data\AZMP_OA\DIS_Data_20190305flags.xlsx')
variables = ['sample_id', 'method', 'data_qc_code']
bio_flags = bio_flags.loc[:,variables]
bio_flags = bio_flags.pivot_table(values='data_qc_code', index='sample_id', columns='method')
bio_flags = bio_flags.rename(columns={'Alkalinity_umol/kg' : 'TAflag'})
bio_flags = bio_flags.rename(columns={'pH_invitro_total_f' : 'pHflag'})
bio_flags = bio_flags.rename(columns={'TIC' : 'TICflag'})
bio_flags.TAflag.update(bio_flags.ALKW_01)
bio_flags.TICflag.update(bio_flags.TICW_01)
bio_flags.pHflag.update(bio_flags.pH_invitro_total_p)
bio_var = ['TAflag', 'TICflag', 'pHflag']
bio_flags = bio_flags.loc[:,bio_var]

df['sample_id'] = df['collector_event_id'] +"_"+ df['StationID'].map(str)
df = df.set_index(df['sample_id'], drop=True)
df = pd.merge(df, bio_flags, how='left', left_index=True, right_index=True)
df = df.reset_index(drop=True)
variables = ['TripID', 'Region', 'StationID', 'timestamp', 'latitude', 'longitude', 'depth', 'NO3', 'O2', 'TA', 'temp_pH', 'pH_25', 'PO4', 'salinity', 'SiO', 'temperature', 'TIC', 'TAflag', 'pHflag', 'TICflag']
df = df.loc[:,variables]
dropGSL = df[(df['Region'] == 'GSL') & (df['timestamp'].dt.year > 2014)].index
df.drop(dropGSL, inplace=True)
df.loc[(df['StationID'].str.contains('CSL', na=False)), 'Region'] = 'GSL'

################remove the Labrador Sea AR7W Line################
df=df[~df.StationID.str.contains('L3')]
###############rename the IML station ID to match rest of IML dataset############
df['StationID'] = df['StationID'].replace('0421', 'TIDM9')
df['StationID'] = df['StationID'].replace('0431', 'TIDM8')
df['StationID'] = df['StationID'].replace('0441', 'TIDM7')
df['StationID'] = df['StationID'].replace('0451', 'TIDM6')
df['StationID'] = df['StationID'].replace('0461', 'TIDM5')
df['StationID'] = df['StationID'].replace('0471', 'TIDM4')
df['StationID'] = df['StationID'].replace('0481', 'TIDM3')
df['StationID'] = df['StationID'].replace('0491', 'TIDM2')
df['StationID'] = df['StationID'].replace('0501', 'TIDM1')
df['StationID'] = df['StationID'].replace('0541', 'TDC1')
df['StationID'] = df['StationID'].replace('0551', 'TDC2')
df['StationID'] = df['StationID'].replace('0561', 'TDC3')
df['StationID'] = df['StationID'].replace('0571', 'TDC4')
df['StationID'] = df['StationID'].replace('0581', 'TDC5')
df['StationID'] = df['StationID'].replace('0591', 'TDC6')
df['StationID'] = df['StationID'].replace('0621', 'IF28')
df['StationID'] = df['StationID'].replace('0631', 'IF27')
df['StationID'] = df['StationID'].replace('0641', 'TCEN1')
df['StationID'] = df['StationID'].replace('0661', 'TCEN3')
df['StationID'] = df['StationID'].replace('0681', 'TCEN5')
df['StationID'] = df['StationID'].replace('0701', 'CH6')
df['StationID'] = df['StationID'].replace('0721', 'TBB1')
df['StationID'] = df['StationID'].replace('0731', 'TBB2')
df['StationID'] = df['StationID'].replace('0741', 'TBB3')
df['StationID'] = df['StationID'].replace('0761', 'TBB5')
df['StationID'] = df['StationID'].replace('0781', 'TBB7')
df['StationID'] = df['StationID'].replace('0821', 'IF10')
df['StationID'] = df['StationID'].replace('0841', 'IF12')
df['StationID'] = df['StationID'].replace('0861', 'CMO1/CH1')
df['StationID'] = df['StationID'].replace('0881', 'CH2')
df['StationID'] = df['StationID'].replace('0911', 'CMO2/CH9')
df['StationID'] = df['StationID'].replace('0921', 'IF1')
df['StationID'] = df['StationID'].replace('0931', 'TASO5')
df['StationID'] = df['StationID'].replace('0951', 'TASO3')
df['StationID'] = df['StationID'].replace('0971', 'TASO1')
df['StationID'] = df['StationID'].replace('0981', 'TSI1')
df['StationID'] = df['StationID'].replace('1001', 'TSI3')
df['StationID'] = df['StationID'].replace('1021', 'TSI5')
df['StationID'] = df['StationID'].replace('1031', 'TSI6')
df['StationID'] = df['StationID'].replace('1061', 'TESL2')
df['StationID'] = df['StationID'].replace('1071', 'TESL3/IML4')
df['StationID'] = df['StationID'].replace('1111', 'TESL6')
df['StationID'] = df['StationID'].replace('1121', 'TESL7')
df['StationID'] = df['StationID'].replace('1151', 'CM03/CH12')
df['StationID'] = df['StationID'].replace('123A1', '14ML')
df['StationID'] = df['StationID'].replace('1401', 'IF37')
df['StationID'] = df['StationID'].replace('1421', 'IF32')
df['StationID'] = df['StationID'].replace('1441', 'IF34')


################### Load NL data into a Pandas DataFrame #####################################
dfNL = pd.read_csv('C:\Users\gibbo\Documents\data\AZMP_OA\AZMP_OA_NL.csv', delimiter=',')
dfNL['timestamp'] = dfNL.Date.astype(str) + ' ' + dfNL.GMT.astype(str)
dfNL['timestamp'] =  pd.to_datetime(dfNL['timestamp'], format='%m/%d/%Y %H:%M:%S')
dfNL.index = dfNL.timestamp
dfNL2018DO = pd.read_excel('C:\Users\gibbo\Documents\data\AZMP_OA\NL_2018_satO2_perc.xlsx')
dfNL2018DO['timestamp'] = dfNL2018DO.Date.astype(str) + ' ' + dfNL2018DO.GMT.astype(str)
dfNL2018DO['timestamp'] =  pd.to_datetime(dfNL2018DO['timestamp'], format='%Y/%m/%d %H:%M:%S')
dfNL2018DO.index = dfNL2018DO.timestamp
dfNL['Dissolved Oxygen (mL/L)'].update(dfNL2018DO['Dissolved Oxygen (mL/L)'])
dfNL = dfNL.reset_index(drop=True)

## Remove space in column names ##
cols = dfNL.columns
cols = cols.map(lambda x: x.replace(' ', '_') if isinstance(x, (str, unicode)) else x)
dfNL.columns = cols

###########################rename columns##########################
dfNL = dfNL.rename(columns={'Ship_Trip' : 'TripID'})
dfNL = dfNL.rename(columns={'Latitude_(Decimal_Degrees)' : 'latitude'})
dfNL = dfNL.rename(columns={'Longiude_(Decimal_Degrees)' : 'longitude'})
dfNL = dfNL.rename(columns={'Station_Name' : 'StationID'})
dfNL = dfNL.rename(columns={'Pressure_(dbar)' : 'depth'})
dfNL = dfNL.rename(columns={'Salinity_(psu)' : 'salinity'})
dfNL = dfNL.rename(columns={'Temperature_(C)' : 'temperature'})
dfNL = dfNL.rename(columns={'Dissolved_Oxygen_(mL/L)' : 'O2'})
dfNL = dfNL.rename(columns={'Calibrated_Oxygen_(mL/L)' : 'oxygen'})
dfNL = dfNL.rename(columns={'Phosphate_Concentration_(mmol_m-3)' : 'PO4'})
dfNL = dfNL.rename(columns={'Silicate_Concentration_(mmol_m-3)' : 'SiO'})
dfNL = dfNL.rename(columns={'Nitrate_Concentration_(mmol_m-3)' : 'NO3'})
dfNL = dfNL.rename(columns={'pH_25_' : 'pH_25'})
dfNL = dfNL.rename(columns={'Discrete_Chlorophyll_a_(mg_m-3)' : 'chla'})
dfNL = dfNL.rename(columns={'Total_P_(mmol_m-3)' : 'PO4'})
dfNL = dfNL.rename(columns={'Total_Si_(mmol_m-3)' : 'SiO'})
dfNL = dfNL.rename(columns={'Total_NO2/3_(mmol_m-3)' : 'NO3'})
dfNL = dfNL.rename(columns={'Measured_TA_(micromol/kg)' : 'TA'})
dfNL = dfNL.rename(columns={'Measure_TCO2_(micromol/kg)' : 'TIC'})
dfNL = dfNL.rename(columns={'Density_(kg_m3)' : 'density'})
dfNL = dfNL.rename(columns={'TA_Data_Flag' : 'TAflag'})
dfNL = dfNL.rename(columns={'TCO2_Data_Flag' : 'TICflag'})
dfNL['Region'] = 'NL'
dfNL['temp_pH']= np.NaN
dfNL['pHflag']= np.NaN

dfNL = dfNL.loc[:,variables]

################################load IML data########################################
df2018 = pd.read_excel('C:\Users\gibbo\Documents\data\AZMP_OA\AZMP_OA_IML2018.xlsx', encoding='utf-8')
df2017 = pd.read_excel('C:\Users\gibbo\Documents\data\AZMP_OA\AZMP_OA_IML2017.xlsx', encoding='utf-8')
df2016 = pd.read_excel('C:\Users\gibbo\Documents\data\AZMP_OA\AZMP_OA_IML2016.xlsx', encoding='utf-8')
df2015 = pd.read_excel('C:\Users\gibbo\Documents\data\AZMP_OA\AZMP_OA_IML2015.xlsx', encoding='utf-8')
df2015['CTD Station (nom/no)'] = df2015['CTD Station (nom/no)'].replace(7, 'TESL2')
df2015['CTD Station (nom/no)'] = df2015['CTD Station (nom/no)'].replace(16, 'TSI3')
df2015['CTD Station (nom/no)'] = df2015['CTD Station (nom/no)'].replace(23, 'TASO3')
df2015['CTD Station (nom/no)'] = df2015['CTD Station (nom/no)'].replace(28, 'CMO2/CH9')
df2015['CTD Station (nom/no)'] = df2015['CTD Station (nom/no)'].replace(39, 'CH2')
df2015['CTD Station (nom/no)'] = df2015['CTD Station (nom/no)'].replace(51, 'TCEN3')
df2015['CTD Station (nom/no)'] = df2015['CTD Station (nom/no)'].replace(56, 'CSL_04')
df2015['CTD Station (nom/no)'] = df2015['CTD Station (nom/no)'].replace(61, 'TIDM3')
df2015['CTD Station (nom/no)'] = df2015['CTD Station (nom/no)'].replace(67, 'TIDM9')
df2015f = pd.read_excel('C:\Users\gibbo\Documents\data\AZMP_OA\AZMP_OA_IML2015f.xlsx', encoding='utf-8')
df2018f = pd.read_excel('C:\Users\gibbo\Documents\data\AZMP_OA\AZMP_OA_IML2018f.xlsx', encoding='utf-8')
df2017f = pd.read_excel('C:\Users\gibbo\Documents\data\AZMP_OA\AZMP_OA_IML2017f.xlsx', encoding='utf-8')
df2016f = pd.read_excel('C:\Users\gibbo\Documents\data\AZMP_OA\AZMP_OA_IML2016f.xlsx', encoding='utf-8')
df2018s = pd.read_excel('C:\Users\gibbo\Documents\data\AZMP_OA\AZMP_OA_IML2018s.xlsx', encoding='utf-8') ###some items in CTD Fichier were changed for an easier 'TripID'
df2017s = pd.read_excel('C:\Users\gibbo\Documents\data\AZMP_OA\AZMP_OA_IML2017s.xlsx', encoding='utf-8') ###some items in CTD Fichier were changed for an easier 'TripID'
dfIML = pd.concat([df2018f, df2018s, df2018, df2017f, df2017s, df2017, df2016f, df2016, df2015, df2015f], axis=0, sort=False)
dfIML = dfIML.reset_index()

################clean up formatting################
cols = dfIML.columns
cols = cols.map(lambda x: x.replace(' ', '_') if isinstance(x, (str, unicode)) else x)
dfIML.columns = cols

dfIML = dfIML.replace(r'^\s+$', np.nan, regex=True)

dfIML = dfIML.rename(columns={'CTD_Date_(jj-mmm-yyyy)' : 'Dates'})
dfIML = dfIML.rename(columns={'CTD_Heure_(GMT)' : 'Times'})
dfIML['timestamp'] = dfIML.Dates.astype(str) + ' ' + dfIML.Times.astype(str)
dfIML['timestamp'] =  pd.to_datetime(dfIML['timestamp'], format='%Y/%m/%d %H:%M:%S')

dfIML = dfIML.rename(columns={'CTD_Fichier_(nom)' : 'TripID'})
dfIML['TripID'] = dfIML['TripID'].replace('-','', regex=True)
dfIML.loc[(dfIML['TripID'].str.contains('iml|TEL|PER', na=False)), 'TripID'] = dfIML['TripID'].str[:10]
#dfIML.to_excel("C:\Users\gibbo\Documents\data\AZMP_OA\AZMP_OA_dfIML.xlsx")

###########################rename columns##########################
dfIML = dfIML.rename(columns={'CTD_Station_(nom/no)' : 'StationID'})
dfIML = dfIML.rename(columns={'CTD_Longitude_(degres)' : 'longitude'})
dfIML = dfIML.rename(columns={'CTD_Latitude_(degres)' : 'latitude'})
dfIML = dfIML.rename(columns={'CTD_PRES_(dbar)' : 'depth'})
dfIML = dfIML.rename(columns={'CTD_TE90_(celsius)' : 'temperature'})
dfIML = dfIML.rename(columns={'CTD_PSAL_(psu)' : 'salinity'})
dfIML = dfIML.rename(columns={'CTD_DOXY_(ml/l)' : 'oxyCTD'})
dfIML = dfIML.rename(columns={'labo_LABT_01_(deg_C)' : 'temp_pH'})
dfIML = dfIML.rename(columns={'Flag_At_' : 'TAflag'})
dfIML = dfIML.rename(columns={'Flag_pH_total_in_situ' : 'pHflag'})

############################rename stationIDs to be consistent##################
dfIML['StationID'] = dfIML['StationID'].replace('CMO3 (CH12)', 'CMO3/CH12')
dfIML['StationID'] = dfIML['StationID'].replace('TESL3 (RIKI)', 'TESL3')
dfIML['StationID'] = dfIML['StationID'].replace('TSI1 (CG)', 'TSI1')
dfIML['StationID'] = dfIML['StationID'].replace('TSI1/Courant_Gaspe', 'TSI1')
dfIML['StationID'] = dfIML['StationID'].replace('TIDM7/6.4', 'TIDM7')
dfIML['StationID'] = dfIML['StationID'].replace('TIDM5(4.5)', 'TIDM5')
dfIML['StationID'] = dfIML['StationID'].replace('TIDM5/4.5', 'TIDM5')
dfIML['StationID'] = dfIML['StationID'].replace('TIDM2(2.2)', 'TIDM2')
dfIML['StationID'] = dfIML['StationID'].replace('CMO1 (CH1)', 'CMO1/CH1')
dfIML['StationID'] = dfIML['StationID'].replace('CMO2 (CH9)', 'CMO2/CH9')
dfIML['StationID'] = dfIML['StationID'].replace('TIDM9 (Shediac)', 'TIDM9')
dfIML['StationID'] = dfIML['StationID'].replace('TIDM9,SHEDIAC', 'TIDM9')
dfIML['StationID'] = dfIML['StationID'].replace('TIDM9/Shediac', 'TIDM9')
dfIML['StationID'] = dfIML['StationID'].replace('TIDM9ShediacValley', 'TIDM9')
dfIML['StationID'] = dfIML['StationID'].replace('CMO1,CH1', 'CMO1/CH1')
dfIML['StationID'] = dfIML['StationID'].replace('TESL3,Riki', 'TESL3')
dfIML['StationID'] = dfIML['StationID'].replace('TESL3/Riki', 'TESL3')
dfIML['StationID'] = dfIML['StationID'].replace('ES2(13ML)', 'ES2/13ML')
dfIML['StationID'] = dfIML['StationID'].replace('If37', 'IF37')

dfIML['Region'] = 'GSL'

#############################calculate mean values of replicate analyses########################
dfIML['NO3'] = dfIML[['labo_NOX_03_(mmol/m**3)', 'labo_NOX_03_(mmol/m**3).1', 'labo_NOX_03_(mmol/m**3).2', 'labo_NOx_03_(mmol/m**3)', 'labo_NOx_03_(mmol/m**3).1', 'labo_NOx_03_(mmol/m**3).2']].mean(axis=1)
dfIML['PO4'] = dfIML[['labo_PO4_03_(mmol/m**3)', 'labo_PO4_03_(mmol/m**3).1', 'labo_PO4_03_(mmol/m**3).2']].mean(axis=1)
dfIML['SiO'] = dfIML[['labo_Si_03_(mmol/m**3)', 'labo_Si_03_(mmol/m**3).1', 'labo_Si_03_(mmol/m**3).2']].mean(axis=1)
dfIML['TA'] = dfIML[['labo_ALKW_01_(umol/kg**1)', 'labo_ALKW_01_(umol/kg**1).1']].mean(axis=1)
dfIML['pH_25'] = dfIML[['labo_LBPHT_02_(Total_scale)', 'labo_LBPHT_02_(Total_scale).1']].mean(axis=1)
dfIML['TIC']= np.NaN
dfIML['O2'] = dfIML[['labo_oxy_02_(ml/l)', 'labo_oxy_02_(ml/l).1', 'labo_OXY_02_(ml/l)', 'labo_OXY_02_(ml/l).1']].mean(axis=1)
dfIML = dfIML.loc[:,variables]


##############Make one df of all three Regions (datasets) of Original Data to plot (includes temperature and salinity data without associated carbonate data#######
df = pd.concat([df, dfNL, dfIML], axis=0)
df = df.reset_index(drop=True)
df['satO2'] = swx.satO2(df['salinity'], df['temperature']) 
df['satO2_perc'] = df['O2']/df['satO2']*100  ######calculate percent oxygen saturation
df['AOU'] = df['O2']-df['satO2']
df['NO3']=df['NO3']/1.025 #####convert nutrient data from uM=mmol/m3/umol/L to umol/kgSW (/1.025)
df['SiO']=df['SiO']/1.025
df['PO4']=df['PO4']/1.025
#df.to_excel("C:\Users\gibbo\Documents\data\AZMP_OA\AZMP_OA_Alldata.xlsx")
dforig = df.copy(deep=True) ###all original data, since some is replaced/modified

#############data clean up for CO2sys = estimate TA from TA-S plot, remove empty rows, remove flagged data, replace nutrient nan with zeros################
#############drop flagged data####################
df.drop(df.loc[df['TAflag']>=3].index, inplace=True)
df.drop(df.loc[df['pHflag']>=3].index, inplace=True)
df.drop(df.loc[df['TICflag']>=3].index, inplace=True)

##################calculate missing TA values based on TA-S relationship = linear regression (GSL tripID  ##########################
##################this is done for GSL spring 2015 (TRIPID IML2015160) using data from spring GSL data from all other years ##############

dflinreg = df.copy(deep=True)
dflinreg.dropna(subset=['TA', 'salinity'], inplace=True)

var_x = 'salinity'
var_y = 'TA'

dflinreg = dflinreg[dflinreg['Region'].str.contains("GSL")]
dflinreg = dflinreg[(dflinreg.timestamp.dt.month>=3) & (dflinreg.timestamp.dt.month<=6)]

linreg = sns.lmplot(var_x, var_y, dflinreg, legend=False, fit_reg=True, ci=95, scatter_kws={'s':20, 'alpha':1}, line_kws={'lw':1}, hue='Region', palette='plasma', height=4, aspect=1.1)
lines = sns.regplot(x=var_x, y=var_y, data=dflinreg, scatter=False, ci=95, color='k', line_kws={'lw': 1})
linreg.set_xlabels(r'Salinity')
linreg.set_ylabels(r'TA $(\rm \mu $mol/kg)')
axes = linreg.axes
axes[0,0].set_xlim(17, 38)
axes[0,0].set_ylim(1700, 2500)
#leg = linreg.axes.flat[0].get_legend()
#new_title = 'Region'
#leg.set_title(new_title)
#new_labels = ['SS', 'GSL', 'NLS']
#for t, l in zip(leg.texts, new_labels): t.set_text(l)
sns.set_style('whitegrid')
sns.despine(top=False, right=False, left=False, bottom=False)   
plt.plot([-2,22], [400,400,], 'black', linewidth=2, linestyle='dashed')

x=dflinreg[var_x].values
y=dflinreg[var_y].values

#xSS=df_SS[var_x].values
#xNLS=df_NLS[var_x].values
xGSL=dflinreg[var_x].values
#ySS=df_SS[var_y].values
#yNLS=df_NLS[var_y].values
yGSL=dflinreg[var_y].values

#labs=['All data', 'SS', 'NLS', 'GSL']
#xs=[x, xSS, xNLS, xGSL]
#ys=[y, ySS, yNLS, yGSL]
#for i,j,k in zip(xs,ys, labs):    
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
T = stats.t.ppf(1-0.025, len(dflinreg.salinity)-2)
ci = std_err*T
#print 'Linear Regression Statistics for ', k
print "linear regression line y=",slope, "x +",intercept
print "r squared: ", r_value**2
print "p_value: ", p_value
print "standard error: ", std_err
print "95% confidence interval: ", ci
print ''

df.loc[((df.timestamp.dt.year==2015) & (df.timestamp.dt.month==6)) & (df['Region'].str.contains('GSL')), 'TA'] = ((slope * df.salinity)+ intercept)

########################################################################

df.dropna(subset=['TA', 'TIC', 'pH_25'], axis=0, thresh=2, inplace=True)

#################### replace NaN in nutrient by zero (does not change anything in Excel CO2SYS)
no_nan = df.SiO[(df.PO4.isna()) | (df.SiO.isna())].size
perc_nan = float(no_nan)/df.shape[0]*100.0
df.PO4 = df.PO4.replace(np.nan, 0)
df.SiO = df.SiO.replace(np.nan, 0)
print(str(no_nan) + ' out of ' + str(df.shape[0]) + ' nutrient data where replaced by zeros (' + str(perc_nan) + ' %)')



##################below is the CO2sys script#############
###############change any parameters here################
df['par1type'] = (1)
df['par2type'] = (2)
df['par3type'] = (3)
df['presin'] = (0)
df['pHscale'] = (1)
df['k1k2c'] = (4)
df['kso4c'] = (1)

####################I separated the data according to sets of carbonate parameters###########
####################for TA and DIC data#######################
df_shelf = df[df['pH_25'].isnull()] ###keeps all rows with a pH value of zero, left with TA and DIC
par1type = df_shelf['par1type'] # The first parameter supplied is of type "1", which is "alkalinity"
par2type = df_shelf['par2type']  # The first parameter supplied is of type "2", which is "DIC"
par3type = df_shelf['par3type']  # The first parameter supplied is of type "3", which is "pH"
presin   = df_shelf['presin']  # Pressure at input conditions
tempout  = df_shelf['temperature']  # Temperature at output conditions.
presout  = df_shelf['depth']  # Pressure    at output conditions.
pHscale  = df_shelf['pHscale']  # pH scale at which the input pH is reported ("1" means "Total Scale")
k1k2c    = df_shelf['k1k2c']  # Choice of H2CO3 and HCO3- dissociation constants K1 and K2 ("4" means "Mehrbach refit")
kso4c    = df_shelf['kso4c']  # Choice of HSO4- dissociation constants KSO4 ("1" means "Dickson")
alk_s  = df_shelf['TA']
dic_s  = df_shelf['TIC']
sal_s  = df_shelf['salinity']
pH_s = df_shelf['pH_25'] # pH measured in the lab (at ~25C)
temp_s = df_shelf['temp_pH'] # Temperature at which the pH was measured (~25C)
si_s   = df_shelf['SiO']
p_s    = df_shelf['PO4']
out_shelf, niceheaders = co2sys.CO2SYS(alk_s, dic_s, par1type, par2type, sal_s, temp_s, tempout, presin, presout, si_s, p_s, pHscale, k1k2c, kso4c)
dfCO2out_shelf = pd.DataFrame.from_dict(out_shelf, orient='index')
dfCO2out_shelf = dfCO2out_shelf.transpose()
dfCO2out_shelf = dfCO2out_shelf[['TAlk', 'TCO2', 'pHoutTOTAL', 'pCO2out', 'OmegaCAout', 'OmegaARout']].copy() 
dfCO2out_shelf = dfCO2out_shelf.set_index(df_shelf.index)
dfCO2shelf = pd.merge(df_shelf, dfCO2out_shelf, how='outer', left_index=True, right_index=True)

##################for TA and pH data############################
df_gulf = df.dropna(subset=['pH_25']) ###drops all rows with a pH value of zero, left with pH and TA
par1type = df_gulf['par1type'] # The first parameter supplied is of type "1", which is "alkalinity"
par2type = df_gulf['par2type']  # The first parameter supplied is of type "2", which is "DIC"
par3type = df_gulf['par3type']  # The first parameter supplied is of type "3", which is "pH"
presin = df_gulf['presin']  # Pressure at input conditions
tempout = df_gulf['temperature']  # Temperature at output conditions.
presout = df_gulf['depth']  # Pressure    at output conditions.
pHscale = df_gulf['pHscale']  # pH scale at which the input pH is reported ("1" means "Total Scale")
k1k2c    = df_gulf['k1k2c']  # Choice of H2CO3 and HCO3- dissociation constants K1 and K2 ("4" means "Mehrbach refit")
kso4c    = df_gulf['kso4c']  # Choice of HSO4- dissociation constants KSO4 ("1" means "Dickson")
alk_s  = df_gulf['TA']
dic_s  = df_gulf['TIC']
sal_s  = df_gulf['salinity']
pH_s = df_gulf['pH_25'] # pH measured in the lab (at ~25C)
temp_s = df_gulf['temp_pH'] # Temperature at which the pH was measured (~25C)
si_s   = df_gulf['SiO']
p_s    = df_gulf['PO4']
out_gulf, niceheaders = co2sys.CO2SYS(alk_s, pH_s, par1type, par3type, sal_s, temp_s, tempout, presin, presout, si_s, p_s, pHscale, k1k2c, kso4c)
dfCO2out_gulf = pd.DataFrame.from_dict(out_gulf, orient='index')
dfCO2out_gulf = dfCO2out_gulf.transpose()
dfCO2out_gulf = dfCO2out_gulf[['TAlk', 'TCO2', 'pHoutTOTAL', 'pCO2out', 'OmegaCAout', 'OmegaARout']].copy() 
dfCO2out_gulf = dfCO2out_gulf.set_index(df_gulf.index)
dfCO2gulf = pd.merge(df_gulf, dfCO2out_gulf, how='outer', left_index=True, right_index=True)


#############combine data and create/cleanup new dataframe for PCA or plotting############
dfCO2sys = pd.concat([dfCO2shelf, dfCO2gulf], axis=0)
variables = ['TripID', 'Region', 'StationID', 'timestamp', 'latitude', 'longitude', 'depth', 'temperature', 'salinity', 'satO2_perc', 'NO3', 'PO4', 'SiO', 'TCO2', 'TAlk', 'pHoutTOTAL', 'pCO2out', 'OmegaCAout', 'OmegaARout']
dfCO2sys = dfCO2sys.loc[:,variables]
dfCO2sys = dfCO2sys.rename(columns={'TAlk' : 'TA'})
dfCO2sys = dfCO2sys.rename(columns={'TCO2' : 'TIC'})
dfCO2sys = dfCO2sys.rename(columns={'pHoutTOTAL' : 'pH'})
dfCO2sys = dfCO2sys.rename(columns={'pCO2out' : 'pCO2'})
dfCO2sys = dfCO2sys.rename(columns={'OmegaCAout' : 'Omega_C'})
dfCO2sys = dfCO2sys.rename(columns={'OmegaARout' : 'Omega_A'})
dfCO2sys = dfCO2sys.dropna(subset=['Omega_C'], axis=0)
dfCO2sys.to_excel("C:\Users\gibbo\Documents\data\AZMP_OA\AZMP_OA_CO2stats.xlsx")

################recombine with original dataset that includes the NaNs for plotting then cleanup##################
var = ['TripID', 'Region', 'StationID', 'timestamp', 'latitude', 'longitude', 'depth', 'temperature', 'salinity', 'satO2_perc', 'NO3', 'PO4', 'SiO']
dforig = dforig.loc[:,var]
dfCO2var = dfCO2sys.loc[:,['TA', 'TIC', 'pH', 'pCO2', 'Omega_A', 'Omega_C']]
dfplot = pd.merge(dforig, dfCO2var, how='left', left_index=True, right_index=True)
dfplot = dfplot.dropna(axis=0, thresh=8)
dfplot.to_excel("C:\Users\gibbo\Documents\data\AZMP_OA\AZMP_OA_plot.xlsx")



##########below is the original CO2SYS function################################
#import pandas as pd
#import co2sys

#par1type = 1  # The first parameter supplied is of type "1", which is "alkalinity"
#par2type = 2  # The first parameter supplied is of type "2", which is "DIC"
#par3type = 3  # The first parameter supplied is of type "3", which is "pH"
#presin   = 0  # Pressure at input conditions
#tempout  = 15.267  # Temperature at output conditions.
#presout  = 3.498  # Pressure    at output conditions.
#pHscale  = 1  # pH scale at which the input pH is reported ("1" means "Total Scale")
#k1k2c    = 4  # Choice of H2CO3 and HCO3- dissociation constants K1 and K2 ("4" means "Mehrbach refit")
#kso4c    = 1  # Choice of HSO4- dissociation constants KSO4 ("1" means "Dickson")
#alk_s  = 2134.7
#dic_s  = 1948.8
#sal_s  = 31.066
#pH_s = 7.892
#temp_s = 25
#si_s   = 1.005
#p_s    = 0.1945
#
#out, niceheaders = co2sys.CO2SYS(alk_s, pH_s, par1type, par3type, sal_s, temp_s, tempout, presin, presout, si_s, p_s, pHscale, k1k2c, kso4c)
#dfCO2out = pd.DataFrame.from_dict(out, orient='index')
#dfCO2sys = dfCO2out.transpose()
#dfCO2sys = dfCO2sys[['TAlk', 'TCO2', 'pHoutTOTAL', 'pCO2out', 'OmegaCAout', 'OmegaARout']].copy() 
#dfCO2sys = dfCO2sys.dropna(axis=0, how='any')
#dfCO2sys.to_excel("C:\Users\gibbo\Documents\data\AZMP_OA\AZMP_OA__NL_CO2sys.xlsx")
