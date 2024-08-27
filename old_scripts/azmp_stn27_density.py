# -*- coding: utf-8 -*-
'''
Station 27 density analysis parameters.

This script will compute water density from temperature and salinity, as well as stratification and Mixed Layer Depth.

** see azmp_stn27.py, azmp_stn27_analysis.py 
for more options and ways to explore the dataset

Frederic.Cyr@dfo-mpo.gc.ca - September 2020

'''

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr
import datetime
import os
import warnings
import matplotlib.dates as mdates
from matplotlib.ticker import NullFormatter
from matplotlib.dates import MonthLocator, DateFormatter
import cmocean
import gsw

## font = {'family' : 'normal',
##         'weight' : 'bold',
##         'size'   : 14}
## plt.rc('font', **font)

## ---- Some custom parameters ---- ##
#s27 = [47.55,-52.59]
#dc = .1
s27 = [47.54667,-52.58667]
dc = .025
#variable = 'salinity'
use_viking = False
QC_SD = 4
# to append a new year to stn27_all_casts.nc
APPEND = False
APPEND_YEAR = 2023
# apply a moving average?
binning=True
move_ave = False
zbin = 5
CASTS_path = '/home/jcoyne/Documents/CASH/Combined_Data/CASTS_new-vertical_v2/'


## ---- Open data and select ---- ##
if os.path.isfile('operation_files/stn27_all_casts.nc'):
    ds = xr.open_dataset('operation_files/stn27_all_casts.nc') 
    if APPEND == True:
        ds2 = xr.open_dataset(CASTS_path + str(APPEND_YEAR) + '.nc')
        ds2 = ds2.where(ds2.instrument_ID!='MEDBA', drop=True) # BATHY GTS message 
        ds2 = ds2.where(ds2.instrument_ID!='MEDTE', drop=True) # TESAC GTS message 
        # Select a depth range
        ds2 = ds2.sel(level=ds2['level']<180)
        ds2 = ds2.sel(level=ds2['level']>0)
        # Select stn27 data according to lat-lon in a box [47.55,-52.59]
        ds2 = ds2.where((ds2.longitude>s27[1]-dc/2) & (ds2.longitude<s27[1]+dc/2), drop=True) # original one
        ds2 = ds2.where((ds2.latitude>s27[0]-dc/2) & (ds2.latitude<s27[0]+dc/2), drop=True)
        ds2 = ds2.sortby('time')
        ds = xr.concat([ds, ds2]) # to be tested!
else:
    ds = xr.open_mfdataset(CASTS_path+'*.nc')
    # Remome GTS datasets
    ds = ds.where(ds.instrument_ID!='MEDBA', drop=True) # BATHY GTS message 
    ds = ds.where(ds.instrument_ID!='MEDTE', drop=True) # TESAC GTS message 
    # Select a depth range
    ds = ds.sel(level=ds['level']<180)
    ds = ds.sel(level=ds['level']>0)
    # Select stn27 data according to lat-lon in a box [47.55,-52.59]
    ds = ds.where((ds.longitude>s27[1]-dc/2) & (ds.longitude<s27[1]+dc/2), drop=True) # original one
    ds = ds.where((ds.latitude>s27[0]-dc/2) & (ds.latitude<s27[0]+dc/2), drop=True)
    ds = ds.sortby('time')
    # Save data
    ds.to_netcdf('stn27_all_casts.nc')

# select variable
daT = ds['temperature']
daS = ds['salinity']
df_hydroT = daT.to_pandas()
df_hydroS = daS.to_pandas()

# Drop when all NaNs... but will still be in netCDF...
# Remove indices when T or S are all NaNs
idx_S = pd.notnull(df_hydroS).any(1).values.nonzero()[0]
df_hydroT = df_hydroT.iloc[idx_S]
df_hydroS = df_hydroS.iloc[idx_S]
idx_T = pd.notnull(df_hydroT).any(1).values.nonzero()[0]
df_hydroT = df_hydroT.iloc[idx_T]
df_hydroS = df_hydroS.iloc[idx_T]

## ---- Open Viking data and concatenate ---- ##
if use_viking:
    # open dataset
    ds_vik = xr.open_mfdataset('/home/cyrf0006/data/dev_database/viking_nc/*_viking.nc')
    # Select a depth range
    ds_vik = ds_vik.sel(level=ds_vik['level']<180)
    ds_vik = ds_vik.sel(level=ds_vik['level']>0)
    # get variable
    da_vikT = ds_vik['temperature']
    da_vikS = ds_vik['salinity']
    df_vikT = da_vikT.to_pandas()
    df_vikS = da_vikS.to_pandas()
    # remove NaNs
    idx_S = pd.notnull(df_vikS).any(1).values.nonzero()[0]
    df_vikT = df_vikT.iloc[idx_S]
    df_vikS = df_vikS.iloc[idx_S]
    idx_T = pd.notnull(df_vikT).any(1).values.nonzero()[0]
    df_vikT = df_vikT.iloc[idx_T]
    df_vikS = df_vikS.iloc[idx_T]

    # concatenate (better if they sahre same vertical resolution)
    df_temp = pd.concat([df_hydroT, df_vikT])
    df_sal = pd.concat([df_hydroS, df_vikS])

    df_temp = df_temp.interpolate(axis=1).where(df_temp.bfill(axis=1).notnull()) # here interpolate and leave NaNs.
    df_sal = df_sal.interpolate(axis=1).where(df_sal.bfill(axis=1).notnull()) # here interpolate and leave NaNs.
else:
    df_temp = df_hydroT.copy()
    df_sal = df_hydroS.copy()

## ---- Vertical binning or moving average to smooth the data ---- ##
if binning:
    old_z = df_temp.columns
    new_z = np.arange((old_z[0]+zbin)/2,old_z[-1],zbin)
    # temperature
    df_temp = df_temp.loc[:,0:new_z[-1]+1] # remove "unnecesary" columns
    dfT = df_temp.groupby(np.arange(len(df_temp.columns))//zbin, axis=1).mean()
    dfT.columns=new_z
    df_temp = dfT
    # salinty
    df_sal = df_sal.loc[:,0:new_z[-1]+1]
    dfS = df_sal.groupby(np.arange(len(df_sal.columns))//zbin, axis=1).mean()
    dfS.columns=new_z
    df_sal = dfS    
elif move_ave:
    df_temp = df_temp.rolling(zrolling, center=True, min_periods=1, axis=1).mean()
    df_sal = df_sal.rolling(zrolling, center=True, min_periods=1, axis=1).mean()

## ---- Compute density ---- ##
Z = df_temp.columns
SP = df_sal.values
PT = df_temp.values
SA = gsw.SA_from_SP(SP, Z, -50, 47)
CT = gsw.CT_from_pt(SA, PT)
RHO = gsw.rho(SA, CT, Z)
SIG0 = gsw.sigma0(SA, CT)
df_SA = pd.DataFrame(SA, index=df_temp.index, columns=df_temp.columns)
df_CT = pd.DataFrame(CT, index=df_temp.index, columns=df_temp.columns)
df_rho = pd.DataFrame(RHO, index=df_temp.index, columns=df_temp.columns)
df_sig = pd.DataFrame(SIG0, index=df_temp.index, columns=df_temp.columns)

## ---- QA/QC check ---- ##
idx_to_remove = []
#plt.close('all')
for i in np.arange(1,13):
    #fig = plt.figure()
    df_tmp = df_rho[df_rho.index.month == i]

    for idx in df_tmp.index:
        if df_tmp.loc[idx].dropna().size == 1: # size 1, remove.
            idx_to_remove.append(idx)
            df_tmp = df_tmp.drop(idx) # remove also from clim for next step
        elif ~np.all(np.diff(df_tmp.loc[idx].dropna().values.squeeze())>0): # check if monotonic
            idx_to_remove.append(idx)
            df_tmp = df_tmp.drop(idx) 
        else: 
            continue

## Remove idx_to_remove.
df_rho.drop(idx_to_remove, inplace=True)
df_sig.drop(idx_to_remove, inplace=True)
df_SA.drop(idx_to_remove, inplace=True)
df_CT.drop(idx_to_remove, inplace=True)
ds = ds.isel(time = ~np.isin(ds.time, np.array(idx_to_remove, dtype='datetime64')))

## ---- Vertically interpolate density ---- ##
df_SA = df_SA.interpolate(method='linear',axis=1).where(df_SA.bfill(axis=1).notnull())
df_CT = df_CT.interpolate(method='linear',axis=1).where(df_CT.bfill(axis=1).notnull())
df_rho = df_rho.interpolate(method='linear',axis=1).where(df_rho.bfill(axis=1).notnull())
df_sig = df_sig.interpolate(method='linear',axis=1).where(df_sig.bfill(axis=1).notnull())

df_rho.to_pickle('S27_rho_raw.pkl')
df_sig.to_pickle('S27_sig_raw.pkl')
df_SA.to_pickle('S27_SA_raw.pkl')
df_CT.to_pickle('S27_CT_raw.pkl')


## ---- Compute N2 and MLD ---- ##
print('Compute N2 and MLD for every index (make take some time...)')
SA = df_SA.values
CT = df_CT.values
N2 = np.full((df_rho.index.size, df_rho.columns.size-1), np.nan)
MLD = np.full((df_rho.index.size), np.nan)
for i,idx in enumerate(df_rho.index):     
        N2_tmp, pmid = gsw.Nsquared(SA[i,:], CT[i,:], Z, 47)
        N2[i,:] =  N2_tmp
        N2_tmp[np.where((pmid<=10) | (pmid>=100))] = np.nan
        if ~np.isnan(N2_tmp).all():
                MLD[i] = pmid[np.nanargmax(N2_tmp)]
print('  Done!')
df_N2 = pd.DataFrame(N2, index=df_CT.index, columns=pmid)
df_N2.to_pickle('S27_N2_raw.pkl')
MLD = pd.Series(MLD, index=df_CT.index)
MLD.to_pickle('S27_MLD_raw.pkl')


## ---- Save monthly averages ---- ##
df_rho.resample('M').mean().to_pickle('S27_rho_monthly.pkl')
MLD.resample('M').mean().to_pickle('S27_MLD_monthly.pkl')

## ---- Compute Stratification (now 5m-150m) ---- ## SINCE 2023 (copy wu)

#Determine density straight from NetCDF
Z = ds.level.values
SP = ds.salinity.values
PT = ds.temperature.values
SA = gsw.SA_from_SP(SP, Z, -50, 47)
CT = gsw.CT_from_pt(SA, PT)
RHO = gsw.rho(SA, CT, Z)

#Remove surface non-bottle measurements
bottles = np.isin(ds.instrument_ID.values, ['BO','FAPBO'])
RHO[~bottles,:2] = np.nan

#Record the relevant data
rho = pd.DataFrame(RHO,index=ds.time, columns=ds.level)

#Define the limits (range or individual level)
lower_limit = [0,5]
upper_limit = [48,53]

#Isolate for the relevant depths
warnings.simplefilter("ignore", category=RuntimeWarning)
indx = [np.nanargmin(np.abs(ds.level - ii)) for ii in lower_limit]
rho_top = np.nanmean(rho.iloc[:,indx[0]:indx[1]], axis=1)
indx = [np.nanargmin(np.abs(ds.level - ii)) for ii in upper_limit]
rho_bottom = np.nanmean(rho.iloc[:,indx[0]:indx[1]], axis=1)

#Remove values where top density is higher than bottom density
unstable_filt = rho_bottom <= rho_top
rho_top[unstable_filt] = np.nan
rho_bottom[unstable_filt] = np.nan

#Calculate the stratification
depth_range = np.mean(upper_limit) - np.mean(lower_limit)
strat = pd.Series((rho_bottom - rho_top)/depth_range, index=ds.time).resample('M').mean()
strat.to_pickle('S27_stratif_0-50_raw.pkl')
strat.resample('M').mean().to_pickle('S27_stratif_0-50_monthly.pkl')

#Define the limits (range or individual level)
lower_limit = [8,13]
upper_limit = [148,153]

#Isolate for the relevant depths
warnings.simplefilter("ignore", category=RuntimeWarning)
indx = [np.nanargmin(np.abs(ds.level - ii)) for ii in lower_limit]
rho_top = np.nanmean(rho.iloc[:,indx[0]:indx[1]], axis=1)
indx = [np.nanargmin(np.abs(ds.level - ii)) for ii in upper_limit]
rho_bottom = np.nanmean(rho.iloc[:,indx[0]:indx[1]], axis=1)

#Remove values where top density is higher than bottom density
unstable_filt = rho_bottom <= rho_top
rho_top[unstable_filt] = np.nan
rho_bottom[unstable_filt] = np.nan

#Calculate the stratification
depth_range = np.mean(upper_limit) - np.mean(lower_limit)
strat = pd.Series((rho_bottom - rho_top)/depth_range, index=ds.time).resample('M').mean()
strat.to_pickle('S27_stratif_10-150_raw.pkl')
strat.resample('M').mean().to_pickle('S27_stratif_10-150_monthly.pkl')

