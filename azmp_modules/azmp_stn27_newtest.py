import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr
import datetime
import os
import matplotlib.dates as mdates
from matplotlib.ticker import NullFormatter
from matplotlib.dates import MonthLocator, DateFormatter
from matplotlib.colors import from_levels_and_colors
import cmocean
import gsw
import unicodedata
import warnings

'''
This is a test script.
The goal is to convert azmp_stn27.py into a function package.
This will make it easier to run and to change from year to year.
'''


'''
# Derived parameter
V = {}
Vanom = {}
Vanom_ticks = {}
CMAP = {}
variable = 'temperature'
V[variable] = np.concatenate([np.arange(-1.75, 0.25, .25), np.arange(1, 14, 1)])  
Vanom[variable] = np.linspace(-5.5, 5.5, 23)
Vanom_ticks[variable] = np.array([-5, -3, -1, 1, 3, 5])
Vanom[variable] = np.delete(Vanom[variable], np.where(Vanom[variable]==0))
CMAP[variable] = cmocean.tools.lighten(cmocean.cm.thermal, .9)
variable = 'salinity'
V[variable] = np.arange(30, 33.5, .25)
Vanom[variable] = np.linspace(-1, 1, 21) 
Vanom[variable] = np.delete(Vanom[variable], np.where(Vanom[variable]==0))
Vanom_ticks[variable] = np.linspace(-1, 1, 11) 
CMAP[variable] = cmocean.cm.haline


#Sample input for stn27_dataisolate
file_location = 'operation_files/stn27_all_casts.nc'
CASTS_path = '/home/jcoyne/Documents/CASH/Combined_Data/CASTS_new-vertical_v2/*.nc'
problem_casts = ['1988_10073164.nc','1990_02183115.nc','1991_11229001.nc']
s27 = [47.54667,-52.58667]
dc = .025
year_clim = [1991, 2020]
binning=True
move_ave = False
current_year = 2023
XLIM = [datetime.date(current_year, 1, 1), datetime.date(current_year, 12, 31)]
'''

def stn27_dataisolate(file_location,CASTS_path,s27_loc,dc,problem_casts):
	#Import CASTS
	if os.path.isfile(os.path.expanduser(file_location)):
		ds = xr.open_dataset(file_location)
	else:
		#If the files not available, make one
		ds = xr.open_mfdataset(CASTS_path)
		#Isolate for the standard depth range
		ds = ds.sel(level=ds['level']<180)
		#Select stn27 data according to lat-lon in a box [47.55,-52.59]
		ds = ds.where((ds.longitude>s27_loc[1]-dc/2).compute() & (ds.longitude<s27_loc[1]+dc/2).compute(), drop=True)
		ds = ds.where((ds.latitude>s27_loc[0]-dc/2).compute() & (ds.latitude<s27_loc[0]+dc/2).compute(), drop=True)
		ds = ds.sortby('time')
		#Remove problem casts
		ds = ds.isel(time = ~np.isin(ds.file_names,problem_casts))
		#Re-state string variables as strings before saving
		for i in ds.variables.keys():
			if ds[i].dtype == 'O':
				ds[i] = ds[i].astype('U40')
		#Save data
		#time_attrs = ds['time'].attrs
		#time_attrs['_FillValue'] = -99.9999
		#ds['time'].attrs = time_attrs
		warnings.simplefilter(action='ignore', category=FutureWarning)
		ds.to_netcdf(file_location)


def stn27_climatology(file_location,year_clim,current_year,zbin=5,binning=True,move_ave=False):
	#Import the isolated data
	ds = xr.open_dataset(file_location)
	#Isolate the temperature and salinity
	da = {'temperature': ds['temperature'], 'salinity': ds['salinity']}
	#Cycle through each variable
	df = {}
	df_monthly = {}
	weekly_clim = {}
	df_year = {}
	df_weekly = {}
	anom = {}
	for variable in da:
		df[variable] = da[variable].to_pandas()
		df[variable].dropna(axis=0,how='all',inplace=True)
		#Cycle through if binning or moving average
		if binning:
			old_z = df[variable].columns
			new_z = np.arange((old_z[0]+zbin)/2,old_z[-1],zbin)
			dfT = df[variable].T.groupby(np.arange(len(df[variable].columns))//zbin).mean().T
			dfT.columns=new_z[:]
			df[variable] = dfT
		elif move_ave:
			df[variable] = df[variable].rolling(zbin,center=True,min_periods=1,axis=1).mean()
		#Determine the monthly average (15th of the month)
		df_monthly[variable] = df[variable].resample('MS').mean()
		df_monthly[variable].index = df_monthly[variable].index + pd.tseries.frequencies.to_offset('14d')
		df_monthly[variable].to_pickle('S27_'+variable+'_monthly.pkl')
		#Determine the climatology
		df_clim_period = df_monthly[variable][(df_monthly[variable].index.year>=year_clim[0]) & (df_monthly[variable].index.year<=year_clim[1])]
		monthly_clim = df_clim_period.groupby(df_clim_period.index.month).mean()
		monthly_clim.index = pd.to_datetime(monthly_clim.index.values, format='%m')
		monthly_clim.to_pickle('S27_' + variable + '_monthly_clim.pkl')
		#Weekly clim (upsample monthly clim to weekly)
		last_row = monthly_clim.iloc[0]
		last_row.name = np.datetime64('1901-01-01')
		monthly_clim = monthly_clim._append(last_row)
		weekly_clim[variable] = monthly_clim.resample('W').mean().interpolate(method='linear')
		weekly_clim[variable] = weekly_clim[variable].iloc[:-1]
		weekly_clim[variable].to_pickle('S27_' + variable + '_weekly_clim.pkl')
		# Update climatology index to current year
		weekly_clim[variable].index = pd.to_datetime(str(current_year) + '-' + weekly_clim[variable].index.month.astype(str) + '-' + weekly_clim[variable].index.day.astype(str))
		#Yearly average
		df_year[variable] = df[variable][df[variable].index.year==current_year]

		#Average to monthly and interpolate temporally, no greater than 2 months
		df_monthly[variable] = df_year[variable].resample('MS').mean()#.interpolate(method='linear')
		s = df_monthly[variable].notnull()
		s = s.ne(s.shift()).cumsum()
		m_df_monthly = pd.DataFrame([df_monthly[variable].groupby([s[i], df_monthly[variable][i].isnull()])[i].transform('size').where(df_monthly[variable][i].isnull()) for i in s.columns]).T
		df_monthly[variable] = df_monthly[variable].interpolate(limit_area='inside', method='linear').mask(m_df_monthly>2)
		#Average to monthly and interpolate temporally, no greater than 2 months (8 weeks)
		df_weekly[variable] = df_year[variable].resample('W-MON').mean()
		s = df_weekly[variable].notnull()
		s = s.ne(s.shift()).cumsum()
		m_df_weekly = pd.DataFrame([df_weekly[variable].groupby([s[i], df_weekly[variable][i].isnull()])[i].transform('size').where(df_weekly[variable][i].isnull()) for i in s.columns]).T
		df_weekly[variable] = df_weekly[variable].interpolate(limit_area='inside', method='linear').mask(m_df_weekly>8)
		weekly = weekly_clim[variable].resample('W-MON').mean().interpolate(method='linear')
		anom[variable] = df_weekly[variable] - weekly
	return df, df_monthly, weekly_clim, df_weekly, df_year, anom


#TEST SCRIPT
#df,df_monthly,weekly_clim,df_weekly,df_year,anom = stn27_climatology(file_location,year_clim,current_year)
#Save the current year monthly average
#for variable in df_monthly:
#	csv_file = 'monthly_' + variable +  '_' + str(current_year) + '.csv'
#	df_monthly[variable].T.to_csv(csv_file, float_format='%.4f')


def stn27_plottinginfo():
	V = {}
	Vanom = {}
	Vanom_ticks = {}
	CMAP = {}
	variable = 'temperature'
	V[variable] = np.concatenate([np.arange(-1.75, 0.25, .25), np.arange(1, 14, 1)])  
	Vanom[variable] = np.linspace(-5.5, 5.5, 23)
	Vanom_ticks[variable] = np.array([-5, -3, -1, 1, 3, 5])
	Vanom[variable] = np.delete(Vanom[variable], np.where(Vanom[variable]==0))
	CMAP[variable] = cmocean.tools.lighten(cmocean.cm.thermal, .9)
	variable = 'salinity'
	V[variable] = np.arange(30, 33.5, .25)
	Vanom[variable] = np.linspace(-1, 1, 21) 
	Vanom[variable] = np.delete(Vanom[variable], np.where(Vanom[variable]==0))
	Vanom_ticks[variable] = np.linspace(-1, 1, 11) 
	CMAP[variable] = cmocean.cm.haline
	return V,Vanom,Vanom_ticks,CMAP


def stn27_climatology_plot(weekly_clim,XLIM):
	#Import the plotting details
	V,Vanom,Vanom_ticks,CMAP = stn27_plottinginfo()
	#Cycle through each variable
	for variable in weekly_clim:
		#Create the figure and plot
		fig, ax = plt.subplots(nrows=1, ncols=1)
		c = plt.contourf(
			weekly_clim[variable].index,
			weekly_clim[variable].columns,
			weekly_clim[variable].values.T,
			V[variable], extend='both', cmap=CMAP[variable])
		plt.ylim([0, 175])
		plt.xlim([XLIM[0], XLIM[1]])
		#Create titles
		plt.ylabel('Depth (m)', fontsize=15, fontweight='bold')
		plt.title('Climatology')
		ax.invert_yaxis()
		cax = fig.add_axes([0.91, .15, 0.01, 0.7])
		cb = plt.colorbar(c, cax=cax, orientation='vertical')
		if variable == 'temperature':
			ccc = ax.contour(
				weekly_clim[variable].index,
				weekly_clim[variable].columns,
				weekly_clim[variable].values.T,
				[0,], colors='k', linewidths=3)
			cb.set_label(r'$\rm T(^{\circ}C)$', fontsize=12, fontweight='normal')
		elif variable == 'salinity':
			cb.set_label(r'S', fontsize=12, fontweight='normal')
		#Format date ticks
		ax.xaxis.set_major_locator(MonthLocator())
		ax.xaxis.set_minor_locator(MonthLocator(bymonthday=15))
		ax.xaxis.set_major_formatter(NullFormatter())
		ax.xaxis.set_minor_formatter(NullFormatter())
		ax.xaxis.label.set_visible(False)
		#Save Figure
		fig.set_size_inches(w=12, h=6)
		outfile_clim = 's27_' + variable + '_clim.png'
		fig.savefig(outfile_clim, dpi=200)
		os.system('convert -trim ' + outfile_clim + ' ' + outfile_clim)
		#Save French Figure
		ax.set_ylabel('Profondeur (m)', fontsize=15, fontweight='bold')
		ax.set_title('Climatologie')
		fig.set_size_inches(w=12, h=6)
		outfile_climFR = 's27_' + variable + '_clim_FR.png'
		fig.savefig(outfile_climFR, dpi=200)
		os.system('convert -trim ' + outfile_climFR + ' ' + outfile_climFR)

#TEST SCRIPT
#stn27_climatology_plot(weekly_clim, V, XLIM, CMAP)


def stn27_currentyear_plot(df_weekly,df_year,current_year,XLIM):
	#Import the plotting details
	V,Vanom,Vanom_ticks,CMAP = stn27_plottinginfo()
	#Cycle through each variable
	for variable in df_weekly:
		#Create figure and plot
		fig, ax = plt.subplots(nrows=1, ncols=1)
		c = plt.contourf(
			df_weekly[variable].index,
			df_weekly[variable].columns,
			df_weekly[variable].values.T,
			V[variable], extend='both', cmap=CMAP[variable])
		plt.plot(
			np.array(df_year[variable].index),
			np.repeat(0, df_year[variable].index.size),
			'|k', markersize=20)
		plt.ylim([0, 175])
		plt.xlim([XLIM[0], XLIM[1]])
		#Create titles
		plt.ylabel('Depth (m)', fontsize=15, fontweight='bold')
		plt.title(str(current_year) + ' observations')
		ax.invert_yaxis()
		cax = fig.add_axes([0.91, .15, 0.01, 0.7])
		cb = plt.colorbar(c, cax=cax, orientation='vertical')
		if variable == 'temperature':
			ccc = ax.contour(
				df_weekly[variable].index,
				df_weekly[variable].columns,
				df_weekly[variable].values.T,
				[0,], colors='k', linewidths=3)
			cb.set_label(r'$\rm T(^{\circ}C)$', fontsize=12, fontweight='normal')
		elif variable == 'salinity':
			cb.set_label(r'S', fontsize=12, fontweight='normal')
		#Format date ticks
		ax.xaxis.set_major_locator(MonthLocator())
		ax.xaxis.set_minor_locator(MonthLocator(bymonthday=15))
		ax.xaxis.set_major_formatter(NullFormatter())
		ax.xaxis.set_minor_formatter(NullFormatter())
		ax.xaxis.label.set_visible(False)
		#Save Figure
		fig.set_size_inches(w=12, h=6)
		outfile_year = 's27_' + variable + '_' + str(current_year) + '.png'
		fig.savefig(outfile_year, dpi=200)
		os.system('convert -trim ' + outfile_year + ' ' + outfile_year)
		#Save French Figure
		ax.set_ylabel('Profondeur (m)', fontsize=15, fontweight='bold')
		ax.set_title('Observations en ' + str(current_year))
		fig.set_size_inches(w=12, h=6)
		outfile_yearFR = 's27_' + variable + '_' + str(current_year) + '_FR.png'
		fig.savefig(outfile_yearFR, dpi=200)
		os.system('convert -trim ' + outfile_yearFR + ' ' + outfile_yearFR)

#TEST SCRIPT
#stn27_currentyear_plot(df_weekly,df_year,current_year,V,XLIM,CMAP)


def stn27_anomaly_plot(anom,current_year,XLIM):
	#Import the plotting details
	V,Vanom,Vanom_ticks,CMAP = stn27_plottinginfo()
	#Cycle through each variable
	for variable in anom:
		#Create figure and plot
		fig, ax = plt.subplots(nrows=1, ncols=1)
		c = plt.contourf(
			anom[variable].index,
			anom[variable].columns,
			anom[variable].values.T,
			Vanom[variable], extend='both', cmap=cmocean.cm.balance)
		plt.ylim([0, 175])
		plt.xlim([XLIM[0], XLIM[1]])
		#Create titles
		plt.ylabel('Depth (m)', fontsize=15, fontweight='bold')
		plt.title(str(current_year) + ' anomaly')
		ax.invert_yaxis()
		cax = fig.add_axes([0.91, .15, 0.01, 0.7])
		cb = plt.colorbar(c, cax=cax, ticks=Vanom_ticks[variable], orientation='vertical')
		if variable == 'temperature':
			cb.set_label(r'$\rm T(^{\circ}C)$', fontsize=12, fontweight='normal')
		elif variable == 'salinity':
			cb.set_label(r'S', fontsize=12, fontweight='normal')
		#Format date ticks
		ax.xaxis.set_major_locator(MonthLocator())
		ax.xaxis.set_minor_locator(MonthLocator(bymonthday=15))
		ax.xaxis.set_major_formatter(NullFormatter())
		ax.xaxis.set_minor_formatter(DateFormatter('%b'))
		#Save Figure
		fig.set_size_inches(w=12, h=6)
		outfile_anom = 's27_' + variable + '_anom_' + str(current_year) + '.png'
		fig.savefig(outfile_anom, dpi=200)
		os.system('convert -trim ' + outfile_anom + ' ' + outfile_anom)
		#Save French Figure
		ax.xaxis.set_minor_formatter(DateFormatter('%m'))
		ax.set_title('Anomalie en ' + str(current_year))
		ax.set_ylabel('Profondeur (m)', fontsize=15, fontweight='bold')
		fig.set_size_inches(w=12, h=6)
		outfile_anomFR = 's27_' + variable + '_anom_' + str(current_year) + '_FR.png'
		fig.savefig(outfile_anomFR, dpi=200)
		os.system('convert -trim ' + outfile_anomFR + ' ' + outfile_anomFR)

#TEST SCRIPT
'''
stn27_anomaly_plot(anom,current_year,V,XLIM,Vanom,Vanom_ticks)
#Covert to subplots and remove individual plots
for variable in anom:
	os.system('montage '+\
		's27_'+variable+'_2023.png '+\
		's27_'+variable+'_clim.png '+\
		's27_'+variable+'_anom_2023.png '+\
		'-tile 1x3 -geometry +10+10  -background white  s27_'+variable+'_subplot_'+str(current_year)+'.png')
	os.system('montage '+\
		's27_'+variable+'_2023_FR.png '+\
		's27_'+variable+'_clim_FR.png '+\
		's27_'+variable+'_anom_2023_FR.png '+\
		'-tile 1x3 -geometry +10+10  -background white  s27_'+variable+'_subplot_'+str(current_year)+'_FR.png')
	os.system('rm '+\
		's27_'+variable+'_2023.png '+\
		's27_'+variable+'_clim.png '+\
		's27_'+variable+'_anom_2023.png')
	os.system('rm '+\
		's27_'+variable+'_2023_FR.png '+\
		's27_'+variable+'_clim_FR.png '+\
		's27_'+variable+'_anom_2023_FR.png')
'''
def station_occupations(file_location,current_year,year_clim,binning=True,move_ave=False,zbin=5):
	#Import the time data
	ds = xr.open_dataset(file_location)
	da_occu = ds['time']
	df_occu = da_occu.to_pandas()
	df_occu = df_occu[df_occu.index.year>=1945]
	#Plot the weekly occupations
	fig, ax = plt.subplots(nrows=1, ncols=1)
	plt.clf()
	plt.plot(df_occu.index.year.values, df_occu.index.isocalendar().week.values, '.k')
	plt.ylabel('Week of year')
	plt.xlabel('Station 27 occupation')
	#Save Figure
	fig.set_size_inches(w=12, h=6)
	outfile_occu = 's27_occupation.png'
	fig.savefig(outfile_occu, dpi=200)
	os.system('convert -trim ' + outfile_occu + ' ' + outfile_occu)
	#Plot weekly occupations plus the no. of weeks per year
	W = df_occu.resample('w').count()
	W = W[W>0]  
	W = W.resample('Y').count()
	fig = plt.figure()
	#First axis
	ax1 = plt.subplot2grid((2, 1), (0, 0))
	plt.bar(W.index.year, W.values)
	ax1.xaxis.label.set_visible(False)
	ax1.tick_params(labelbottom='off')
	ax1.set_ylabel('No. of weekly occupations')
	plt.title('Station 27 occupation')
	#Second axis
	ax2 = plt.subplot2grid((2, 1), (1, 0))
	plt.plot(df_occu.index.year.values, df_occu.index.isocalendar().week.values, '.k')
	ax2.set_ylabel('week of year')
	# Save Figure
	fig.set_size_inches(w=12, h=12)
	outfile_occu = 's27_occupation_stats.png'
	fig.savefig(outfile_occu, dpi=200)
	os.system('convert -trim ' + outfile_occu + ' ' + outfile_occu)
	#Sampling this year
	A = df_occu.resample('M').count()
	print('Stn27 sampling in ' + str(current_year) + ':')
	print(A[A.index.year==current_year])
	#Plot a simple annual time series for temperature
	df = ds['temperature'].to_pandas()
	df.dropna(axis=0,how='all',inplace=True)
	#Cycle through if binning or moving average
	if binning:
		old_z = df.columns
		new_z = np.arange((old_z[0]+zbin)/2,old_z[-1],zbin)
		dfT = df.T.groupby(np.arange(len(df.columns))//zbin).mean().T
		dfT.columns=new_z[:]
		df = dfT
	elif move_ave:
		df = df.rolling(zbin,center=True,min_periods=1,axis=1).mean()
	fig = plt.figure(1)
	fig.clf()
	df_tmp = df.mean(axis=1).resample('As').mean()
	df_tmp = df_tmp[df_tmp.index.year>=1950]
	df_clim = df_tmp[(df_tmp.index.year>=year_clim[0]) & (df_tmp.index.year<=year_clim[1])]
	mean = df_clim.mean()
	std = df_clim.std()
	df_tmp.plot(linewidth=2)
	plt.fill_between([df_tmp.index[0], df_tmp.index[-1]], [mean+.5*std, mean+.5*std], [mean-.5*std, mean-.5*std], facecolor='gray', alpha=.2)
	plt.ylabel(r'$\rm T(^{\circ}C)$', fontsize=14)
	plt.xlabel(' ')
	plt.grid()
	plt.title(r'Station 27 - Average temperature (0-176m)', fontsize=14)
	fig.set_size_inches(w=12,h=6)
	outfile_meanT = 's27_meanT.png'
	fig.savefig(outfile_meanT, dpi=300)
	os.system('convert -trim ' + outfile_occu + ' ' + outfile_occu)



#########################################################################
#THIS IS THE START OF azmp_stn27_density.py associated FUNCTIONS

def density_calculator(df,file_location):
	#Define each of the necessary variables
	df_hydroT = df['temperature']
	df_hydroS = df['salinity']
	#Ensure that the measurement for each are present
	df_hydroS = df_hydroS.iloc[np.unique(df_hydroS.index.values,return_index=True)[1]]
	df_hydroT = df_hydroT.iloc[np.unique(df_hydroT.index.values,return_index=True)[1]]
	df_hydroS = df_hydroS.iloc[np.isin(df_hydroS.index.values,df_hydroT.index)]
	df_hydroT = df_hydroT.iloc[np.isin(df_hydroT.index.values,df_hydroS.index)]
	#Remove nans from each
	idx_S = pd.notnull(df_hydroS).any(axis=1).values.nonzero()[0]
	df_hydroT = df_hydroT.iloc[idx_S]
	df_hydroS = df_hydroS.iloc[idx_S]
	idx_T = pd.notnull(df_hydroT).any(axis=1).values.nonzero()[0]
	df_hydroT = df_hydroT.iloc[idx_T]
	df_hydroS = df_hydroS.iloc[idx_T]
	df_temp = df_hydroT.copy()
	df_sal = df_hydroS.copy()
	#Compute density
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
	#Quick QA/QC check
	idx_to_remove = []
	idx_to_sort = []
	for i in np.arange(1,13):
		#Isolate the densities for each month
		df_tmp = df_rho[df_rho.index.month == i]
		for idx in df_tmp.index:
			#Remove casts with one measurement
			if df_tmp.loc[idx].dropna().size == 1:
				idx_to_remove.append(idx)
				df_tmp = df_tmp.drop(idx)
			#Check to see if the cast is monotonic
			elif np.sum(~(np.diff(df_tmp.loc[idx].dropna().values.squeeze())>0)) > 5:
				idx_to_remove.append(idx)
				df_tmp = df_tmp.drop(idx)
			elif np.sum(~(np.diff(df_tmp.loc[idx].dropna().values.squeeze())>0)) > 0:
				idx_to_sort.append(idx)
				df_tmp = df_tmp.drop(idx)
			else:
				continue

	#Remove idx_to_remove
	df_rho.drop(idx_to_remove, inplace=True)
	df_sig.drop(idx_to_remove, inplace=True)
	df_SA.drop(idx_to_remove, inplace=True)
	df_CT.drop(idx_to_remove, inplace=True)

	#Update df_rho and df_sig with sorted values
	for i,value in enumerate(idx_to_sort):
		raw_data = df_rho.loc[value].values
		values_sorted = np.sort(raw_data[~np.isnan(raw_data)])
		raw_data[~np.isnan(raw_data)] = values_sorted
		df_rho.loc[value] = raw_data
		raw_data = df_sig.loc[value].values
		values_sorted = np.sort(raw_data[~np.isnan(raw_data)])
		raw_data[~np.isnan(raw_data)] = values_sorted
		df_sig.loc[value] = raw_data

	ds = xr.open_dataset(file_location)
	ds = ds.isel(time = ~np.isin(ds.time, np.array(idx_to_remove, dtype='datetime64')))
	#Vertically interpolate density
	df_SA = df_SA.interpolate(method='linear',axis=1).where(df_SA.bfill(axis=1).notnull())
	df_CT = df_CT.interpolate(method='linear',axis=1).where(df_CT.bfill(axis=1).notnull())
	df_rho = df_rho.interpolate(method='linear',axis=1).where(df_rho.bfill(axis=1).notnull())
	df_sig = df_sig.interpolate(method='linear',axis=1).where(df_sig.bfill(axis=1).notnull())
	#Save each
	df_rho.to_pickle('S27_rho_raw.pkl')
	df_sig.to_pickle('S27_sig_raw.pkl')
	df_SA.to_pickle('S27_SA_raw.pkl')
	df_CT.to_pickle('S27_CT_raw.pkl')
	return df_rho,df_sig,df_SA,df_CT,ds,Z

def MLD_calculator(df_SA,df_CT,df_rho,Z):
	#Define each of the necessary variables
	SA = df_SA.values
	CT = df_CT.values
	N2 = np.full((df_rho.index.size, df_rho.columns.size-1), np.nan)
	MLD = np.full((df_rho.index.size), np.nan)
	#Cycle through each date
	for i,idx in enumerate(df_rho.index):
		N2_tmp, pmid = gsw.Nsquared(SA[i,:], CT[i,:], Z, 47)
		N2[i,:] =  N2_tmp
		N2_tmp[np.where((pmid<=10) | (pmid>=100))] = np.nan
		if ~np.isnan(N2_tmp).all():
			MLD[i] = pmid[np.nanargmax(N2_tmp)]
	#Save the MLD and N2
	df_N2 = pd.DataFrame(N2, index=df_CT.index, columns=pmid)
	df_N2.to_pickle('S27_N2_raw.pkl')
	MLD = pd.Series(MLD, index=df_CT.index)
	MLD.to_pickle('S27_MLD_raw.pkl')
	#Also save the monthly averages
	df_rho.resample('M').mean().to_pickle('S27_rho_monthly.pkl')
	MLD.resample('M').mean().to_pickle('S27_MLD_monthly.pkl')

def stratification_calculator(file_location,lower_limit={'0-50':[0,5],'10-150':[8,13]},upper_limit={'0-50':[48,53],'10-150':[148,153]}):
	#Here, we need to determine density straight from the netcdf
	ds = xr.open_dataset(file_location)
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
	#Isolate for the relevant depths
	for i in lower_limit:
		warnings.simplefilter("ignore", category=RuntimeWarning)
		indx = [np.nanargmin(np.abs(ds.level - ii)) for ii in lower_limit[i]]
		rho_top = np.nanmean(rho.iloc[:,indx[0]:indx[1]], axis=1)
		indx = [np.nanargmin(np.abs(ds.level - ii)) for ii in upper_limit[i]]
		rho_bottom = np.nanmean(rho.iloc[:,indx[0]:indx[1]], axis=1)
		#Remove values where top density is higher than bottom density
		unstable_filt = rho_bottom <= rho_top
		rho_top[unstable_filt] = np.nan
		rho_bottom[unstable_filt] = np.nan
		#Calculate the stratification
		depth_range = np.mean(upper_limit[i]) - np.mean(lower_limit[i])
		strat = pd.Series((rho_bottom - rho_top)/depth_range, index=ds.time).resample('M').mean()
		strat.to_pickle('S27_stratif_'+i+'_raw.pkl')
		strat.resample('M').mean().to_pickle('S27_stratif_'+i+'_monthly.pkl')



#########################################################################
#THIS IS THE START OF azmp_stn27_analysis.py associated FUNCTIONS


#Some parameters that might be needed
#year_clim = [1991, 2020]
#current_year = 2023
#XLIM = [datetime.date(1945, 1, 1), datetime.date(2023, 12, 31)]
#french_months = ['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D']
#years_flag = [1950,1980]
#df = pd.read_pickle('S27_temperature_monthly.pkl')

def anomaly_calculator(df,years_flag,current_year,year_clim):
	#Flag the bad years
	for y in years_flag:
		df[df.index.year == y] = np.nan
	#Limit to the 1950s and after
	df = df[df.index.year>=1950]
	df = df[df.index.year<current_year+1]
	#Ensure that there is data for the whole year
	if df.index.size%12 != 0:
		last_month = df.index.month[-1]+1
		missing_months = 12 - df.index.size%12
		#Cycle through each of the missing months and set to nans
		for i in np.arange(last_month,last_month+missing_months):
			pd.options.mode.chained_assignment = None
			df.loc[pd.to_datetime(str(current_year)+'-'+"%.2d" % i+'-15')] = df.loc[str(current_year)+'-01-15']*np.nan
	#Determine the vertically averaged temperature
	ts_stack = df.groupby([(df.index.year),(df.index.month)]).mean()
	ts_stack.index = ts_stack.index.set_names(['year', 'month'])
	#Isolate, calculate the period of climatology
	df_clim_period = df[(df.index.year>=year_clim[0]) & (df.index.year<=year_clim[1])]
	monthly_clim_mean = df_clim_period.groupby(df_clim_period.index.month).mean()
	monthly_clim_stdv = df_clim_period.groupby(df_clim_period.index.month).std()
	ts_monthly_clim = monthly_clim_mean.mean(axis=1)
	ts_monthly_std = monthly_clim_stdv.mean(axis=1)
	#Tile the climatology for however many years are present
	years = len(df.index.year.unique())
	monthly_clim_mean = pd.concat([monthly_clim_mean] * years)
	monthly_clim_stdv = pd.concat([monthly_clim_stdv] * years)
	monthly_clim_mean.set_index(ts_stack.index, inplace=True)
	monthly_clim_stdv.set_index(ts_stack.index, inplace=True)
	#Calculate the anomalies
	anom_mi = ts_stack-monthly_clim_mean
	std_anom_mi = (ts_stack-monthly_clim_mean) / monthly_clim_stdv
	anom_monthly = anom_mi.mean(axis=1)
	std_anom_monthly = std_anom_mi.mean(axis=1)
	monthly_stdanom = std_anom_monthly.unstack()
	monthly_anom = anom_monthly.unstack()
	#Calculate the annual anomalies
	anom_std = monthly_stdanom.mean(axis=1)
	anom_std.index = pd.to_datetime(anom_std.index, format='%Y')
	anom = monthly_anom.mean(axis=1) 
	anom.index = pd.to_datetime(anom.index, format='%Y')
	#Annual mean is given by annual anomaly + monthly clim
	annual_mean = anom + ts_monthly_clim.mean()
	return df,anom_std,anom,annual_mean,ts_monthly_clim


def anomaly_plotter(anom_std,variable,XLIM,YLIM=[-3,3]):
	#Plot the vertically averaged anomaly
	df1 = anom_std[anom_std>0]
	df2 = anom_std[anom_std<0]
	fig = plt.figure(1)
	fig.clf()
	width = 200
	p1 = plt.bar(df1.index, np.squeeze(df1.values), width, alpha=0.8, color='indianred', zorder=10)
	p2 = plt.bar(df2.index, np.squeeze(df2.values), width, bottom=0, alpha=0.8, color='steelblue', zorder=10)
	plt.fill_between([anom_std.index[0], anom_std.index[-1]+np.timedelta64(365,'D')], [-.5, -.5], [.5, .5], facecolor='gray', alpha=.2)
	plt.ylabel('Normalized anomaly')
	plt.title('Station 27 - Average '+variable+' (0-176m)')
	plt.ylim(YLIM)
	plt.xlim(XLIM)
	plt.grid()
	# Save Figure
	fig.set_size_inches(w=12,h=6)
	if variable == 'salinity':
		fig_name = 's27_vert_sal_anomaly.png'
	elif variable == 'temperature':
		fig_name = 's27_vert_temp_anomaly.png'
	fig.savefig(fig_name, dpi=300)
	os.system('convert -trim ' + fig_name + ' ' + fig_name)
	# Save French Figure
	plt.ylabel(u'Anomalie standardizée')
	if variable == 'salinity':
		plt.title(u'Station 27 - Salinité moyenne (0-176m)')
		fig_name = 's27_vert_sal_anomalyFR.png'
	elif variable == 'temperature':
		plt.title(u'Station 27 - Température moyenne (0-176m)')
		fig_name = 's27_vert_temp_anomalyFR.png'
	fig.savefig(fig_name, dpi=300)
	os.system('convert -trim ' + fig_name + ' ' + fig_name)


def climatology_plotter(ts_monthly_clim,annual_mean,variable):
	#Create a vertical average climatology figure
	fig = plt.figure(1)
	fig.clf()
	annual_mean.plot()
	plt.plot([annual_mean.index[0]-np.timedelta64(365,'D'),annual_mean.index[-1]+np.timedelta64(365,'D')],[ts_monthly_clim.mean(), ts_monthly_clim.mean()], '--k', linewidth=3)
	plt.fill_between([annual_mean.index[0]-np.timedelta64(365,'D'),annual_mean.index[-1]+np.timedelta64(365,'D')], [ts_monthly_clim.mean()+annual_mean.std(), ts_monthly_clim.mean()+annual_mean.std()], [ts_monthly_clim.mean()-annual_mean.std(), ts_monthly_clim.mean()-annual_mean.std()], facecolor='gray', alpha=.2)
	plt.ylabel(r'Mean '+variable+' ($^\circ$C)')
	plt.title('Station 27 - Average '+variable+' (0-176m)')
	#plt.xlim(XLIM)
	plt.ylim([-.5, 1.75])
	plt.grid()
	# Save Figure
	fig.set_size_inches(w=12,h=6)
	if variable == 'temperature':
		fig_name = 's27_vert_temp_annual_mean.png'
	elif variable == 'salinity':
		fig_name = 's27_vert_sal_annual_mean.png'
	fig.savefig(fig_name, dpi=300)
	os.system('convert -trim ' + fig_name + ' ' + fig_name)


def CIL_calculator(df):
	#Ensure that df is for temperature
	#Select the summer
	df_temp_summer = df[(df.index.month>=6) & (df.index.month<=8)]
	df_temp_summer = df_temp_summer.resample('As').mean()
	df_temp_summer = df_temp_summer[df_temp_summer.columns[(df_temp_summer.columns>=10)]] # remove top 10m to be sure
	df_temp_summer = df_temp_summer.interpolate(method='linear',axis=1).where(df_temp_summer.bfill(axis=1).notnull())
	#Set up empty variables
	cil_temp = np.full(df_temp_summer.index.shape, np.nan)
	cil_core = np.full(df_temp_summer.index.shape, np.nan)
	cil_coredepth = np.full(df_temp_summer.index.shape, np.nan)
	cil_thickness = np.full(df_temp_summer.index.shape, np.nan)
	zbin = df_temp_summer.columns[2] - df_temp_summer.columns[1]
	#Cycle through each year
	for idx,YEAR in enumerate(df_temp_summer.index.year):
		#Get single year
		tmp = df_temp_summer.iloc[idx]
		#CIL core
		cil_core[idx] = np.nanmin(tmp.values)
		#CIL core depth
		z_idx = np.where(tmp==np.nanmin(tmp.values))
		cil_coredepth[idx] = np.mean(tmp.index[z_idx])
		#CIL thickness
		CIL_idxs = np.squeeze(np.where(tmp<=0))
		cil_thickness[idx] =  np.size(CIL_idxs)*zbin
		#CIL temp (found a bug Jan 2023. loc was used instead of iloc, introducing a 10m shift)
		cil_temp[idx] =  tmp.iloc[CIL_idxs].mean()
	#Convert to pandas series    
	cil_temp = pd.Series(cil_temp, index=df_temp_summer.index)
	cil_core = pd.Series(cil_core, index=df_temp_summer.index)
	cil_coredepth = pd.Series(cil_coredepth, index=df_temp_summer.index)
	cil_thickness = pd.Series(cil_thickness, index=df_temp_summer.index)
	cil_summer_stats = pd.concat([cil_temp, cil_core, cil_coredepth, cil_thickness], axis=1, keys=['CIL temp', 'CIL core T', 'CIL core depth', 'CIL thickness'])
	cil_summer_stats.to_pickle('S27_CIL_summer_stats.pkl')
	return cil_temp,cil_core,cil_coredepth,cil_thickness

def CIL_plotter(cil_stat,title,title_FR,save_title,XLIM,year_clim=[1991,2020],YLIM=[-5,5]):
	#Define the variable and climatology
	my_ts = cil_stat
	my_clim = my_ts[(my_ts.index.year>=year_clim[0]) & (my_ts.index.year<=year_clim[1])]
	clim_mean = my_clim.mean()
	clim_std = my_clim.std()
	#Calculate the anomaly
	anom = (my_ts - clim_mean)
	anom_std = anom / clim_std
	df1 = anom_std[anom_std>0]
	df2 = anom_std[anom_std<0]
	#Plot the figure
	fig = plt.figure(1)
	fig.clf()
	width = 200
	p1 = plt.bar(df1.index, np.squeeze(df1.values), width, alpha=0.8, color='indianred', zorder=10)
	p2 = plt.bar(df2.index, np.squeeze(df2.values), width, bottom=0, alpha=0.8, color='steelblue', zorder=10)
	plt.fill_between([anom_std.index[0], anom_std.index[-1]+np.timedelta64(365,'D')], [-.5, -.5], [.5, .5], facecolor='gray', alpha=.2)
	plt.ylabel('Normalized anomaly')
	plt.title('Station 27 - '+title)
	plt.ylim(YLIM)
	plt.xlim(XLIM)
	plt.grid()
	#Save Figure
	fig.set_size_inches(w=12,h=6)
	fig_name = save_title+'.png'
	fig.savefig(fig_name, dpi=300)
	os.system('convert -trim ' + fig_name + ' ' + fig_name)    
	#Save French Figure
	plt.ylabel(u'Anomalie standardizée')
	plt.title(u'Station 27 - '+title_FR)
	fig_name = save_title+'FR.png'
	fig.savefig(fig_name, dpi=300)
	os.system('convert -trim ' + fig_name + ' ' + fig_name)

#TEST SCRIPT
#CIL_plotter(cil_temp,'CIL mean temperature','Température moyenne de la CIF','s27_CILtemp_anomaly')
#CIL_plotter(cil_core,'CIL core temperature','Température du coeur de la CIF','s27_CILcore_anomaly')
#CIL_plotter(cil_coredepth,'CIL core depth','Profondeur du coeur de la CIF','s27_CILcoredepth_anomaly')
#CIL_plotter(cil_thickness,'CIL thickness','Épaisseur de la CIF','s27_CILthickness_anomaly')





#strat_shallow_path = 'S27_stratif_0-50_monthly.pkl'
#strat_deep_path = 'S27_stratif_10-150_monthly.pkl'

def stratification_plotter(strat_shallow_path,strat_deep_path,years_flag,current_year,year_clim):
	#Import the stratification data
	strat_monthly_shallow = pd.read_pickle(strat_shallow_path)
	strat_monthly_deep = pd.read_pickle(strat_deep_path)
	#Flag the bad years
	for y in years_flag:
		strat_monthly_shallow[strat_monthly_shallow.index.year == y] = np.nan
		strat_monthly_deep[strat_monthly_deep.index.year == y] = np.nan
	#Limit to the 1950s and after
	strat_monthly_shallow = strat_monthly_shallow[strat_monthly_shallow.index.year>=1950]
	strat_monthly_shallow = strat_monthly_shallow[strat_monthly_shallow.index.year<current_year+1]
	strat_monthly_deep = strat_monthly_deep[strat_monthly_deep.index.year>=1950]
	strat_monthly_deep = strat_monthly_deep[strat_monthly_deep.index.year<current_year+1]
	#Determine the anomalies
	anom = {}
	anom_std = {}
	strat_monthly_clim = {}
	for i in ['shallow','deep']:
		if i == 'shallow':
			strat = strat_monthly_shallow
		elif i == 'deep':
			strat = strat_monthly_deep
		strat_stack = strat.groupby([(strat.index.year),(strat.index.month)]).mean()
		strat_unstack = strat_stack.unstack()
		#Compute clim
		strat_clim_period = strat[(strat.index.year>=year_clim[0]) & (strat.index.year<=year_clim[1])]
		strat_monthly_stack = strat_clim_period.groupby([(strat_clim_period.index.year),(strat_clim_period.index.month)]).mean()
		strat_monthly_clim[i] = strat_monthly_stack.groupby(level=1).mean()
		strat_monthly_std = strat_monthly_stack.groupby(level=1).std()
		monthly_anom = strat_unstack - strat_monthly_clim[i]
		monthly_stdanom = (strat_unstack - strat_monthly_clim[i]) /  strat_monthly_std
		#See if each year has enough data present
		nom = np.sum(~np.isnan(monthly_anom.values),axis=1)
		monthly_anom.iloc[nom<3] = np.nan
		monthly_stdanom.iloc[nom<3] = np.nan
		anom_std[i] = monthly_stdanom.mean(axis=1)
		anom_std[i].index = pd.to_datetime(anom_std[i].index, format='%Y')
		anom[i] = monthly_anom.mean(axis=1)*1000 
		anom[i].index = pd.to_datetime(anom[i].index, format='%Y')
	return strat_monthly_shallow,strat_monthly_deep,anom,anom_std,strat_monthly_clim


def stratification_barplot(anom_std):
	#Cycle through shallow and deep
	for i in ['shallow','deep']:
		df1 = anom_std[i][anom_std[i]>0]
		df2 = anom_std[i][anom_std[i]<0]
		fig = plt.figure(1)
		fig.clf()
		width = 200
		p1 = plt.bar(df1.index, np.squeeze(df1.values), width, alpha=0.8, color='indianred', zorder=10)
		p2 = plt.bar(df2.index, np.squeeze(df2.values), width, bottom=0, alpha=0.8, color='steelblue', zorder=10)
		plt.ylabel('Normalized anomaly')
		if i == 'shallow':
		    plt.title('Station 27 - Stratification (5-50m)')
		elif i == 'deep':
		    plt.title('Station 27 - Stratification (10-150m)')
		plt.ylim([-1.5, 1.5])
		plt.grid()
		# Save Figure
		fig.set_size_inches(w=12,h=6)
		fig_name = 's27_stratif_bar_'+i+'.png'
		fig.savefig(fig_name, dpi=300)
		os.system('convert -trim ' + fig_name + ' ' + fig_name)
		# Save French Figure
		plt.ylabel(u'Anomalie standardizée')
		fig_name = 's27_stratif_bar_'+i+'_FR.png'
		fig.savefig(fig_name, dpi=300)
		os.system('convert -trim ' + fig_name + ' ' + fig_name)

def stratification_timeseries(anom):
	#Cycle through shallow and deep
	for i in ['shallow','deep']:
		fig = plt.figure(1)
		plt.clf()
		anom[i].plot(color='gray')
		anom[i].rolling(5, min_periods=3).mean().plot(color='k', linewidth=3, linestyle='--')
		plt.grid()
		plt.ylabel(r'Stratification anomaly $\rm (g\,m^{-4})$')
		plt.xlabel(' ')
		#Save Figure
		fig.set_size_inches(w=12,h=6)
		fig_name = 's27_stratif_plot_'+i+'.png'
		fig.savefig(fig_name, dpi=300)
		os.system('convert -trim ' + fig_name + ' ' + fig_name)
		#Same in French
		plt.ylabel(r'Anomalie de stratification $\rm (g\,m^{-4})$')
		fig_name = 's27_stratif_plot_'+i+'_FR.png'
		fig.savefig(fig_name, dpi=300)
		os.system('convert -trim ' + fig_name + ' ' + fig_name)

def stratification_timeseries_mean(
	anom,
	strat_monthly_clim,
	year_start = [1950,1990],
	year_end = [1990,2022],
	markers = {'shallow': ['dotted','x'],'deep': ['--','.']},
	titles = {'shallow': '0-50m','deep': '10-150m'},
	y_range = {'shallow': [15,30], 'deep': [13,18]}):
	#New stratification figure with means
	strat_virtual = {}
	for i in ['shallow','deep']:
		strat_virtual[i] = anom[i] + strat_monthly_clim[i].mean()*1000
	#Plot the figure
	fig =plt.figure(5)
	plt.clf()
	x = 1
	for i in ['shallow','deep']:
		plt.subplot(2,1,x)
		strat_virtual[i].plot(color='gray', linestyle=' ', marker=markers[i][1], label=titles[i])
		strat_virtual[i].rolling(5, min_periods=3).mean().plot(color='k', linewidth=3, linestyle=markers[i][0], label=titles[i])
		for ii in np.arange(np.size(year_start)):
			mean = strat_virtual[i][(strat_virtual[i].index.year>=year_start[ii]) & (strat_virtual[i].index.year<=year_end[ii])].mean()
			stdv = strat_virtual[i][(strat_virtual[i].index.year>=year_start[ii]) & (strat_virtual[i].index.year<=year_end[ii])].std()
			plt.hlines(mean,
				xmin=np.datetime64(str(year_start[ii])+'-01','M'),
				xmax=np.datetime64(str(year_end[ii])+'-01','M'),
				color='tab:red', linewidth=1.5)
			plt.text(np.datetime64(str(np.mean([year_start[ii],year_end[ii]]).astype(int))+'-01','M'),
				mean+1,
				str(np.round(mean,2))+'$\pm$'+str(np.round(stdv,2))+' (gm$^{-4}$)',
				color='tab:red', fontweight='bold', horizontalalignment='center')
		plt.ylabel(r'$\frac{\Delta \sigma}{\Delta z}$ $\rm (g\,m^{-4})$')
		plt.ylim(y_range[i])
		plt.grid()
		plt.legend(loc=3, ncol=4)
		x += 1
	# Save Figure
	fig.set_size_inches(w=8,h=6)
	fig_name = 's27_stratif_plot_means_withlines.png'
	fig.savefig(fig_name, dpi=300)
	os.system('convert -trim ' + fig_name + ' ' + fig_name)


def stratification_currentyear_barplot(strat_monthly_shallow,strat_monthly_deep,current_year,year_clim=[1991,2020]):
	#Cycle through shallow and deep
	for i in ['shallow','deep']:
		if i == 'shallow':
			strat = strat_monthly_shallow
		elif i == 'deep':
			strat = strat_monthly_deep
		#Determine the climatology
		strat_clim_period = strat[(strat.index.year>=year_clim[0]) & (strat.index.year<=year_clim[1])]
		strat_monthly_stack = strat_clim_period.groupby([(strat_clim_period.index.year),(strat_clim_period.index.month)]).mean()
		strat_monthly_clim = strat_monthly_stack.groupby(level=1).mean()
		strat_monthly_std = strat_monthly_stack.groupby(level=1).std()
		#Fill in the current year if not complete
		year_index = pd.date_range(str(current_year)+"-01-31", periods=12, freq="M")
		strat_current_year = pd.Series(np.nan, index = year_index)
		strat_tmp = strat[strat.index.year==current_year]
		strat_current_year.loc[strat_tmp.index] = strat_tmp
		#Sort out the index
		strat_current_year.index=strat_current_year.index.month
		monthly_anom = strat_current_year - strat_monthly_clim
		monthly_std_anom = monthly_anom/strat_monthly_std
		monthly_std_anom.index = year_index.strftime('%b')
		monthly_anom.index = year_index.strftime('%b')
		#Plot the monthly anomaly
		ind = np.arange(len(monthly_anom.keys()))
		width = 0.35
		fig, ax = plt.subplots()
		rects1 = ax.bar(ind - width/2, strat_monthly_clim.values*1000, width, yerr=strat_monthly_std.values*.5*1000, label='1991-2020', color='lightblue', zorder=1)
		rects2 = ax.bar(ind + width/2, np.squeeze(strat_current_year.values*1000), width, yerr=None, label=str(current_year), color='steelblue', zorder=1)
		#Add some text for labels, title and custom x-axis tick labels, etc.
		ax.set_ylabel(r'$\rm \frac{\Delta \rho}{\Delta z} (g\,m^{-4})$')
		ax.set_title('Station 27 - Stratification')
		ax.set_xticks(ind)
		ax.set_xticklabels(monthly_anom.index)
		ax.legend()
		ax.yaxis.grid() # horizontal lines
		#Save Figure
		fig.set_size_inches(w=12,h=6)
		fig_name = 's27_stratif_monthly_'+i+'.png'
		fig.savefig(fig_name, dpi=300)
		os.system('convert -trim ' + fig_name + ' ' + fig_name)
		#Save French Figure
		french_months = ['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D']
		ax.set_xticklabels(french_months)
		fig_name = 's27_stratif_monthly_'+i+'_FR.png'
		fig.savefig(fig_name, dpi=300)
		os.system('convert -trim ' + fig_name + ' ' + fig_name)

#MLD_path = 'S27_MLD_monthly.pkl'

def MLD_processor(MLD_path,years_flag,year_clim,current_year):
	#Import the MLD data
	mld = pd.read_pickle(MLD_path)
	#Flag the bad years
	for y in years_flag:
		mld[mld.index.year == y] = np.nan
	#Limit to the 1950s and after
	mld = mld[mld.index.year>=1950]
	mld = mld[mld.index.year<current_year+1]
	#Stack the months
	mld_stack = mld.groupby([(mld.index.year),(mld.index.month)]).mean()
	mld_unstack = mld_stack.unstack()
	#Compute the climatology
	mld_clim_period = mld[(mld.index.year>=year_clim[0]) & (mld.index.year<=year_clim[1])]
	mld_monthly_stack = mld_clim_period.groupby([(mld_clim_period.index.year),(mld_clim_period.index.month)]).mean()
	mld_monthly_clim = mld_monthly_stack.groupby(level=1).mean()
	mld_monthly_std = mld_monthly_stack.groupby(level=1).std()
	#Compute the anomaly
	monthly_anom = mld_unstack - mld_monthly_clim 
	monthly_stdanom = (mld_unstack - mld_monthly_clim)/mld_monthly_std
	anom_std = monthly_stdanom.mean(axis=1) 
	anom_std.index = pd.to_datetime(anom_std.index, format='%Y')
	anom = monthly_anom.mean(axis=1)
	anom.index = pd.to_datetime(anom.index, format='%Y')
	#Save the MLD anomaly
	anom_std_mld = anom_std.copy()
	anom_std_mld.index = anom_std_mld.index.year
	anom_std_mld.to_csv('MLD_std_anom.csv', sep=',', float_format='%0.3f')
	anom_std_mld.rolling(5, min_periods=3).mean().to_csv('MLD_std_anom_5Y.csv', sep=',', float_format='%0.3f')
	return mld,anom,anom_std


def MLD_barplot(anom_std):
	df1 = anom_std[anom_std>0]
	df2 = anom_std[anom_std<0]
	fig = plt.figure(1)
	fig.clf()
	width = 200
	p1 = plt.bar(df1.index, np.squeeze(df1.values), width, alpha=0.8, color='indianred', zorder=10)
	p2 = plt.bar(df2.index, np.squeeze(df2.values), width, bottom=0, alpha=0.8, color='steelblue', zorder=10)
	plt.ylabel(r'Normalized anomaly')
	plt.title('Station 27 - Mixed layer depth')
	#plt.xlim(XLIM)
	plt.ylim([-3, 3])
	plt.grid()
	fig.set_size_inches(w=12,h=6)
	fig_name = 's27_mld_bar.png'
	fig.savefig(fig_name, dpi=300)
	os.system('convert -trim ' + fig_name + ' ' + fig_name)
	# French Figure
	plt.ylabel(u'Anomalie standardizée')
	plt.title(u'Station 27 - Profondeur couche de mélange')
	fig_name = 's27_mld_barFR.png'
	fig.savefig(fig_name, dpi=300)
	os.system('convert -trim ' + fig_name + ' ' + fig_name)

def MLD_timeseries(anom):
	fig = plt.figure(1)
	plt.clf()
	anom.plot(color='gray')
	anom.rolling(5, min_periods=3).mean().plot(color='k', linewidth=3, linestyle='--')
	plt.grid()
	plt.ylabel(r'MLD anomaly (m)')
	plt.xlabel(' ')
	fig.set_size_inches(w=12,h=6)
	fig_name = 's27_mld_plot.png'
	fig.savefig(fig_name, dpi=200)
	os.system('convert -trim ' + fig_name + ' ' + fig_name)

	# French Figure
	plt.ylabel(r'Anomalie de la PCM (m)')
	fig_name = 's27_mld_plotFR.png'
	fig.savefig(fig_name, dpi=200)
	os.system('convert -trim ' + fig_name + ' ' + fig_name)

def MLD_currentyear_barplot(mld,current_year,year_clim=[1991,2020]):
	#Determine the climatology
	mld_clim_period = mld[(mld.index.year>=year_clim[0]) & (mld.index.year<=year_clim[1])]
	mld_monthly_stack = mld_clim_period.groupby([(mld_clim_period.index.year),(mld_clim_period.index.month)]).mean()
	mld_monthly_clim = mld_monthly_stack.groupby(level=1).mean()
	mld_monthly_std = mld_monthly_stack.groupby(level=1).std()
	mld_current_year = mld[mld.index.year==current_year]
	#Fill in the current year if not complete
	year_index = pd.date_range(str(current_year)+"-01-31", periods=12, freq="M")
	mld_current_year = pd.Series(np.nan, index = year_index)
	mld_tmp = mld[mld.index.year==current_year]
	mld_current_year.loc[mld_tmp.index] = mld_tmp
	#Sort out the index
	mld_current_year.index=mld_monthly_std.index
	monthly_anom = mld_current_year - mld_monthly_clim
	monthly_std_anom = monthly_anom/mld_monthly_std
	monthly_std_anom.index = year_index.strftime('%b')
	monthly_anom.index = year_index.strftime('%b')
	#Plot the monthly anomaly
	ind = np.arange(len(monthly_anom.keys()))  # the x locations for the groups
	width = 0.35  # the width of the bars
	fig, ax = plt.subplots()
	rects1 = ax.bar(ind - width/2, mld_monthly_clim.values, width, yerr=mld_monthly_std.values*.5,
		label='1991-2020', color='lightblue', zorder=1)
	rects2 = ax.bar(ind + width/2, np.squeeze(mld_current_year.values), width, yerr=None,
		label=str(current_year), color='steelblue', zorder=1)
	#Add some text for labels, title and custom x-axis tick labels, etc.
	ax.set_ylabel(r'MLD (m)')
	ax.set_title('Station 27 - Mixed layer depth')
	ax.set_xticks(ind)
	ax.set_xticklabels(monthly_anom.index)
	ax.legend()
	ax.yaxis.grid() # horizontal lines
	fig.set_size_inches(w=12,h=6)
	fig_name = 's27_mld_monthly.png'
	fig.savefig(fig_name, dpi=300)
	os.system('convert -trim ' + fig_name + ' ' + fig_name)
	#French Figure
	ax.set_ylabel(r'PCM (m)')
	ax.set_title(u'Station 27 - Profondeur couche de mélange')
	french_months = ['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D']
	ax.set_xticklabels(french_months)
	fig_name = 's27_mld_monthly_FR.png'
	fig.savefig(fig_name, dpi=300)
	os.system('convert -trim ' + fig_name + ' ' + fig_name)


#########################################################################
#THIS IS THE START OF azmp_stn27_scorecards.py associated FUNCTIONS

def scorecard_TS_processor(file_location,year_clim,year_plot,variable):
	#Import the data and isolate to post-1950
	df_var = pd.read_pickle(file_location)
	df_var = df_var[df_var.index.year>=1950]
	#Calculate the anomaly
	#Isolate variable by year, month
	var_stack = df_var.groupby([(df_var.index.year),(df_var.index.month)]).mean()
	var_stack.index = var_stack.index.set_names(['year','month'])
	#Isolate the climatology period
	var_clim_period = df_var[(df_var.index.year>=year_clim[0]) & (df_var.index.year<=year_clim[1])]
	monthly_clim_mean = var_clim_period.groupby(var_clim_period.index.month).mean()
	monthly_clim_stdv = var_clim_period.groupby(var_clim_period.index.month).std()
	#Tile the climatology for the number of years
	var_years = len(df_var.index.year.unique())
	monthly_clim_mean = pd.concat([monthly_clim_mean]*var_years)
	monthly_clim_stdv = pd.concat([monthly_clim_stdv]*var_years)
	monthly_clim_mean = monthly_clim_mean.iloc[:var_stack.index.size,:]
	monthly_clim_stdv = monthly_clim_stdv.iloc[:var_stack.index.size,:]
	#Set multi-index to clim (using ts_stack index)
	monthly_clim_mean.set_index(var_stack.index, inplace=True)
	monthly_clim_stdv.set_index(var_stack.index, inplace=True)
	#Calculate anomalies (mi = miltiindex)
	anom_mi = var_stack-monthly_clim_mean
	std_anom_mi = (var_stack-monthly_clim_mean)/monthly_clim_stdv
	#Vertically-averaged anomalies to three different groups
	var_0_btm_anom_monthly = anom_mi.mean(axis=1)
	var_0_btm_std_anom_monthly = std_anom_mi.mean(axis=1)
	var_0_50_anom_monthly = anom_mi[anom_mi.columns[(anom_mi.columns<=50)]]
	var_0_50_std_anom_monthly = std_anom_mi[std_anom_mi.columns[(std_anom_mi.columns<=50)]]
	var_150_btm_anom_monthly = anom_mi[anom_mi.columns[(anom_mi.columns>=150)]]
	var_150_btm_std_anom_monthly = std_anom_mi[std_anom_mi.columns[(std_anom_mi.columns>=150)]]
	#Convert to yearly
	ave_anom = var_0_btm_anom_monthly.unstack().mean(axis=1)
	ave_anom_std = var_0_btm_std_anom_monthly.unstack().mean(axis=1)
	surf_anom = var_0_50_anom_monthly.unstack().mean(axis=1)
	surf_anom_std = var_0_50_std_anom_monthly.unstack().mean(axis=1)
	bot_anom = var_150_btm_anom_monthly.unstack().mean(axis=1)
	bot_anom_std = var_150_btm_std_anom_monthly.unstack().mean(axis=1)
	ave_anom.index = pd.to_datetime(ave_anom.index, format='%Y')
	surf_anom.index = pd.to_datetime(surf_anom.index, format='%Y')
	bot_anom.index = pd.to_datetime(bot_anom.index, format='%Y')
	ave_anom_std.index = pd.to_datetime(ave_anom_std.index, format='%Y')
	surf_anom_std.index = pd.to_datetime(surf_anom_std.index, format='%Y')
	bot_anom_std.index = pd.to_datetime(bot_anom_std.index, format='%Y')
	#Merge the anomalies
	anom = pd.concat([ave_anom,surf_anom,bot_anom], axis=1, keys=[variable+' 0-176m',variable+' 0-50m',variable+' 150-176m'])
	anom_std = pd.concat([ave_anom_std,surf_anom_std,bot_anom_std], axis=1, keys=[variable+' 0-176m',variable+' 0-50m',variable+' 150-176m'])
	anom = anom[anom.index.year>=1947]
	anom_std = anom_std[anom_std.index.year>=1947]
	#Vertically average plot years
	var_0_btm_clim_monthly = var_stack.mean(axis=1).unstack()
	var_0_50_clim_monthly = var_stack[var_stack.columns[(var_stack.columns<=50)]].unstack()
	var_150_btm_clim_monthly = var_stack[var_stack.columns[(var_stack.columns>=150)]].unstack()
	#Convert to yearly
	var_0_btm_clim = var_0_btm_clim_monthly[var_0_btm_clim_monthly.index>=year_plot[0]].mean(axis=1)
	var_0_50_clim = var_0_50_clim_monthly[var_0_50_clim_monthly.index>=year_plot[0]].mean(axis=1)
	var_150_btm_clim = var_150_btm_clim_monthly[var_150_btm_clim_monthly.index>=year_plot[0]].mean(axis=1)
	#Save pkl for climate indices (whole timeseries)
	anom_std.to_pickle('s27_'+variable.lower()+'_std_anom.pkl')
	#Keep only relevant window
	anom_std = anom_std[anom_std.index.year>=year_plot[0]]
	#Annual clims (for inclusion)
	clim_period_annual = pd.concat([var_0_btm_clim,var_0_50_clim,var_150_btm_clim], axis=1, keys=[variable+' 0-176m',variable+' 0-50m',variable+' 150-176m'])
	return anom,anom_std,clim_period_annual

def scorecard_CIL_processor(file_location,year_clim,year_plot):
	#Load the saved data
	df_CIL = pd.read_pickle(file_location)
	df_CIL = df_CIL[df_CIL.index.year>=year_plot[0]]
	#Compute anomalies 
	CIL_clim_period = df_CIL[(df_CIL.index.year>=year_clim[0]) & (df_CIL.index.year<=year_clim[1])]
	#Anomaly
	CIL_anom = (df_CIL - CIL_clim_period.mean())
	CIL_anom_std = CIL_anom / CIL_clim_period.std()
	return CIL_anom,CIL_anom_std,CIL_clim_period

def scorecard_MLD_processor(file_location,year_clim,year_plot):
	#Load pickled data
	mld_monthly = pd.read_pickle(file_location)
	mld = mld_monthly
	#Flag some of the data
	mld[mld.index=='2019-03-15']=np.nan
	mld[mld.index=='1980-03-15']=np.nan
	mld = mld[mld.index.year>=year_plot[0]]
	#Stack the months
	mld_stack = mld.groupby([(mld.index.year),(mld.index.month)]).mean()
	mld_unstack = mld_stack.unstack()
	#Compute the climatology
	mld_clim_period = mld[(mld.index.year>=year_clim[0]) & (mld.index.year<=year_clim[1])]
	mld_monthly_stack = mld_clim_period.groupby([(mld_clim_period.index.year),(mld_clim_period.index.month)]).mean()
	mld_monthly_clim = mld_monthly_stack.groupby(level=1).mean()
	mld_monthly_std = mld_monthly_stack.groupby(level=1).std()
	monthly_anom = mld_unstack - mld_monthly_clim 
	monthly_stdanom = (mld_unstack - mld_monthly_clim)/mld_monthly_std
	#Seasonal and annual anomalies
	mld_winter_anom = monthly_stdanom[[1,2,3]].mean(axis=1)
	mld_spring_anom = monthly_stdanom[[4,5,6]].mean(axis=1)
	mld_summer_anom = monthly_stdanom[[7,8,9]].mean(axis=1)
	mld_fall_anom = monthly_stdanom[[10,11,12]].mean(axis=1)
	mld_annual_anom = monthly_stdanom.mean(axis=1)
	#Mean values for climatology period (for last columns)
	mld_winter_clim = mld_monthly_clim[[1,2,3]]
	mld_spring_clim = mld_monthly_clim[[4,5,6]]
	mld_summer_clim = mld_monthly_clim[[7,8,9]]
	mld_fall_clim = mld_monthly_clim[[10,11,12]]
	mld_annual_clim = mld_monthly_clim
	#Concat clim and anomaly
	mld_clim = pd.concat([mld_winter_clim, mld_spring_clim, mld_summer_clim, mld_fall_clim, mld_annual_clim], axis=1, keys=['MLD winter', 'MLD spring', 'MLD summer', 'MLD fall', 'MLD annual'])
	mld_anom_std = pd.concat([mld_winter_anom, mld_spring_anom, mld_summer_anom, mld_fall_anom, mld_annual_anom], axis=1, keys=['MLD winter', 'MLD spring', 'MLD summer', 'MLD fall', 'MLD annual'])
	return mld_clim,mld_anom_std

def scorecard_strat_processor(file_location,year_clim,year_plot):
	#Load the pickled data
	strat_monthly = pd.read_pickle(file_location)
	strat = strat = strat_monthly
	#Flag some data
	strat[strat.index=='2019-03-15']=np.nan
	strat[strat.index=='1980-03-15']=np.nan
	strat = strat[strat.index.year>=year_plot[0]]
	#Stack months
	strat_stack = strat.groupby([(strat.index.year),(strat.index.month)]).mean()
	strat_unstack = strat_stack.unstack()
	#Compute clim
	strat_clim_period = strat[(strat.index.year>=year_clim[0]) & (strat.index.year<=year_clim[1])]
	strat_monthly_stack = strat_clim_period.groupby([(strat_clim_period.index.year),(strat_clim_period.index.month)]).mean()
	strat_monthly_clim = strat_monthly_stack.groupby(level=1).mean()
	strat_monthly_std = strat_monthly_stack.groupby(level=1).std()
	monthly_anom = strat_unstack - strat_monthly_clim
	monthly_stdanom = (strat_unstack - strat_monthly_clim)/strat_monthly_std
	nom = np.sum(~np.isnan(monthly_anom.values),axis=1)
	#Seasonal and annual anomalies
	strat_winter_anom = monthly_stdanom[[1,2,3]].mean(axis=1)
	strat_spring_anom = monthly_stdanom[[4,5,6]].mean(axis=1)
	strat_summer_anom = monthly_stdanom[[7,8,9]].mean(axis=1)
	strat_fall_anom = monthly_stdanom[[10,11,12]].mean(axis=1)
	monthly_stdanom.iloc[nom < 3] = np.nan
	strat_annual_anom = monthly_stdanom.mean(axis=1)
	#Mean values for climatology period (for last columns)
	strat_winter_clim = strat_monthly_clim[[1,2,3]]
	strat_spring_clim = strat_monthly_clim[[4,5,6]]
	strat_summer_clim = strat_monthly_clim[[7,8,9]]
	strat_fall_clim = strat_monthly_clim[[10,11,12]]
	strat_annual_clim = strat_monthly_clim
	#Concat clim and anomaly
	strat_clim = pd.concat([strat_winter_clim, strat_spring_clim, strat_summer_clim, strat_fall_clim, strat_annual_clim], axis=1, keys=['strat winter', 'strat spring', 'strat summer', 'strat fall', 'strat annual'])
	strat_anom_std = pd.concat([strat_winter_anom, strat_spring_anom, strat_summer_anom, strat_fall_anom, strat_annual_anom], axis=1, keys=['strat winter', 'strat spring', 'strat summer', 'strat fall', 'strat annual'])
	return strat_clim,strat_anom_std

def scorecard_plotter(var_data,var_data_clim,collabel,collabel_FR,savename,years_present=True):
	#Build the colormap
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
	cmap_r, norm_r = from_levels_and_colors(levels, np.flipud(colors), extend='both')
	#Common parameters
	hcell,wcell = 0.5, 0.7
	hpad,wpad = .8, 0

	#Define the variable
	if savename=='CIL':
		my_df = var_data.drop('CIL core depth', axis=1).T
	else:
		my_df = var_data.T
	my_df['MEAN'] = var_data_clim.mean()
	my_df['SD'] = var_data_clim.std()
	#Get text values + cell color
	vals = np.around(my_df.values,1)
	vals[vals==-0.] = 0.
	vals_color = vals.copy()
	vals_color[:,-1] = 0 #No color to last two columns (mean and STD)
	vals_color[:,-2] = 0
	#Custom colour flips
	if savename=='CIL':
		vals_color[2,:] = vals_color[2,:]*-1
	if savename.startswith('strat'):
		for i in np.arange(0,5):
			vals[i,-1] = np.around(my_df.values[i,-1],3)
			vals[i,-2] = np.around(my_df.values[i,-2],3)
	#Set up the rows and columns (+1 for years)
	if years_present:
		#Preamble, determine the years
		year_list = var_data.index.year.astype('str')
		year_list = [i[2:4] for i in year_list] #2-digit year
		year_list.append(r'$\rm \overline{x}$') #add 2 extra columns
		year_list.append(r'sd')
		year_list_FR = year_list[:]
		year_list_FR[-1] = u'ET'
		nrows, ncols = my_df.index.size+1, my_df.columns.size
	else:
		nrows, ncols = my_df.index.size, my_df.columns.size
	fig=plt.figure(figsize=(ncols*wcell+wpad, nrows*hcell+hpad)) # size+1 because of year's row
	ax = fig.add_subplot(111)
	ax.axis('off')
	#Plot the table
	header = ax.table(cellText=[['']],
		colLabels=[collabel],
		loc='center'
		)
	header.set_fontsize(13)
	if years_present:
		the_table=ax.table(cellText=vals, rowLabels=my_df.index, colLabels=year_list,
			loc='center', cellColours=cmap(norm(vals_color)), cellLoc='center',
			bbox=[0, 0, 1, 0.5]
			)
	else:
		the_table=ax.table(cellText=vals, rowLabels=my_df.index, colLabels=None,
			loc='center', cellColours=cmap(norm(vals_color)), cellLoc='center',
			bbox=[0, 0, 1, 0.5]
			)
	#Change font color to white where needed:
	the_table.auto_set_font_size(False)
	the_table.set_fontsize(13)
	table_props = the_table.properties()
	last_columns = np.arange(vals.shape[1]-2, vals.shape[1])
	if years_present:
		for key, cell in the_table.get_celld().items():
			cell_text = cell.get_text().get_text() 
			if is_number(cell_text) == False:
				pass
			elif key[0] == 0: #year's row = no color
				pass
			elif key[1] in last_columns:
				cell._text.set_color('darkslategray')
			elif (float(cell_text) <= -2) | (float(cell_text) >= 2) :
				cell._text.set_color('white')
			elif (cell_text=='nan'):
				cell._set_facecolor('darkgray')
				cell._text.set_color('darkgray')
	else:
		for key, cell in the_table.get_celld().items():
			cell_text = cell.get_text().get_text() 
			if is_number(cell_text) == False:
				pass
			elif key[1] in last_columns:
				cell._text.set_color('darkslategray')
			elif (float(cell_text) <= -2) | (float(cell_text) >= 2) :
				cell._text.set_color('white')
			elif (cell_text=='nan'):
				cell._set_facecolor('darkgray')
				cell._text.set_color('darkgray')

	#Save the figure, english
	plt.savefig('scorecards_s27_'+savename+'.png', dpi=300)
	os.system('convert -trim scorecards_s27_'+savename+'.png scorecards_s27_'+savename+'.png')

	#French table
	header = ax.table(cellText=[['']],
		colLabels=[collabel_FR],
		loc='center'
		)
	header.set_fontsize(13)
	if years_present:
		the_table=ax.table(cellText=vals, rowLabels=my_df.index, colLabels=year_list_FR,
			loc='center', cellColours=cmap(norm(vals_color)), cellLoc='center',
			bbox=[0, 0, 1, 0.5]
			)
	else:
		the_table=ax.table(cellText=vals, rowLabels=my_df.index, colLabels=None,
			loc='center', cellColours=cmap(norm(vals_color)), cellLoc='center',
			bbox=[0, 0, 1, 0.5]
			)
	#Change font color to white where needed:
	the_table.auto_set_font_size(False)
	the_table.set_fontsize(13)
	table_props = the_table.properties()
	last_columns = np.arange(vals.shape[1]-2, vals.shape[1])
	if years_present:
		for key, cell in the_table.get_celld().items():
			cell_text = cell.get_text().get_text() 
			if is_number(cell_text) == False:
				pass
			elif key[0] == 0: #year's row = no color
				pass
			elif key[1] in last_columns:
				cell._text.set_color('darkslategray')
			elif (float(cell_text) <= -2) | (float(cell_text) >= 2) :
				cell._text.set_color('white')
			elif (cell_text=='nan'):
				cell._set_facecolor('darkgray')
				cell._text.set_color('darkgray')
	else:
		for key, cell in the_table.get_celld().items():
			cell_text = cell.get_text().get_text() 
			if is_number(cell_text) == False:
				pass
			elif key[1] in last_columns:
				cell._text.set_color('darkslategray')
			elif (float(cell_text) <= -2) | (float(cell_text) >= 2) :
				cell._text.set_color('white')
			elif (cell_text=='nan'):
				cell._set_facecolor('darkgray')
				cell._text.set_color('darkgray')

	#Save the figure, french
	plt.savefig('scorecards_s27_'+savename+'_FR.png', dpi=300)
	os.system('convert -trim scorecards_s27_'+savename+'_FR.png scorecards_s27_'+savename+'_FR.png')


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


