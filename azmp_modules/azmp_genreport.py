'''
Module azmp_genreport.py 

A series of function to generate figures for AZMP ResDoc and presentation on annual conditions.

Contains following functions:
- nao(YEAR
- ao(YEAR)
- amo(YEAR)

-

** For an example of call of the different functions, check azmp_genreportYYYY.py

----------

Atlantic Zone Monitoring Program @NAFC:
https://azmp-nl.github.io/

'''

__author__ = 'Frederic.Cyr@dfo-mpo.gc.ca'
__version__ = '0.1'

import os
import re
import numpy as np
import pandas as pd
import sys
import netCDF4
from sys import version_info
import xarray as xr
import h5py
import matplotlib.pyplot as plt

#from scipy.interpolate import griddata
#from scipy.interpolate import interp1d  # to remove NaNs in profiles
#from shapely.geometry import Point
#from shapely.geometry.polygon import Polygon
#from shapely.ops import cascaded_union
#from area import area # external fns to compute surface area
#from seawater import extras as swx
# maps
os.environ['PROJ_LIB'] = '/home/cyrf0006/anaconda3/share/proj'
#from mpl_toolkits.basemap import Basemap

    
def nao(YEAR,nao_file_loc,url_loc='https://www.cpc.ncep.noaa.gov/products/precip/CWlink/pna/norm.nao.monthly.b5001.current.ascii.table'):
    '''
    Update the winter North Atlantic Oscillation
    run in /home/cyrf0006/AZMP/state_reports/airTemp

    Generate data:
    NAO_annual.pkl
    NAO_winter.pkl
    NAO_summer.pkl

    Generate Figures:
    NAO_winter_1950-YYYY.png'
    NAO_winter_1950-YYYY_FR.png'

    Adapted from nao.py
    '''

    # Adjust fontsize/weight
    font = {'family' : 'normal',
            'weight' : 'bold',
            'size'   : 18}
    plt.rc('font', **font)

    # Download and save up-to-date  NAO index from NOAA (data.csv) if needed
    #url = 'https://www.ncdc.noaa.gov/teleconnections/nao/data.csv' (until 2020...)
    url = url_loc
    nao_file = nao_file_loc
    import urllib3
    http = urllib3.PoolManager()
    r = http.request('GET', url)

    #Split the data
    r_data = str(r.data)[2:-1].split('\\n')
    months = re.split(r'\s{2,}', r_data[0])
    data = np.array([re.split(r'\s{2,}', i) for i in r_data[1:-1]])
    years = data[:,0]
    data = data[:,1:].astype(float)
    pd_data = pd.DataFrame(data, columns=months[1:], index=years)
    pd_data.to_csv(nao_file_loc)

    # Reload using pandas
    col_names = ["Year", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
    df = pd.read_csv(nao_file, header=0)
    df.set_axis(col_names, axis=1, inplace=True)

    # Set index
    df = df.set_index('Year')
    df = df.stack()
    df.index = pd.to_datetime(df.index.get_level_values(1) + '-' + df.index.get_level_values(0).astype('str'))

    ## ----  plot Winter NAO ---- ####
    # Select only DJF
    df_winter_djf = df[(df.index.month==12) | (df.index.month==1) | (df.index.month==2)]
    df_winter = df[(df.index.month==12) | (df.index.month==1) | (df.index.month==2) | (df.index.month==3)]
    df_summer = df[(df.index.month>=6) & (df.index.month<=9)]

    # Start Dec-1950
    df_winter_djf = df_winter_djf[df_winter_djf.index>pd.to_datetime('1950-10-01')]
    df_winter = df_winter[df_winter.index>pd.to_datetime('1950-10-01')]
    df_summer = df_summer[df_summer.index>pd.to_datetime('1950-01-01')]

    # Average 3 consecutive values (DJF average); We loose index.
    df_winter_djf = df_winter_djf.groupby(np.arange(len(df_winter_djf))//3).mean()
    df_winter = df_winter.groupby(np.arange(len(df_winter))//4).mean()
    df_summer = df_summer.groupby(np.arange(len(df_summer))//3).mean()

    # Reset index using years only
    year_unique = pd.unique(df.index.year)[1:,]
    df_winter = df_winter.iloc[np.arange(0, year_unique.size)] 
    df_winter_djf = df_winter_djf.iloc[np.arange(0, year_unique.size)] 
    df_winter.index = year_unique
    df_winter_djf.index = year_unique

    df_summer = df_summer.iloc[np.arange(0, year_unique.size)] 
    df_summer.index = year_unique

    # pickle DataFrame for scorecards:
    df.to_pickle('NAO_monthly.pkl')
    df.to_csv('NAO_monthly.csv')
    df_annual = df.resample('As').mean()
    df_annual.index = df_annual.index.year
    df_annual.to_pickle('NAO_annual.pkl')
    df_winter.to_pickle('NAO_winter.pkl')
    df_summer.to_pickle('NAO_summer.pkl')


    #nao_clim = df_winter[(df_winter.index>=1981) & (df_winter.index<=2010)]
    #df_winter = (df_winter - df_winter.mean()) / nao_clim.std()


    ## ---- plot winter NAO bar plots ---- ##
    #df_winter[df_winter.index==2021]=np.nan # Remove 2021 for 2020 ResDoc
    df1 = df_winter[df_winter>0]
    df2 = df_winter[df_winter<0]

    ## ---- plot winter NAO bar plots #2 ---- ##
    fig = plt.figure(1)
    fig.clf()
    width = .9
    p1 = plt.bar(df1.index, np.squeeze(df1.values), width, alpha=0.8, color='steelblue')
    p2 = plt.bar(df2.index, np.squeeze(df2.values), width, bottom=0, alpha=0.8, color='indianred')
    #p1 = plt.bar(df1.index[-1], np.squeeze(df1.values[-1]), width, alpha=.5, color='white')
    plt.ylabel('NAO subindex')
    plt.title('Winter NAO average (DJFM)')
    ticks = plt.gca().xaxis.get_ticklocs()
    plt.fill_between([ticks[0]-1, ticks[-1]+1], [-.5, -.5], [.5, .5], facecolor='gray', alpha=.2)
    plt.xlim([1950, YEAR+1])
    plt.grid()
    fig.set_size_inches(w=15,h=9)
    fig_name = 'NAO_winter_1950-' + str(YEAR) + '.png'
    #plt.annotate('data source: www.ncdc.noaa.gov/teleconnections/', xy=(.58, .01), xycoords='figure fraction', annotation_clip=False, FontSize=12)
    fig.savefig(fig_name, dpi=300)
    os.system('convert -trim -bordercolor White -border 10x10 ' + fig_name + ' ' + fig_name)

    # French Figure
    fig = plt.figure(1)
    fig.clf()
    width = .9
    p1 = plt.bar(df1.index, np.squeeze(df1.values), width, alpha=0.8, color='steelblue')
    p2 = plt.bar(df2.index, np.squeeze(df2.values), width, bottom=0, alpha=0.8, color='indianred')
    #p1 = plt.bar(df1.index[-1], np.squeeze(df1.values[-1]), width, alpha=.5, color='white')
    plt.ylabel('indice ONA')
    plt.title('Oscillation Nord-Atlantique hivernale (DJFM)')
    plt.grid()
    fig.set_size_inches(w=15,h=9)
    fig_name = 'NAO_winter_1950-' + str(YEAR) + '_FR.png'
    plt.annotate('source données: www.ncdc.noaa.gov/teleconnections/', xy=(.58, .01), xycoords='figure fraction', annotation_clip=False, fontsize=12)
    fig.savefig(fig_name, dpi=300)
    os.system('convert -trim ' + fig_name + ' ' + fig_name)

    return None

def ao(YEAR,ao_file_loc,url_loc='https://www.cpc.ncep.noaa.gov/products/precip/CWlink/daily_ao_index/monthly.ao.index.b50.current.ascii.table'):
    '''
    Update the Arctic Oscillation
    run in /home/cyrf0006/AZMP/state_reports/airTemp

    Generate data:
    AO_annual.pkl

    Generate Figures:
    AO_winter_1950-YYYY.png'

    Adapted from ao.py
    '''
    # Adjust fontsize/weight
    font = {'family' : 'normal',
            'weight' : 'bold',
            'size'   : 18}
    plt.rc('font', **font)

    # Download and save up-to-date  AO index from NOAA (data.csv) if needed
    # url = 'https://www.ncdc.noaa.gov/teleconnections/ao/data.csv'
    url = url_loc
    ao_file = ao_file_loc
    if os.path.exists(ao_file):
        py3 = version_info[0] > 2 #creates boolean value for test that Python major version > 2        
        response_isnt_good = True
        while response_isnt_good:
            if py3:
                response = input('Do you what to update '  + ao_file + '? [y/n]')
            else:
                response = raw_input('Do you what to update '  + ao_file + '? [y/n]')

            if response == 'y':
                import urllib3
                http = urllib3.PoolManager()
                r = http.request('GET', url)
                open(ao_file_loc, 'wb').write(r.data)
                response_isnt_good = False
            elif response == 'n':
                response_isnt_good = False
            else:
                print(' -> Please answer "y" or "n"')

    # Reload using pandas
    df = pd.read_csv(ao_file, header=0, delimiter=r"\s+", error_bad_lines=False)

    # Set index
    df.index.name = 'Year'
    df = df.stack()
    df.index = pd.to_datetime(df.index.get_level_values(1) + '-' + df.index.get_level_values(0).astype('str'))

    # Resample
    # pickle DataFrame for scorecards:
    #df.to_pickle('AO_monthly.pkl')
    df_annual = df.resample('As').mean()
    df_annual.index = df_annual.index.year
    df_annual.to_pickle('AO_annual.pkl')

    ## ## ---- plot winter AO bar plots ---- ##
    df1 = df_annual[df_annual>0]
    df2 = df_annual[df_annual<0]

    ## ---- plot winter AO bar plots #2 ---- ##
    fig = plt.figure(4)
    fig.clf()
    width = .9
    p1 = plt.bar(df1.index, np.squeeze(df1.values), width, alpha=0.8, color='steelblue')
    p2 = plt.bar(df2.index, np.squeeze(df2.values), width, bottom=0, alpha=0.8, color='indianred')
    p1 = plt.bar(df1.index[-1], np.squeeze(df1.values[-1]), width, alpha=.3, color='black')
    plt.ylabel('AO index')
    #plt.xlabel('Year')
    plt.title('Annual AO average')
    plt.grid()
    fig.set_size_inches(w=15,h=9)
    fig_name = 'AO_1950-' + str(2021) + '.png'
    plt.annotate('data source: www.ncdc.noaa.gov/teleconnections/', xy=(.58, .01), xycoords='figure fraction', annotation_clip=False, FontSize=12)
    fig.savefig(fig_name, dpi=300)
    os.system('convert -trim ' + fig_name + ' ' + fig_name)

    return None





def amo(YEAR,amo_file_loc,url_loc='https://www.esrl.noaa.gov/psd/data/correlation/amon.us.data'):
    '''
    Update the Atlantic Multidecadal Oscillation
    run in /home/cyrf0006/AZMP/state_reports/airTemp

    Generate data:
    AMO_monthly.pkl
    AMO_annual.pkl

    Generate Figures:
    AMO_monthly_1950-YYYY.png'
    AMO_bar_1950-YYYY.png'

    Adapted from amo.py
    '''

    # Adjust fontsize/weight
    font = {'family' : 'normal',
            'weight' : 'bold',
            'size'   : 18}
    plt.rc('font', **font)

    # Download and save up-to-date  AMO index from NOAA (data.csv) if needed
    url = url_loc
    amo_file = amo_file_loc
    if os.path.exists(amo_file):
        py3 = version_info[0] > 2 #creates boolean value for test that Python major version > 2        
        response_isnt_good = True
        while response_isnt_good:
            if py3:
                response = input('Do you what to update '  + amo_file + '? [y/n]')
            else:
                response = raw_input('Do you what to update '  + amo_file + '? [y/n]')

            if response == 'y':
                import urllib3
                http = urllib3.PoolManager()
                r = http.request('GET', url)
                open(amo_file_loc, 'wb').write(r.data)
                response_isnt_good = False
            elif response == 'n':
                response_isnt_good = False
            else:
                print(' -> Please answer "y" or "n"')


    # Reload using pandas
    col_names = ["Year", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
    df = pd.read_csv(amo_file, header=1, delimiter=r"\s+", error_bad_lines=False, names=col_names)
    # Remove last rows
    df = df.iloc[0:-4]

    # Set index
    df = df.set_index('Year')
    df = df.stack()      
    df.index = pd.to_datetime(df.index.get_level_values(1) + '-' + df.index.get_level_values(0))

    # Remove wrong values
    df = df.astype(float)
    df = df[df>-99]

    ## ----  plot Monthly AMO + 5year running mean ---- ##
    fig = plt.figure(1)
    fig.clf()
    plt.plot(df)
    plt.plot(df.rolling(window=60, center=True).mean(), linewidth=3)
    plt.ylabel('AMO')
    plt.xlabel('Year')
    plt.grid()
    plt.title('AMO unsmoothed, detrended from the Kaplan SST V2')
    plt.legend(['monthly', '5-year smooth'])
    fig_name = 'AMO_monthly_1950-' + str(YEAR) + '.png'
    fig.set_size_inches(w=12,h=9)
    fig.savefig(fig_name, dpi=300)

    # pickle DataFrame for scorecards:
    df.to_pickle('AMO_monthly.pkl')
    df_annual = df.resample('As').mean()
    df_annual.index = df_annual.index.year
    df_annual.to_pickle('AMO_annual.pkl')

    ## ## ---- plot annual AMO bar plots ---- ##
    df1 = df_annual[df_annual>0]
    df2 = df_annual[df_annual<0]
    
    fig = plt.figure(4)
    fig.clf()
    width = .9
    p1 = plt.bar(df1.index, np.squeeze(df1.values), width, alpha=0.8, color='steelblue')
    p2 = plt.bar(df2.index, np.squeeze(df2.values), width, bottom=0, alpha=0.8, color='indianred')
    p1 = plt.bar(df1.index[-1], np.squeeze(df1.values[-1]), width, alpha=.3, color='black')
    plt.ylabel('AMO index')
    #plt.xlabel('Year')
    plt.title('AMO unsmoothed, detrended from the Kaplan SST V2')
    plt.grid()
    fig.set_size_inches(w=15,h=9)
    fig_name = 'AMO_bar_1950-' + str(YEAR) + '.png'
    plt.annotate('data source: http://www.esrl.noaa.gov/psd/data/timeseries/AMO/', xy=(.58, .01), xycoords='figure fraction', annotation_clip=False, FontSize=12)
    fig.savefig(fig_name, dpi=300)
    os.system('convert -trim ' + fig_name + ' ' + fig_name)

    return None
