'''
AZMP script to extract SST data in standard reporting boxes from the BIO's remote sensing server:
ftp://ftp.dfo-mpo.gc.ca/bometrics
ftp://ftp.dfo-mpo.gc.ca/bometrics/noaa/stats/boxes/

For oofline use, you can download all data with:
wget -m ftp://ftp.dfo-mpo.gc.ca/bometrics/noaa/stats/boxes/*.stat

check in : /home/cyrf0006/AZMP/annual_meetings/2019
http://www.bio.gc.ca/science/data-donnees/base/data-donnees/sst-en.php
Frederic.Cyr@dfo-mpo.gc.ca - February 2019

'''

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import datetime
import matplotlib.dates as mdates
import pandas as pd
import os
from sys import version_info
import re

## ----  Parameters to edit (getting data infile or offline) ---- ##
#prefix = 'ftp://ftp.dfo-mpo.gc.ca/bometrics/noaa/stats/boxes/'
prefix = '/home/cyrf0006/data/BIO_remote/bometrics/noaa/stats/boxes/'
# Get region names
df_box = pd.read_excel('/home/cyrf0006/github/AZMP-NL/utils/SST_boxes.xslx')

## ---- Loop on NL regions and store in a dataFrame ---- ##
dfs = []
df_labels = []
for box in df_box[(df_box.region=='NL') | (df_box.region=='NS') | (df_box.region=='NS')].box_name.values:
    df = pd.read_csv(prefix + box +'_sst.stat', delimiter='\s+')
    df = df.rename(columns={'date-id':'date'})
    # Set index (need first to swap (a,b) by (7,15))
    date_tmp = df.date.str[-1].apply(lambda x: re.sub('a','07',x))
    date_day = date_tmp.str[-1].apply(lambda x: re.sub('b','21',x)) 
    date = df.date.map(lambda x: str(x)[:-1]) # strip last character
    df['date'] = date.astype(str) + date_day
    df = df.set_index('date')
    df.index = pd.to_datetime(df.index, format='%Y%b%d')

    dfs.append(df)
    df_labels.append(box)

# convert the list of DataFrames into a single multiindex DataFrame
df_all = pd.concat(dfs, keys=df_labels, axis=0)    


## ---- Just mean SST now ---- ##
df_sst = df_all.mean_sst 
df_sst = df_sst.unstack(level=0)
df_sst = df_sst.replace(-999.00000, np.NaN)
df_sst = df_sst[df_sst.index.year<=2018]
df_sst = df_sst.resample('As').mean()

df_sst.to_pickle('SSTs_bometrics_annual.pkl')
