'''
Plot G. Han's Lab Current transport

Frederic.Cyr@dfo-mpo.gc.ca - February 2019

'''

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
#import datetime
#import matplotlib.dates as mdates
import pandas as pd
import os
#from sys import version_info
#import re

# Adjust fontsize/weight
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 12}
plt.rc('font', **font)

## ----  Read data ---- ##
df = pd.read_csv('/home/cyrf0006/AZMP/altimetry/lc_index_1993_2018.txt', header=None, delimiter='\s+')
df = df.set_index(0)
df = df.rename(columns={1 : 'NL'})
df = df.rename(columns={2 : 'SS'})

ave_ss = .6
std_ss = .3
ave_nl = 13
std_nl = 1.4

df.NL = df.NL*std_nl + ave_nl
df.SS = df.SS*std_ss + ave_ss

XLIMS = [1992, 2019]

## ---- plot ---- ##
fig = plt.figure(1)
plt.clf()

# AX1 - LC
ax1 = plt.subplot2grid((2, 1), (0, 0), rowspan=1, colspan=1)
plt.plot(df.index, df.NL)
plt.fill_between(df.index, np.repeat(df.NL.mean(),df.index.size)-df.NL.std()/2,np.repeat(df.NL.mean(),df.index.size)+df.NL.std()/2 , facecolor='steelblue', interpolate=True , alpha=.3)
plt.legend(['NL', r'$\rm \overline{x}\pm\,0.5 sd$'])
plt.xlim([df.index.min(), df.index.max()])
ax1.tick_params(labelbottom='off')
ax1.xaxis.label.set_visible(False)
ax1.set_ylabel('Transport (Sv)', fontsize=14, fontweight='bold')

# AX2 - current year
ax2 = plt.subplot2grid((2, 1), (1, 0), rowspan=1, colspan=1)
plt.plot(df.index, df.SS)
plt.fill_between(df.index, np.repeat(df.SS.mean(),df.index.size)-df.SS.std()/2,np.repeat(df.SS.mean(),df.index.size)+df.SS.std()/2 , facecolor='steelblue', interpolate=True , alpha=.3)
plt.legend(['SS', r'$\rm \overline{x}\pm\,0.5 sd$'])
plt.xlim([df.index.min(), df.index.max()])
ax2.set_ylabel('Transport (Sv)', fontsize=14, fontweight='bold')


# Save fig
fig.set_size_inches(w=12, h=9)
fig_name = 'LC_transport.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)


### ---- Anomalies plot ---- ##
df = pd.read_csv('/home/cyrf0006/AZMP/altimetry/lc_index_1993_2018.txt', header=None, delimiter='\s+')
df = df.set_index(0)
df = df.rename(columns={1 : 'NL'})
df = df.rename(columns={2 : 'SS'})
df_SS = df.SS
df_NL = df.NL


# plot
fig = plt.figure(4)
fig.clf()
ax = plt.subplot2grid((2, 1), (0, 0))
df1 = df_NL[df_NL>0]
df2 = df_NL[df_NL<0]
width = .8
p1 = plt.bar(df1.index, np.squeeze(df1.values), width, alpha=0.8, color='indianred')
p2 = plt.bar(df2.index, np.squeeze(df2.values), width, bottom=0, alpha=0.8, color='steelblue')
plt.fill_between([1990, 2020], [-.5, -.5], [.5, .5], facecolor='gray', alpha=.2)
#plt.xticks(df.index[::1], rotation='vertical')
plt.ylabel('LC index')
plt.xlabel('Year')
plt.title('NL slope')
plt.xlim(XLIMS)
plt.ylim([-2,2])
plt.text(XLIMS[0], -1.6,  r'$\,\rm \overline{x} = $' + str(ave_nl) + r' $\pm$ ' +  str(std_nl) + ' Sv')
plt.grid()
ax.xaxis.label.set_visible(False)
ax.tick_params(labelbottom='off')

ax2 = plt.subplot2grid((2, 1), (1, 0))
df3 = df_SS[df_SS>0]
df4 = df_SS[df_SS<0]
p1 = plt.bar(df3.index, np.squeeze(df3.values), width, alpha=0.8, color='indianred')
p2 = plt.bar(df4.index, np.squeeze(df4.values), width, bottom=0, alpha=0.8, color='steelblue')
plt.fill_between([1990, 2020], [-.5, -.5], [.5, .5], facecolor='gray', alpha=.2)
#plt.xticks(df.index[::1], rotation='vertical')
plt.title('SS slope')
plt.ylabel('LC index')
plt.xlim(XLIMS)
plt.ylim([-2,2])
plt.text(XLIMS[0], 1.5,  r'$\,\rm \overline{x} = $' + str(ave_ss) + r' $\pm$ ' +  str(std_ss) + ' Sv')
plt.grid()

fig.set_size_inches(w=7,h=5)
fig_name = 'LC_index.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)


# Save French Figure
