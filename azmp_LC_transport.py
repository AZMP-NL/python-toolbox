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


INFILE = '/home/cyrf0006/github/AZMP-NL/external_data/LC_transport/lc_index.txt'

# Adjust fontsize/weight
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 12}
plt.rc('font', **font)

## ----  Read data ---- ## (I just drop column that was there)
df = pd.read_csv(INFILE, header=None, delimiter='\s+', names=['year', 'NL', 'SS'])
df = df.set_index('year')
df.to_csv('LC_index.csv', float_format='%.2f')

ave_ss = .6
std_ss = .3
ave_nl = 13
std_nl = 1.4

df.NL = df.NL*std_nl + ave_nl
df.SS = df.SS*std_ss + ave_ss

df.to_csv('LC_transport.csv', float_format='%.2f')

XLIMS = [1990, 2024]

## ---- plot ---- ##
fig = plt.figure(1)
plt.clf()

# AX1 - LC
ax1 = plt.subplot2grid((2, 1), (0, 0), rowspan=1, colspan=1)
plt.plot(df.index.values, df.NL.values)
plt.fill_between(XLIMS, np.repeat(df.NL.mean(),2)-df.NL.std()/2,np.repeat(df.NL.mean(),2)+df.NL.std()/2 , facecolor='steelblue', interpolate=True , alpha=.3)
plt.legend(['NL', r'$\rm \overline{x}\pm\,0.5 sd$'])
plt.xlim([df.index.min(), df.index.max()])
ax1.tick_params(labelbottom='off')
ax1.xaxis.label.set_visible(False)
ax1.set_ylabel('Transport (Sv)', fontsize=14, fontweight='bold')

# AX2 - current year
ax2 = plt.subplot2grid((2, 1), (1, 0), rowspan=1, colspan=1)
plt.plot(df.index.values, df.SS.values)
plt.fill_between(XLIMS, np.repeat(df.SS.mean(),2)-df.SS.std()/2,np.repeat(df.SS.mean(),2)+df.SS.std()/2 , facecolor='steelblue', interpolate=True , alpha=.3)
plt.legend(['SS', r'$\rm \overline{x}\pm\,0.5 sd$'])
plt.xlim([df.index.min(), df.index.max()])
ax2.set_ylabel('Transport (Sv)', fontsize=14, fontweight='bold')


# Save fig
fig.set_size_inches(w=12, h=9)
fig_name = 'LC_transport.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)


### ---- Anomalies plot ---- ##
df = pd.read_csv(INFILE, header=None, delimiter='\s+', names=['year', 'NL', 'SS'])
df = df.set_index('year')
df_SS = df.SS
df_NL = df.NL



# plot
fig = plt.figure(4)
fig.clf()
ax = plt.subplot2grid((2, 1), (0, 0))
df1 = df_NL[df_NL>0]
df2 = df_NL[df_NL<0]
width = .5
p1 = plt.bar(df1.index, np.squeeze(df1.values), width, alpha=1, color='indianred', zorder=10)
p2 = plt.bar(df2.index, np.squeeze(df2.values), width, bottom=0, alpha=1, color='steelblue', zorder=10)
plt.fill_between(XLIMS, [-.5, -.5], [.5, .5], facecolor='gray', alpha=.2)
#plt.xticks(df.index[::1], rotation='vertical')
plt.ylabel('LC index')
plt.title('NL slope')
plt.xlim(XLIMS)
plt.ylim([-2,2])
plt.text(XLIMS[0], -1.8,  r'$\,\rm \overline{x} = $' + str(ave_nl) + r' $\pm$ ' +  str(std_nl) + ' Sv')
plt.grid()
ax.tick_params(labelbottom=False) # can pass a series of params here

ax2 = plt.subplot2grid((2, 1), (1, 0))
df3 = df_SS[df_SS>0]
df4 = df_SS[df_SS<0]
p1 = plt.bar(df3.index, np.squeeze(df3.values), width, alpha=1, color='indianred', zorder=10)
p2 = plt.bar(df4.index, np.squeeze(df4.values), width, bottom=0, alpha=1, color='steelblue', zorder=10)
plt.fill_between(XLIMS, [-.5, -.5], [.5, .5], facecolor='gray', alpha=.2)
#plt.xticks(df.index[::1], rotation='vertical')
plt.title('SS slope')
plt.ylabel('LC index')
plt.xlim(XLIMS)
plt.ylim([-2,2])
plt.text(XLIMS[0], -1.8,  r'$\,\rm \overline{x} = $' + str(ave_ss) + r' $\pm$ ' +  str(std_ss) + ' Sv')
plt.grid()

fig.set_size_inches(w=7,h=5)
fig_name = 'LC_index.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)


# Save French Figure
fig = plt.figure(5)
fig.clf()
ax = plt.subplot2grid((2, 1), (0, 0))
df1 = df_NL[df_NL>0]
df2 = df_NL[df_NL<0]
width = .5
p1 = plt.bar(df1.index, np.squeeze(df1.values), width, alpha=1, color='indianred', zorder=10)
p2 = plt.bar(df2.index, np.squeeze(df2.values), width, bottom=0, alpha=1, color='steelblue', zorder=10)
plt.fill_between(XLIMS, [-.5, -.5], [.5, .5], facecolor='gray', alpha=.2)
#plt.xticks(df.index[::1], rotation='vertical')
plt.ylabel('Index CL')
plt.title('Talus du Labrador')
plt.xlim(XLIMS)
plt.ylim([-2,2])
plt.text(XLIMS[0], -1.6,  r'$\,\rm \overline{x} = $' + str(ave_nl) + r' $\pm$ ' +  str(std_nl) + ' Sv')
plt.grid()
ax.tick_params(labelbottom=False) # can pass a series of params here

ax2 = plt.subplot2grid((2, 1), (1, 0))
df3 = df_SS[df_SS>0]
df4 = df_SS[df_SS<0]
p1 = plt.bar(df3.index, np.squeeze(df3.values), width, alpha=1, color='indianred', zorder=10)
p2 = plt.bar(df4.index, np.squeeze(df4.values), width, bottom=0, alpha=1, color='steelblue', zorder=10)
plt.fill_between(XLIMS, [-.5, -.5], [.5, .5], facecolor='gray', alpha=.2)
#plt.xticks(df.index[::1], rotation='vertical')
plt.title('Talus Néo-Écossais')
plt.ylabel('Index CL')
plt.xlim(XLIMS)
plt.ylim([-2,2])
plt.text(XLIMS[0], 1.5,  r'$\,\rm \overline{x} = $' + str(ave_ss) + r' $\pm$ ' +  str(std_ss) + ' Sv')
plt.grid()

fig.set_size_inches(w=7,h=5)
fig_name = 'LC_index_FR.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)
