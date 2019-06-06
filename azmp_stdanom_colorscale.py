# -*- coding: utf-8 -*-
'''
Standard colorscale for standardized anomalies

'''

import numpy as  np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.colors import from_levels_and_colors
import pandas as pd
import os
import unicodedata

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

# Build the dataFrame
data = np.linspace(-3.25, 3.25, 14)
data_lims = np.linspace(-3.5, 3.5, 15)

text = []
for idx, txt in enumerate(data_lims):
    if idx<data_lims.size-1:
        if txt<0:
            text.append(str(txt+.1) + ' to ' + str(data_lims[idx+1]))
        elif txt>=0:
            text.append(str(txt) + ' to ' + str(data_lims[idx+1]-.1))
            
text[0] = '< ' + str(data_lims[1])
text[-1] = '> ' + str(data_lims[-2])

df = pd.DataFrame([data], columns=text)
vals = df.values

# Build the colormap
## vmin = -3.5
## vmax = 3.5
## midpoint = 0
## levels = np.linspace(vmin, vmax, 15)
## midp = np.mean(np.c_[levels[:-1], levels[1:]], axis=1)
## #colvals = np.interp(midp, [vmin, midpoint, vmax], [-1, 0., 1])
## normal = plt.Normalize(-3.5, 3.5)

## reds = plt.cm.Reds(np.linspace(0,1, num=7))
## blues = plt.cm.Blues_r(np.linspace(0,1, num=7))
## whites = [(1,1,1,1)]*1
## colors = np.vstack((blues[0:-1,:], whites, reds[1:,:]))
## colors = np.concatenate([[colors[0,:]], colors, [colors[-1,:]]], 0)
## cmap, norm = from_levels_and_colors(midp, colors, extend='both')
## #cmap, norm = from_levels_and_colors(colvals, colors)

# Build the colormap
vmin = -3.51
vmax = 3.51
midpoint = 0
levels = np.linspace(vmin, vmax, 15)
midp = np.mean(np.c_[levels[:-1], levels[1:]], axis=1)
colvals = np.interp(midp, [vmin, midpoint, vmax], [-1, 0., 1])
normal = plt.Normalize(-3.51, 3.51)
reds = plt.cm.Reds(np.linspace(0,1, num=7))
blues = plt.cm.Blues_r(np.linspace(0,1, num=7))
whites = [(1,1,1,1)]*2
colors = np.vstack((blues[0:-1,:], whites, reds[1:,:]))
colors = np.concatenate([[colors[0,:]], colors, [colors[-1,:]]], 0)
cmap, norm = from_levels_and_colors(levels, colors, extend='both')

# build the table
nrows, ncols = vals.shape
hcell, wcell = 3, 1.8
hpad, wpad = 0, 0
fig=plt.figure(figsize=(ncols*wcell+wpad, nrows*hcell+hpad))
ax = fig.add_subplot(111)
ax.axis('off')
the_table=ax.table(cellText=None, rowLabels=None, colLabels=df.columns,
                    loc='center', cellColours=cmap(normal(vals)), cellLoc='center',
                    bbox=[0, 0, 1, 0.5]
                    )
# change font color to white where needed:
the_table.auto_set_font_size(False)
the_table.set_fontsize(12.5)
for key, cell in the_table.get_celld().items():
    cell._text.set_fontsize(15)
    cell._text.set_weight('bold')

# Save Figure
plt.savefig("anomaly_colors.png", dpi=300)
os.system('convert -trim anomaly_colors.png anomaly_colors.png')

