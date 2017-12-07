# A first test to read Excel nutrient file and export to Pandas.

# Check in:
#  /home/cyrf0006/research/AZMP_database/biochem

import pandas
import matplotlib.pyplot as plt
import numpy as np

df = pandas.read_excel('~/github/AZMP-NL/data/AZMP_Nutrients_1999_2016_Good_Flags_Only.xlsx')

#print the column names
print df.columns

#get the values for a given column
#values = df['sas_date'].values[1:10]

#get a data frame with selected columns
#FORMAT = ['Col_1', 'Col_2', 'Col_3']
#df_selected = df[FORMAT]


plt.plot(df['salinity'], df['NO3'], '.k')
plt.plot(df['temp'], df['NO3'], '.k')
plt.plot(df['PO4'], df['NO3'], '.k')
plt.show()

plt.hist(df['PO4'], 50)
plt.show()


PO4 = df['PO4'].values
NO3 = df['NO3'].values
idx = [i for i,v in enumerate(NO3) if v > 1e-4]
ratio =  PO4[idx]/NO3[idx]
    
plt.hist(ratio, np.linspace(0,10,1000))
plt.show()
plt.plot(df['NO3'][idx], ratio, '.k')
plt.show()


# A lot to play with!

# Test xarray
import xarray as xr
xr.Dataset.from_dataframe(df)
