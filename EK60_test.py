#!/usr/bin/env python3

import matplotlib.pyplot as plt
from echonix import ek60, echogram, imaging, raw

# Displays an echogram of an interesting krill swarm at 120kHz

#filename = r'/home/cyrf0006/research/AZMP_surveys/TEL176/EK60Data/TEL176-D20170726-T062139.raw'
filename = r'/home/cyrf0006/research/AZMP_surveys/TEL176/EK60Data/TEL176-D20170276-T095507.raw'
#filename = r'/home/cyrf0006/research/AZMP_surveys/TEL176/EK60Data/TEL176-D20170726-T132835.raw'
#### ---- To print ome info ---- ###
#d = raw.load_raw(filename)
#print(d[0].configurationheader.surveyname)

#### ---- Single frequency ---- ####
## frequency = 120000
## Sv, r = ek60.raw_to_sv(filename, frequency)
## echogram.egshow(Sv, max = -50, min = -95, range=r)

## cbar = plt.colorbar()
## cbar.set_label('Sv (dB)', rotation=90)

## plt.xlabel('Sample')
## plt.ylabel('Range /m')

## plt.title('Test figure')

## plt.show()

#### ---- Composite ---- ####

Sv38, r = ek60.raw_to_sv(filename, 38000)
Sv120, r = ek60.raw_to_sv(filename, 120000)
Sv200, r = ek60.raw_to_sv(filename, 200000)

im = imaging.composite(Sv38, Sv120, Sv200, min = -95, max = -50)

echogram.imshow(im, range=r)

plt.xlabel('Sample')
plt.ylabel('Range /m')

plt.title('Test plot EK60')
plt.colorbar()

plt.show()
