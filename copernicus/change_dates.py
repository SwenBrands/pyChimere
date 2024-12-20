#!/usr/bin/env python

from netCDF4 import Dataset
import xarray as xr
import numpy as np
import pandas as pd
import netCDF4
import numpy.matlib
from netCDF4 import num2date, date2num

#set input variables
import sys
tardate=str(sys.argv[1])
taryear=str(sys.argv[2])
tarmonth=str(sys.argv[3])

#tardate='20171014'
#taryear='2017'
#tarmonth='10'

print(tardate)
print(taryear)
print(tarmonth)

filename = '/mnt/lustre/scratch//home/cesga/orballo/OP/DATOS/CICC/chimere2017r4/MACC/'+tardate+'/gasmet_'+taryear+tarmonth+'.nc'
print(filename)

## EXECUTE #############################################################
nc = xr.open_dataset(filename)

#dates = nc.variables['time']
#dates = str(dates.values)
#dates = dates.replace("-", "")
#dates = dates.replace("\n", "")
#dates = dates.replace("\n", "")
#dates = dates.replace("T", ".")
#dates = dates.replace(":00:00.000000000","")
#dates = dates.replace(" ",",")
#dates = dates.replace("[","")
#dates = dates.replace("]","")
#dates = dates.replace("'","")
#dates = dates.split(',')

#for i in xrange(len(dates)):
    #dates[i] = dates[i][0:8]

#dates=[float(i) for i in dates]
#for i in xrange(len(dates)):
    #dates[i] = dates[i]+0.125*i

##overwrite dates with netCDF4 package
#dataset = Dataset(filename,'r+')
#dataset.variables['time'][:] = dates

#convert from kg/kg to ug/m3
dust1 = nc.variables['dust1']
#dust2 = nc.variables['dust2']
#dust3 = nc.variables['dust3']

#air density as retrieved from https://www.ecmwf.int/en/forecasts/documentation-and-support/60-model-levels
density = np.array((0.000150,0.000402,0.000672,0.001025,0.001485,0.002087,0.002814,0.003688,0.004740,0.006018,0.007599,0.009582,0.012082,0.015130,0.018863,0.023518,0.029322,0.036557,0.045579,0.056826,0.070849,0.087826,0.107127,0.129274,0.154328,0.182375,0.213467,0.247616,0.284796,0.324943,0.367189,0.403822,0.441886,0.481280,0.521847,0.563389,0.605674,0.648445,0.691426,0.734331,0.776863,0.818729,0.859634,0.899295,0.937436,0.973801,1.008150,1.040270,1.069971,1.097099,1.121533,1.143192,1.162039,1.178087,1.191403,1.202112,1.210406,1.216546,1.220871,1.223803))

nc.close()

dataset.close()
#dset['var'][:][dset['var'][:] < 0] = -1
