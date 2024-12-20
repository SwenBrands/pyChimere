#!/usr/bin/env python

#this scripts loads the CHIMERE output files and searches out the data at the grid-boxes nearest to those defined in <tarlat> and <tarlon>, the output is then save in nc format.

import xarray as xr
import numpy as np
import pandas as pd
import sys
import dask
import os
import numpy.matlib
import string
from datetime import datetime
from dateutil.parser import parse
from matplotlib import pyplot as plt
import scipy.stats

expname=str(sys.argv[1]) #name of the experiment

homedir = os.getenv("HOME")
lustre = os.getenv("LUSTRE")

## DEFINE USER OPTIONS
startdate = '2018-04-01 04:00:00'
enddate = '2018-04-30 23:00:00'
variables = ['PM10','O3'] #which variables are to be loaded
domain = 'fine' #target domain, MAKE SURE THAT THIS IS THE RIGHT DOMAIN
selmodel = 1
nanlimit = 60 #in percent
titlesize = 8
xticksize = 5
yticksize = 7
obslabel = 'OBS'
oplabel = 'OP'
ylabel = 'ug/m3'
minval = 1 #minimum value allowed in observations, set to nan if lower

## EXECUTE #############################################################

## load metadata from metadata.py
tarlat = [] #target latitudes, will be loaded by metadata.py
tarlon = [] #target longitudes, will be loaded by metadata.py
station_names_meta = [] #station names, will be loaded by metadata.py

execfile("metadata.py")

#for each day, load the model data at the gridboxes nearest to the station locations defined in tarlat y tarlon. Finally, concatenate the time series
src_dirs = [lustre+'/SWEN/DATA/OBSDATA',lustre+'/SWEN/DATA/MODELDATA',lustre+'/OP/DATOS/RESULTADOS/chimere2017r4/'+expname+'/validation']
#src_dirs = [lustre+'/SWEN/DATA/OBSDATA', lustre+'/SWEN/DATA/MODELDATA', lustre+'/OP/DATOS/RESULTADOS/chimere2017r4/exp0/validation', lustre+'/OP/DATOS/RESULTADOS/chimere2017r4/exp0a/validation', lustre+'/OP/DATOS/RESULTADOS/chimere2017r4/exp0a/validation', lustre+'/OP/DATOS/RESULTADOS/chimere2017r4/exp1/validation', lustre+'/OP/DATOS/RESULTADOS/chimere2017r4/exp2/validation', lustre+'/OP/DATOS/RESULTADOS/chimere2017r4/exp5/validation']
savedir = lustre+'/SWEN/val_chimere/series'

#load example file to check the size of the time dim
srcfile=src_dirs[0]+'/'+variables[0]+'_*.nc'
nc=xr.open_mfdataset(srcfile)
cstring='nc.'+variables[0]
dummy = eval(cstring)
dummy = dummy.sel(time=slice(startdate,enddate))
#init output array
nc.close()

#load the data
ARRAY = np.zeros((len(src_dirs),len(variables),len(dummy),dummy.shape[1])) 
for ii in xrange(len(src_dirs)):
    for vv in xrange(len(variables)):    
        #load the data
        srcfile=src_dirs[ii]+'/'+variables[vv]+'_*.nc'
        print('INFO: loading '+srcfile+'...')
        nc=xr.open_mfdataset(srcfile)
        cstring='nc.'+variables[vv]
        data = eval(cstring)
        data = data.sel(time=slice(startdate,enddate))                
        lat = nc.id.lat
        lon = nc.id.lon
        names = nc.id.names
        dates = nc.time
        dates = dates.sel(time=slice(startdate,enddate))
        ARRAY[ii,vv,:,:] = np.expand_dims(data.values,0)          
        #valind = np.where(ARRAY[0,:,1]>0)

#convert dates to pandas datetimeindex
dates = pd.date_range(start=dates.values[0], end=dates.values[-1], freq='H')
#find indexes for the start of the day (used when results are plotted below)
dates_uniq = np.unique(dates.date).astype(str)
#find index of hour 12 for each unique day
xtick = [np.where(dates == dates_uniq[dd]+' 12:00:00')[0][0] for dd in range(len(dates_uniq))]#
dates_uniq = dates_uniq.tolist()
#to save space in the figure, only the days of the months are going to be plotted.
days_uniq = [dates_uniq[ii][-2:] for ii in xrange(len(dates_uniq))]

#start the verfication
OBSARR = ARRAY[0,] #ARRAY containing observations (3D)
MODARR = ARRAY[1:,] #ARRAY containing model output (4D since various model versions are considered
BIAS = np.zeros((len(MODARR),len(variables),ARRAY.shape[-1]))
RHO = np.zeros((len(MODARR),len(variables),ARRAY.shape[-1]))
RATIO = np.zeros((len(MODARR),len(variables),ARRAY.shape[-1]))
MAE = np.zeros((len(MODARR),len(variables),ARRAY.shape[-1]))
for ii in xrange(len(MODARR)):
    MODSTEP = MODARR[ii,]
    for vv in xrange(len(variables)):
        for st in xrange(ARRAY.shape[-1]):
            #calc nanmean and nanstd for outlier definition
            meanval = np.nanmean(OBSARR,axis=1)
            stdval = np.nanstd(OBSARR,axis=1)
            nanind = np.where((np.isnan(OBSARR[vv,:,st])) | (OBSARR[vv,:,st] < minval) | (OBSARR[vv,:,st] > meanval[vv,st]+5*stdval[vv,st]))
            #nanind = np.where((np.isnan(OBSARR[vv,:,st])) | (OBSARR[vv,:,st] < minval))
            OBSARR[vv,nanind,st] = np.nan #set values below minval to nan
            nrobs = len(OBSARR[vv,:,st])
            nanfreq = len(nanind[0]) / nrobs *100
            print('NaN percentage is:'+str(nanfreq))
            if nanfreq > nanlimit:
               print (str(st)+' only has nan entries for '+variables[vv])
               BIAS[ii,vv,st] = np.nan
               RHO[ii,vv,st] = np.nan
               RATIO[ii,vv,st] = np.nan
               MAE[ii,vv,st] = np.nan               
            else: # in case there are sufficient data
               valind = np.where(~np.isnan(OBSARR[vv,:,st]))
               moddata = MODSTEP[vv,valind,st]
               obsdata = OBSARR[vv,valind,st]
               
               #BIAS[ii,vv,st] = (np.mean(moddata,axis=1)-np.mean(obsdata,axis=1))/np.mean(obsdata,axis=1)*100
               BIAS[ii,vv,st] = np.mean(moddata,axis=1)-np.mean(obsdata,axis=1)
               stepmat = np.concatenate((moddata,obsdata))
               stepmat = np.corrcoef(stepmat)
               RHO[ii,vv,st] = stepmat[0,1]
               del stepmat
               #RHO[ii,vv,st] = scipy.stats.spearmanr(MODSTEP[vv,valind,st],OBSARR[vv,valind,st],axis=1).correlation
               #RATIO[ii,vv,st] = (np.std(MODSTEP[vv,valind,st],axis=1)-np.std(OBSARR[vv,valind,st],axis=1))/np.std(OBSARR[vv,valind,st],axis=1)*100
               pctmod = np.percentile(np.squeeze(moddata), [75, 25])
               iqrmod = pctmod[0]-pctmod[1]
               pctobs = np.percentile(np.squeeze(obsdata), [75, 25])
               iqrobs = pctobs[0]-pctobs[1]
               RATIO[ii,vv,st] = iqrmod/iqrobs*100
               MAE[ii,vv,st] = np.mean(np.absolute(moddata-obsdata))
               
               #plot the time series
               series1 = np.transpose(np.stack((OBSARR[vv,:,st],MODARR[0,vv,:,st])))
               series2 = np.transpose(np.stack((OBSARR[vv,:,st],MODARR[selmodel,vv,:,st])))               
               fig, axes = plt.subplots(nrows=2, ncols=1, sharey=True)
               fig.suptitle(station_names_csv[st]+' '+variables[vv], fontsize=titlesize, fontweight='bold')               
               #set parameters for all axes
               plt.setp(axes, xticks=xtick, xticklabels=days_uniq)
               #plot first pair of time series
               plt.sca(axes[0])
               plt.title('MAE:'+str(np.around(MAE[0,vv,st], decimals=2))+' RHO:'+str(np.around(RHO[0,vv,st], decimals=2))+' BIAS:'+str(np.around(BIAS[0,vv,st],decimals=2)), fontsize=xticksize)
               plt.xticks(fontsize=xticksize)
               plt.yticks(fontsize=yticksize)
               plt.ylabel(ylabel, fontsize=yticksize) 
               plt.plot(series1)
               plt.legend([obslabel,oplabel], fontsize=xticksize)
               #plot second pair of time series
               plt.sca(axes[1])
               plt.title('MAE:'+str(np.around(MAE[1,vv,st], decimals=2))+' RHO:'+str(np.around(RHO[1,vv,st], decimals=2))+' BIAS:'+str(np.around(BIAS[1,vv,st],decimals=2)), fontsize=xticksize)
               plt.xticks(fontsize=xticksize)
               plt.yticks(fontsize=yticksize)
               plt.ylabel(ylabel, fontsize=yticksize)
               plt.plot(series2)
               plt.legend([obslabel,expname], fontsize=xticksize)
               #save figure
               fig.savefig(src_dirs[2]+'/series/'+variables[vv]+'/'+variables[vv]+'_'+str(st)+'_'+domain+'.png', dpi=300)
               plt.close(fig)
               plt.close()
               del(fig,axes)

    
    


