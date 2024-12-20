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
homedir = os.getenv("HOME")
lustre = os.getenv("LUSTRE")

## DEFINE USER OPTIONS
#startdate = '2018-01-17T03:00:00.000000000'
#enddate = '2018-01-31T23:00:00.000000000'
startdate = '2018-01-17 03:00:00'
#enddate = '2018-02-15 23:00:00'
enddate = '2018-01-30 23:00:00'
variables = ['PM10','O3','NO2','SO2'] #which variables are to be loaded
domain = 'fine' #target domain, MAKE SURE THAT THIS IS THE RIGHT DOMAIN
stationind = 0
expid = 1
#xlabel = ['6x8km_TROPOS', '4x4xkm_GS1_HTAP', '4x4km_GS2_HTAP', '6x8km_GS2_IG', '4x4km_GS2_HTAP_OFF', '6x8km_GS2_HTAP_ON']
xlabel = ['6x8_TROP_997', '4x4_GS1_HTAP_999', '4x4_GS2_HTAP_999', '6x8_GS1_HTAP_997']

## EXECUTE #############################################################

## load metadata from metadata.py
tarlat = [] #target latitudes, will be loaded by metadata.py
tarlon = [] #target longitudes, will be loaded by metadata.py
station_names_meta = [] #station names, will be loaded by metadata.py

execfile("metadata.py")

#for each day, load the model data at the gridboxes nearest to the station locations defined in tarlat y tarlon. Finally, concatenate the time series
#src_dirs = [lustre+'/SWEN/DATA/OBSDATA',lustre+'/SWEN/DATA/MODELDATA',lustre+'/OP/DATOS/RESULTADOS/chimere2017r4/exp1/validation',lustre+'/OP/DATOS/RESULTADOS/chimere2017r4/exp2/validation',lustre+'/OP/DATOS/RESULTADOS/chimere2017r4/exp3/validation',lustre+'/OP/DATOS/RESULTADOS/chimere2017r4/exp4/validation',lustre+'/OP/DATOS/RESULTADOS/chimere2017r4/exp5/validation']
src_dirs = [lustre+'/SWEN/DATA/OBSDATA',lustre+'/SWEN/DATA/MODELDATA', lustre+'/OP/DATOS/RESULTADOS/chimere2017r4/exp1/validation', lustre+'/OP/DATOS/RESULTADOS/chimere2017r4/exp2/validation', lustre+'/OP/DATOS/RESULTADOS/chimere2017r4/exp5/validation']
savedir = lustre+'/SWEN/val_chimere'

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
            valind=np.where(~np.isnan(OBSARR[vv,:,st]))
            BIAS[ii,vv,st] = (np.mean(MODSTEP[vv,valind,st],axis=1)-np.mean(OBSARR[vv,valind,st],axis=1))/np.mean(OBSARR[vv,valind,st],axis=1)*100
            stepmat = np.concatenate((MODSTEP[vv,valind,st],OBSARR[vv,valind,st]))
            stepmat = np.corrcoef(stepmat)
            RHO[ii,vv,st] = stepmat[0,1]
            del stepmat
            #RHO[ii,vv,st] = scipy.stats.spearmanr(MODSTEP[vv,valind,st],OBSARR[vv,valind,st],axis=1).correlation
            #RATIO[ii,vv,st] = (np.std(MODSTEP[vv,valind,st],axis=1)-np.std(OBSARR[vv,valind,st],axis=1))/np.std(OBSARR[vv,valind,st],axis=1)*100
            RATIO[ii,vv,st] = (np.subtract(*np.percentile(MODSTEP[vv,valind,st], [75, 25])) / np.subtract(*np.percentile(OBSARR[vv,valind,st], [75, 25])))*100
            MAE[ii,vv,st] = np.mean(np.absolute(MODSTEP[vv,valind,st]-OBSARR[vv,valind,st]))

#and plot the results
for vv in xrange(len(variables)):
    fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=False, sharey=False)
    xtick = range(len(src_dirs))
    xtick = xtick[1:]    
    fig.suptitle('Validation results for '+variables[vv], fontsize=10, fontweight='bold')
    stepmat = BIAS[:,vv,:]
    mask = ~np.isnan(stepmat)
    stepmat = [d[m] for d, m in zip(stepmat, mask)]
    ax1.boxplot(np.transpose(stepmat))    
    ax1.set_ylabel('BIAS (%)', fontsize=10, color='blue')
    ax1.grid(color='gray', linestyle='--', linewidth=0.5)
    stepmat = RHO[:,vv,:]
    mask = ~np.isnan(stepmat)
    stepmat = [d[m] for d, m in zip(stepmat, mask)]
    ax2.boxplot(np.transpose(stepmat))
    ax2.set_ylabel('CORR. COEFF.', fontsize=10, color='green')
    ax2.grid(color='gray', linestyle='--', linewidth=0.5)
    stepmat = MAE[:,vv,:]
    mask = ~np.isnan(stepmat)
    stepmat = [d[m] for d, m in zip(stepmat, mask)]
    ax3.boxplot(np.transpose(stepmat))
    ax3.set_ylabel('MAE', fontsize=10, color='red')
    ax3.set_xticks(xtick)
    ax3.set_xticklabels(xlabel, size=8)
    ax3.grid(color='gray', linestyle='--', linewidth=0.5)
    fig.savefig('boxplot_'+domain+'_'+variables[vv]+'.png', dpi=300)
    plt.close(fig)
    del(fig,ax1,ax2,ax3)    
    
    fig, (ax1, ax2) = plt.subplots(2, sharex=False, sharey=True)
    fig.suptitle('Validation results for '+variables[vv], fontsize=10, fontweight='bold')
    series1 = np.transpose(np.stack((OBSARR[vv,:,stationind],MODARR[0,vv,:,stationind])))
    series2 = np.transpose(np.stack((OBSARR[vv,:,stationind],MODARR[expid,vv,:,stationind])))
    ax1.plot(series1)
    ax2.plot(series2)
    fig.savefig('series_'+domain+'_'+variables[vv]+'.png', dpi=300)
    plt.close(fig)
    del(fig,ax1,ax2)
    
    

#plt.plot(np.transpose(ARRAY[(0,1),:,1]))
#plt.savefig('newop.png')
#plt.close()

#plt.plot(np.transpose(ARRAY[(0,2),:,1]))
#plt.savefig('oldop.png')
#plt.close()

    
    


