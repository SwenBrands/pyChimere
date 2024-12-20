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

valmode=str(sys.argv[1]) #temporal or spatial

homedir = os.getenv("HOME")
lustre = os.getenv("LUSTRE")

## DEFINE USER OPTIONS
startdate = '2018-04-01 04:00:00'
enddate = '2018-04-25 23:00:00'
variables = ['PM10','O3','NO2'] #which variables are to be loaded, there is a problem with outliers in SO2, correct this in future versions
domain = 'fine' #target domain, MAKE SURE THAT THIS IS THE RIGHT DOMAIN
stationind = 0
#xlabel = ['OP', '5343', '5564', '2314', '2315', '10344c', '11344', '10384']
xlabel = ['OP', '5564', '2314', '2315', '10344c', '11344', '10384', '10345']
nanlimit = 60 #in percent
titlesize = 5
labelsize = 4
minval = 1 #minimum value allowed in observations, set to nan if lower
RESpath = lustre+'/OP/DATOS/RESULTADOS/chimere2017r4'

## EXECUTE #############################################################

## load metadata from metadata.py
tarlat = [] #target latitudes, will be loaded by metadata.py
tarlon = [] #target longitudes, will be loaded by metadata.py
station_names_csv = [] #station names, will be loaded by metadata.py
execfile("metadata.py")

##for each day, load the model data at the gridboxes nearest to the station locations defined in tarlat y tarlon. Finally, concatenate the time series
#src_dirs = [lustre+'/SWEN/DATA/OBSDATA', lustre+'/SWEN/DATA/MODELDATA', RESpath+'/h5v3e4p3/validation', RESpath+'/h5v5e6p4/validation',
#            RESpath+'/h2v3e1p4/validation', RESpath+'/h2v3e1p5/validation', RESpath+'/h10v3e4p4c/validation', RESpath+'/h11v3e4p4/validation',
#            RESpath+'/h10v3e8p4/validation']
src_dirs = [lustre+'/SWEN/DATA/OBSDATA', lustre+'/SWEN/DATA/MODELDATA', RESpath+'/h5v5e6p4/validation',
            RESpath+'/h2v3e1p4/validation', RESpath+'/h2v3e1p5/validation', RESpath+'/h10v3e4p4c/validation', RESpath+'/h11v3e4p4/validation',
            RESpath+'/h10v3e8p4/validation', RESpath+'/h10v3e4p5/validation']

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
BIAS90 = np.zeros((len(MODARR),len(variables),ARRAY.shape[-1]))
RHO = np.zeros((len(MODARR),len(variables),ARRAY.shape[-1]))
RATIO = np.zeros((len(MODARR),len(variables),ARRAY.shape[-1]))
MAE = np.zeros((len(MODARR),len(variables),ARRAY.shape[-1]))
PCT5 = np.zeros((len(MODARR),len(variables),ARRAY.shape[-1]))
PCT90 = np.zeros((len(MODARR),len(variables),ARRAY.shape[-1]))
PCT95 = np.zeros((len(MODARR),len(variables),ARRAY.shape[-1]))
BIAS90 = np.zeros((len(MODARR),len(variables),ARRAY.shape[-1]))
for ii in xrange(len(MODARR)):
    MODSTEP = MODARR[ii,]
    for vv in xrange(len(variables)):
        for st in xrange(ARRAY.shape[-1]):
            #check validation mode
            if valmode == 'temporal': #check whether to validate time series or spatial fields
                print('temporal validation mode ...')                  
            elif valmode == 'spatial':
                print('spatial validation mode ....')
                moddata = np.transpose(MODARR)
                obsdata = np.transpose(OBSARR)
            else:
                raise Exception('check entry for <valmode>!')
            
            #calc nanmean and nanstd for outlier definition
            meanval = np.nanmean(OBSARR,axis=1)
            stdval = np.nanstd(OBSARR,axis=1)
            #nanind = np.where((np.isnan(OBSARR[vv,:,st])) | (OBSARR[vv,:,st] < minval) | (OBSARR[vv,:,st] > meanval[vv,st]+5*stdval[vv,st]))[0]
            nanind = np.where((np.isnan(OBSARR[vv,:,st])) | (OBSARR[vv,:,st] < minval))[0]
            OBSARR[vv,nanind,st] = np.nan #set values below minval to nan
            nrobs = len(OBSARR[vv,:,st])
            nanfreq = len(nanind) / nrobs *100
            print('NaN percentage is:'+str(nanfreq))
            if nanfreq > nanlimit:
               print (str(st)+' only has too many nan entries for '+variables[vv])
               BIAS[ii,vv,st] = np.nan
               RHO[ii,vv,st] = np.nan
               RATIO[ii,vv,st] = np.nan
               MAE[ii,vv,st] = np.nan
               PCT5[ii,vv,st] = np.nan
               PCT90[ii,vv,st] = np.nan
            else: # in case there are sufficient data
               moddata = np.squeeze(MODSTEP[vv,:,st])
               obsdata = np.squeeze(OBSARR[vv,:,st])
               valind = np.where(~np.isnan(obsdata))[0]
               moddata = moddata[valind]
               obsdata = obsdata[valind]               
               BIAS[ii,vv,st] = np.mean(moddata,axis=0)-np.mean(obsdata,axis=0)
               RHO[ii,vv,st] = scipy.stats.pearsonr(moddata, obsdata)[0]              
               pctmod = np.percentile(moddata, [75, 25])
               iqrmod = pctmod[0]-pctmod[1]
               pctobs = np.percentile(obsdata, [75, 25])
               iqrobs = pctobs[0]-pctobs[1]
               RATIO[ii,vv,st] = iqrmod/iqrobs-1
               MAE[ii,vv,st] = np.mean(np.absolute(moddata-obsdata))
               PCT5[ii,vv,st] = np.percentile(moddata,5)-np.percentile(obsdata,5)
               PCT5obs = np.percentile(obsdata,5)
               PCT90obs = np.percentile(obsdata,90)
               PCT95obs = np.percentile(obsdata,95)
               PCT5mod = np.percentile(moddata,5)
               PCT90mod = np.percentile(moddata,90)
               PCT95mod = np.percentile(moddata,95)
               PCT5[ii,vv,st] = PCT5mod-PCT5obs
               PCT90[ii,vv,st] = PCT90mod-PCT90obs
               PCT95[ii,vv,st] = PCT95mod-PCT95obs
               ind90obs = np.where(obsdata>PCT90obs)[0]
               ind90mod = np.where(moddata>PCT90mod)[0]
               BIAS90[ii,vv,st] = np.mean(moddata[ind90mod],axis=0)-np.mean(obsdata[ind90obs],axis=0)

xticks = range(len(src_dirs))
xticks = xticks[1:]
#and plot the results
for vv in xrange(len(variables)):
    fig = plt.plot([1,2,3,4])
    plt.title('Validation results for '+variables[vv], fontsize=8, fontweight='bold')
    
    plt.subplot(231)
    xtick = range(len(src_dirs))
    xtick = xtick[1:]   
    stepmat = BIAS[:,vv,:]
    mask = ~np.isnan(stepmat)
    stepmat = [d[m] for d, m in zip(stepmat, mask)]
    plt.boxplot(np.transpose(stepmat))    
    plt.title('BIAS', fontsize=titlesize, color='blue')
    plt.grid(color='gray', linestyle='--', linewidth=0.5)    
    plt.xticks(xtick, fontsize=labelsize)
    plt.yticks(fontsize=labelsize)
    locs, labels = plt.xticks()
    plt.xticks(locs, xlabel)
        
    plt.subplot(232)
    xtick = range(len(src_dirs))
    xtick = xtick[1:]    
    stepmat = RHO[:,vv,:]
    mask = ~np.isnan(stepmat)
    stepmat = [d[m] for d, m in zip(stepmat, mask)]
    plt.boxplot(np.transpose(stepmat))
    plt.title('CORR. COEFF.', fontsize=titlesize, color='blue')
    plt.grid(color='gray', linestyle='--', linewidth=0.5)   
    plt.xticks(xtick, fontsize=labelsize)
    plt.yticks(fontsize=labelsize)
    locs, labels = plt.xticks()
    plt.xticks(locs, xlabel)
            
    plt.subplot(233)
    xtick = range(len(src_dirs))
    xtick = xtick[1:]    
    stepmat = MAE[:,vv,:]
    mask = ~np.isnan(stepmat)
    stepmat = [d[m] for d, m in zip(stepmat, mask)]
    plt.boxplot(np.transpose(stepmat))
    plt.title('MAE', fontsize=titlesize, color='blue')
    plt.grid(color='gray', linestyle='--', linewidth=0.5)
    plt.xticks(xtick, fontsize=labelsize)
    plt.yticks(fontsize=labelsize)
    locs, labels = plt.xticks()
    plt.xticks(locs, xlabel)
        
    plt.subplot(234)
    xtick = range(len(src_dirs))
    xtick = xtick[1:]    
    stepmat = PCT5[:,vv,:]
    mask = ~np.isnan(stepmat)
    stepmat = [d[m] for d, m in zip(stepmat, mask)]
    plt.boxplot(np.transpose(stepmat))
    plt.title('PCT5', fontsize=titlesize, color='blue')
    plt.grid(color='gray', linestyle='--', linewidth=0.5)
    plt.xticks(xtick, fontsize=labelsize)
    plt.yticks(fontsize=labelsize)
    locs, labels = plt.xticks()
    plt.xticks(locs, xlabel)
    
    plt.subplot(235)
    xtick = range(len(src_dirs))
    xtick = xtick[1:]    
    stepmat = PCT90[:,vv,:]
    mask = ~np.isnan(stepmat)
    stepmat = [d[m] for d, m in zip(stepmat, mask)]
    plt.boxplot(np.transpose(stepmat))
    plt.title('PCT90', fontsize=titlesize, color='blue')
    plt.grid(color='gray', linestyle='--', linewidth=0.5)
    plt.xticks(xtick, fontsize=labelsize)
    plt.yticks(fontsize=labelsize)
    locs, labels = plt.xticks()
    plt.xticks(locs, xlabel)
    
    plt.subplot(236)
    xtick = range(len(src_dirs))
    xtick = xtick[1:]    
    stepmat = BIAS90[:,vv,:]
    mask = ~np.isnan(stepmat)
    stepmat = [d[m] for d, m in zip(stepmat, mask)]
    plt.boxplot(np.transpose(stepmat))
    plt.title('BIAS90', fontsize=titlesize, color='blue')
    plt.grid(color='gray', linestyle='--', linewidth=0.5)
    plt.xticks(xtick, fontsize=labelsize)
    plt.yticks(fontsize=labelsize)
    locs, labels = plt.xticks()
    plt.xticks(locs, xlabel)
    
    plt.savefig('./boxplots/boxplot_'+valmode+'_'+domain+'_'+variables[vv]+'.png', dpi=300)
    plt.close()
    del(fig)
