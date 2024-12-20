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
import scipy.stats
import itertools
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

expname=str(sys.argv[1]) #name of the experiment

homedir = os.getenv("HOME")
lustre = os.getenv("LUSTRE")

## DEFINE USER OPTIONS
#startdate = '2018-06-20 04:00:00'
#enddate = '2018-08-31 23:00:00'
startdate = '2018-06-21 04:00:00'
enddate = '2018-08-21 23:00:00'
variables = ['PM10','PM25','O3','NO2','SO2'] #which variables are to be loaded
domain = 'fine' #target domain, MAKE SURE THAT THIS IS THE RIGHT DOMAIN
selmodel = 1
nanlimit = 60 #in percent
titlesize = 8
xticksize = 3
yticksize = 7
obslabel = 'OBS'
oplabel = 'OP'
ylabel = 'ug/m3'
minval = 1 #minimum value allowed in observations, set to nan if lower
percentiles = [5, 95, 95]
linewidth = 1

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
ARRAY = np.zeros((len(variables),len(src_dirs),len(dummy),dummy.shape[1])) 
for vv in range(len(variables)):
    for ii in range(len(src_dirs)):        
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
        ARRAY[vv,ii,:,:] = np.expand_dims(data.values,0)
        
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
OBSARR = ARRAY[:,0] #ARRAY containing observations (3D)
MODARR = ARRAY[:,1:] #ARRAY containing model output (4D since various model versions are considered
for vv in range(len(variables)):
    obsdata = OBSARR[vv,]
    #find outliers and set to nan    
    maxumb = np.nanmean(obsdata, axis = 0) + 5*np.nanstd(obsdata, axis = 0)
    maxumb = np.tile(maxumb,[obsdata.shape[0],1])
    nanind = (obsdata < minval) | (obsdata > maxumb) | (np.isnan(obsdata))
    obsdata[nanind] = np.nan
    #delete stations with too many outliers
    nanperc = np.sum(nanind, axis=0)/obsdata.shape[0]*100
    retainind = np.where(nanperc < nanlimit)[0]
    #filter station names
    station_names = [station_names_csv[ii] for ii in retainind]
    obsdata = obsdata[:,retainind]
    BIAS = np.zeros((MODARR.shape[1],len(retainind)))
    MAE = np.zeros((MODARR.shape[1],len(retainind)))
    HIBIAS = np.zeros((MODARR.shape[1],len(retainind)))
    LOBIAS = np.zeros((MODARR.shape[1],len(retainind)))
    RHO = np.zeros((MODARR.shape[1],len(retainind)))
    PCTdiff = np.zeros((MODARR.shape[1],len(percentiles),len(retainind)))
    for ii in range(MODARR.shape[1]):
        moddata = MODARR[vv,ii,]
        moddata = moddata[:,retainind]            
        BIAS[ii,] = np.nanmean(moddata, axis=0) - np.nanmean(obsdata, axis=0)
        MAE[ii,] = np.nanmean(np.absolute(moddata-obsdata),axis=0)
        PCTobs = np.nanpercentile(obsdata,percentiles,axis=0)
        PCTmod = np.nanpercentile(moddata,percentiles,axis=0)
        PCTdiff[ii,] = PCTmod - PCTobs        
        hi_index = np.where(obsdata>PCTobs[2,])[0]
        lo_index = np.where(obsdata<PCTobs[0,])[0]        
        #indexmod = np.where(moddata>PCTmod[2,])[0]        
        HIBIAS[ii,] = np.nanmean(moddata[hi_index,], axis=0) - np.nanmean(obsdata[hi_index,], axis=0)
        LOBIAS[ii,] = np.nanmean(moddata[lo_index,], axis=0) - np.nanmean(obsdata[lo_index,], axis=0)
        #vectorize RHO calculation in future versions!
        for st in range(obsdata.shape[1]):
            nanind = np.where(~np.isnan(obsdata[:,st]))[0]
            RHO[ii,st] = scipy.stats.pearsonr(moddata[nanind,st], obsdata[nanind,st])[0]
               
    #plot the time series
    MODARRstep = MODARR[vv,:,:,retainind]
    for st in range(obsdata.shape[1]):
        series1 = np.transpose(np.stack((obsdata[:,st],MODARRstep[st,0,:])))
        series2 = np.transpose(np.stack((obsdata[:,st],MODARRstep[st,selmodel,:])))               
        fig, axes = plt.subplots(nrows=2, ncols=1, sharey=True)
        fig.suptitle(station_names[st]+' '+variables[vv], fontsize=titlesize, fontweight='bold')               
        #set parameters for all axes
        plt.setp(axes, xticks=xtick, xticklabels=days_uniq)
        #plot first pair of time series
        plt.sca(axes[0])
        plt.title('MAE:'+str(np.around(MAE[0,st], decimals=2))+'  RHO:'+str(np.around(RHO[0,st], decimals=2))+'  BIAS:'+str(np.around(BIAS[0,st],decimals=2)) +'  LOBIAS:'+str(np.around(LOBIAS[0,st],decimals=2)) +'  HIBIAS:'+str(np.around(HIBIAS[0,st],decimals=2)), fontsize=xticksize)
        plt.xticks(fontsize=xticksize)
        plt.yticks(fontsize=yticksize)
        plt.ylabel(ylabel, fontsize=yticksize) 
        plt.plot(series1,lw=linewidth)
        plt.legend([obslabel,oplabel], fontsize=xticksize)
        plt.xlim((xtick[0], xtick[-1]))
        #plot second pair of time series
        plt.sca(axes[1])
        plt.title('MAE:'+str(np.around(MAE[1,st], decimals=2))+'  RHO:'+str(np.around(RHO[1,st], decimals=2))+'  BIAS:'+str(np.around(BIAS[1,st],decimals=2)) +'  LOBIAS:'+str(np.around(LOBIAS[1,st],decimals=2)) +'  HIBIAS:'+str(np.around(HIBIAS[1,st],decimals=2)), fontsize=xticksize)
        plt.xticks(fontsize=xticksize)
        plt.yticks(fontsize=yticksize)
        plt.ylabel(ylabel, fontsize=yticksize)
        plt.plot(series2,lw=linewidth)
        plt.legend([obslabel,expname], fontsize=xticksize)
        plt.xlim((xtick[0], xtick[-1]))
        #save figure
        fig.savefig(src_dirs[2]+'/series/'+variables[vv]+'/'+variables[vv]+'_'+station_names[st]+'_'+domain+'.png', dpi=300)
        plt.close(fig)
        plt.close()
        del(fig,axes)

    
    


