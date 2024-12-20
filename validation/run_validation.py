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
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import scipy.stats

valmode=str(sys.argv[1]) #temporal or spatial
aggval=str(sys.argv[2]) #mean, max, min, or instan
stationtype=str(sys.argv[3]) #background-rural, background-suburban, industry-urban, industry-suburban, industry-rural, traffic-urban, traffic-suburban, sea-side, all

homedir = os.getenv("HOME")
lustre = os.getenv("LUSTRE")

## DEFINE USER OPTIONS
startdate = '2018-06-20 04:00:00'
enddate = '2018-08-31 23:00:00'
variables = ['PM10','PM25','O3','NO2','SO2'] #which variables are to be loaded, there is a problem with outliers in SO2, correct this in future versions
domain = 'fine' #target domain, MAKE SURE THAT THIS IS THE RIGHT DOMAIN
#xlabel = ['TR', '2008', 'mx3', 'smx3', 'smx3ompi', 'm', 'msd', 'mhr', 'mhr10', 'mhrsd', 'mhr10sd', 'shr']
xlabel = ['CM10', 'CM20', 'CS10', 'CS20', 'FM10', 'FM20', 'FS10', 'FS20', 'FM10SD', 'FM20SD']
nanlimit = 60 #in percent
titlesize = 5
labelsize = 4
minval = 1 #minimum value allowed in observations, set to nan if lower
RESpath = lustre+'/OP/DATOS/RESULTADOS/chimere2017r4'
percentiles = [5, 95, 95]
fliers = {'markersize': 10, 'markeredgecolor': 'None', 'marker': '.'}
markersymbol = 'o'
nrstations = 4

## for each day, load the model data at the gridboxes nearest to the station locations defined in tarlat y tarlon. Finally, concatenate the time series
#src_dirs = [lustre+'/SWEN/DATA/OBSDATA', lustre+'/SWEN/DATA/MODELDATA', RESpath+'/jja_03/validation',
            #RESpath+'/jja_02/validation', RESpath+'/jja_05/validation', RESpath+'/jja_05_ompi/validation',
            #RESpath+'/jja_06/validation', RESpath+'/jja_07/validation', RESpath+'/jja_08/validation', RESpath+'/jja_10/validation',
            #RESpath+'/jja_09/validation', RESpath+'/jja_11/validation', RESpath+'/jja_12/validation']
src_dirs = [lustre+'/SWEN/DATA/OBSDATA', RESpath+'/jja_10_gal15r/validation', RESpath+'/jja_02_gal15r/validation',RESpath+'/jja_13_gal15r/validation', \
            RESpath+'/jja_12_gal15r/validation', RESpath+'/jja_10/validation', RESpath+'/jja_08/validation',  RESpath+'/jja_13/validation', \
            RESpath+'/jja_12/validation', RESpath+'/jja_11/validation', RESpath+'/jja_09/validation']


## EXECUTE #############################################################

## load metadata from metadata.py
tarlat = [] #target latitudes, will be loaded by metadata.py
tarlon = [] #target longitudes, will be loaded by metadata.py
station_names_csv = [] #station names, will be loaded by metadata.py
execfile("metadata.py")

#optionally choose subset of stations
if stationtype == 'background-rural':
    stationfilter = [3,10]
elif stationtype == 'background-suburban':
    stationfilter = [0,6,9]
elif stationtype == 'background':
    stationfilter = [0,3,6,9,10]
elif stationtype == 'traffic-urban':
    stationfilter = [1,2,4,5,7,8,43,44]
elif stationtype == 'traffic-suburban':
    stationfilter = [45,46]
elif stationtype == 'industry-urban':
    stationfilter = [18,32,33,36,48]
elif stationtype == 'industry-suburban':
    stationfilter = [11,12,16,19,23,26,38,42,47]
elif stationtype == 'industry-rural':
    stationfilter = [13,14,15,20,21,22,24,25,27,28,29,30,31,34,35,37,39,40]
elif stationtype == 'sea-side':
    stationfilter = [9,15]
elif stationtype == 'all':
    stationfilter = range(len(station_names_csv))
else:
    raise Exception('check entry for <stationtype>!')

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
del(vv,ii)        

#filter stations
ARRAY = ARRAY[:,:,:,stationfilter]
station_names_csv = [station_names_csv[ff] for ff in stationfilter]

#opionally calculate aggregated values
if aggval in ('min','max','mean'):
    print('INFO: daily '+aggval+' data are used...')
    ARR_agg = xr.DataArray(ARRAY, coords=[variables, src_dirs, dates, station_names_csv], dims=['variable', 'experiment', 'time', 'loc'], name='aggdata')
    #get daily dates 
    dates_day = pd.date_range(dates.time.values[0],dates.time.values[-1],freq='D').date.astype('str')

    #then find max values for each day
    ARRAY = np.zeros((ARRAY.shape[0],ARRAY.shape[1],len(dates_day),ARRAY.shape[3]))
    for tt in range(len(dates_day)):
        if aggval == 'max':
            ARRAY[:,:,tt,:] = np.max(ARR_agg.sel(time=dates_day[tt]),axis=2)
        elif aggval == 'min':
            ARRAY[:,:,tt,:] = np.min(ARR_agg.sel(time=dates_day[tt]),axis=2)
        elif aggval == 'mean':
            ARRAY[:,:,tt,:] = np.mean(ARR_agg.sel(time=dates_day[tt]),axis=2)
        else:
            raise Exception('check entry for <aggval>!')
elif aggval == 'hourly':
    print('INFO: hourly data are used...')
else:
    raise Exception('check entry for <aggval>!')

#start the verfication
OBSARR = ARRAY[:,0] #ARRAY containing observations (3D)
MODARR = ARRAY[:,1:] #ARRAY containing model output (4D since various model versions are considered

##check validation mode
if valmode == 'spatial':
    print('INFO: spatial validation mode, arrays are transposed!')
    MODARR = np.swapaxes(MODARR,3,2)
    OBSARR = np.swapaxes(OBSARR,1,2)
elif valmode == 'temporal':
    print('INFO: temporal validation mode, arrays are not transposed!')
else:
    raise Exception('check entry for <valmode>!')

for vv in range(len(variables)):
    obsdata = OBSARR[vv,]
    #find out columns with too many nans or outlier values
    maxumb = np.nanmean(obsdata, axis = 0) + 5*np.nanstd(obsdata, axis = 0)
    maxumb = np.tile(maxumb,[obsdata.shape[0],1])
    outlierind_obs = (obsdata < minval) | (obsdata > maxumb) | (np.isnan(obsdata))
    obsdata[outlierind_obs] = np.nan
    #delete stations with too many outliers
    nanperc_obs = np.sum(outlierind_obs, axis=0)/obsdata.shape[0]*100
    retainind_obs = np.where(nanperc_obs < nanlimit)[0]
    obsdata = obsdata[:,retainind_obs]
    
    BIAS = np.zeros((MODARR.shape[1],len(retainind_obs)))
    MAE = np.zeros((MODARR.shape[1],len(retainind_obs)))
    HIBIAS = np.zeros((MODARR.shape[1],len(retainind_obs)))
    LOBIAS = np.zeros((MODARR.shape[1],len(retainind_obs)))
    RHO = np.zeros((MODARR.shape[1],len(retainind_obs)))
    PCTDIFF = np.zeros((MODARR.shape[1],len(percentiles),len(retainind_obs)))    
    for ii in range(MODARR.shape[1]):
        moddata = MODARR[vv,ii,]
        moddata[outlierind_obs] = np.nan
        moddata = moddata[:,retainind_obs]           
        BIAS[ii,] = np.nanmean(moddata, axis=0) - np.nanmean(obsdata, axis=0)
        MAE[ii,] = np.nanmean(np.absolute(moddata-obsdata),axis=0)
        PCTobs = np.nanpercentile(obsdata,percentiles,axis=0)
        PCTmod = np.nanpercentile(moddata,percentiles,axis=0)
        PCTDIFF[ii,] = PCTmod - PCTobs        
        hi_index = np.where(obsdata>PCTobs[2,])[0]
        lo_index = np.where(obsdata<PCTobs[0,])[0]        
        #indexmod = np.where(moddata>PCTmod[2,])[0]        
        HIBIAS[ii,] = np.nanmean(moddata[hi_index,], axis=0) - np.nanmean(obsdata[hi_index,], axis=0)
        LOBIAS[ii,] = np.nanmean(moddata[lo_index,], axis=0) - np.nanmean(obsdata[lo_index,], axis=0)
        #vectorize RHO calculation in future versions!
        for st in range(obsdata.shape[1]):
            nanind = np.where(~np.isnan(obsdata[:,st]))[0]
            RHO[ii,st] = scipy.stats.pearsonr(moddata[nanind,st], obsdata[nanind,st])[0]    
    
    #plot the results
    fig = plt.figure()
    fig.suptitle('Verfication for '+variables[vv]+' '+aggval+' '+stationtype, fontsize=8, fontweight='bold')
    plt.subplot(221)
    xtick = range(len(src_dirs))
    xtick = xtick[1:]
    retainind=np.where(~np.any(np.isnan(BIAS),axis=0))
    BIAS = BIAS[:,retainind[0]]
    plt.plot([0.5,MODARR.shape[1]+0.5],[0,0],'r')
    if len(station_names_csv) < nrstations:
        X = np.transpose(np.tile(range(MODARR.shape[1]),(len(station_names_csv),1)))
        plt.plot(X+1, BIAS, linestyle='None', marker=markersymbol)
        #plt.legend(station_names_csv)
    else:
        plt.boxplot(np.transpose(BIAS), flierprops=fliers)
        plt.ylim([np.min(BIAS), np.max(BIAS)])
        
    plt.title('BIAS', fontsize=titlesize, color='blue')
    plt.grid(color='gray', linestyle='--', linewidth=0.5)    
    plt.xticks(xtick, fontsize=labelsize)
    plt.yticks(fontsize=labelsize)
    locs, labels = plt.xticks()
    plt.xticks(locs, xlabel)
    plt.xlim([0.5, MODARR.shape[1]+0.5])       
        
    plt.subplot(222)
    xtick = range(len(src_dirs))
    xtick = xtick[1:]
    retainind=np.where(~np.any(np.isnan(RHO),axis=0))
    RHO = RHO[:,retainind[0]]
    plt.plot([0.5,MODARR.shape[1]+0.5],[0,0],'r')
    if len(station_names_csv) < nrstations:
        X = np.transpose(np.tile(range(MODARR.shape[1]),(len(station_names_csv),1)))
        plt.plot(X+1, RHO, linestyle='None', marker=markersymbol)
        #plt.legend(station_names_csv)
    else:
        plt.boxplot(np.transpose(RHO), flierprops=fliers)
        plt.ylim([np.min(RHO), np.max(RHO)])
    
    plt.title('CORR. COEFF.', fontsize=titlesize, color='blue')
    plt.grid(color='gray', linestyle='--', linewidth=0.5)   
    plt.xticks(xtick, fontsize=labelsize)
    plt.yticks(fontsize=labelsize)
    locs, labels = plt.xticks()
    plt.xticks(locs, xlabel)
    plt.xlim([0.5, MODARR.shape[1]+0.5])
            
    plt.subplot(223)
    xtick = range(len(src_dirs))
    xtick = xtick[1:]
    retainind=np.where(~np.any(np.isnan(MAE),axis=0))
    MAE = MAE[:,retainind[0]]    
    if len(station_names_csv) < nrstations:
        X = np.transpose(np.tile(range(MODARR.shape[1]),(len(station_names_csv),1)))
        plt.plot(X+1, MAE, linestyle='None', marker=markersymbol)
        #plt.legend(station_names_csv)
    else:
        plt.boxplot(np.transpose(MAE), flierprops=fliers)
        plt.ylim([np.min(MAE), np.max(MAE)])
    
    plt.title('MAE', fontsize=titlesize, color='blue')
    plt.grid(color='gray', linestyle='--', linewidth=0.5)
    plt.xticks(xtick, fontsize=labelsize)
    plt.yticks(fontsize=labelsize)
    locs, labels = plt.xticks()
    plt.xticks(locs, xlabel)
    plt.xlim([0.5, MODARR.shape[1]+0.5])
    
    ##plot BIAS for low observed values
    #plt.subplot(223)
    #xtick = range(len(src_dirs))
    #xtick = xtick[1:]
    #retainind=np.where(~np.any(np.isnan(LOBIAS),axis=0))
    #LOBIAS = LOBIAS[:,retainind[0]]
    #plt.boxplot(np.transpose(LOBIAS), flierprops=fliers)
    #plt.title('BIAS'+str(percentiles[0]), fontsize=titlesize, color='blue')
    #plt.grid(color='gray', linestyle='--', linewidth=0.5)
    #plt.xticks(xtick, fontsize=labelsize)
    #plt.yticks(fontsize=labelsize)
    #locs, labels = plt.xticks()
    #plt.xticks(locs, xlabel)
    
    ##plot differences in the FIRST percentile
    #plt.subplot(234)
    #xtick = range(len(src_dirs))
    #xtick = xtick[1:]
    #plt.boxplot(np.transpose(PCTDIFF[:,0,:]), flierprops=fliers)
    #plt.title('PCT'+str(percentiles[0]), fontsize=titlesize, color='blue')
    #plt.grid(color='gray', linestyle='--', linewidth=0.5)
    #plt.xticks(xtick, fontsize=labelsize)
    #plt.yticks(fontsize=labelsize)
    #locs, labels = plt.xticks()
    #plt.xticks(locs, xlabel)
    
    ##plot differences in the SECOND percentile
    #plt.subplot(235)
    #xtick = range(len(src_dirs))
    #xtick = xtick[1:]    
    #plt.boxplot(np.transpose(PCTDIFF[:,1,:]), flierprops=fliers)
    #plt.title('PCT'+str(percentiles[1]), fontsize=titlesize, color='blue')
    #plt.grid(color='gray', linestyle='--', linewidth=0.5)
    #plt.xticks(xtick, fontsize=labelsize)
    #plt.yticks(fontsize=labelsize)
    #locs, labels = plt.xticks()
    #plt.xticks(locs, xlabel)
    
    #plot BIAS for high observed values    
    plt.subplot(224)
    xtick = range(len(src_dirs))
    xtick = xtick[1:]
    retainind=np.where(~np.any(np.isnan(HIBIAS),axis=0))
    HIBIAS = HIBIAS[:,retainind[0]]
    plt.plot([0.5,MODARR.shape[1]+0.5],[0,0],'r')
    if len(station_names_csv) < nrstations:
        X = np.transpose(np.tile(range(MODARR.shape[1]),(len(station_names_csv),1)))
        plt.plot(X+1, HIBIAS, linestyle='None', marker=markersymbol)
        #plt.legend(station_names_csv)
    else:
        plt.boxplot(np.transpose(HIBIAS), flierprops=fliers)
        plt.ylim([np.min(HIBIAS), np.max(HIBIAS)])
    
    plt.title('BIAS'+str(percentiles[2]), fontsize=titlesize, color='blue')
    #plt.title('PCT'+str(percentiles[2]), fontsize=titlesize, color='blue')
    plt.grid(color='gray', linestyle='--', linewidth=0.5)
    plt.xticks(xtick, fontsize=labelsize)
    plt.yticks(fontsize=labelsize)
    locs, labels = plt.xticks()
    plt.xticks(locs, xlabel)
    plt.xlim([0.5, MODARR.shape[1]+0.5])
    
    plt.savefig('./boxplots/'+stationtype+'/boxplot_'+valmode+'_'+domain+'_'+variables[vv]+'_'+aggval+'_'+stationtype+'.png', dpi=300)
    plt.close(fig)
    del(fig)
