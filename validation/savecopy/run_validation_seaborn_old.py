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
import seaborn as sns

valmode=str(sys.argv[1]) #temporal or spatial
aggval=str(sys.argv[2]) #mean, max, min, or instan
stationtype=str(sys.argv[3]) #background-rural, background-suburban, industry-urban, industry-suburban, industry-rural, traffic-urban, traffic-suburban, sea-side, traffic, all

homedir = os.getenv("HOME")
lustre = os.getenv("LUSTRE")

## DEFINE USER OPTIONS
startdate = '2018-06-20 04:00:00'
enddate = '2018-08-31 23:00:00'
#variables = ['PM25','O3','NO2','SO2','PM10'] #which variables are to be loaded, there is a problem with outliers in SO2, correct this in future versions
variables = ['PM25']
domain = 'fine' #target domain, MAKE SURE THAT THIS IS THE RIGHT DOMAIN
#xlabel = ['TR', '2008', 'mx3', 'smx3', 'smx3ompi', 'm', 'msd', 'mhr', 'mhr10', 'mhrsd', 'mhr10sd', 'shr']
#xlabel = ['CM10', 'CM20', 'CS10', 'CS20', 'FM10', 'FM20', 'FS10', 'FS20', 'FM10D', 'FM20D']
xlabel = ['CM10', 'CM20', 'CS10', 'CS20', 'FM10', 'FM20', 'FS10', 'FS20']
hores = [xlabel[ii][0] for ii in range(len(xlabel))]
layers = [xlabel[ii][2:4] for ii in range(len(xlabel))] 
nanlimit = 20 #in percent, 60 in former versions
titlesize = 5
labelsize = 4
minval = 0 #minimum value allowed in observations, set to nan if lower
RESpath = lustre+'/OP/DATOS/RESULTADOS/chimere2017r4'
percentiles = [5, 95, 95]
fliers = {'markersize': 10, 'markeredgecolor': 'None', 'marker': '.'}
markersymbol = 'o'
nrstations = 5
groupby = 'hores'
palette = ["m", "g"]
scalefactor = 0.5
gridstyle = 'whitegrid'
legend = False
kind = 'box'
nr_rows = 3
nr_cols = 2
extension = '.pdf'
fliersize = 2
plottitle = 'yes'

## for each day, load the model data at the gridboxes nearest to the station locations defined in tarlat y tarlon. Finally, concatenate the time series
#src_dirs = [lustre+'/SWEN/DATA/OBSDATA', lustre+'/SWEN/DATA/MODELDATA', RESpath+'/jja_03/validation',
            #RESpath+'/jja_02/validation', RESpath+'/jja_05/validation', RESpath+'/jja_05_ompi/validation',
            #RESpath+'/jja_06/validation', RESpath+'/jja_07/validation', RESpath+'/jja_08/validation', RESpath+'/jja_10/validation',
            #RESpath+'/jja_09/validation', RESpath+'/jja_11/validation', RESpath+'/jja_12/validation']
#src_dirs = [lustre+'/SWEN/DATA/OBSDATA', RESpath+'/jja_10_gal15r/validation', RESpath+'/jja_02_gal15r/validation',RESpath+'/jja_13_gal15r/validation', \
            #RESpath+'/jja_12_gal15r/validation', RESpath+'/jja_10/validation', RESpath+'/jja_08/validation',  RESpath+'/jja_13/validation', \
            #RESpath+'/jja_12/validation', RESpath+'/jja_11/validation', RESpath+'/jja_09/validation']
src_dirs = [lustre+'/SWEN/DATA/OBSDATA', RESpath+'/jja_10_gal15r/validation', RESpath+'/jja_02_gal15r/validation',RESpath+'/jja_13_gal15r/validation', \
            RESpath+'/jja_12_gal15r/validation', RESpath+'/jja_10/validation', RESpath+'/jja_08/validation',  RESpath+'/jja_13/validation', \
            RESpath+'/jja_12/validation']


## EXECUTE #############################################################
sns.set(context='paper', font_scale = scalefactor, style=gridstyle)
#sns.set_palette(palette, n_colors=None, desat=1, color_codes=False)
sns.set_palette(palette)
sns.set_style("ticks", {"xtick.major.size": 4, "ytick.major.size": 4}) #seems to have no effect

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
elif stationtype == 'traffic-urban':
    stationfilter = [1,2,4,5,7,8,43,44]
elif stationtype == 'traffic-suburban':
    stationfilter = [45,46]
elif stationtype == 'traffic':
    stationfilter = [1,2,4,5,7,8,43,44,45,46]
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
            ARRAY[:,:,tt,:] = np.amax(ARR_agg.sel(time=dates_day[tt]),axis=2)
        elif aggval == 'min':
            ARRAY[:,:,tt,:] = np.amin(ARR_agg.sel(time=dates_day[tt]),axis=2)
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
    station_names = [station_names_csv[ii] for ii in retainind_obs]
    
    BIAS = np.zeros((MODARR.shape[1],len(retainind_obs)))
    MAE = np.zeros((MODARR.shape[1],len(retainind_obs)))
    RATIO = np.zeros((MODARR.shape[1],len(retainind_obs)))
    HIBIAS_obs = np.zeros((MODARR.shape[1],len(retainind_obs)))
    LOBIAS_obs = np.zeros((MODARR.shape[1],len(retainind_obs)))
    HIBIAS_mod = np.zeros((MODARR.shape[1],len(retainind_obs)))
    LOBIAS_mod = np.zeros((MODARR.shape[1],len(retainind_obs)))
    RHO = np.zeros((MODARR.shape[1],len(retainind_obs)))
    PCTDIFF = np.zeros((len(percentiles),MODARR.shape[1],len(retainind_obs)))
    for ii in range(MODARR.shape[1]):
        moddata = MODARR[vv,ii,]
        moddata[outlierind_obs] = np.nan
        moddata = moddata[:,retainind_obs]           
        BIAS[ii,:] = np.nanmean(moddata, axis=0) - np.nanmean(obsdata, axis=0)
        RATIO[ii,:] = np.nanstd(moddata, axis=0) / np.nanstd(obsdata, axis=0)
        MAE[ii,:] = np.nanmean(np.absolute(moddata-obsdata),axis=0)
        PCTobs = np.nanpercentile(obsdata,percentiles,axis=0)
        PCTmod = np.nanpercentile(moddata,percentiles,axis=0)
        PCTDIFF[:,ii,:] = PCTmod - PCTobs
        hi_index_obs = np.where(obsdata>PCTobs[2,])[0]
        lo_index_obs = np.where(obsdata<PCTobs[0,])[0]
        hi_index_mod = np.where(moddata>PCTmod[2,])[0]
        lo_index_mod = np.where(moddata<PCTmod[0,])[0] 
        #indexmod = np.where(moddata>PCTmod[2,])[0]        
        HIBIAS_obs[ii,:] = np.nanmean(moddata[hi_index_obs,], axis=0) - np.nanmean(obsdata[hi_index_obs,], axis=0)
        LOBIAS_obs[ii,:] = np.nanmean(moddata[lo_index_obs,], axis=0) - np.nanmean(obsdata[lo_index_obs,], axis=0)
        HIBIAS_mod[ii,:] = np.nanmean(moddata[hi_index_mod,], axis=0) - np.nanmean(obsdata[hi_index_mod,], axis=0)
        LOBIAS_mod[ii,:] = np.nanmean(moddata[lo_index_mod,], axis=0) - np.nanmean(obsdata[lo_index_mod,], axis=0)
        ##vectorize RHO calculation in future versions!
        #corrme = pd.DataFrame(np.concatenate((obsdata,moddata),axis=1))        
        for st in range(obsdata.shape[1]):
            nanind = np.where(~np.isnan(obsdata[:,st]))[0]
            RHO[ii,st] = scipy.stats.pearsonr(moddata[nanind,st], obsdata[nanind,st])[0]
    
    #plot the results
    fig, axs = plt.subplots(nr_rows, nr_cols, squeeze=False)
    if plottitle == 'yes':
        fig.suptitle('Verfication for '+variables[vv]+' '+aggval+' '+stationtype, fontsize=8, fontweight='bold')
    
    retainind=np.where(~np.any(np.isnan(BIAS),axis=0))[0]
    BIAS = BIAS[:,retainind]
    legend_bias = [station_names[ii] for ii in retainind]
    if len(station_names) < nrstations:
        axs[0][0].plot([0.5,MODARR.shape[1]+0.5],[0,0],'r')
        X = np.transpose(np.tile(range(MODARR.shape[1]),(len(station_names),1)))
        axs[0][0].plot(X+1, BIAS, linestyle='None', marker=markersymbol)
        axs[0][0].legend(legend_bias)
    else:        
        explabels = np.repeat(xlabel,len(retainind))
        horeslabels = np.repeat(hores,len(retainind)) 
        layerslabels = np.repeat(layers,len(retainind)) 
        BIASpd = pd.DataFrame(dict(bias=BIAS.ravel(),experiment=explabels,hores=horeslabels,layers=layerslabels))
        sns.catplot(kind=kind, y='experiment', x='bias', data=BIASpd, hue=groupby, orient='h', legend=legend, fliersize=fliersize, ax = axs[0][0])
        axs[0][0].plot([0,0],[0,MODARR.shape[1]+2],'r')
        #sns.despine(trim=True)
        
    retainind=np.where(~np.any(np.isnan(RHO),axis=0))[0]
    RHO = RHO[:,retainind]
    legend_rho = [station_names[ii] for ii in retainind]
    if len(station_names) < nrstations:
        #axs[0][1].plot([0.5,MODARR.shape[1]+0.5],[0,0],'r')
        X = np.transpose(np.tile(range(MODARR.shape[1]),(len(station_names),1)))
        axs[0][1].plot(X+1, RHO, linestyle='None', marker=markersymbol)
    else:
        explabels = np.repeat(xlabel,len(retainind))
        horeslabels = np.repeat(hores,len(retainind)) 
        layerslabels = np.repeat(layers,len(retainind)) 
        RHOpd = pd.DataFrame(dict(r=RHO.ravel(),experiment=explabels,hores=horeslabels,layers=layerslabels))
        sns.catplot(kind=kind, y='experiment', x='r', data=RHOpd, hue=groupby, orient='h', legend=legend, fliersize=fliersize, ax = axs[0][1])
        #axs[0][1].plot([0,0],[0,MODARR.shape[1]+2],'r')

    retainind=np.where(~np.any(np.isnan(RATIO),axis=0))[0]
    RATIO = RATIO[:,retainind]
    legend_ratio = [station_names[ii] for ii in retainind]
    if len(station_names) < nrstations:
        #axs[1][0].plot([0.5,MODARR.shape[1]+0.5],[0,0],'r')
        X = np.transpose(np.tile(range(MODARR.shape[1]),(len(station_names),1)))
        axs[1][0].plot(X+1, RATIO, linestyle='None', marker=markersymbol)        
    else:
        explabels = np.repeat(xlabel,len(retainind))
        horeslabels = np.repeat(hores,len(retainind)) 
        layerslabels = np.repeat(layers,len(retainind)) 
        RATIOpd = pd.DataFrame(dict(ratio=RATIO.ravel(),experiment=explabels,hores=horeslabels,layers=layerslabels))
        sns.catplot(kind=kind, y='experiment', x='ratio', data=RATIOpd, hue=groupby, orient='h', legend=legend, fliersize=fliersize, ax = axs[1][0])
        axs[1][0].plot([1,1],[0,MODARR.shape[1]+2],'r')
            
    retainind=np.where(~np.any(np.isnan(MAE),axis=0))[0]
    MAE = MAE[:,retainind]
    legend_mae = [station_names[ii] for ii in retainind]
    if len(station_names) < nrstations:
        #axs[1][1].plot([0.5,MODARR.shape[1]+0.5],[0,0],'r')
        X = np.transpose(np.tile(range(MODARR.shape[1]),(len(station_names),1)))
        axs[1][1].plot(X+1, MAE, linestyle='None', marker=markersymbol)        
    else:
        explabels = np.repeat(xlabel,len(retainind))
        horeslabels = np.repeat(hores,len(retainind)) 
        layerslabels = np.repeat(layers,len(retainind)) 
        MAEpd = pd.DataFrame(dict(mae=MAE.ravel(),experiment=explabels,hores=horeslabels,layers=layerslabels))
        sns.catplot(kind=kind, y='experiment', x='mae', data=MAEpd, hue=groupby, orient='h', legend=legend, fliersize=fliersize, ax = axs[1][1])
        #axs[1][0].plot([0,0],[0,MODARR.shape[1]+2],'r')
    
    if aggval in ('mean','max'):
        retainind=np.where(~np.any(np.isnan(HIBIAS_obs),axis=0))[0]
        HIBIAS_obs = HIBIAS_obs[:,retainind]
        legend_hibias_obs = [station_names[ii] for ii in retainind]
        if len(station_names) < nrstations:
            axs[2][0].plot([0.5,MODARR.shape[1]+0.5],[0,0],'r')
            X = np.transpose(np.tile(range(MODARR.shape[1]),(len(station_names),1)))
            axs[2][0].plot(X+1, HIBIAS_obs, linestyle='None', marker=markersymbol)        
        else:
            explabels = np.repeat(xlabel,len(retainind))
            horeslabels = np.repeat(hores,len(retainind)) 
            layerslabels = np.repeat(layers,len(retainind)) 
            HIBIAS_obs_pd = pd.DataFrame(dict(hibias_obs=HIBIAS_obs.ravel(),experiment=explabels,hores=horeslabels,layers=layerslabels))
            sns.catplot(kind=kind, y='experiment', x='hibias_obs', data=HIBIAS_obs_pd, hue=groupby, orient='h', legend=legend, fliersize=fliersize, ax = axs[2][0])
            axs[2][0].plot([0,0],[0,MODARR.shape[1]+2],'r')
        
        retainind=np.where(~np.any(np.isnan(HIBIAS_mod),axis=0))[0]
        HIBIAS_mod = HIBIAS_mod[:,retainind]
        legend_hibias_mod = [station_names[ii] for ii in retainind]
        if len(station_names) < nrstations:
            axs[2][1].plot([0.5,MODARR.shape[1]+0.5],[0,0],'r')
            X = np.transpose(np.tile(range(MODARR.shape[1]),(len(station_names),1)))
            axs[2][1].plot(X+1, HIBIAS_mod, linestyle='None', marker=markersymbol)
        else:
            explabels = np.repeat(xlabel,len(retainind))
            horeslabels = np.repeat(hores,len(retainind)) 
            layerslabels = np.repeat(layers,len(retainind)) 
            HIBIAS_mod_pd = pd.DataFrame(dict(hibias_mod=HIBIAS_mod.ravel(),experiment=explabels,hores=horeslabels,layers=layerslabels))
            sns.catplot(kind=kind, y='experiment', x='hibias_mod', data=HIBIAS_mod_pd, hue=groupby, orient='h', legend=legend, fliersize=fliersize, ax = axs[2][1])
            axs[2][1].plot([0,0],[0,MODARR.shape[1]+2],'r')
    elif aggval == 'min':
        retainind=np.where(~np.any(np.isnan(LOBIAS_obs),axis=0))[0]
        LOBIAS_obs = LOBIAS_obs[:,retainind]
        legend_lobias_obs = [station_names[ii] for ii in retainind]
        if len(station_names) < nrstations:
            plt.plot([0.5,MODARR.shape[1]+0.5],[0,0],'r')
            X = np.transpose(np.tile(range(MODARR.shape[1]),(len(station_names),1)))
            plt.plot(X+1, LOBIAS_obs, linestyle='None', marker=markersymbol)        
        else:
            explabels = np.repeat(xlabel,len(retainind))
            horeslabels = np.repeat(hores,len(retainind)) 
            layerslabels = np.repeat(layers,len(retainind)) 
            LOBIAS_obs_pd = pd.DataFrame(dict(lobias_obs=LOBIAS_obs.ravel(),experiment=explabels,hores=horeslabels,layers=layerslabels))
            sns.catplot(kind=kind, y='experiment', x='lobias_obs', data=LOBIAS_obs_pd, hue=groupby, orient='h', legend=legend, fliersize=fliersize, ax = axs[2][0])
            axs[2][0].plot([0,0],[0,MODARR.shape[1]+2],'r')
        
        retainind=np.where(~np.any(np.isnan(LOBIAS_mod),axis=0))[0]
        LOBIAS_mod = LOBIAS_mod[:,retainind]
        legend_lobias_mod = [station_names[ii] for ii in retainind]
        if len(station_names) < nrstations:
            plt.plot([0.5,MODARR.shape[1]+0.5],[0,0],'r')
            X = np.transpose(np.tile(range(MODARR.shape[1]),(len(station_names),1)))
            plt.plot(X+1, LOBIAS_mod, linestyle='None', marker=markersymbol)        
        else:
            explabels = np.repeat(xlabel,len(retainind))
            horeslabels = np.repeat(hores,len(retainind)) 
            layerslabels = np.repeat(layers,len(retainind)) 
            LOBIAS_mod_pd = pd.DataFrame(dict(lobias_mod=LOBIAS_mod.ravel(),experiment=explabels,hores=horeslabels,layers=layerslabels))
            sns.catplot(kind=kind, y='experiment', x='lobias_mod', data=LOBIAS_mod_pd, hue=groupby, orient='h', legend=legend, fliersize=fliersize, ax = axs[2][1])
            axs[2][1].plot([0,0],[0,MODARR.shape[1]+2],'r')
    else:
        raise Exception('check entry for <aggval>!')
    
    fig.savefig('./boxplots/'+stationtype+'/boxplot_'+valmode+'_'+domain+'_'+variables[vv]+'_'+aggval+'_'+stationtype+extension, dpi=300)
    plt.close(fig)

#calculate performance gain for last variable
horup_vert10 = (np.mean(np.median(MAE[[4,6],:],axis=1)) - np.mean(np.median(MAE[[0,2],:],axis=1)))/ np.mean(np.median(MAE[[0,2],:],axis=1))*100
horup_vert20 = (np.mean(np.median(MAE[[5,7],:],axis=1)) - np.mean(np.median(MAE[[1,3],:],axis=1)))/ np.mean(np.median(MAE[[1,3],:],axis=1))*100

vertup_horc = (np.mean(np.median(MAE[[1,3],:],axis=1)) - np.mean(np.median(MAE[[0,2],:],axis=1)))/ np.mean(np.median(MAE[[0,2],:],axis=1))*100
vertup_horf = (np.mean(np.median(MAE[[5,7],:],axis=1)) - np.mean(np.median(MAE[[4,6],:],axis=1)))/ np.mean(np.median(MAE[[4,6],:],axis=1))*100

tosaprc_vert10 = (np.mean(np.median(MAE[[2,6],:],axis=1)) - np.mean(np.median(MAE[[0,4],:],axis=1)))/ np.mean(np.median(MAE[[0,4],:],axis=1))*100
tosaprc_vert20 = (np.mean(np.median(MAE[[3,7],:],axis=1)) - np.mean(np.median(MAE[[1,5],:],axis=1)))/ np.mean(np.median(MAE[[1,5],:],axis=1))*100

tosaprc_horc = (np.mean(np.median(MAE[[2,3],:],axis=1)) - np.mean(np.median(MAE[[0,1],:],axis=1)))/ np.mean(np.median(MAE[[0,1],:],axis=1))*100
tosaprc_horf = (np.mean(np.median(MAE[[6,7],:],axis=1)) - np.mean(np.median(MAE[[5,6],:],axis=1)))/ np.mean(np.median(MAE[[5,6],:],axis=1))*100

MAEmed = np.median(MAE,axis=1)
perdev = (MAEmed[1:]-MAEmed[0]) / MAEmed[0]*100
print(perdev)

