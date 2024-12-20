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
aggval=str(sys.argv[2]) #mean, max, min, or hourly
stationtype=str(sys.argv[3]) #A for all, I for industry, B for background, T for traffic, S for sulphur dioxide

homedir = os.getenv("HOME")
lustre = os.getenv("LUSTRE")

## DEFINE USER OPTIONS
startdate = '2018-06-20 04:00:00'
enddate = '2018-08-31 23:00:00'
#variables = ['PM25','O3','NO2','PM10','SO2'] #which variables are to be loaded, there is a problem with outliers in SO2, correct this in future versions
variables = ['NO2','PM10','NO2','PM25']
#xlabel = ['TR', '2008', 'mx3', 'smx3', 'smx3ompi', 'm', 'msd', 'mhr', 'mhr10', 'mhrsd', 'mhr10sd', 'shr']
#xlabel = ['CM10', 'CM20', 'CS10', 'CS20', 'FM10', 'FM20', 'FS10', 'FS20', 'FM10D', 'FM20D']
xlabel = ['CS10', 'CM10', 'CS20', 'CM20', 'FS10', 'FM10', 'FS20', 'FM20']
hores = [xlabel[ii][0] for ii in range(len(xlabel))]
layers = [xlabel[ii][2:4] for ii in range(len(xlabel))] 
nanlimit = 20 #in percent, 60 in former versions
titlesize = 5
labelsize = 4
minval = 0 #minimum value allowed in observations, set to nan if lower
RESpath = lustre+'/OP/DATOS/RESULTADOS/chimere2017r4'
percentiles = [20, 80, 80]
fliers = {'markersize': 10, 'markeredgecolor': 'None', 'marker': '.'}
markersymbol = 'o'
nrstations = 4
groupby = 'hores' #hores, layers
palette = ["m", "g"]
scalefactor = 1
gridstyle = 'whitegrid'
legend = False
kind = 'box'
nr_rows = 2
nr_cols = 2
extension = '.pdf'
fliersize = 2
plottitle = 'yes'
colormap = 'PiYG'
maeskill_min = 0 #not implemented for the moment
maeskill_max = 10 # dito
xlim_bias = [-40,320] #in percent, -80,40 for max O3, PM / -100,120 for max NO2, -110,300 for min NO2, O3, -60,100 for min PM10
xlim_rho = [-0.15,0.7] # -0.1, 1  for max O3, PM / -0.45,0.7 for max NO2, -0.15,0.7 for min NO2, O3, -0.15,0.7 for min PM10
xlim_ratio = [0.6,3] # 0.1 1 for max O3, PM / 0,2.2 for max NO2, 0,2 for min NO2, O3, 0,2 for min PM10
xlim_skill = [-35,25] #-35, 25 for max O3, PM / -70,50 for max NO2, -47,52 for min NO2, O3, -35,25 for min PM10
xlim_skill_ext = [-60,60] #-60, 60 for max O3, PM -60,60 for max NO2, -60,60 for min NO2, O3, -60 60 for min PM10
box_linewidth = 1
showcaps=False
boxprops = dict(linestyle='None')
averval = 'median'

## for each day, load the model data at the gridboxes nearest to the station locations defined in tarlat y tarlon. Finally, concatenate the time series
#src_dirs = [lustre+'/SWEN/DATA/OBSDATA', lustre+'/SWEN/DATA/MODELDATA', RESpath+'/jja_03/validation',
            #RESpath+'/jja_02/validation', RESpath+'/jja_05/validation', RESpath+'/jja_05_ompi/validation',
            #RESpath+'/jja_06/validation', RESpath+'/jja_07/validation', RESpath+'/jja_08/validation', RESpath+'/jja_10/validation',
            #RESpath+'/jja_09/validation', RESpath+'/jja_11/validation', RESpath+'/jja_12/validation']
#src_dirs = [lustre+'/SWEN/DATA/OBSDATA', RESpath+'/jja_10_gal15r/validation', RESpath+'/jja_02_gal15r/validation',RESpath+'/jja_13_gal15r/validation', \
            #RESpath+'/jja_12_gal15r/validation', RESpath+'/jja_10/validation', RESpath+'/jja_08/validation',  RESpath+'/jja_13/validation', \
            #RESpath+'/jja_12/validation', RESpath+'/jja_11/validation', RESpath+'/jja_09/validation']
src_dirs = [lustre+'/SWEN/DATA/OBSDATA', RESpath+'/jja_13_gal15r/validation', RESpath+'/jja_10_gal15r/validation', RESpath+'/jja_12_gal15r/validation', \
            RESpath+'/jja_02_gal15r/validation', RESpath+'/jja_13/validation', RESpath+'/jja_10/validation', RESpath+'/jja_12/validation', RESpath+'/jja_08/validation']


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

emclass = pd.Index(emclass)
if stationtype == 'A':
    print('INFO: All stations are used for verification...')
    #stationfilter = range(len(station_names_csv))
    stationfilterT = list(np.where(emclass.get_loc('T'))[0])
    stationfilterB = list(np.where(emclass.get_loc('B'))[0])
    stationfilterI = list(np.where(emclass.get_loc('I'))[0])
    stationfilterIS = list(np.where(emclass.get_loc('IS'))[0])
    stationfilterTS = list(np.where(emclass.get_loc('TS'))[0])
    stationfilter = stationfilterT+stationfilterB+stationfilterI+stationfilterIS+stationfilterTS
elif stationtype in ('B'):
    stationfilter = list(np.where(emclass.get_loc(stationtype))[0])
elif stationtype in ('T'):
    stationfilterT = list(np.where(emclass.get_loc('T'))[0])
    stationfilterTS = list(np.where(emclass.get_loc('TS'))[0])
    stationfilter = stationfilterT+stationfilterTS
elif stationtype in ('I'):
    stationfilterI = list(np.where(emclass.get_loc('I'))[0])
    stationfilterIS = list(np.where(emclass.get_loc('IS'))[0])
    stationfilter = stationfilterI+stationfilterIS
elif stationtype in ('S'):
    stationfilterIS = list(np.where(emclass.get_loc('IS'))[0])
    stationfilterTS = list(np.where(emclass.get_loc('TS'))[0])
    stationfilter = stationfilterIS+stationfilterTS
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

MAE_SKILL_MN = np.zeros((len(src_dirs)-1,len(variables)))
HIMAE_obs_SKILL_MN = np.zeros((len(src_dirs)-1,len(variables)))
LOMAE_obs_SKILL_MN = np.zeros((len(src_dirs)-1,len(variables)))
HIMAE_mod_SKILL_MN = np.zeros((len(src_dirs)-1,len(variables)))
LOMAE_mod_SKILL_MN = np.zeros((len(src_dirs)-1,len(variables)))
for vv in range(len(variables)):
    obsdata = OBSARR[vv,]
    #find out columns with too many nans or outlier values
    maxumb = np.nanmean(obsdata, axis = 0) + 10*np.nanstd(obsdata, axis = 0)
    maxumb = np.tile(maxumb,[obsdata.shape[0],1])
    negativeind_obs = (obsdata < minval)
    outlierind_obs = (obsdata > maxumb) | (np.isnan(obsdata))
    obsdata[negativeind_obs] = 0
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
    HIMAE_obs = np.zeros((MODARR.shape[1],len(retainind_obs)))
    LOMAE_obs = np.zeros((MODARR.shape[1],len(retainind_obs)))
    HIMAE_mod = np.zeros((MODARR.shape[1],len(retainind_obs)))
    LOMAE_mod = np.zeros((MODARR.shape[1],len(retainind_obs)))
    RHO = np.zeros((MODARR.shape[1],len(retainind_obs)))
    PCTDIFF = np.zeros((len(percentiles),MODARR.shape[1],len(retainind_obs)))
    for ii in range(MODARR.shape[1]):
        moddata = MODARR[vv,ii,]
        moddata[outlierind_obs] = np.nan
        moddata = moddata[:,retainind_obs]           
        BIAS[ii,:] = (np.nanmean(moddata, axis=0) - np.nanmean(obsdata, axis=0)) / np.nanmean(obsdata, axis=0)*100
        RATIO[ii,:] = np.nanstd(moddata, axis=0) / np.nanstd(obsdata, axis=0)
        MAE[ii,:] = np.nanmean(np.absolute(moddata-obsdata),axis=0)
        PCTobs = np.nanpercentile(obsdata,percentiles,axis=0)
        PCTmod = np.nanpercentile(moddata,percentiles,axis=0)
        PCTDIFF[:,ii,:] = PCTmod - PCTobs

        ##conditional biases
        obs_large_obs = np.zeros(obsdata.shape)
        obs_large_obs.fill(np.nan)
        obs_large_mod = np.copy(obs_large_obs)
        obs_small_obs = np.copy(obs_large_obs)
        obs_small_mod = np.copy(obs_large_obs)
        mod_large_obs = np.copy(obs_large_obs)
        mod_large_mod = np.copy(obs_large_obs)
        mod_small_obs = np.copy(obs_large_obs)
        mod_small_mod = np.copy(obs_large_obs)
        
        obsdata_pd = pd.DataFrame(obsdata)
        moddata_pd = pd.DataFrame(moddata)
        for pp in range(obsdata.shape[1]):
            #select values smaller than PCT5
            getind_obs = list(np.where(obsdata_pd[pp] < PCTobs[0][pp]))            
            getind_mod = list(np.where(moddata_pd[pp] < PCTmod[0][pp]))
            obs_small_obs[getind_obs,pp] = obsdata[getind_obs,pp]
            obs_small_mod[getind_mod,pp] = obsdata[getind_mod,pp]
            mod_small_obs[getind_obs,pp] = moddata[getind_obs,pp]
            mod_small_mod[getind_mod,pp] = moddata[getind_mod,pp]
            del(getind_obs,getind_mod)
            #select values larger than PCT95
            getind_obs = list(np.where(obsdata_pd[pp] > PCTobs[2][pp]))            
            getind_mod = list(np.where(moddata_pd[pp] > PCTmod[2][pp]))
            obs_large_obs[getind_obs,pp] = obsdata[getind_obs,pp]
            obs_large_mod[getind_mod,pp] = obsdata[getind_mod,pp]
            mod_large_obs[getind_obs,pp] = moddata[getind_obs,pp]
            mod_large_mod[getind_mod,pp] = moddata[getind_mod,pp]
            del(getind_obs,getind_mod)

        #then calculate conditional biases        
        HIBIAS_obs[ii,:] = (np.nanmean(mod_large_obs, axis=0) - np.nanmean(obs_large_obs, axis=0)) / np.nanmean(obs_large_obs, axis=0)*100
        LOBIAS_obs[ii,:] = (np.nanmean(mod_small_obs, axis=0) - np.nanmean(obs_small_obs, axis=0)) / np.nanmean(obs_small_obs, axis=0)*100
        HIBIAS_mod[ii,:] = (np.nanmean(mod_large_mod, axis=0) - np.nanmean(obs_large_mod, axis=0)) / np.nanmean(obs_large_mod, axis=0)*100
        LOBIAS_mod[ii,:] = (np.nanmean(mod_small_mod, axis=0) - np.nanmean(obs_small_mod, axis=0)) / np.nanmean(obs_small_mod, axis=0)*100
        HIMAE_obs[ii,:] = np.nanmean(np.absolute(mod_large_obs-obs_large_obs),axis=0)
        LOMAE_obs[ii,:] = np.nanmean(np.absolute(mod_small_obs-obs_small_obs),axis=0)
        HIMAE_mod[ii,:] = np.nanmean(np.absolute(mod_large_mod-obs_large_mod),axis=0)
        LOMAE_mod[ii,:] = np.nanmean(np.absolute(mod_small_mod-obs_small_mod),axis=0)
        
        ##vectorize RHO calculation in future versions!
        #corrme = pd.DataFrame(np.concatenate((obsdata,moddata),axis=1))        
        for st in range(obsdata.shape[1]):
            getind = np.where((~np.isnan(obsdata[:,st])) & (~np.isnan(moddata[:,st])))[0] #index to discard outlier data in obs or mod
            RHO[ii,st] = scipy.stats.pearsonr(moddata[getind,st], obsdata[getind,st])[0] #retain only the valid data
    
    #calculate skill scores for MAE and extreme MAES for each station, average them over all stations then refill MAE_SKILL_MN
    MAE_REF = np.tile(MAE[0,],(MAE.shape[0],1))
    MAE_SKILL = (1-(MAE/MAE_REF))*100
    MAE_SKILL_MN[:,vv] = eval('np.nan'+averval+'(MAE_SKILL,axis=1)')
    
    HIMAE_obs_REF=np.tile(HIMAE_obs[0,],(HIMAE_obs.shape[0],1))
    HIMAE_obs_SKILL = (1-(HIMAE_obs/HIMAE_obs_REF))*100
    HIMAE_obs_SKILL_MN[:,vv] = eval('np.nan'+averval+'(HIMAE_obs_SKILL,axis=1)')

    LOMAE_obs_REF=np.tile(LOMAE_obs[0,],(LOMAE_obs.shape[0],1))
    LOMAE_obs_SKILL = (1-(LOMAE_obs/LOMAE_obs_REF))*100
    LOMAE_obs_SKILL_MN[:,vv] = eval('np.nan'+averval+'(LOMAE_obs_SKILL,axis=1)')
    
    HIMAE_mod_REF=np.tile(HIMAE_mod[0,],(HIMAE_mod.shape[0],1))
    HIMAE_mod_SKILL = (1-(HIMAE_mod/HIMAE_mod_REF))*100
    HIMAE_mod_SKILL_MN[:,vv] = eval('np.nan'+averval+'(HIMAE_mod_SKILL,axis=1)')
    
    LOMAE_mod_REF=np.tile(LOMAE_mod[0,],(LOMAE_mod.shape[0],1))
    LOMAE_mod_SKILL = (1-(LOMAE_mod/LOMAE_mod_REF))*100
    LOMAE_mod_SKILL_MN[:,vv] = eval('np.nan'+averval+'(LOMAE_mod_SKILL,axis=1)')
    
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
        sns.catplot(kind=kind, y='experiment', x='bias', data=BIASpd, hue=groupby, orient='h', legend=legend, fliersize=fliersize, ax = axs[0][0], showcaps=showcaps)
        axs[0][0].plot([0,0],[0,MODARR.shape[1]+2],'r')
        axs[0][0].set_xlim(xlim_bias)
        
    retainind=np.where(~np.any(np.isnan(RHO),axis=0))[0]
    RHO = RHO[:,retainind]
    legend_rho = [station_names[ii] for ii in retainind]
    if len(station_names) < nrstations:
        axs[0][1].plot([0.5,MODARR.shape[1]+0.5],[0,0],'r')
        X = np.transpose(np.tile(range(MODARR.shape[1]),(len(station_names),1)))
        axs[0][1].plot(X+1, RHO, linestyle='None', marker=markersymbol)
    else:
        explabels = np.repeat(xlabel,len(retainind))
        horeslabels = np.repeat(hores,len(retainind)) 
        layerslabels = np.repeat(layers,len(retainind)) 
        RHOpd = pd.DataFrame(dict(r=RHO.ravel(),experiment=explabels,hores=horeslabels,layers=layerslabels))
        sns.catplot(kind=kind, y='experiment', x='r', data=RHOpd, hue=groupby, orient='h', legend=legend, fliersize=fliersize, ax = axs[0][1], showcaps=showcaps)
        axs[0][1].set_xlim(xlim_rho)
        axs[0][1].plot([0,0],[0,MODARR.shape[1]+2],'r')        

    retainind=np.where(~np.any(np.isnan(RATIO),axis=0))[0]
    RATIO = RATIO[:,retainind]
    legend_ratio = [station_names[ii] for ii in retainind]
    if len(station_names) < nrstations:
        axs[1][0].plot([0.5,MODARR.shape[1]+0.5],[0,0],'r')
        X = np.transpose(np.tile(range(MODARR.shape[1]),(len(station_names),1)))
        axs[1][0].plot(X+1, RATIO, linestyle='None', marker=markersymbol)        
    else:
        explabels = np.repeat(xlabel,len(retainind))
        horeslabels = np.repeat(hores,len(retainind)) 
        layerslabels = np.repeat(layers,len(retainind)) 
        RATIOpd = pd.DataFrame(dict(ratio=RATIO.ravel(),experiment=explabels,hores=horeslabels,layers=layerslabels))
        sns.catplot(kind=kind, y='experiment', x='ratio', data=RATIOpd, hue=groupby, orient='h', legend=legend, fliersize=fliersize, ax = axs[1][0], showcaps=showcaps)
        axs[1][0].set_xlim(xlim_ratio)
        axs[1][0].plot([1,1],[0,MODARR.shape[1]+2],'r')
            
    retainind=np.where(~np.any(np.isnan(MAE_SKILL),axis=0))[0]
    MAE_SKILL = MAE_SKILL[:,retainind]
    legend_mae = [station_names[ii] for ii in retainind]
    if len(station_names) < nrstations:
        axs[1][1].plot([0.5,MODARR.shape[1]+0.5],[0,0],'r')
        X = np.transpose(np.tile(range(MODARR.shape[1]),(len(station_names),1)))
        axs[1][1].plot(X+1, MAE_SKILL, linestyle='None', marker=markersymbol)        
    else:
        explabels = np.repeat(xlabel,len(retainind))
        horeslabels = np.repeat(hores,len(retainind)) 
        layerslabels = np.repeat(layers,len(retainind)) 
        MAE_SKILLpd = pd.DataFrame(dict(mae_skill=MAE_SKILL.ravel(),experiment=explabels,hores=horeslabels,layers=layerslabels))
        sns.catplot(kind=kind, y='experiment', x='mae_skill', data=MAE_SKILLpd, hue=groupby, orient='h', legend=legend, fliersize=fliersize, ax = axs[1][1], showcaps=showcaps)
        axs[1][1].set_xlim(xlim_skill)
        axs[1][1].plot([0,0],[0,MODARR.shape[1]+2],'r')
    
    fig.savefig('./boxplots/'+stationtype+'/boxplot_'+valmode+'_'+variables[vv]+'_'+aggval+'_'+stationtype+extension, dpi=300)
    plt.close(fig)

##generate a pcolor plot and save
MAE_SKILL_MN = np.flipud(MAE_SKILL_MN)
xlabel = np.flipud(xlabel) #invert experiment labels
xlabel = xlabel.astype(str).tolist() #then turn to list containing strings

#get correct colorbar
cblim_mae = 25 #25 for nanmedian
cblim_himae_obs = 25 #25 for nanmediann
cblim_himae_mod = 25 #25 for nanmedian
#cblim_mae = np.max(np.abs(MAE_SKILL_MN))
#cblim_himae_obs = np.max(np.abs(HIMAE_obs_SKILL_MN))
#cblim_himae_mod = np.max(np.abs(HIMAE_mod_SKILL_MN))

if aggval in ('hourly','max','mean'):
    printscores = ['MAE','HIMAE_obs','LOMAE_obs','HIMAE_mod','LOMAE_mod']
elif aggval == 'min':
    printscores = ['MAE','HIMAE_obs','HIMAE_mod']
else:
    raise Exception('check entry for <aggval>!')
for iii in ['MAE','HIMAE_obs','LOMAE_obs','HIMAE_mod','LOMAE_mod']:
    fig = plt.figure()
    plt.pcolor(eval(iii+'_SKILL_MN'), cmap=colormap, vmin=cblim_himae_mod*-1, vmax=cblim_himae_mod)
    plt.xticks(np.array(range(len(variables)))+0.5,variables, fontsize = labelsize+12)
    plt.yticks(np.array(range(len(xlabel)))+0.5,xlabel, fontsize = labelsize+12)
    cb = plt.colorbar()
    cb.set_label(label=iii+' skill score (%)',weight='bold',fontsize=labelsize+12)
    cb.ax.tick_params(labelsize=labelsize+8)
    fig.savefig('./boxplots/'+stationtype+'/pcolor_'+iii+'_'+valmode+'_'+aggval+'_'+stationtype+extension, dpi=300)
    plt.close(fig)
