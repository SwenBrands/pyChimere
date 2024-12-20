# -*- coding: utf-8 -*-

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
from mpl_toolkits.basemap import Basemap
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import scipy.stats
import seaborn as sns
from utiles import crea_cmap

aggval=str(sys.argv[1]) #mean, max, min, or hourly
stationtype=str(sys.argv[2]) #A for all, I for industry, B for background, T for traffic, S for sulphur dioxide

homedir = os.getenv("HOME")
lustre = os.getenv("STORE2")

## DEFINE USER OPTIONS
startdate = '2018-06-21 04:00:00'
enddate = '2018-08-21 23:00:00'
#variables = ['O3','NO2','PM25','PM10','SO2'] #which variables are to be loaded, there is a problem with outliers in SO2, correct this in future versions
variables = ['PM25']
#xlabel = ['FM20Hclim', 'H5', 'CS10H', 'CM10H', 'CS20H', 'CM20H', 'FS10H', 'FM10H', 'FM10HSD', 'FS20H', 'FM20H', 'FM20HSD', 'FM20HCM', 'FM20HNBE', 'FM20HGAL3', 'FM20HGAL3c', 'OP', \
#            'CS10E', 'CM10E', 'CS20E', 'CM20E', 'FS10E', 'FM10E', 'FS20E', 'FM20E'] #CM = "coarse meteorology", NBE = no biogenic emissions, SD = statisical downscaling, c=corrected
xlabel = ['CS10E', 'CS20E', 'FS10E', 'FS20E', 'FM20H'] #CM = "coarse meteorology", NBE = no biogenic emissions, SD = statisical downscaling, c=corrected

hores = [xlabel[ii][0] for ii in range(len(xlabel))]
layers = [xlabel[ii][2:4] for ii in range(len(xlabel))] 
nanlimit = 20 #in percent, 60 in former versions
titlesize = 5
labelsize = 4
minval_obs = 0 #minimum value allowed in observations, set to nan if lower
RESpath = lustre+'/OP/DATOS/RESULTADOS/chimere2017r4' #path to the model and obs data
savepath = lustre+'/SWEN/publications/2018_Brands_chimere/figs/overlay_maps' #path to save the figures
percentile = 90
markersymbol = 'o'
colormap = 'hot_r' #'terrain', 'Greys'
markersize = 40
transparency = 1
extension = '.pdf'
dpival = 200
plottitle = 'yes'
maeskill_min = 0 #not implemented for the moment
maeskill_max = 10 # dito
xlim_bias = [-40,320] #in percent, -80,40 for max O3, PM / -100,120 for max NO2, -110,300 for min NO2, O3, -60,100 for min PM10
xlim_rho = [-0.15,0.7] # -0.1, 1  for max O3, PM / -0.45,0.7 for max NO2, -0.15,0.7 for min NO2, O3, -0.15,0.7 for min PM10
xlim_ratio = [0.6,3] # 0.1 1 for max O3, PM / 0,2.2 for max NO2, 0,2 for min NO2, O3, 0,2 for min PM10
xlim_skill = [-35,25] #-35, 25 for max O3, PM / -70,50 for max NO2, -47,52 for min NO2, O3, -35,25 for min PM10
xlim_skill_ext = [-60,60] #-60, 60 for max O3, PM -60,60 for max NO2, -60,60 for min NO2, O3, -60 60 for min PM10
color_points = 'grey' #red
bias_threshold = 4.0 #threshold for spatial mean bias to decide on color of the coastline

## for each day, load the model data at the gridboxes nearest to the station locations defined in tarlat y tarlon. Finally, concatenate the time series
src_dir_obs = lustre+'/SWEN/DATA/OBSDATA'
#experiments = ['jja_16', 'hindcast_05', 'jja_13_gal15r', 'jja_10_gal15r', 'jja_12_gal15r', 'jja_02_gal15r', 'jja_13', 'jja_10', 'jja_11', 'jja_12', \
#            'jja_08', 'jja_09', 'jja_14', 'jja_15', 'jja_06', 'jja_06_corr_gal3', 'jja_05', 'jja_24', 'jja_23', 'jja_22', 'jja_21', 'jja_20', 'jja_19', 'jja_18', 'jja_17']

experiments = ['jja_24', 'jja_22', 'jja_20', 'jja_18', 'jja_12']

#domains = ['gal0504r', 'gal3','gal15r','gal15r','gal15r','gal15r','gal0504r','gal0504r','gal0504r','gal0504r','gal0504r','gal0504r','gal0504r','gal0504r', \
#            'gal3','gal3','gal3', 'gal15r', 'gal15r', 'gal15r', 'gal15r', 'gal0504r', 'gal0504r', 'gal0504r', 'gal0504r']

domains = ['gal15r','gal15r','gal0504r','gal0504r','gal0504r']


## EXECUTE #############################################################
# Creamos el custom colormap:
refdates = pd.date_range(startdate,enddate,freq='H')

## load metadata from metadata.py
tarlat = [] #target latitudes, will be loaded by metadata.py
tarlon = [] #target longitudes, will be loaded by metadata.py
station_names_csv = [] #station names, will be loaded by metadata.py
execfile('metadata_rev1.py')
execfile(homedir+'/OP/LANZAR/chimere2017r4/figures/get_rgbs_new.py')

emclass = pd.Index(emclass)
if stationtype == 'A':
    print('INFO: All stations are used for verification...')
    #stationfilter = range(len(station_names_csv))
    stationfilterT = list(np.where(emclass.get_loc('T'))[0])
    stationfilterB = list(np.where(emclass.get_loc('B'))[0])
    stationfilterBC = list(np.where(emclass.get_loc('BC'))[0])
    stationfilterI = list(np.where(emclass.get_loc('I'))[0])
    #stationfilterIC = list(np.where(emclass.get_loc('IC'))[0])
    stationfilterIS = list(np.where(emclass.get_loc('IS'))[0])
    stationfilterISC = list(np.where(emclass.get_loc('ISC'))[0])
    stationfilterTS = list(np.where(emclass.get_loc('TS'))[0])
    stationfilter = stationfilterT+stationfilterB+stationfilterBC+stationfilterI+stationfilterIS+stationfilterISC+stationfilterTS
elif stationtype in ('B'):
    stationfilterB = list(np.where(emclass.get_loc('B'))[0])
    stationfilterBC = list(np.where(emclass.get_loc('BC'))[0])
    stationfilter = stationfilterB+stationfilterBC
elif stationtype in ('T'):
    stationfilterT = list(np.where(emclass.get_loc('T'))[0])
    stationfilterTS = list(np.where(emclass.get_loc('TS'))[0])
    stationfilter = stationfilterT+stationfilterTS
elif stationtype in ('I'):
    stationfilterI = list(np.where(emclass.get_loc('I'))[0])
    stationfilterIC = list(np.where(emclass.get_loc('IC'))[0])
    stationfilterIS = list(np.where(emclass.get_loc('IS'))[0])
    stationfilterISC = list(np.where(emclass.get_loc('ISC'))[0])
    stationfilter = stationfilterI+stationfilterIC+stationfilterIS+stationfilterISC
elif stationtype in ('S'):
    stationfilterIS = list(np.where(emclass.get_loc('IS'))[0])
    stationfilterTS = list(np.where(emclass.get_loc('TS'))[0])
    stationfilter = stationfilterIS+stationfilterTS
elif stationtype in ('C'):
    stationfilterBC = list(np.where(emclass.get_loc('BC'))[0])
    stationfilterTC = list(np.where(emclass.get_loc('TC'))[0])
    stationfilterISC = list(np.where(emclass.get_loc('ISC'))[0])
    stationfilter = stationfilterBC+stationfilterTC+stationfilterISC
    #stationfilter = stationfilterBC+stationfilterTC
elif stationtype in ('NOTR'): #excludes traffic stations
    stationfilterB = list(np.where(emclass.get_loc('B'))[0])
    stationfilterI = list(np.where(emclass.get_loc('I'))[0])
    stationfilterIS = list(np.where(emclass.get_loc('IS'))[0])
    stationfilterBC = list(np.where(emclass.get_loc('BC'))[0])
    stationfilterISC = list(np.where(emclass.get_loc('ISC'))[0])
    stationfilter = stationfilterB+stationfilterI+stationfilterIS+stationfilterBC
else:
    raise Exception('check entry for <stationtype>!')

#load example file to check the size of the time dim
srcfile=src_dir_obs+'/'+variables[0]+'_*.nc'
nc=xr.open_mfdataset(srcfile)
cstring='nc.'+variables[0]
dummy = eval(cstring)
dummy = dummy.sel(time=slice(startdate,enddate))
#init output array
nc.close()

#init for errors in the mean
pattern_rho_mean = np.zeros((len(variables),len(experiments)))
pattern_pval_mean = np.copy(pattern_rho_mean)
pattern_std_ratio_mean = np.copy(pattern_rho_mean)
pattern_bias_mean = np.copy(pattern_rho_mean)
pattern_mae_mean = np.copy(pattern_rho_mean)
#init for errors in sigma
pattern_rho_std = np.copy(pattern_rho_mean)
pattern_pval_std = np.copy(pattern_rho_mean)
pattern_std_ratio_std = np.copy(pattern_rho_mean)
pattern_bias_std = np.copy(pattern_rho_mean)
pattern_mae_std = np.copy(pattern_rho_mean)
for vv in range(len(variables)):
    
    #define coastline color
    if variables[vv] == 'O3':
        color_coast_mean = 'cyan'
    elif variables[vv] in ('PM10','PM25') and aggval == 'min':
        color_coast_mean = 'cyan'
    else:
        color_coast_mean = 'blue'
    
    color_coast_std = 'blue'
    
    #load the observations for this variable
    obsfile=src_dir_obs+'/'+variables[vv]+'_20180930.nc'
    print('INFO: loading '+obsfile+'...')
    nc=xr.open_dataset(obsfile)
    cstring='nc.'+variables[vv]
    data = eval(cstring)
    obsarray = data.sel(time=slice(startdate,enddate)).values    
    lat = nc.id.lat
    lon = nc.id.lon
    names = nc.id.names
    dates_obs = nc.time
    dates_obs = dates_obs.sel(time=slice(startdate,enddate))
    #obsarray = data.values
    #close nc file
    nc.close()
    #filter stations
    obsarray = obsarray[:,stationfilter]
    lat = lat[stationfilter]
    lon = lon[stationfilter]
    #filter out the target stations
    station_names = [station_names_csv[ff] for ff in stationfilter]
    
    ##find out columns with too many nans or outlier values
    maxumb = np.nanmean(obsarray, axis = 0) + 10*np.nanstd(obsarray, axis = 0)
    maxumb = np.tile(maxumb,[obsarray.shape[0],1])
    negativeind_obs = (obsarray < minval_obs)
    outlierind_obs = (obsarray > maxumb) | (np.isnan(obsarray))
    obsarray[negativeind_obs] = 0
    obsarray[outlierind_obs] = np.nan
    #delete stations with too many outliers
    nanperc_obs = np.sum(outlierind_obs, axis=0)/obsarray.shape[0]*100
    retainind_obs = np.where(nanperc_obs < nanlimit)[0]
    obsarray = obsarray[:,retainind_obs]
    lat = lat[retainind_obs]
    lon = lon[retainind_obs]
    station_names = [station_names[ii] for ii in retainind_obs]
        
    #get list of chimere folders, one folder per date
    chim_folders = [str(dates_obs.values[ii])[0:10].replace('-','') for ii in range(len(dates_obs))]
    chim_folders = np.unique(np.array(chim_folders)).tolist()
    #load the full grid of the model experiments for each day, then load and stack the data

    for mm in range(len(experiments)):
        src_dir_mod = RESpath+'/'+experiments[mm]
        #rename for use within in the loop
        obsarray_step = obsarray
        dates_obs_step = dates_obs
        
        #check if save directory already exists and create if it does not
        if os.path.isdir(savepath+'/'+stationtype):
            print(savepath+'/'+stationtype+' already exists..')
        else:
            os.mkdir(savepath+'/'+stationtype)
        
        #check if save directory already exists and create if it does not
        if os.path.isdir(savepath+'/'+stationtype+'/'+xlabel[mm]):
            print(savepath+'/'+stationtype+'/'+xlabel[mm]+' already exists..')
        else:
            os.mkdir(savepath+'/'+stationtype+'/'+xlabel[mm])
            
        print('loading model data for experiment '+src_dir_mod)
        for dd in range(len(chim_folders)):
            print(chim_folders[dd])
            modfile = src_dir_mod+'/'+chim_folders[dd]+'/*out.'+chim_folders[dd]+'03_*'+domains[mm]+'.nc'
            print('processing '+modfile)
            if dd == 0:
                print('searching for nearest neighbours in '+modfile)
                #load model lats and lons for use within nearest neighbour search
                nc = xr.open_mfdataset(modfile)
                lat_mod_step = nc.lat.values[:,0]
                lon_mod_step = nc.lon.values[0,:]
                nc.close()
                #find index values of the nearest neighbours                
                ind_lat = np.zeros(len(lat)).astype('int')
                ind_lon = np.zeros(len(lat)).astype('int')
                nn_lat = np.zeros(len(lat))
                nn_lon = np.zeros(len(lat))
                for nn in xrange(len(lat)):
                    #latmat = np.matlib.repmat(lat[nn],lat_mod_step.shape[0],lat_mod_step.shape[1])
                    #lonmat = np.matlib.repmat(lon[nn],lon_mod_step.shape[0],lon_mod_step.shape[1])
                    latvec = np.matlib.repeat(lat[nn],len(lat_mod_step))
                    lonvec = np.matlib.repeat(lon[nn],len(lon_mod_step))
                    dlat = abs(lat_mod_step-latvec)
                    dlon = abs(lon_mod_step-lonvec) 
                    ind_lat[nn] = int(np.argmin(dlat))
                    ind_lon[nn] = int(np.argmin(dlon))
                    nn_lat[nn] = lat_mod_step[int(ind_lat[nn])]
                    nn_lon[nn] = lon_mod_step[int(ind_lon[nn])]
                    print('target latitude is '+str(lat[nn]))
                    print('nn is '+str(nn_lat[nn]))
                    print('target longitude is '+str(lon[nn]))
                    print('nn is '+str(nn_lon[nn]))            
            
            nc = xr.open_mfdataset(modfile)
            #get the dimensions of the data array, 4 dims in case of original chimere file, 3 dims in case of the postprocessed file called "reduced...nc"
            nrdims = len(nc.variables[variables[vv]].shape)
            south_north = nc.south_north.values
            west_east = nc.west_east.values
            lat_mod = nc.lat.values
            lon_mod = nc.lon.values
            dates_mod_step = nc.Times.values[0:-1] #don't load the last date (3 UTC) to avoid repetitions
            if nrdims == 4:
                data_step = nc.variables[variables[vv]].values[0:-1,0,:,:]  #don't load the last date (3 UTC) to avoid repetitions, filter first level
            elif nrdims == 3:
                data_step = nc.variables[variables[vv]].values[0:-1,:,:]  #don't load the last date (3 UTC) to avoid repetitions, in this case the first level has been already filtered in previous postprocessing step
            else:
                raise Exception('unknown number of dimensions for input data array!')
                
            if dd == 0:
                dates_mod = dates_mod_step
                modarray = data_step
            else:
                dates_mod = np.append(dates_mod,dates_mod_step)
                modarray = np.append(modarray,data_step,axis=0)
            del data_step
            nc.close()
        
        #convert format of the model dates and select dates defined in refdates
        dates_mod = list(dates_mod)
        dates_mod = [dates_mod[ii].replace('_',' ') for ii in range(len(dates_mod))]
        dates_mod = pd.date_range(dates_mod[0],dates_mod[-1],freq='H') #bring dates_mod to the same format as refdates
        #brind to xarray dataarray format
        arr_mod = xr.DataArray(modarray, coords=[dates_mod,south_north,west_east], dims=['time','south_north','west_east'], name=variables[vv])
        #filter target dates as defined in refdays
        arr_mod = arr_mod.sel(time=slice(refdates[0], refdates[-1]))
        #obtain pandas datetime format
        dates_mod = pd.date_range(arr_mod.time.values[0],arr_mod.time.values[-1],freq='H')
        
        #repeat for observations
        dates_obs_step = dates_obs_step.values
        #dates_obs_step = [dates_obs_step[ii].astype('string').replace('T',' ') for ii in range(len(dates_obs_step))]
        dates_obs_step = pd.date_range(dates_obs_step[0],dates_obs_step[-1],freq='H') #also bring dates_obs_step to pandas datetime format
        arr_obs = xr.DataArray(obsarray_step, coords=[dates_obs_step, station_names], dims=['time', 'loc'], name=variables[vv])
        
        #from hereon, dates, refdates and dates_mod must cover the same period and have the same format
        
        #optionally calculate aggregated values for the observations
        if aggval in ('min','max','mean'):
            print('INFO: daily '+aggval+' are calculated..')            
            #get daily dates 
            dates_day = pd.date_range(dates_obs_step[0],dates_obs_step[-1],freq='D').date.astype('str')
            #then find max values for each day
            obsarray_step = np.zeros((len(dates_day),obsarray_step.shape[1])) #here, obsarray_step is overwritten!
            modarray = np.zeros((len(dates_day),arr_mod.shape[1],arr_mod.shape[2])) #here, modarray is overwritten
            for tt in range(len(dates_day)):
                if aggval == 'max':
                    obsarray_step[tt,:] = np.amax(arr_obs.sel(time=dates_day[tt]),axis=0)
                    modarray[tt,:,:] = np.amax(arr_mod.sel(time=dates_day[tt]),axis=0)
                elif aggval == 'min':
                    obsarray_step[tt,:] = np.amin(arr_obs.sel(time=dates_day[tt]),axis=0)
                    modarray[tt,:,:] = np.amin(arr_mod.sel(time=dates_day[tt]),axis=0)
                elif aggval == 'mean':
                    obsarray_step[tt,:] = np.mean(arr_obs.sel(time=dates_day[tt]),axis=0)
                    modarray[tt,:,:] = np.mean(arr_mod.sel(time=dates_day[tt]),axis=0)
                else:
                    raise Exception('check entry for <aggval>!')
        elif aggval == 'hourly':
            print('INFO: hourly data are calculated...')
        else:
            raise Exception('check entry for <aggval>!')
        
        ##initialize results arrays for the model and observations
        mean_mod = np.nanmean(modarray,axis=0)
        std_mod = np.nanstd(modarray,axis=0)
        pct_mod = np.nanpercentile(modarray,percentile,axis=0)
        # MEAN_mod_hiobs = np.zeros((MODARR.shape[1],len(retainind_obs)))
        # MEAN_mod_loobs = np.zeros((MODARR.shape[1],len(retainind_obs)))
        # MEAN_mod_himod = np.zeros((MODARR.shape[1],len(retainind_obs)))
        # MEAN_mod_ lomod = np.zeros((MODARR.shape[1],len(retainind_obs)))
    
        mean_obs = np.nanmean(obsarray_step,axis=0)
        std_obs = np.nanstd(obsarray_step,axis=0)
        pct_obs = np.nanpercentile(obsarray_step,percentile,axis=0)
        # MEAN_obs_hiobs = np.copy(MEAN_mod_hiobs)
        # MEAN_obs_loobs = np.copy(MEAN_mod_loobs)
        # MEAN_obs_himod = np.copy(MEAN_mod_himod)
        # MEAN_obs_lomod = np.copy(MEAN_mod_lomod)
        
        #calculate spatial error for the model mean (bias, rho sigmaratio, mae)
        bias_mean = mean_mod[ind_lat,ind_lon] - mean_obs
        pattern_bias_mean[vv,mm] = np.mean(bias_mean)
        pattern_rho_step_mean = scipy.stats.spearmanr(mean_mod[ind_lat,ind_lon],mean_obs)
        pattern_rho_mean[vv,mm] = pattern_rho_step_mean[0]
        pattern_pval_mean[vv,mm] = pattern_rho_step_mean[1]
        pattern_std_ratio_mean[vv,mm] = np.std(mean_mod[ind_lat,ind_lon]) / np.std(mean_obs)
        pattern_mae_mean[vv,mm] = np.mean(abs(bias_mean))
        
        #calculate spatial error for the model standard deviation
        bias_std = std_mod[ind_lat,ind_lon] - std_obs
        pattern_bias_std[vv,mm] = np.mean(bias_std)
        pattern_rho_step_std = scipy.stats.spearmanr(std_mod[ind_lat,ind_lon],std_obs)
        pattern_rho_std[vv,mm] = pattern_rho_step_std[0]
        pattern_pval_std[vv,mm] = pattern_rho_step_std[1]
        pattern_std_ratio_std[vv,mm] = np.std(std_mod[ind_lat,ind_lon]) / np.std(std_obs)
        pattern_mae_std[vv,mm] = np.mean(abs(bias_std))
        
        del(obsarray_step)
        # #define color for the coastlines
        # if pattern_bias[vv,mm] > bias_threshold:
            # color_coast = color_coast_pos
        # else:
            # color_coast = color_coast_neg
        
        #plot model results overlain with station results        
        minlat = 41.7
        minlon = -9.5
        maxlat = 43.9
        maxlon = -6.7
        lat_0 = 43 #42.4
        lon_0 = -9
        #get parallels and meridians
        parallels = np.arange(np.ceil(minlat),np.floor(maxlat+1))
        meridians = np.arange(np.ceil(minlon),np.floor(maxlon+1))
        #minval = np.nanmin(np.concatenate((mean_mod.flatten(),mean_obs)))
        #maxval = np.nanmax(np.concatenate((mean_mod.flatten(),mean_obs)))
        if variables[vv] == 'PM25' and aggval == 'max':
            minval_mean = np.nanmin(mean_obs)
            maxval_mean = np.sort(mean_obs)[-2]
            minval_std = np.nanmin(std_obs)
            maxval_std = np.nanmax(std_obs)
        elif variables[vv] == 'PM10' and aggval == 'max':
            minval_mean = np.nanmin(mean_obs)
            maxval_mean = np.sort(mean_obs)[-4]
            minval_std = np.nanmin(std_obs)
            maxval_std = np.nanmax(std_obs)
        else:
            minval_mean = np.nanmin(mean_obs)
            maxval_mean = np.nanmax(mean_obs)
            minval_std = np.nanmin(std_obs)
            maxval_std = np.nanmax(std_obs)
        
        ## optionally build colormap
        #cbounds = np.linspace(minval,maxval,len(hex_no2))
        #colormap = crea_cmap(cbounds, hex_no2, hex_no2[0], hex_no2[-1])
        
        #build figure for the mean values
        fig = plt.figure()
        mymap = Basemap(projection='cea',llcrnrlon=minlon,llcrnrlat=minlat,urcrnrlon=maxlon,urcrnrlat=maxlat,lat_ts=lat_0)#lat_ts=43
        #mymap.drawcoastlines(color='black')
        #mymap.drawcountries(color='black')
        #mymap.readshapefile(homedir+'/OP/LANZAR/chimere2017r4/figures/shapes/municipios', 'municipios',linewidth=0.35, color='k', antialiased=1)          
        mymap.readshapefile(homedir+'/OP/LANZAR/chimere2017r4/figures/shapes/espana', 'espana',linewidth=0.35, color=color_coast_mean, antialiased=1)
        mymap.readshapefile(homedir+'/OP/LANZAR/chimere2017r4/figures/shapes/portugal', 'portugal',linewidth=0.35, color=color_coast_mean, antialiased=1)
        X, Y = mymap(lon_mod,lat_mod)
        lon1, lat2 = mymap(lon,lat)
        #mymap.pcolormesh(lon_mod, lat_mod, mean_mod, cmap=colormap, vmin=minval, vmax=maxval)
        #mymap.scatter(lon, lat, s=None, c=mean_obs, marker='o', edgecolors='grey', cmap=colormap, norm=None, vmin=minval, vmax=maxval)
        mymap.pcolormesh(X, Y, mean_mod, cmap=colormap, vmin=minval_mean, vmax=maxval_mean)
        mymap.scatter(lon1, lat2, s=markersize, c=mean_obs, marker='o', edgecolors=color_points, cmap=colormap, norm=None, vmin=minval_mean, vmax=maxval_mean, alpha=transparency, zorder=5)

        #ax.set_title(u'Concentracións de '+species+' no '+dates[tt])
        mymap.drawparallels(parallels,labels=[True,False,False,True], color='None', size=10)
        mymap.drawmeridians(meridians,labels=[True,False,False,True], color='None', size=10)
        plt.title('Bias: '+str(round(pattern_bias_mean[vv,mm],1))+' R: '+str(round(pattern_rho_mean[vv,mm],2))+' Ratio: '+str(round(pattern_std_ratio_mean[vv,mm],2))+' MAE: '+str(round(pattern_mae_mean[vv,mm],1)))
        plt.colorbar()
        savefile = savepath+'/'+stationtype+'/'+xlabel[mm]+'/'+variables[vv]+'_mean_'+xlabel[mm]+'_'+domains[mm]+'_'+aggval+'_'+stationtype+extension
        fig.savefig(savefile, dpi=dpival)
        plt.close('all')
        
        #build figure for the std values
        fig = plt.figure()
        mymap = Basemap(projection='cea',llcrnrlon=minlon,llcrnrlat=minlat,urcrnrlon=maxlon,urcrnrlat=maxlat,lat_ts=lat_0)#lat_ts=43
        mymap.readshapefile(homedir+'/OP/LANZAR/chimere2017r4/figures/shapes/espana', 'espana',linewidth=0.35, color=color_coast_std, antialiased=1)
        mymap.readshapefile(homedir+'/OP/LANZAR/chimere2017r4/figures/shapes/portugal', 'portugal',linewidth=0.35, color=color_coast_std, antialiased=1)
        X, Y = mymap(lon_mod,lat_mod)
        lon1, lat2 = mymap(lon,lat)
        mymap.pcolormesh(X, Y, std_mod, cmap=colormap, vmin=minval_std, vmax=maxval_std)
        mymap.scatter(lon1, lat2, s=markersize, c=std_obs, marker='o', edgecolors=color_points, cmap=colormap, norm=None, vmin=minval_std, vmax=maxval_std, alpha=transparency, zorder=5)

        #ax.set_title(u'Concentracións de '+species+' no '+dates[tt])
        mymap.drawparallels(parallels,labels=[True,False,False,True], color='None', size=10)
        mymap.drawmeridians(meridians,labels=[True,False,False,True], color='None', size=10)
        plt.title('SBIAS: '+str(round(pattern_bias_std[vv,mm],1))+' SR: '+str(round(pattern_rho_std[vv,mm],2))+' SRATIO: '+str(round(pattern_std_ratio_std[vv,mm],2))+' SMAE: '+str(round(pattern_mae_std[vv,mm],1)))
        plt.colorbar()
        savefile = savepath+'/'+stationtype+'/'+xlabel[mm]+'/'+variables[vv]+'_std_'+xlabel[mm]+'_'+domains[mm]+'_'+aggval+'_'+stationtype+extension
        fig.savefig(savefile, dpi=dpival)
        plt.close('all')



# #start the verfication
# OBSARR = ARRAY[:,0] #ARRAY containing observations (3D)
# MODARR = ARRAY[:,1:] #ARRAY containing model output (4D since various model versions are considered

# ##check validation mode
# if valmode == 'spatial':
    # print('INFO: spatial validation mode, arrays are transposed!')
    # MODARR = np.swapaxes(MODARR,3,2)
    # OBSARR = np.swapaxes(OBSARR,1,2)
# elif valmode == 'temporal':
    # print('INFO: temporal validation mode, arrays are not transposed!')
# else:
    # raise Exception('check entry for <valmode>!')

# for vv in range(len(variables)):
    # obsdata = OBSARR[vv,]
    # #find out columns with too many nans or outlier values
    # maxumb = np.nanmean(obsdata, axis = 0) + 10*np.nanstd(obsdata, axis = 0)
    # maxumb = np.tile(maxumb,[obsdata.shape[0],1])
    # negativeind_obs = (obsdata < minval)
    # outlierind_obs = (obsdata > maxumb) | (np.isnan(obsdata))
    # obsdata[negativeind_obs] = 0
    # obsdata[outlierind_obs] = np.nan
    # #delete stations with too many outliers
    # nanperc_obs = np.sum(outlierind_obs, axis=0)/obsdata.shape[0]*100
    # retainind_obs = np.where(nanperc_obs < nanlimit)[0]
    # obsdata = obsdata[:,retainind_obs]
    # station_names = [station_names_csv[ii] for ii in retainind_obs]
    

    
    # for ii in range(MODARR.shape[1]):
        # moddata = MODARR[vv,ii,]
        # moddata[outlierind_obs] = np.nan
        # moddata = moddata[:,retainind_obs]           
        
        # #get results for the model
        # MEAN_mod[ii,:] = np.nanmean(moddata, axis=0)
        # STD_mod[ii,:] = np.nanstd(moddata, axis=0)
        # PCT_mod = np.nanpercentile(moddata,percentiles,axis=0)
        
        # #get results for the observations
        # MEAN_obs[ii,:] = np.nanmean(obsdata, axis=0)
        # STD_obs[ii,:] = np.nanstd(obsdata, axis=0)
        # PCT_obs = np.nanpercentile(obsdata,percentiles,axis=0)
        
        # ##conditional biases
        # obs_large_obs = np.zeros(obsdata.shape)
        # obs_large_obs.fill(np.nan)
        # obs_large_mod = np.copy(obs_large_obs)
        # obs_small_obs = np.copy(obs_large_obs)
        # obs_small_mod = np.copy(obs_large_obs)
        # mod_large_obs = np.copy(obs_large_obs)
        # mod_large_mod = np.copy(obs_large_obs)
        # mod_small_obs = np.copy(obs_large_obs)
        # mod_small_mod = np.copy(obs_large_obs)
        
        # obsdata_pd = pd.DataFrame(obsdata)
        # moddata_pd = pd.DataFrame(moddata)
        # for pp in range(obsdata.shape[1]):
            # #select values smaller than PCT5
            # getind_obs = list(np.where(obsdata_pd[pp] < PCTobs[0][pp]))            
            # getind_mod = list(np.where(moddata_pd[pp] < PCTmod[0][pp]))
            # obs_small_obs[getind_obs,pp] = obsdata[getind_obs,pp]
            # obs_small_mod[getind_mod,pp] = obsdata[getind_mod,pp]
            # mod_small_obs[getind_obs,pp] = moddata[getind_obs,pp]
            # mod_small_mod[getind_mod,pp] = moddata[getind_mod,pp]
            # del(getind_obs,getind_mod)
            # #select values larger than PCT95
            # getind_obs = list(np.where(obsdata_pd[pp] > PCTobs[2][pp]))            
            # getind_mod = list(np.where(moddata_pd[pp] > PCTmod[2][pp]))
            # obs_large_obs[getind_obs,pp] = obsdata[getind_obs,pp]
            # obs_large_mod[getind_mod,pp] = obsdata[getind_mod,pp]
            # mod_large_obs[getind_obs,pp] = moddata[getind_obs,pp]
            # mod_large_mod[getind_mod,pp] = moddata[getind_mod,pp]
            # del(getind_obs,getind_mod)

        # #then calculate conditional mean values      
        # MEAN_mod_hiobs[ii,:] = np.nanmean(mod_large_obs, axis=0)
        # MEAN_mod_loobs[ii,:] = np.nanmean(mod_small_obs, axis=0)
        # MEAN_obs_hiobs[ii,:] = np.nanmean(obs_large_obs, axis=0)
        # MEAN_obs_loobs[ii,:] = np.nanmean(obs_small_obs, axis=0)
        
        # MEAN_mod_himod[ii,:] = np.nanmean(mod_large_mod, axis=0)
        # MEAN_mod_lomod[ii,:] = np.nanmean(mod_small_mod, axis=0)
        # MEAN_obs_himod[ii,:] = np.nanmean(obs_large_mod, axis=0)
        # MEAN_obs_lomod[ii,:] = np.nanmean(obs_small_mod, axis=0)

        # #bis hierhin
    
    # #plot the results
    # fig, axs = plt.subplots(nr_rows, nr_cols, squeeze=False)
    # if plottitle == 'yes':
        # fig.suptitle('Verfication for '+variables[vv]+' '+aggval+' '+stationtype, fontsize=8, fontweight='bold')
    
    # retainind=np.where(~np.any(np.isnan(BIAS),axis=0))[0]
    # BIAS = BIAS[:,retainind]
    # legend_bias = [station_names[ii] for ii in retainind]
    # if len(station_names) < nrstations:
        # axs[0][0].plot([0.5,MODARR.shape[1]+0.5],[0,0],'r')
        # X = np.transpose(np.tile(range(MODARR.shape[1]),(len(station_names),1)))
        # axs[0][0].plot(X+1, BIAS, linestyle='None', marker=markersymbol)
        # axs[0][0].legend(legend_bias)
    # else:        
        # explabels = np.repeat(xlabel,len(retainind))
        # horeslabels = np.repeat(hores,len(retainind)) 
        # layerslabels = np.repeat(layers,len(retainind)) 
        # BIASpd = pd.DataFrame(dict(bias=BIAS.ravel(),experiment=explabels,hores=horeslabels,layers=layerslabels))
        # sns.catplot(kind=kind, y='experiment', x='bias', data=BIASpd, hue=groupby, orient='h', legend=legend, fliersize=fliersize, ax = axs[0][0], showcaps=showcaps)
        # #axs[0][0].set_xlim(xlim_bias)
        # axs[0][0].plot([0,0],[0,MODARR.shape[1]+2],'r')
        
    # retainind=np.where(~np.any(np.isnan(RHO),axis=0))[0]
    # RHO = RHO[:,retainind]
    # legend_rho = [station_names[ii] for ii in retainind]
    # if len(station_names) < nrstations:
        # axs[0][1].plot([0.5,MODARR.shape[1]+0.5],[0,0],'r')
        # X = np.transpose(np.tile(range(MODARR.shape[1]),(len(station_names),1)))
        # axs[0][1].plot(X+1, RHO, linestyle='None', marker=markersymbol)
    # else:
        # explabels = np.repeat(xlabel,len(retainind))
        # horeslabels = np.repeat(hores,len(retainind)) 
        # layerslabels = np.repeat(layers,len(retainind)) 
        # RHOpd = pd.DataFrame(dict(r=RHO.ravel(),experiment=explabels,hores=horeslabels,layers=layerslabels))
        # sns.catplot(kind=kind, y='experiment', x='r', data=RHOpd, hue=groupby, orient='h', legend=legend, fliersize=fliersize, ax = axs[0][1], showcaps=showcaps)
        # #axs[0][1].set_xlim(xlim_rho)
        # axs[0][1].plot([0,0],[0,MODARR.shape[1]+2],'r')        

    # retainind=np.where(~np.any(np.isnan(RATIO),axis=0))[0]
    # RATIO = RATIO[:,retainind]
    # legend_ratio = [station_names[ii] for ii in retainind]
    # if len(station_names) < nrstations:
        # axs[1][0].plot([0.5,MODARR.shape[1]+0.5],[0,0],'r')
        # X = np.transpose(np.tile(range(MODARR.shape[1]),(len(station_names),1)))
        # axs[1][0].plot(X+1, RATIO, linestyle='None', marker=markersymbol)        
    # else:
        # explabels = np.repeat(xlabel,len(retainind))
        # horeslabels = np.repeat(hores,len(retainind)) 
        # layerslabels = np.repeat(layers,len(retainind)) 
        # RATIOpd = pd.DataFrame(dict(ratio=RATIO.ravel(),experiment=explabels,hores=horeslabels,layers=layerslabels))
        # sns.catplot(kind=kind, y='experiment', x='ratio', data=RATIOpd, hue=groupby, orient='h', legend=legend, fliersize=fliersize, ax = axs[1][0], showcaps=showcaps)
        # #axs[1][0].set_xlim(xlim_ratio)
        # axs[1][0].plot([1,1],[0,MODARR.shape[1]+2],'r')
            
    # retainind=np.where(~np.any(np.isnan(MAE_SKILL),axis=0))[0]
    # MAE_SKILL = MAE_SKILL[:,retainind]
    # legend_mae = [station_names[ii] for ii in retainind]
    # if len(station_names) < nrstations:
        # axs[1][1].plot([0.5,MODARR.shape[1]+0.5],[0,0],'r')
        # X = np.transpose(np.tile(range(MODARR.shape[1]),(len(station_names),1)))
        # axs[1][1].plot(X+1, MAE_SKILL, linestyle='None', marker=markersymbol)        
    # else:
        # explabels = np.repeat(xlabel,len(retainind))
        # horeslabels = np.repeat(hores,len(retainind)) 
        # layerslabels = np.repeat(layers,len(retainind)) 
        # MAE_SKILLpd = pd.DataFrame(dict(mae_skill=MAE_SKILL.ravel(),experiment=explabels,hores=horeslabels,layers=layerslabels))
        # sns.catplot(kind=kind, y='experiment', x='mae_skill', data=MAE_SKILLpd, hue=groupby, orient='h', legend=legend, fliersize=fliersize, ax = axs[1][1], showcaps=showcaps)
        # #axs[1][1].set_xlim(xlim_skill)
        # axs[1][1].plot([0,0],[0,MODARR.shape[1]+2],'r')
    
    # fig.savefig('./boxplots/'+stationtype+'/boxplot_'+valmode+'_'+variables[vv]+'_'+aggval+'_'+stationtype+extension, dpi=300)
    # plt.close(fig)

# ##generate a pcolor plot and save
# MAE_SKILL_MN = np.flipud(MAE_SKILL_MN)
# xlabel = np.flipud(xlabel) #invert experiment labels
# xlabel = xlabel.astype(str).tolist() #then turn to list containing strings

# #get correct colorbar
# cblim_mae = 25 #25 for nanmedian
# cblim_himae_obs = 25 #25 for nanmediann
# cblim_himae_mod = 25 #25 for nanmedian
# #cblim_mae = np.max(np.abs(MAE_SKILL_MN))
# #cblim_himae_obs = np.max(np.abs(HIMAE_obs_SKILL_MN))
# #cblim_himae_mod = np.max(np.abs(HIMAE_mod_SKILL_MN))

# if aggval in ('hourly','max','mean'):
    # printscores = ['MAE','HIMAE_obs','LOMAE_obs','HIMAE_mod','LOMAE_mod']
# elif aggval == 'min':
    # printscores = ['MAE','HIMAE_obs','HIMAE_mod']
# else:
    # raise Exception('check entry for <aggval>!')
# for iii in ['MAE','HIMAE_obs','LOMAE_obs','HIMAE_mod','LOMAE_mod']:
    # fig = plt.figure()
    # plt.pcolor(eval(iii+'_SKILL_MN'), cmap=colormap, vmin=cblim_himae_mod*-1, vmax=cblim_himae_mod)
    # plt.xticks(np.array(range(len(variables)))+0.5,variables, fontsize = labelsize+12)
    # plt.yticks(np.array(range(len(xlabel)))+0.5,xlabel, fontsize = labelsize+12)
    # cb = plt.colorbar()
    # cb.set_label(label=iii+' skill score (%)',weight='bold',fontsize=labelsize+12)
    # cb.ax.tick_params(labelsize=labelsize+8)
    # fig.savefig('./boxplots/'+stationtype+'/pcolor_'+iii+'_'+valmode+'_'+aggval+'_'+stationtype+extension, dpi=300)
    # plt.close(fig)
