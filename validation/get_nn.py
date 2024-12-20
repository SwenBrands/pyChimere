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

expname=str(sys.argv[1]) #name of the experiment
domain=str(sys.argv[2]) #name of the domain in lower case letters
reduced=str(sys.argv[3]) #yes or no

homedir = os.getenv("HOME")
lustre = os.getenv("STORE2")

## DEFINE USER OPTIONS
#startdates = (20180401,20180402,20180403,20180404,20180405,20180406,20180407,20180408,20180409,20180410,20180411,20180412,20180413,20180414,20180415,20180416,20180417,20180418,20180419,20180420,20180421,20180422,20180423,20180424,20180425,20180426,20180427,20180428,20180429,20180430) #which startdates are to be loaded
#enddates = (20180402,20180403,20180404,20180405,20180406,20180407,20180408,20180409,20180410,20180411,20180412,20180413,20180414,20180415,20180416,20180417,20180418,20180419,20180420,20180421,20180422,20180423,20180424,20180425,20180426,20180427,20180428,20180429,20180430,20180501) #needed to load the data
#startdates=(20180620,20180621,20180622,20180623,20180624,20180625,20180626,20180627,20180628,20180629,20180630,20180701,20180702,20180703,20180704,20180705,20180706,20180707,20180708,20180709,20180710,20180711,20180712,20180713,20180714,20180715,20180716,20180717,20180718,20180719,20180720,20180721,20180722,20180723,20180724,20180725,20180726,20180727,20180728,20180729,20180730,20180731,20180801,20180802,20180803,20180804,20180805,20180806,20180807,20180808,20180809,20180810,20180811,20180812,20180813,20180814,20180815,20180816,20180817,20180818,20180819,20180820,20180821,20180822,20180823,20180824,20180825,20180826,20180827,20180828,20180829,20180830,20180831)
#enddates=(20180621,20180622,20180623,20180624,20180625,20180626,20180627,20180628,20180629,20180630,20180701,20180702,20180703,20180704,20180705,20180706,20180707,20180708,20180709,20180710,20180711,20180712,20180713,20180714,20180715,20180716,20180717,20180718,20180719,20180720,20180721,20180722,20180723,20180724,20180725,20180726,20180727,20180728,20180729,20180730,20180731,20180801,20180802,20180803,20180804,20180805,20180806,20180807,20180808,20180809,20180810,20180811,20180812,20180813,20180814,20180815,20180816,20180817,20180818,20180819,20180820,20180821,20180822,20180823,20180824,20180825,20180826,20180827,20180828,20180829,20180830,20180831,20180901)
startdates=(20180620,20180621,20180622,20180623,20180624,20180625,20180626,20180627,20180628,20180629,20180630,20180701,20180702,20180703,20180704,20180705,20180706,20180707,20180708,20180709,20180710,20180711,20180712,20180713,20180714,20180715,20180716,20180717,20180718,20180719,20180720,20180721,20180722,20180723,20180724,20180725,20180726,20180727,20180728,20180729,20180730,20180731,20180801,20180802,20180803,20180804,20180805,20180806,20180807,20180808,20180809,20180810,20180811,20180812,20180813,20180814,20180815,20180816,20180817,20180818,20180819,20180820,20180821)
enddates=(20180621,20180622,20180623,20180624,20180625,20180626,20180627,20180628,20180629,20180630,20180701,20180702,20180703,20180704,20180705,20180706,20180707,20180708,20180709,20180710,20180711,20180712,20180713,20180714,20180715,20180716,20180717,20180718,20180719,20180720,20180721,20180722,20180723,20180724,20180725,20180726,20180727,20180728,20180729,20180730,20180731,20180801,20180802,20180803,20180804,20180805,20180806,20180807,20180808,20180809,20180810,20180811,20180812,20180813,20180814,20180815,20180816,20180817,20180818,20180819,20180820,20180821,20180822)


variables = ['O3','PM10','PM25','NO2','SO2'] #which variables are to be loaded
tarlat = [] #target latitudes, will be loaded by metadata.py
tarlon = [] #target longitudes, will be loaded by metadata.py
station_names_meta = [] #station names, will be loaded by metadata.py
level = 0 #target model level

## load metadata from metadata.py
execfile("metadata.py")

#for each day, load the model data at the gridboxes nearest to the station locations defined in tarlat y tarlon. Finally, concatenate the time series
RES=lustre+'/OP/DATOS/RESULTADOS/chimere2017r4/'+expname
savedir=RES+'/validation'

if not os.path.exists(savedir):
    os.mkdir(savedir)

#DATAEND=np.empty((len(variables),0,len(tarlat)))
#DATESEND=[]
for ii in xrange(len(startdates)):
    #load the data
    if reduced == 'yes':
        srcfile=RES+'/'+str(startdates[ii])+'/reduced_out.'+str(startdates[ii])+'03_'+str(enddates[ii])+'03_'+domain+'.nc'
    else:
        srcfile=RES+'/'+str(startdates[ii])+'/out.'+str(startdates[ii])+'03_'+str(enddates[ii])+'03_'+domain+'.nc'
        
    print('INFO: loading '+srcfile+'...')
    nc = xr.open_dataset(srcfile)
    lat = nc.lat.values
    lon = nc.lon.values
    #get dimensions of the input file
    nrdims = len(nc.O3.values.shape)
    
    #define number of timesteps
    nrhours = nc.Times.shape[0]-1
    dates = nc.Times[0:nrhours]
    dates = dates.values
    dates = [string.replace(dates[dd],'-','/') for dd in range(len(dates))]
    dates = [string.replace(dates[dd],'_',' ') for dd in range(len(dates))]
    dates = [parse(dates[dd]) for dd in range(len(dates))]
    #generate pandas datetime index
    dates = pd.date_range(dates[0],dates[-1], freq='H')
    
    #find index values of the nearest neighbours
    ind_lon = np.zeros(len(tarlat))
    ind_lat = np.zeros(len(tarlat))
    nn_lat = np.zeros(len(tarlat))
    nn_lon = np.zeros(len(tarlat))
    for nn in xrange(len(tarlat)):
        latmat = np.matlib.repmat(tarlat[nn],lat.shape[0],lat.shape[1])
        lonmat = np.matlib.repmat(tarlon[nn],lon.shape[0],lon.shape[1])
        dlat = abs(lat-latmat)
        dlon = abs(lon-lonmat) 
        ind_lat_step = zip(*np.where(dlat == np.min(dlat)))
        ind_lat[nn] = ind_lat_step[1][0]
        ind_lon_step = zip(*np.where(dlon == np.min(dlon)))
        ind_lon[nn] = ind_lon_step[0][1]
        nn_lat[nn] = lat[int(ind_lat[nn]),int(ind_lon[nn])]
        nn_lon[nn] = lon[int(ind_lat[nn]),int(ind_lon[nn])]
        print('target latitude is '+str(tarlat[nn]))
        print('nn is '+str(nn_lat[nn]))
        print('target longitude is '+str(tarlon[nn]))
        print('nn is '+str(nn_lon[nn]))        
        
    #load values for each variable
    data = np.zeros((len(variables),nrhours,len(ind_lat)))
    for vv in xrange(len(variables)):
        cstring='nc.'+variables[vv]
        print(cstring)
        data_step = eval(cstring)
        if nrdims == 4:
            data = data_step.values[0:nrhours,level,ind_lat.astype('int'),ind_lon.astype('int')]
        elif nrdims == 3:
            data = data_step.values[0:nrhours,ind_lat.astype('int'),ind_lon.astype('int')]
        else:
            raise Exception('check dimensions of the input nc file!')
    
        #DATAEND=np.append(DATAEND,data,axis=1)  
        #DATESEND.append(dates)
            
        #create an xarray dataset with dates, add lats and lons, then save to netCDF
        output = xr.DataArray(data, coords=[dates, range(len(nn_lon))], dims=['time', 'id'], name = variables[vv])
        output.id.attrs['lat']=nn_lat
        output.id.attrs['lon']=nn_lon
        output.id.attrs['names']=station_names_meta
            
        #save one nc file per date
        savename=RES+'/validation/'+variables[vv]+'_'+str(startdates[ii])+'.nc'
        print(savename)
        output.to_netcdf(savename)
        output.close()

    #close source nc file
    nc.close()
    
    


