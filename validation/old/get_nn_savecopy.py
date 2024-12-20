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
homedir = os.getenv("HOME")
lustre = os.getenv("LUSTRE")

## DEFINE USER OPTIONS
startdates = (20180117,20180118,20180119,20180120,20180121,20180122,20180123,20180124,20180125,20180126,20180127,20180128,20180129,20180130) #which startdates are to be loaded
enddates = (20180118,20180119,20180120,20180121,20180122,20180123,20180124,20180125,20180126,20180127,20180128,20180129,20180130,20180131) #needed to load the data
variables = ['O3'] #which variables are to be loaded
tarlat = [] #target latitudes, will be loaded by metadata.py
tarlon = [] #target longitudes, will be loaded by metadata.py
station_names_meta = [] #station names, will be loaded by metadata.py
domain = 'gal0504r' #target domain
level = 0 #target model level

## load metadata from metadata.py
execfile("metadata.py")

#for each day, load the model data at the gridboxes nearest to the station locations defined in tarlat y tarlon. Finally, concatenate the time series
RES=lustre+'/OP/DATOS/RESULTADOS/chimere2017r4/exp1'
savedir=RES+'/validation'

if not os.path.exists(savedir):
    os.mkdir(savedir)

#DATAEND=np.empty((len(variables),0,len(tarlat)))
#DATESEND=[]
for ii in xrange(len(startdates)):
    #load the data
    srcfile=RES+'/'+str(startdates[ii])+'/out.'+str(startdates[ii])+'03_'+str(enddates[ii])+'03_'+domain+'.nc'
    print('INFO: loading '+srcfile+'...')
    nc=xr.open_dataset(srcfile)
    lat=nc.lat.values
    lon=nc.lon.values
    #define number of timesteps
    nrhours = nc.Time.shape[0]-1
    dates = nc.Times[0:nrhours]
    dates = dates.values
    dates = [string.replace(dates[dd],'-','/') for dd in range(len(dates))]
    dates = [string.replace(dates[dd],'_',' ') for dd in range(len(dates))]
    dates = [parse(dates[dd]) for dd in range(len(dates))]
    
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
        nn_lat[nn] = lat[ind_lat[nn],ind_lon[nn]]
        nn_lon[nn] = lon[ind_lat[nn],ind_lon[nn]]
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
        data[vv,:,:] = data_step.values[0:nrhours,level,ind_lat.astype('int'),ind_lon.astype('int')]
    
    #DATAEND=np.append(DATAEND,data,axis=1)  
    #DATESEND.append(dates)
        
    #create an xarray dataset with dates, add lats and lons, then save to netCDF
    output = xr.DataArray(data, coords=[variables, dates, range(len(nn_lon))], dims=['variable', 'time', 'id'])
    output.id.attrs['lat']=nn_lat
    output.id.attrs['lon']=nn_lon
    output.id.attrs['names']=station_names_meta
        
    #save one nc file per date
    savename=RES+'/validation/data_'+str(startdates[ii])+'.nc'
    print(savename)
    output.to_netcdf(savename)
    nc.close()
    
    


