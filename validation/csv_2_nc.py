#!/usr/bin/env python

#this scripts loads the csv files for the observations and old model data, formats it to be comparable to the new model data and then saves to ncCHIMERE output files and searches out the data at the grid-boxes nearest to those defined in <tarlat> and <tarlon>, the output is then save in nc format.

## WICHTIG: IN DIESER VERSION von csv_2_nc.py WERDEN DIE LETZTEN BEIDEN STATIONEN / SPALTEN GELOESCHT, SIEHE LINIEN 105-106. Pruefe ob metadata.py tarlat, tarlon und station_names_csv entsprechend definiert!

import xarray as xr
import numpy as np
import pandas as pd
import sys
import os
import numpy.matlib
import csv
from datetime import datetime
from dateutil.parser import parse

#get environment
homedir = os.getenv("HOME")
lustre = os.getenv("LUSTRE")

#defines source and output directories
src_obs=lustre+'/SWEN/DATA/OBSDATA'
src_mod=lustre+'/SWEN/DATA/MODELDATA'

#define input parameters
variables = ['O3','PM10','PM25','NO2','SO2'] #which variables are to be loaded
tarlat = [] #target latitudes, will be loaded by metadata.py, ASSURE THAT THE LATS ARE IN THE SAME ORDER THEN THE DATA ROWS IN THE CSV FILES
tarlon = [] #target longitudes, will be loaded by metadata.py, DITO!
station_names_csv = [] #station names, will be loaded by metadata.py
savelabel = '20180930' #label used to save the files
changehour = '2018-03-25 02:00:00' #'2018-03-24 02:00:00', time slice that is missing in Anthony's csv files, will be fileld with nans

## EXECUTE #############################################################
## load metadata from metadata.py
execfile("metadata.py")

for vv in xrange(len(variables)):
    obsfile=src_obs+'/MED_'+variables[vv]+'.csv'
    modfile=src_mod+'/MOD_'+variables[vv]+'.csv'
    
    #observation file
    ifile = open(obsfile, "rb")
    reader = csv.reader(ifile, delimiter=';')    
    #create a list containing the csv's rows 
    obs=[]
    for row in reader:
        #print row
        obs.append(row)
    ifile.close()
    del(row,reader,ifile)
    
    #model file
    ifile = open(modfile, "rb")
    reader = csv.reader(ifile, delimiter=';')
    #create a list containing the csv's rows 
    mod=[]
    for row in reader:
        #print row
        mod.append(row)            
    ifile.close()
    del(row,reader,ifile)

    #get the dates
    dates = [None]*(len(obs)-1)    
    for rr in range(len(obs)-1): #first row is not read
        obsstep = obs[rr+1]
        dates[rr] = obsstep[0]
    
    #dates = [parse(dates[dd]) for dd in range(len(dates))]
    dates_csv = [datetime.strptime(dates[ii], '%d/%m/%Y %H:%M:%S') for ii in range(len(dates))]

    #generate pandas datetime index
    dates = pd.date_range(dates_csv[0],dates_csv[-1], freq='H')
    
    #then get the data and fill in with NaNs at the time jump winter -> summer
    obsdata = np.zeros((len(dates),len(obs[0])-2)) #first row is not read, last one is empty
    moddata = np.zeros((len(dates),len(obs[0])-2))
    
    hourind = dates.get_loc(changehour)
    for rr in range(len(dates)): #first row is not read
        if rr == hourind:
            print('change in hour is on: '+str(dates[hourind]))
            obsstep = obs[rr+1]
            modstep = mod[rr+1]
            obsdata[rr,:] = np.repeat(np.nan,len(obsstep)-2)
            moddata[rr,:] = np.repeat(np.nan,len(modstep)-2)
        elif rr < hourind:
            obsstep = obs[rr+1]
            modstep = mod[rr+1]
            obsdata[rr,:] = [np.float(obsstep[cc+1]) for cc in range(len(obsstep)-2)]
            moddata[rr,:] = [np.float(modstep[cc+1]) for cc in range(len(modstep)-2)]
        else:
            obsstep = obs[rr]
            modstep = mod[rr]
            obsdata[rr,:] = [np.float(obsstep[cc+1]) for cc in range(len(obsstep)-2)]
            moddata[rr,:] = [np.float(modstep[cc+1]) for cc in range(len(modstep)-2)]

    ##filter out obsdata, moddata, GETIND MUST COINCIDE WITH the index of <station_names> in the MED and MOD csv files !!
    #obsdata = obsdata[:,tarind]
    #moddata = moddata[:,tarind]
    
    #set nans
    obsdata[np.where(obsdata==-999)] = np.nan
    moddata[np.where(moddata==-999)] = np.nan
    
    ##through out last two stations to compensate Anthony's error, make sure that this is consistent with output from metadata.py
    #obsdata = obsdata[:,0:-2]
    #moddata = moddata[:,0:-2]
    
    ##expand dims and append
    #obsdata = np.expand_dims(obsdata,0)
    #moddata = np.expand_dims(moddata,0)
    #OBSDATA_end = np.append(OBSDATA_end,obsdata,axis=0)
    #MODDATA_end = np.append(MODDATA_end,moddata,axis=0)

    #create an xarray dataframe
    obs_output = xr.DataArray(obsdata, coords=[dates, range(len(tarlon))], dims=['time', 'id'], name=variables[vv])
    obs_output.id.attrs['lat']=tarlat
    obs_output.id.attrs['lon']=tarlon
    obs_output.id.attrs['names']=station_names_csv

    mod_output = xr.DataArray(moddata, coords=[dates, range(len(tarlon))], dims=['time', 'id'], name = variables[vv])
    mod_output.id.attrs['lat']=tarlat
    mod_output.id.attrs['lon']=tarlon
    mod_output.id.attrs['names']=station_names_csv

    #and finally save to nc format
    savename_obs=src_obs+'/'+variables[vv]+'_'+savelabel+'.nc'
    savename_mod=src_mod+'/'+variables[vv]+'_'+savelabel+'.nc'
    print(savename_obs)
    obs_output.to_netcdf(savename_obs)
    obs_output.close()
    print(savename_mod)
    mod_output.to_netcdf(savename_mod)    
    mod_output.close()
    
