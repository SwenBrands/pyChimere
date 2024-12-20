# -*- coding: utf-8 -*-

#load python packages
import xarray as xr
import numpy as np
import os
from matplotlib import pyplot as plt

species = 'NO2'
domain = 'GAL3'
emis = 'mix2'
months = ['01','02','03','04','05','06','07','08','09','10','11','12']
#months = ['01']

#define the paths
homedir = os.getenv("HOME")
lustre = os.getenv("LUSTRE")
srcpath = lustre+'/OP/PRED/chimere2017r4'

#load land sea mask of the correpsonding domain
maskfile = homedir+'/OP/LANZAR/chimere2017r4/domains/'+domain+'/LANDUSE_USGS_'+domain+'.nc'
nc = xr.open_dataset(maskfile)
seafrac = nc.Ocean.values
nc.close()

#replace the data for each month
for month in months:
    #load auxiliary file
    auxfile = lustre+'/OP/DATOS/CICC/chimere2017r4/EMISSIONS/htap/melchior/usgs/'+domain+'/NOSD/EMI_2010/EMIS.'+domain+'.'+month+'.'+species+'.s.nc'
    nc = xr.open_dataset(auxfile)
    dataaux = nc.variables[species].values
    nc.close()
    
    #load file to be modified    
    #srcfile = lustre+'/OP/DATOS/CICC/chimere2017r4/EMISSIONS/'+emis+'/EMIS.'+domain+'.'+month+'.'+species+'.s.nc'
    srcfile = lustre+'/OP/DATOS/CICC/chimere2017r4/EMISSIONS/'+emis+'/melchior/usgs/'+domain+'/NOSD/EMI_2010/EMIS.'+domain+'.'+month+'.'+species+'.s.nc'
    nc = xr.open_dataset(srcfile)
    data = nc.variables[species].values    
    
    #replicate seafrac and search sea points
    seafrac_step = np.tile(seafrac,(data.shape[0],data.shape[1],data.shape[2],1,1))
    getind = np.where(seafrac_step == 1)
    
    #replace data
    data[getind] = dataaux[getind]
    nc.variables[species].values = data
    
    #then save the modified data
    outfile = lustre+'/OP/DATOS/CICC/chimere2017r4/EMISSIONS/'+emis+'/melchior/usgs/'+domain+'/NOSD/EMI_2010/new_EMIS.'+domain+'.'+month+'.'+species+'.s.nc'
    nc.to_netcdf(outfile)
    nc.close()
