# -*- coding: utf-8 -*-

#load python packages
import numpy as np
import xarray as xr
from mpl_toolkits.basemap import Basemap
import pandas as pd
import os
import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import matplotlib.colors as colors

domain = 'GAL3'
species = ['APINEN','BaP_fin','BbF_fin','BCAR_coa','BCAR_fin','BkF_fin','C2H4','C2H5OH','C2H6', \
            'C3H6','C5H8','CH3CHO','CH3COE','CH3OH','CH4','CO','H2SO4_fin','HCHO','HONO', 'NC4H10', \
            'NH3', 'NO2', 'NO', 'OCAR_coa', 'OCAR_fin', 'OXYL', 'PPM_big', 'PPM_coa', 'PPM_fin', \
            'SO2', 'TMB', 'TOL', 'TPPM_coa', 'TPPM_fin']
            
tarmonths = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']

lustre = os.getenv("LUSTRE")
srcpath = lustre+'/OP/DATOS/CICC/chimere2017r4/EMISSIONS/htap/melchior/usgs/'+domain+'/NOSD_corrected/EMI_2010/test'
savepath = lustre+'/OP/DATOS/CICC/chimere2017r4/EMISSIONS/htap/melchior/usgs/'+domain+'/NOSD_corrected/EMI_2010'

## EXECUTE #############################################################
for sp in range(len(species)):
    for mm in range(len(tarmonths)):
        filename = srcpath+'/EMIS.'+domain+'.'+str(tarmonths[mm])+'.'+str(species[sp])+'.s.nc'
        #outfilename = savepath+'/new_EMIS.'+domain+'.'+str(tarmonths[mm])+'.'+str(species[sp])+'.s.nc'
        outfilename = savepath+'/EMIS.'+domain+'.'+str(tarmonths[mm])+'.'+str(species[sp])+'.s.nc' 
        nc = xr.open_dataset(filename)
        lats = nc.lat
        lons = nc.lon
        latsmask = lats[1:,0:-1]
        lonsmask = lons[1:,0:-1]       
        wholefile = nc.load()
        data = wholefile.variables[species[sp]].values
        mask = data[:,:,:,1:,0:-1]
        data[:,:,:,0:-1,1:] = mask
        wholefile.variables[species[sp]].values = data
        #outfile = xr.DataArray(data)
        wholefile.to_netcdf(outfilename)
        wholefile.close()
        nc.close()

