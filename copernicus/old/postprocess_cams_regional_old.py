#!/usr/bin/env python

#from netCDF4 import Dataset
import xarray as xr
import numpy as np
import pandas as pd
import sys
import dask
import os
import numpy.matlib
#from matplotlib import pyplot as plt

#call with postprocess_cams_regional.py 20180120 00

##set input variables
tardate=str(sys.argv[1])
tarhour=str(sys.argv[2])
leadtime=int(sys.argv[3]) #leadtime in hours

#tardate='20180520'
#tarhour='00'
#leadtime=72

taryear=tardate[0:4]
tarmonth=tardate[4:6]

#define home
homedir = os.getenv("HOME")
lustre = os.getenv("LUSTRE")

#dimensions of the CAMS input files are (latitude: 451, level: 60, longitude: 900, time: 11)
source=lustre+'/OP/DATOS/CICC/CAMS/'+tardate
source_clim=lustre+'/OP/DATOS/CICC/chimere2017r4/MACC/myclim_regional'

variables_in = ('c2h6','hcho','ch4','co','hno3','c5h8','no2','go3','pan','so2','aermr01','aermr02','aermr03','aermr04','aermr05','aermr06','aermr08','aermr10','aermr11','lnsp','aermr11','aermr11','aermr11','aermr11','aermr11','aermr11','aermr11','aermr11','aermr11') #set variable name as defined in CAMS
variables_out = ('C2H6','CH2O','CH4','CO','HNO3','ISOP','NO2','O3','PAN','SO2','SS1','SS2','SS3','DUST1','DUST2','DUST3','OM','BC','SO4','PSFC','BIGALK','BIGENE','CH3CHO','GLYOXAL','H2O2','NH3','TOLUENE','C10H16','TEMP') #set variable name as required by CHIMERE
units = ('ppb','ppb','ppb','ppb','ppb','ppb','ppb','ppb','ppb','ppb','mmr','mmr','mmr','ug/m3','ug/m3','ug/m3','ug/m3','ug/m3','ug/m3','Pa','ppb','ppb','ppb','ppb','ppb','ppb','ppb','ppb','K')
aero_type = ('gas','gas','gas','gas','gas','gas','gas','gas','gas','gas','aer','aer','aer','dust','dust','dust','aer','aer','aer','none','gas','gas','gas','gas','gas','gas','gas','gas','none') #set aerosol type as required by CHIMERE
aero_low = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.03, 0.5, 5.0, 0.03, 0.55, 0.9, 0.1, 0.1, 0.1, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0] #refers to the aerosol size in um, required by chimere for sea-salt and dust aerosol
aero_high = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.5, 5.0 ,20.0 ,0.55, 0.9, 20.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0] # "
mmass = [30.0, 30.0, 16.0, 28.0, 63.0, 68.0, 46.0, 48.0, 121.0, 64.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 96.0, 1.0, 72.0, 56.0, 44.0, 28.0, 38.0, 17.0, 92.0, 136.0, 1.0] #set the molar mass, required by CHIMERE
density = np.array((0.000150,0.000402,0.000672,0.001025,0.001485,0.002087,0.002814,0.003688,0.004740,0.006018,0.007599,0.009582,0.012082,0.015130,0.018863,0.023518,0.029322,0.036557,0.045579,0.056826,0.070849,0.087826,0.107127,0.129274,0.154328,0.182375,0.213467,0.247616,0.284796,0.324943,0.367189,0.403822,0.441886,0.481280,0.521847,0.563389,0.605674,0.648445,0.691426,0.734331,0.776863,0.818729,0.859634,0.899295,0.937436,0.973801,1.008150,1.040270,1.069971,1.097099,1.121533,1.143192,1.162039,1.178087,1.191403,1.202112,1.210406,1.216546,1.220871,1.223803))

#variables_in = ('aermr01','aermr02')
#variables_out = ('SS1','SS2')
#aero_type = ('aer','aer')
#aero_low = [0.03, 0.55]
#aero_high = [0.55, 0.9]
#mmass = [100.0, 100.0]

#for unit conversion recall:
# 1 ppb = 1 ug/kg
# 1 ppm = 1 mg/kg

### EXECUTE #############################################################
#get alpha and beta coefficients for vertical grid
alpha = np.array((20.000000,38.425343,63.647804,95.636963,134.483307,180.584351,234.779053,298.495789,373.971924,464.618134,575.651001,713.218079,883.660522,1094.834717,1356.474609,1680.640259,2082.273926,2579.888672,3196.421631,3960.291504,4906.708496,6018.019531,7306.631348,8765.053711,10376.126953,12077.446289,13775.325195,15379.805664,16819.474609,18045.183594,19027.695313,19755.109375,20222.205078,20429.863281,20384.480469,20097.402344,19584.330078,18864.750000,17961.357422,16899.468750,15706.447266,14411.124023,13043.218750,11632.758789,10209.500977,8802.356445,7438.803223,6144.314941,4941.778320,3850.913330,2887.696533,2063.779785,1385.912598,855.361755,467.333588,210.393890,65.889244,7.367743,0.000000,0.000000))
beta = np.array((0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000076,0.000461,0.001815,0.005081,0.011143,0.020678,0.034121,0.051690,0.073534,0.099675,0.130023,0.164384,0.202476,0.243933,0.288323,0.335155,0.383892,0.433963,0.484772,0.535710,0.586168,0.635547,0.683269,0.728786,0.771597,0.811253,0.847375,0.879657,0.907884,0.931940,0.951822,0.967645,0.979663,0.988270,0.994019,0.997630,1.000000))
hoursperday=[0, 1, 2, 3, 4, 5, 6, 7]
hoursperday=numpy.matlib.repmat(hoursperday,1,leadtime/24+1)
hoursperday=hoursperday[0]
hoursperday=hoursperday[1:-6]
outputdir=source+'/aux'

for ii in xrange(len(variables_in)):
    #load file in cache
    srcvar = variables_in[ii]
    tarvar = variables_out[ii]
    srcfile = source+'/z_cams_c_ecmf_'+tardate+tarhour+'0000_prod_fc_ml_*_'+srcvar+'.nc'
    outputfile = outputdir+'/'+variables_out[ii]+'.nc' 
    print('loading files '+srcfile)        
    nc = xr.open_mfdataset(srcfile, concat_dim='time', engine='netcdf4')
    
    #load and process lat and lon
    lats = nc.variables['latitude'][:]
    lons = nc.variables['longitude'][:]
    lats = lats.values
    lons = lons.values
    #get indices for target region
    latsind=np.squeeze(np.array(np.where((lats > 30) & (lats < 60))))
    lonsind_west = np.array(np.where(lons > 340))
    lonsind_east = np.array(np.where(lons < 10))
    lonsind = np.squeeze(np.concatenate((lonsind_east,lonsind_west),axis=1)) 
    lats = lats[latsind]
    lons = lons[lonsind]
        
    #load and process the time variable
    dates = nc.variables['time']    
    dates = str(dates.values)
    dates = dates.replace("-", "")
    dates = dates.replace("\n", "")
    dates = dates.replace("\n", "")
    dates = dates.replace("T", ".")
    dates = dates.replace(":00:00.000000000","")
    dates = dates.replace(" ",",")
    dates = dates.replace("[","")
    dates = dates.replace("]","")
    dates = dates.replace("'","")
    dates = dates.split(',')    
    
    for i in xrange(len(dates)):
        dates[i] = dates[i][0:8]
        
    dates=[float(i) for i in dates]    
    
    for i in xrange(len(hoursperday)):             
        dates[i] = dates[i]+0.125*hoursperday[i]
    
    ##overwrite lats, lons in case the climatological mean values are available only and use dates from dummy variables defined in variables_in
    #if tarvar in ('BIGALK','BIGENE','CH3CHO','GLYOXAL','H2O2','NH3','TOLUENE','C10H16','TEMP','SS1','SS2','SS3'):
    if tarvar in ('BIGALK','BIGENE','CH3CHO','GLYOXAL','H2O2','NH3','TOLUENE','C10H16','TEMP'): #for these variables, the climatological mean values are loaded    
        #add path to the climatological files
        srcfile_clim = source_clim+'/clim_'+tarvar+'_'+tarmonth+'.nc'        
        print('loading climatological values from '+srcfile_clim)        
        nc_clim = xr.open_dataset(srcfile_clim, engine='netcdf4')
        
        #load and process lat and lon
        lats_clim = nc_clim.variables['lat'][:]
        lons_clim = nc_clim.variables['lon'][:]
        lats_clim = lats_clim.values
        lons_clim = lons_clim.values
        #get indices for target region
        latsind_clim=np.squeeze(np.array(np.where((lats_clim > 30) & (lats_clim < 60))))
        lonsind_clim=np.squeeze(np.array(np.where((lons_clim > 340) | (lons_clim < 10))))
                
        ##recall to add the following two lines in case xr.open_dataset is used and variables are reloaded (which is not the case in the current version of this script)
        #nc_clim.close()
        #nc_clim = xr.open_dataset(srcfile_clim, engine='netcdf4')
        
        #overwrite lats and lons
        lats = lats_clim[latsind_clim]
        lons = lons_clim[lonsind_clim]
        
        data = nc_clim.variables[tarvar][:,latsind_clim,lonsind_clim]
        
        if tarvar in ('BIGALK','BIGENE','CH3CHO','GLYOXAL','H2O2','NH3','TOLUENE','C10H16'):
           #convert kg/kg to ppb (both mass mixing ratios, recall 1 ppb = 1 ug/kg)
           data = data*10**9
           #then convert ppb (mass missing ratio) to ppb (mixing ratio by volume)
           data = data*(29./mmass[ii]) #29 is the molar mass of air in g/mole, see http://paos.colorado.edu/~toohey/5710faq.html
        elif tarvar in ('SS1','SS2','SS3'):
           ## 1 ppb = 0.001 mg/kg
           #print('converting kg/kg to ug/m3 for '+tarvar)
           #data = data*10**9 #convert kg/kg to ppb
           ##then weigh by air density
           #dens3d = np.tile(density,(data.shape[1],data.shape[2],1))
           #dens3d = np.moveaxis(dens3d,-1,0)
           #data = data*dens3d
           print('no conversion needed for '+tarvar)
        elif tarvar in ('TEMP'):
           print('no need to convert to ppb for '+tarvar)
        
        levs = nc.variables['level'][:] #load the 60 model levels from the dummy variable
                
        ##append climatological mean value 25 times (3 hour interval, leadtime = 72 hours) to mimic a time series, required by option 2 for the lateral boundary conditions in the CHIMERE namelist (chimere.par) 
        data = data.values
        data = np.tile(data,(len(dates),1,1,1))
        output = xr.DataArray(data, coords=[dates, levs, np.double(lats), np.double(lons)], dims=['time', 'lev', 'lat', 'lon'], name=variables_out[ii])
        nc_clim.close()        
        
        #define variable attributes
        if tarvar in ('BIGALK','BIGENE','CH3CHO','GLYOXAL','H2O2','NH3','TOLUENE','C10H16'):
           output.attrs['refpres'] = 1
           output.attrs['typ'] = str(aero_type[ii])
           output.attrs['molarmass'] = float(mmass[ii])
        elif tarvar in ('SS1','SS2','SS3'):
           output.attrs['refpres'] = 1
           output.attrs['typ'] = str(aero_type[ii])
           output.attrs['molarmass'] = float(mmass[ii])
           #set aerosol diamter
           output.attrs['bin_low'] = float(aero_low[ii])
           output.attrs['bin_high'] = float(aero_high[ii])
        
    ##load and process the data    
    #elif srcvar in ('aermr01','aermr02','aermr03','aermr04','aermr05','aermr06','aermr08','aermr10','aermr11'):
    elif tarvar in ('DUST1','DUST2','DUST3','OM','BC','SO4','SS1','SS2','SS3'):
        #this a 4D variable
        data = nc.variables[srcvar][:,:,latsind,lonsind]
        levs = nc.variables['level'][:] #load 3d levels
        if tarvar in ('DUST1','DUST2','DUST3','OM','BC','SO4'):
           ## 1 ppb = 0.001 mg/kg
           print('converting kg/kg to ug/m3 for '+tarvar)
           data = data*10**9 #convert kg/kg to ppb
           ##then weigh by air density
           dens4d = np.tile(density,(data.shape[0],data.shape[2],data.shape[3],1))        
           dens4d = np.moveaxis(dens4d,-1,1)
           data = data*dens4d
        elif tarvar in ('SS1','SS2','SS3'):
           print('no conversion needed for '+tarvar)
                           
        output = xr.DataArray(data, coords=[dates, levs, np.double(lats), np.double(lons)], dims=['time', 'lev', 'lat', 'lon'], name=variables_out[ii])
        output.attrs['refpres'] = 1
        output.attrs['typ'] = str(aero_type[ii])
        output.attrs['molarmass'] = float(mmass[ii])
        #set aerosol diamter
        output.attrs['bin_low'] = float(aero_low[ii])
        output.attrs['bin_high'] = float(aero_high[ii])
    
    #elif srcvar in ('c2h6','hcho','ch4','co','hno3','c5h8','hno3','no2','go3','pan','so2'):
    elif tarvar in ('C2H6','CH2O','CH4','CO','HNO3','ISOP','NO2','O3','PAN','SO2'):
        print('converting kg/kg to ppb for '+tarvar)
        data = nc.variables[srcvar][:,:,latsind,lonsind]
        
        #convert kg/kg to ppb (both mass mixing ratios, recall 1 ppb = 1 ug/kg)
        data = data*10**9
        #then convert ppb (mass missing ratio) to ppb (mixing ratio by volume)
        data = data*(29./mmass[ii]) #29 is the molar mass of air in g/mole, see http://paos.colorado.edu/~toohey/5710faq.html
       
        levs = nc.variables['level'][:] #load 3d levels        
        output = xr.DataArray(data, coords=[dates, levs, np.double(lats), np.double(lons)], dims=['time', 'lev', 'lat', 'lon'], name=variables_out[ii])
        output.attrs['refpres'] = 1
        output.attrs['typ'] = str(aero_type[ii])
        output.attrs['molarmass'] = float(mmass[ii])
                
    #elif srcvar == 'lnsp':
    elif tarvar == 'PSFC':
        #this is a 3D variable
        print('INFO: '+tarvar+' normal, i.e. non logarithmic, pressure has been already obtained by get_cams.sh!')        
        data = nc.variables[srcvar][:,latsind,lonsind]
        ##convert logarithmic surface pressure to Pa
        #print('converting log(Pa) to Pa for '+tarvar)
        #data = np.exp(data)
        output = xr.DataArray(data, coords=[dates, np.double(lats), np.double(lons)], dims=['time', 'lat', 'lon'], name=variables_out[ii])        
    else:        
        raise Exception('CAMS input variable not found!')  
    
    ##then assign units and save to netCDF
    print('INFO: the unit for '+variables_out[ii]+' is '+units[ii])
    print('INFO: the species class for '+variables_out[ii]+' is '+aero_type[ii])
    print('INFO: the molar mass for '+variables_out[ii]+' is '+str(mmass[ii]))
    print('INFO: the lower aerosol size for '+variables_out[ii]+' is '+str(aero_low[ii]))
    print('INFO: the upper aerosol size for '+variables_out[ii]+' is '+str(aero_high[ii]))    
    output.attrs['units'] = units[ii]
    output.to_netcdf(outputfile,mode='w',format='NETCDF3_64BIT') 
    output.close()
    del output
    del data
    del lons
    del lats
    nc.close()
    print('INFO: postprocess_cams_regional.py has been run successfully!')
