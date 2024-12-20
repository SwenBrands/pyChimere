#!/usr/bin/env python

#from netCDF4 import Dataset
import xarray as xr
import numpy as np
import pandas as pd
import sys
import os
import dask
#import numpy.matlib

#call with postprocess_cams_regional_137.py 20180120 00

#set input variables
tardate=str(sys.argv[1])
tarhour=str(sys.argv[2])

#tardate='20181101'
#tarhour='00'

taryear=tardate[0:4]
tarmonth=tardate[4:6]

#define home
homedir = os.getenv("HOME")
lustre = os.getenv("OP_STORE")
bcdir = os.getenv("BCDIR")
cdir = os.getenv("CDIR")
auxdir = os.getenv("AUXDIR")
op_chimere = os.getenv("OP_CHIMERE")

#dimensions of the CAMS input files are (latitude: 451, level: 60, longitude: 900, time: 11)
outformat = 'float64'
source = bcdir+'/'+tardate+'/'+tarhour
fctype = 'prod' #data type, 'prod' or 'test'

variables_in = ('c2h6','hcho','ch4_c','co','hno3','c5h8','no2','go3','pan','so2','aermr01','aermr02','aermr03','aermr04','aermr05','aermr06','aermr08','aermr10','aermr11','lnsp','aermr11','aermr11','aermr11','aermr11','aermr11','aermr11','aermr11','aermr11','aermr11') #set variable name as defined in CAMS
variables_out = ('C2H6','CH2O','CH4','CO','HNO3','ISOP','NO2','O3','PAN','SO2','SS1','SS2','SS3','DUST1','DUST2','DUST3','OM','BC','SO4','PSFC','BIGALK','BIGENE','CH3CHO','GLYOXAL','H2O2','NH3','TOLUENE','C10H16','TEMP') #set variable name as required by CHIMERE
units = ('kg','kg','kg','kg','kg','kg','kg','kg','kg','kg','mmr','mmr','mmr','kg','kg','kg','kg','kg','kg','Pa','ppb','ppb','ppb','ppb','ppb','ppb','ppb','ppb','K')
aero_type = ('gas','gas','gas','gas','gas','gas','gas','gas','gas','gas','aer','aer','aer','dust','dust','dust','aer','aer','aer','none','gas','gas','gas','gas','gas','gas','gas','gas','none') #set aerosol type as required by CHIMERE
aero_low = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.03, 0.5, 5.0, 0.03, 0.55, 0.9, 0.1, 0.1, 0.1, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0] #refers to the aerosol size in um, required by chimere for sea-salt and dust aerosol
aero_high = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.5, 5.0 ,20.0 ,0.55, 0.9, 20.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0] # "
mmass = [30.0, 30.0, 16.0, 28.0, 63.0, 68.0, 46.0, 48.0, 121.0, 64.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 96.0, 1.0, 72.0, 56.0, 44.0, 28.0, 38.0, 17.0, 92.0, 136.0, 1.0] #set the molar mass, required by CHIMERE
critval = np.array(1*10**-12) #works with np.array(1*10**-15) and np.array(1*10**-8)
minval = np.array(1*10**-12) #works with np.array(1*10**-15) and np.array(1*10**-8)

#for unit conversion recall:
# 1 ppb = 1 ug/kg
# 1 ppm = 1 mg/kg

### EXECUTE #############################################################
#get alpha and beta coefficients for the vertical grid of the L60 version of CAMS, BOTH ARE CURRENTLY NOT USED IN THE SCRIPT !!
alpha = np.array((20.000000,38.425343,63.647804,95.636963,134.483307,180.584351,234.779053,298.495789,373.971924,464.618134,575.651001,713.218079,883.660522,1094.834717,1356.474609,1680.640259,2082.273926,2579.888672,3196.421631,3960.291504,4906.708496,6018.019531,7306.631348,8765.053711,10376.126953,12077.446289,13775.325195,15379.805664,16819.474609,18045.183594,19027.695313,19755.109375,20222.205078,20429.863281,20384.480469,20097.402344,19584.330078,18864.750000,17961.357422,16899.468750,15706.447266,14411.124023,13043.218750,11632.758789,10209.500977,8802.356445,7438.803223,6144.314941,4941.778320,3850.913330,2887.696533,2063.779785,1385.912598,855.361755,467.333588,210.393890,65.889244,7.367743,0.000000,0.000000))
beta = np.array((0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000076,0.000461,0.001815,0.005081,0.011143,0.020678,0.034121,0.051690,0.073534,0.099675,0.130023,0.164384,0.202476,0.243933,0.288323,0.335155,0.383892,0.433963,0.484772,0.535710,0.586168,0.635547,0.683269,0.728786,0.771597,0.811253,0.847375,0.879657,0.907884,0.931940,0.951822,0.967645,0.979663,0.988270,0.994019,0.997630,1.000000))

#options for generation of date vector in CHIMERE format
if tarhour == '00':
    if int(tardate) > int('20180806'):
        hoursperday=[1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6]
    else:        
        hoursperday=[1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7, 0, 1]
elif tarhour == '12': #12 UTC available from 20180807 onwards
    hoursperday=[5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2]
else:
    raise Exception('check entry for <tarhour>!')

outputdir=auxdir

#check the number of vertical layers
srcfile_test = source+'/z_cams_c_ecmf_'+tardate+tarhour+'0000_'+fctype+'_fc_ml_003_go3.nc' #load ozone values at forecast hour 03 just to check the number of vertical levels, not used for any other purpose 
nc_test = xr.open_dataset(srcfile_test)
camslevels = nc_test.level.values.shape[0]
print('INFO: CAMS data comes on '+str(camslevels)+' model levels...')
nc_test.close()

#load densities needes for variable conversion and define level-specific variables
source_clim = cdir+'/myclim_regional' #this is the path to the L60 climatological files previously processed with prepare_clim_regional.py
path_modeldens = op_chimere+'/copernicus/model_densities.py' #path to the model densities
execfile(path_modeldens) #load L60 and L137 densities
density = density_60 #L60 densities are used in any case since 60 levels out of 137 levels are cut out in case 137 levels are available (see variable <correspondence>)
if camslevels == 137:
    #index of L137 level nearest ot each L60 level, see https://confluence.ecmwf.int/display/COPSRV/Correspondence+between+60+and+137+model+level+definitions
    correspondence = np.array([6, 9, 11, 13, 15, 16, 18, 19, 21, 22, 24, 25, 27, 29, 31, 33, 35, 38, 40, 43, 46, 50, 53, 56, 59, 63, 66, 69, 72, 75, 77, 80, 82, 85, 87, 89, 92, 94, 96, 98, 99, 101, 103, 105, 107, 109, 111, 113, 115, 117, 119, 121, 123, 126, 128, 130, 132, 134, 136, 137])-1
elif camslevels == 60:
    print('Info: no correspondence needed for 60 model levels...')
else:
    raise Exception('unknown entry for the number of CAMS levles!')

for ii in xrange(len(variables_in)):
    #define variable name in nc file and within the file name
    srcvar = variables_in[ii]
    filenamevar = variables_in[ii]
    
    #rename filenamevar and srcvar in case of ch4_c. Prior to 20180710 it was called 'ch4' in the both the filename and necdf, between 20180710 and 20181120
    #filename was ch4_c and netcdf variable name ch4, from 20181121 both where changed to 'ch4_c'
    if srcvar == 'ch4_c':
        if (int(tardate) > int('20180709')) & (int(tardate) < int('20181121')):         
            srcvar = 'ch4'
            print('INFO: changing srcvar to ch4!')
        elif int(tardate) >= int('20181121'):
            print('INFO: no need to change srcvar nor filenamevar for ch4_c!')
        elif int(tardate) <= int('20180709'):
            srcvar = 'ch4'
            filenamevar = 'ch4'
            print('INFO: changing both, srcvar and filenamevar, to ch4!')        
        else:
            raise Exception('unknown time range!')
    else:
        print('INFO: no need to convert '+srcvar)
    
    tarvar = variables_out[ii]     
    srcfile = source+'/z_cams_c_ecmf_'+tardate+tarhour+'0000_'+fctype+'_fc_ml_*_'+filenamevar+'.nc'
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
    if tarvar in ('BIGALK','BIGENE','CH3CHO','GLYOXAL','H2O2','NH3','TOLUENE','C10H16','TEMP','SS1','SS2','SS3'): #for these variables, the climatological mean values are loaded    
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
        
        data = nc_clim.variables[tarvar][:,latsind_clim,lonsind_clim].astype(outformat).values
        
        if tarvar in ('BIGALK','BIGENE','CH3CHO','GLYOXAL','H2O2','NH3','TOLUENE','C10H16'):
           ##convert kg/kg to ppb (both mass mixing ratios, recall 1 ppb = 1 ug/kg)
           data = data*10**9
           ##then convert ppb (mass missing ratio) to ppb (mixing ratio by volume)
           data = data*(29./mmass[ii]) #29 is the molar mass of air in g/mole, see http://paos.colorado.edu/~toohey/5710faq.html
           print('no conversion needed for '+tarvar)
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
           print('no conversion needed for '+tarvar)
        else:
           raise Exception('check entry for <tarvar>!')
        
        levs = nc_clim.variables['lev'][:].values #load the 60 model levels from the climatological file
                
        ##append climatological x times (where x is the lenght of dates) to mimic a time series, required by option 2 for the lateral boundary conditions in the CHIMERE namelist (chimere.par) 
        #data = data.values
        data = np.tile(data,(len(dates),1,1,1))
        
        #generate xarray data array
        minind = np.where(data<critval)
        data[minind] = minval
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
        
    elif tarvar in ('C2H6','CH2O','CH4','CO','HNO3','ISOP','NO2','O3','PAN','SO2','DUST1','DUST2','DUST3','OM','BC','SO4'):
        #this a 4D variable
        if camslevels == 60:
            data = nc.variables[srcvar][:,:,latsind,lonsind].astype(outformat).values
            levs = nc.variables['level'][:].values #load 3d levels
        elif camslevels == 137:
            data = nc.variables[srcvar][:,correspondence,latsind,lonsind].astype(outformat).values
            levs = np.array(range(1,61)).astype('int32')
        else:
            raise Exception('unknown entry for the number of CAMS levles!')

        #generate xarray data array
        minind = np.where(data<critval)
        data[minind] = minval
        output = xr.DataArray(data, coords=[dates, levs, np.double(lats), np.double(lons)], dims=['time', 'lev', 'lat', 'lon'], name=variables_out[ii])
        output.attrs['refpres'] = 1
        output.attrs['typ'] = str(aero_type[ii])
        output.attrs['molarmass'] = float(mmass[ii])        
        #optionally set the aerosol diameter
        if tarvar in ('DUST1','DUST2','DUST3','OM','BC','SO4'):
            print('INFO: Setting the aerosol diameter for '+tarvar+'...')
            output.attrs['bin_low'] = float(aero_low[ii])
            output.attrs['bin_high'] = float(aero_high[ii])

    #elif srcvar == 'lnsp':
    elif tarvar == 'PSFC':
        #this is a 3D variable
        print('INFO: '+tarvar+' normal, i.e. non logarithmic, pressure has been already obtained by get_cams.sh!')        
        data = nc.variables[srcvar][:,latsind,lonsind].astype(outformat).values
        ##convert logarithmic surface pressure to Pa
        #print('converting log(Pa) to Pa for '+tarvar)
        #data = np.exp(data)
        
        #generate xarray data array
        minind = np.where(data<critval)
        data[minind] = minval
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

#load the single netCDF files, add global attribute <refpres> and merge into one
print('INFO: loading the individual boundary conditions files into a single one...')
loadlist = [auxdir+'/C2H6.nc', auxdir+'/CH2O.nc', auxdir+'/CH4.nc', auxdir+'/CO.nc', auxdir+'/HNO3.nc', auxdir+'/ISOP.nc', auxdir+'/NO2.nc', auxdir+'/O3.nc', auxdir+'/PAN.nc', auxdir+'/SO2.nc', auxdir+'/SS1.nc', auxdir+'/SS2.nc', auxdir+'/SS3.nc', auxdir+'/DUST1.nc', auxdir+'/DUST2.nc', auxdir+'/DUST3.nc', auxdir+'/OM.nc', auxdir+'/BC.nc', auxdir+'/SO4.nc', auxdir+'/PSFC.nc', auxdir+'/BIGALK.nc', auxdir+'/BIGENE.nc', auxdir+'/CH3CHO.nc', auxdir+'/GLYOXAL.nc', auxdir+'/NH3.nc', auxdir+'/H2O2.nc', auxdir+'/TOLUENE.nc', auxdir+'/C10H16.nc', auxdir+'/TEMP.nc', auxdir+'/coeffs.nc']
joint_outputfile =  auxdir+'/gasmet_'+taryear+tarmonth+'.nc'
nc = xr.open_mfdataset(loadlist,engine='netcdf4')
#nc = nc.transpose("time", "lon", "lat", "lev")
nc.attrs['refpres'] = 1.
print('INFO: saving the joint boundary conditions file to NETCDF4 format')
#nc.to_netcdf(joint_outputfile,mode='w',format='NETCDF3_64BIT') #is unaccpetably slow
#nc.load().to_netcdf(joint_outputfile,mode='w',format='NETCDF3_64BIT') is unacceptably slow
nc.to_netcdf(joint_outputfile,mode='w',format='NETCDF4') #this is the fastest option, but requires transformation from netcdf4 to netcdf3 64bit-offset in op_nueva.sh
nc.close()

print('INFO: postprocess_cams_regional_137.py has been run successfully!')
