# -*- coding: utf-8 -*-

##this script loads all CHIMERE outputfile in a given directory to a single xarray dataArray, removes the double entries at 03 UTC and then saves to a single nc file

#load python packages
import numpy as np
import os
import xarray as xr
import os
import sys

##set input variables
homedir = os.getenv("HOME")
lustre = os.getenv("STORE_METEO")
chemtype = 'melchior_low_chemonly_reduced'
srcfile = lustre+'/OP/DATOS/RESULTADOS/CHIMERE_NOVO/2021/reduced*nc'
savefile = lustre+'/OP/DATOS/RESULTADOS/CHIMERE_NOVO/2021/chimere_gal3_hourly_surface_prediction_2021010103_2022010102.nc'

#srcfile = lustre+'/OP/DATOS/RESULTADOS/chimere2017r4/hindcast_05_gal3/reduced*nc'
#savefile = lustre+'/OP/DATOS/RESULTADOS/chimere2017r4/hindcast_05_gal3/chimere_gal3_hourly_surface_2018010103_2019010102.nc'

#level = 0

##set 3d variables to be retained
#melchior, low output detail
if chemtype == 'melchior_low':
    outvars = ['O3','NO2','NO','PAN','HNO3','H2O2','HONO','SO2','CO','OH','NOX','TOTPAN','NOY','HOX','ROX','ANT','BIO','COVS','ANTS','OX','pDUST','pBCAR','pOCAR', \
            'pPPM','pSALT','pAnA1D','pAnBmP','pBiA1D','pBiBmP','pISOPA1','pH2SO4','pHNO3','pNH3','pWATER','PM25','PM25bio','PM25ant','PM10','PM10bio','PM10ant','hlay','thlay', \
            'airm','relh','temp','winz','winm','winw','sphu','kzzz','clwc','cliq','cice','pres','dpeu','dped','dpdu','dpdd','flxu','flxd','jO3','jNO2','hght']
elif chemtype == 'melchior_low_chemonly':
    outvars = ['O3','NO2','NO','PAN','HNO3','H2O2','HONO','SO2','CO','OH','NOX','TOTPAN','NOY','HOX','ROX','ANT','BIO','COVS','ANTS','OX','pDUST','pBCAR','pOCAR', \
            'pPPM','pSALT','pAnA1D','pAnBmP','pBiA1D','pBiBmP','pISOPA1','pH2SO4','pHNO3','pNH3','pWATER','PM25','PM25bio','PM25ant','PM10','PM10bio','PM10ant','hght']
elif chemtype == 'melchior_low_chemonly_reduced':
    outvars = ['O3','NO2','NO','PAN','HNO3','H2O2','HONO','SO2','CO','NOX','TOTPAN','pDUST','pBCAR','pOCAR', \
            'pPPM','pSALT','pISOPA1','pH2SO4','pHNO3','pNH3','PM25','PM25bio','PM25ant','PM10','PM10bio','PM10ant','hght']
#melchior, high output detail
elif chemtype == 'melchior_full':
    outvars = ['O3','NO2','NO','PAN','HNO3','H2O2','HONO','SO2','CO','CH4','C2H6','NC4H10','C2H4','C3H6','OXYL','C5H8','HCHO','CH3CHO','GLYOX','MGLYOX','CH3COE','NH3', \
            'APINEN','BPINEN','LIMONE','TERPEN','OCIMEN','HUMULE','TOL','TMB','AnA1D','AnBmP','BiA1D','BiBmP','ISOPA1','PPMAQ','H2SO4AQ','HNO3AQ','NH3AQ','AnA1DAQ', \
            'AnBmPAQ','BiA1DAQ','BiBmPAQ','ISOPA1AQ','SALTAQ','DUSTAQ','OCARAQ','BCARAQ','TOTPAN','NOX','OX','HOX','NOY','ROOH','HCNM','RNO3','pBCAR','pDUST','pOCAR', \
            'pPPM','pSALT','pAnA1D','pAnBmP','pBiA1D','pBiBmP','pISOPA1','pH2SO4','pHNO3','pNH3','pWATER','PM25','PM25bio','PM25ant','PM10','PM10bio','PM10ant','hlay', \
            'thlay','airm','relh','temp','winz','winm','winw','sphu','kzzz','clwc','cliq','cice','pres','dpeu','dped','dpdu','dpdd','flxu','flxd','jO3','jNO2','hght']
#saprc, high output detail
elif chemtype == 'saprc_full':
    outvars = ['O3','NO2','AACD','ACET','ACYE','AFG1','AFG2','AFG3','ALK1','ALK2','ALK3','ALK4','ALK5','APINEN','ARO1','ARO2','BACL','BALD','BENZ','BPINEN','BZCO3','BZO','C2H4', \
            'C5H8','CCHO','CH4','CO','COOH','CRES','FACD','GLY','H2','H2O2','HCHO','HNO3','HNO4','HO2','HONO','HUMULE','IPRD','LIMONE','MACO3','MACR','MAPAN','MECO3','MEK', \
            'MEO2','MEOH','MGLY','MVK','N2O5','NO','NO3','NPHE','OCIMEN','OH','OLE1','OLE2','PACD','PAN','PAN2','PBZN','PRD2','R2O2','RCHO','RCO3','RNO3','RO2N','RO2R','SO2', \
            'H2SO4','TBUO','TERPEN','XOOH','TOL','TMB','AnA1D','AnBmP','BiA1D','BiBmP','ISOPA1','PPMAQ','H2SO4AQ','HNO3AQ','NH3AQ','AnA1DAQ','AnBmPAQ','BiA1DAQ','BiBmPAQ', \
            'ISOPA1AQ','SALTAQ','DUSTAQ','OCARAQ','BCARAQ','NOX','TOTPAN','NOY','HOX','ROX','ANT','BIO','COVS','ANTS','OX','pBCAR','pDUST','pOCAR','pPPM','pSALT','pAnA1D', \
            'pAnBmP','pBiA1D','pBiBmP','pISOPA1','pH2SO4','pHNO3','pNH3','pWATER','PM25','PM25bio','PM25ant','PM10','PM10bio','PM10ant','hlay','thlay','airm','relh','temp', \
            'winz','winm','winw','sphu','kzzz','clwc','cliq','cice','pres','dpeu','dped','dpdu','dpdd','flxu','flxd','jO3','jNO2','hght']
else:
    raise Exception('check entry for <chemtype>!')

##EXECUTE ##############################################################
nc = xr.open_mfdataset(srcfile,concat_dim='Times')
#extract common variables
times = nc.Times.values
#times = np.delete(times, np.arange(24, times.shape[0], 25)) #delete the double 03 UTC entries
lats = nc.lat.values[0,] #open_mfdataset erroneusly stacks the lats and lons as many times as there are entries in 'Times', this is corrected here
lons = nc.lon.values[0,]
#load variables in a loop
for vv in range(len(outvars)):
    print('INFO: processing '+outvars[vv])
    #variable = nc.variables[outvars[vv]][:,level,:,:]
    variable = nc.variables[outvars[vv]]
    data = variable.values
    #data = np.delete(data, np.arange(24, data.shape[0], 25), axis=0) #delete the double 03 UTC entries
    xrarr = xr.DataArray(data, coords=[times, np.unique(lats), np.unique(lons)], dims=['Times', 'south_north', 'west_east'], attrs=variable.attrs, name=outvars[vv])
    #create a dataset of assign variable to this dataset    
    if vv == 0:
        ds = xrarr.to_dataset(name = outvars[vv])
    else:
        ds[outvars[vv]] = xrarr        
    del(variable,data,xrarr)
#add latitude and longitude as 2d variables
xrarr = xr.DataArray(lats, coords=[np.unique(lats), np.unique(lons)], dims=['south_north', 'west_east'], name='lat')
ds['lat'] = xrarr
xrarr = xr.DataArray(lons, coords=[np.unique(lats), np.unique(lons)], dims=['south_north', 'west_east'], name='lon')
ds['lon'] = xrarr 
nc.close()

#save the output in nc format
ds.to_netcdf(savefile)
ds.close()

    
