#!/usr/bin/env python
# -*- coding: utf-8 -*-

# import xarray as xr
# import numpy as np
# import pandas as pd
# import sys
# import dask
# import os
# import numpy.matlib
# import matplotlib
# matplotlib.use('Agg')
# from matplotlib import pyplot as plt
# import matplotlib.colors as colors
# from mpl_toolkits.basemap import Basemap, addcyclic

import numpy as np
import os
import xarray as xr
from mpl_toolkits.basemap import Basemap, addcyclic
import pandas as pd
import os
import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import matplotlib.colors as colors
import Ngl
from utiles import crea_cmap

tardate=str(sys.argv[1])
tarhour=str(sys.argv[2])
variable=str(sys.argv[3])
level=int(sys.argv[4])

#tardate='20180916'
#tarhour='12'
#variable = ('aermr06') #set variable name as defined in CAMS
#level = 59 #model level

colormap = 'gist_ncar'
leadtime=87
mmass = 100.0
filetype = 'test' #'prod' or 'test'

##set input variables
#variable = ('c2h6','hcho','ch4','co','hno3','c5h8','no2','go3','pan','so2','aermr01','aermr02','aermr03','aermr04','aermr05','aermr06','aermr08','aermr10','aermr11','lnsp','aermr11','aermr11','aermr11','aermr11','aermr11','aermr11','aermr11','aermr11','aermr11') #set variable name as defined in CAMS
#mmass = [30.0, 30.0, 16.0, 28.0, 63.0, 68.0, 46.0, 48.0, 121.0, 64.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 96.0, 1.0, 72.0, 56.0, 44.0, 28.0, 38.0, 17.0, 92.0, 136.0, 1.0] #set the molar mass, required by CHIMERE

#define home
homedir = os.getenv("HOME")
lustre = os.getenv("STORE_METEO")

#dimensions of the CAMS input files are (latitude: 451, level: 60, longitude: 900, time: 11)
srcdir = lustre+'/OP/DATOS/CICC/CAMS_137/'+tardate+'/'+tarhour

#for unit conversion recall:
# 1 ppb = 1 ug/kg
# 1 ppm = 1 mg/kg

### EXECUTE ############################################################
taryear=tardate[0:4]
tarmonth=tardate[4:6]
density = np.array((0.000150,0.000402,0.000672,0.001025,0.001485,0.002087,0.002814,0.003688,0.004740,0.006018,0.007599,0.009582,0.012082,0.015130,0.018863,0.023518,0.029322,0.036557,0.045579,0.056826,0.070849,0.087826,0.107127,0.129274,0.154328,0.182375,0.213467,0.247616,0.284796,0.324943,0.367189,0.403822,0.441886,0.481280,0.521847,0.563389,0.605674,0.648445,0.691426,0.734331,0.776863,0.818729,0.859634,0.899295,0.937436,0.973801,1.008150,1.040270,1.069971,1.097099,1.121533,1.143192,1.162039,1.178087,1.191403,1.202112,1.210406,1.216546,1.220871,1.223803))

if variable == "CO":
    cbounds = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400]
elif variable == "O3":
    cbounds = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200]
elif variable in ('aermr01'):
    cbounds = np.array([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200]) *10**-10
elif variable in ('aermr02'):
    cbounds = np.array([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200]) *10**-9
elif variable in ('aermr03'):
    cbounds = np.array([0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400])*10**-8
elif variable in ('aermr04','aermr05','aermr06'):
    cbounds = np.array([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200]) *10**-10
elif variable == "pBCAR":
    cbounds = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
elif variable in ("pNA","pH2SO4","pHCL","pWATER"):    
    cbounds = [0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 47.5, 50]
else:
     print('ATTENTION: unknown entry for <species>')

#options for generation of date vector in CHIMERE format
if tarhour == '00':
    hoursperday=[1, 2, 3, 4, 5, 6, 7, 0]
elif tarhour == '12':
    hoursperday=[5, 6, 7, 0, 1, 2, 3, 4]
else:
    print('ATTENTION: unknown entry for <species>')
hoursperday=numpy.matlib.repmat(hoursperday,1,leadtime/24+1)
hoursperday=hoursperday[0]
hoursperday=hoursperday[0:-2]

#load nc files
srcfile = srcdir+'/z_cams_c_ecmf_'+tardate+tarhour+'0000_'+filetype+'_fc_ml_*_'+variable+'.nc'
print('loading files '+srcfile)        
nc = xr.open_mfdataset(srcfile, concat_dim='time', engine='netcdf4')
#load and process lat and lon
lats = nc.latitude.values
lons = nc.longitude.values
mask = np.where(lons>180)[0]
lons[mask] = lons[mask] - 360 
XX,YY = np.meshgrid(lons,lats)
        
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
        
#variable in ('SS1','SS2','SS3','DUST1','DUST2','DUST3','OM','BC','SO4'):
if variable in ('aermr01','aermr02','aermr03','aermr04','aermr05','aermr06','aermr08','aermr10','aermr11'):
   #this a 4D variable
    data = eval('nc.'+variable)
    
    #print('converting kg/kg to ppb for '+variable)
    ## 1 ppb = 0.001 mg/kg
    #print('converting kg/kg to ug/m3 for '+variable)
    #data = data*10**9 #convert kg/kg to ppb
    ##then weigh by air density
    #dens4d = np.tile(density,(data.shape[0],data.shape[2],data.shape[3],1))    
    #dens4d = np.moveaxis(dens4d,-1,1)
    #data = data*dens4d               
    
    levs = nc.variables['level'][:] #load 3d levels
elif variable in ('C2H6','CH2O','CH4','CO','HNO3','ISOP','NO2','O3','PAN','SO2'):
    print('converting kg/kg to ppb for '+variable)        
    data = nc.variables[variable][:,:,latsind,lonsind]
    #convert kg/kg to ppb (both mass mixing ratios, recall 1 ppb = 1 ug/kg)
    data = data*10**9
    #then convert ppb (mass missing ratio) to ppb (mixing ratio by volume)
    data = data*(29./mmass[ii]) #29 is the molar mass of air in g/mole, see http://paos.colorado.edu/~toohey/5710faq.html
    levs = nc.variables['level'][:] #load 3d levels        
             
elif srcvar == 'lnsp':
    #this is a 3D variable
    print('INFO: '+variable+' has been already postprocessed by get_cams.sh!')        
    data = nc.variables[variable][:,latsind,lonsind]
    ##convert logarithmic surface pressure to Pa
    #print('converting log(Pa) to Pa for '+variable)
    #data = np.exp(data)        
else:        
    raise Exception('CAMS input variable not found!')  
    
##then get the model level and plot the data
print('INFO: retain model level: '+str(level))
indlev = np.where(levs==level)[0][0]
data = data[:,indlev,:,:].values

if not os.path.exists(srcdir+'/figs'):
   os.makedirs(srcdir+'/figs')

##plot the variable for each time step
#indlon = np.where(lons > 180)
#lons[indlon] = lons[indlon]-360
X,Y = np.meshgrid(lons,lats)
for tt in range(len(dates)):
    fig1 = plt.figure()
    mymap = Basemap(projection='cyl', resolution='l', llcrnrlon=np.min(lons),llcrnrlat=np.min(lats),urcrnrlon=np.max(lons),urcrnrlat=np.max(lats))    
    #mymap.pcolormesh(lons,lats, data[tt,:,:], latlon=True, cmap=colormap, vmin=cbounds[0], vmax=cbounds[-1])
    mymap.contourf(XX,YY, data[tt,:,:], cbounds, cmap=colormap, vmin=cbounds[0], vmax=cbounds[-1])
    mymap.drawcoastlines()
    plt.colorbar()
    plt.title(str(dates[tt]))
    #plt.contourf(data[tt,:,:])
    #ax = plt.axes(projection=crs.PlateCarree())
    #ax.coastlines(lw=1) 
    #data[tt,:,:].plot()
    savename=srcdir+'/figs/'+str(variable)+'_level'+str(level)+'_'+str(tt)+'.png'
    print(savename)
    fig1.savefig(savename, dpi=300)
    plt.close(fig1)
    
nc.close()    




