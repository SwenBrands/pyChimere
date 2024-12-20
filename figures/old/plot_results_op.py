# -*- coding: utf-8 -*-

#load python packages
from netCDF4 import Dataset
import numpy as np
import os
#import xarray as xr
#import cartopy.crs as crs
from matplotlib import pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.basemap import Basemap
from utiles import crea_cmap
import pandas as pd
import os
import sys

#set input variables
startdate=[sys.argv[1]] # defines first input argument, i.e. the start date of the simulation
enddate=[sys.argv[2]] # defines first input argument, i.e. the start date of the simulation
domain=str(sys.argv[3]) # defines the domain to be mapped
species=str(sys.argv[4]) # defines first input argument, i.e. the chemical species

print(startdate)
print(enddate)
print(domain)
print(species)

#aux=startdate[0]
#aux=aux[0:8]

#define the paths
homedir = os.getenv("HOME")
lustre = os.getenv("LUSTRE")
srcpath = lustre+'/OP/PRED/chimere2017r4'
#srcpath = lustre+'/OP/DATOS/RESULTADOS/chimere2017r4/exp0i/'+aux
rundir = homedir+'/OP/LANZAR/chimere2017r4'
scriptpath = homedir+'/OP/LANZAR/chimere2017r4/figures'
savepath = srcpath+'/figs' 
hlev = [0] #(0, 7)

leadtime = 25 #leadtime of the forecast in hours, 25 or 73
starthour = 03 #init of the forecasts

## EXECUTE #############################################################
os.chdir(scriptpath)

BLUE = '#6699cc'
GRAY = '#999999'
under = '#320032'
over  = '#fe00c8'

rgbs = [
 '#960096', #violeta
 '#fa00fe', #violeta
 '#9600fe', #violeta
 '#3200fe', #azul
 '#0032fe', #azul
 '#0096fe', #azul
 '#00e6fe', #azul
 '#0ef2ee', #azul
 '#00e677', #verde
 '#00e650', #verde 
 '#00fa00', #verde
 '#c5ed12', #amarillo
 '#fee100', #amarillo
 '#feae00', #naranja
 '#e67d00', #naranja
 '#e66400', #naranja
 '#c8321d', #marron
 '#aa001d', #marron
 '#c80064', #violeta
 '#f20c86', #violeta
 '#f40e29']
 #'#e0ef0e', 
 
if species == "CO":
    cbounds = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400]
elif species == "O3":
    cbounds = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200]
elif species in ("PM10","PM10bio","PM10ant","PM25","PM25bio","PM25ant","NO2","SO2","pSALT","pDUST"):
    cbounds = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200]
    #cbounds = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100]
    ##cbounds = list(np.array(cbounds)*10)
elif species == "pBCAR":
    cbounds = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
elif species in ("pNA","pH2SO4","pHCL","pWATER"):    
    cbounds = [0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 47.5, 50]
else:
     print('ATTENTION: unknown entry for <species>')

#cbounds = np.arange(0, 200, 10)
#cnorm = colors.BoundaryNorm(boundaries=cbounds, ncolors=len(cbounds-1))

# Creamos el custom colormap:
colormap = crea_cmap(cbounds, rgbs, under, over)
#colormap = 'gist_rainbow_r'

#create a DatetimeIndex
stepdate = str(startdate[0])
dstart = stepdate[6:8]+'/'+stepdate[4:6]+'/'+stepdate[0:4]+' 03:00:00'
rng = pd.date_range(dstart, periods=leadtime, freq='H')

for dd in np.arange(len(startdate)):    
    #filename = srcpath+'/out.'+str(startdate[dd])+'03_'+str(enddate[dd])+'03_'+domain+'.nc'
    filename = srcpath+'/out.'+str(startdate[dd])+'_'+str(enddate[dd])+'_'+domain+'.nc'
    print(filename)
   
    nc = Dataset(filename, mode='r')
    lons = nc.variables['lon'][:]
    lats = nc.variables['lat'][:]
    pres = nc.variables['pres'][:]
    tarvar = nc.variables[species][:]
    tarvar_units = nc.variables[species].units
    dates = nc.variables['Times'][:]
    nc.close()
    #v = np.linspace(-.1, 20, 200, endpoint=True)

    for tt in np.arange(tarvar.shape[0]):
        for hh in np.array(hlev):
            fig1 = plt.figure()
            ax = fig1.add_axes([0.1,0.1,0.8,0.8])
            
            # map setup
            if domain in {'gal05r', 'gal0504r','gal3'}:
                minlat=40.8
                maxlat=44.5
                minlon=-11
                maxlon=-5.4
                lon_0=-9; lat_0=0
                mymap = Basemap(projection='tmerc',lon_0=lon_0,lat_0=lat_0, #tmerc
                k_0=0.9996,rsphere=(6378137.00,6356752.314245179),
                llcrnrlon=minlon,llcrnrlat=minlat,urcrnrlon=maxlon,urcrnrlat=maxlat,resolution='c')
                X, Y = mymap(lons,lats)
                mymap.contourf(X,Y, tarvar[tt,hh,], cbounds, cmap=colormap)
                #mymap.pcolormesh(lons,lats, tarvar[tt,hh,], cmap=colormap, latlon=True, vmin=cbounds[0], vmax=cbounds[-1])
                mymap.readshapefile(scriptpath+'/shapes/municipios', 'municipios',linewidth=0.5, color='k', antialiased=1)          
                mymap.readshapefile(scriptpath+'/shapes/espana', 'espana',linewidth=0.5, color='k', antialiased=1)
                mymap.readshapefile(scriptpath+'/shapes/portugal', 'portugal',linewidth=0.5, color='k', antialiased=1)    
            elif domain in {'pib27','ib15r', 'ib16r', 'km12', 'ib1914r', 'km36', 'gal1511r'}:
                mymap = Basemap(projection='merc', resolution='i', llcrnrlon=np.min(lons),llcrnrlat=np.min(lats),urcrnrlon=np.max(lons),urcrnrlat=np.max(lats))              
                X, Y = mymap(lons,lats) #activate if contourf is used
                mymap.contourf(X,Y, tarvar[tt,hh,], cbounds, cmap=colormap)
                #mymap.pcolormesh(lons,lats, tarvar[tt,hh,], cmap=colormap, latlon=True, vmin=cbounds[0], vmax=cbounds[-1])
                mymap.drawcoastlines()
               # mymap.drawparallels(np.arange(40,48,0.5),labels=[1,0,0,0],linewidth=0.5,color='gray')
               # mymap.drawmeridians(np.arange(-12,-6,0.5),labels=[0,0,0,1],linewidth=0.5,color='gray')
            else:
                print("ATTENTION: Check entry for <domain>")            
           
            #map.contourf(lons, lats, tarvar[tt,hh,], cbounds, cmap=colormap)
                       
            #mymap.ax = ax
            cbar = plt.colorbar()
            cbar.set_label('ug/m3')
            
            ##figure finetuning
            ax.set_title(u'Concentracións de '+species+' no '+str(rng[tt]))
                        
            #if not os.path.exists(savepath+'/'+str(startdate[dd])+'/lvl'+str(hh)):
            #   os.makedirs(savepath+'/'+str(startdate[dd])+'/lvl'+str(hh))
               
            #savename=savepath+'/'+str(startdate[dd])+'/lvl'+str(hh)+'/'+species+'_'+str(startdate[dd])+'_'+str(int(tt)+3)
            
            if not os.path.exists(savepath+'/'+str(domain)+'/'+str(species)):
               os.makedirs(savepath+'/'+str(domain)+'/'+str(species))
               
            savename=savepath+'/'+str(domain)+'/'+str(species)+'/'+str(domain)+'_'+str(species)+'_'+str(startdate[dd])+'_'+str(int(tt)+3)+'.png'
            
            print(savename)
            fig1.savefig(savename, dpi=300)
            plt.close(fig1)
            
#write logfile
logfile=rundir+'/FLAG/map_'+str(startdate[0])+'_'+species+'.flag'
file = open(logfile,'w')
file.write('figures for '+species+' have been generated')
file.close()

