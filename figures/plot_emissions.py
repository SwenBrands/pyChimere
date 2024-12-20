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

#set input variables
import sys
domain=str(sys.argv[1]) # defines the domain to be mapped, CAPITAL letters
species=str(sys.argv[2]) # defines first input argument, i.e. the chemical species
tarmonth=str(sys.argv[3]) # defines the month to be mapped, TWO numbers, e.g. 02
emistype=str(sys.argv[4]) # defines the emission type, either area (s) o pointwise (p)

print(domain)
print(species)
print(tarmonth)
print(emistype)

lustre = os.getenv("LUSTRE")
srcpath = lustre+'/OP/DATOS/CICC/chimere2017r4/EMISSIONS/EMEP_GAL3_pop1r1ag1sh0lsp0prtr1_globcover_melchior/EMI_2018'
#srcpath = lustre+'/OP/DATOS/CICC/chimere2017r4/EMISSIONS/myEMEP01x01_GAL3_pop1r1ag1sh0lsp0_globcover_melchior/EMI_2018'

hourvec = [0, 5, 10, 15, 20] #daytime hours to be plotted
#dayvec = [0, 3, 6] #days of the week to be plotted
dayvec = [3]
colormap = 'PiYG_r'
markersize = 20
minval = 1

filename = srcpath+'/EMIS.'+str(domain)+'.'+str(tarmonth)+'.'+str(species)+'.'+emistype+'.nc'    
print(filename)
nc = xr.open_dataset(filename)
lons = nc.variables['lon'].values
lats = nc.variables['lat'].values

if emistype == 's':
    height = nc.variables['EMEP_levels'].values
elif emistype == 'p':
    height = [20]
else:
    raise Exception('check entry for <emistype>!')

print('heigth levels in meters are: '+str(height))
VAR = nc.variables[species].values
VAR_units = str(nc.variables[species].attrs['units'])
nc.close()

##ideal colorbar for SO2
#cbounds = [0, 10**9, 10*10**9, 15*10**9, 20*10**9, 25*10**9, 30*10**9, 35*10**9, 40*10**9, 45*10**9, 50*10**9, 55*10**9, 60*10**9, 65*10**9, 70*10**9, 75*10**9, 80*10**9]

if species in ('NO2','PPM_coa','PPM_fin'):
    cbmax = 0.5*10**11
elif species == 'APINEN':
    cbmax = 1*10**8
elif species == 'NH3':
    cbmax = 2*10**12
elif species == 'SO2':
    cbmax = 0.6*10**13
else:
    cbmax = np.max(VAR)

cbounds = np.linspace(0,cbmax,31)

minlon = np.min(lons)
minlat = np.min(lats)
maxlon = np.max(lons)
maxlat = np.max(lats)
halfres_lat = np.round(np.max(np.gradient(lats))/2,4)
halfres_lon = np.round(np.max(np.gradient(lons))/2,4)

for tt in xrange(len(hourvec)): #time of the day 1-24
    for dd in xrange(len(dayvec)): # day of the week 1-7
        for hh in xrange(len(height)): #EMEP height level
            if emistype == 's':
                tarvar = VAR[hourvec[tt],dayvec[dd],hh,:,:]
                savepath = srcpath+'/figs_area'
            elif emistype == 'p':
                tarvar = VAR[hourvec[tt],dayvec[dd],:]
                lons_point = nc.Sources_lon.values
                lats_point = nc.Sources_lat.values
                savepath = srcpath+'/figs_point'
                ##exclude very small emission points
                outind = np.where(tarvar < minval)[0]
                tarvar = np.delete(tarvar, outind)
                lons_point = np.delete(lons_point, outind)
                lats_point = np.delete(lats_point, outind) 
            else:
                raise Exception('check entry for <emistype>!')
            
            ##get the colorbar intervals
            #cbounds = np.linspace(np.min(tarvar),np.percentile(tarvar,90),20)#
            
            fig1 = plt.figure()
            #ax = fig1.add_axes([0.1,0.1,0.8,0.8])
            # map setup
            mymap = Basemap(projection='cyl', resolution='i', llcrnrlon=minlon,llcrnrlat=minlat,urcrnrlon=maxlon,urcrnrlat=maxlat)   
            if emistype == 's':
                X, Y = mymap(lons,lats)
                mymap.pcolormesh(np.round(lons-halfres_lon,4),np.round(lats-halfres_lat,4), tarvar, cmap=colormap, latlon=True, vmin=np.array(cbounds[0]), vmax=np.array(cbounds[-1]))
                if domain in {'GAL3','GAL0504R','GAL0504S'}:
                   mymap.readshapefile('./shapes/municipios', 'municipios',linewidth=0.5, color='k', antialiased=1)
                   mymap.readshapefile('./shapes/espana', 'espana',linewidth=0.5, color='k', antialiased=1)
                   mymap.readshapefile('./shapes/portugal', 'portugal',linewidth=0.5, color='k', antialiased=1)
                elif domain in {'IB16R','PIB27','GAL15R'}:
                   mymap.readshapefile('./shapes/espana', 'espana',linewidth=0.5, color='k', antialiased=1)
                   mymap.readshapefile('./shapes/portugal', 'portugal',linewidth=0.5, color='k', antialiased=1)  
                else:
                  print("ATTENTION: Check entry for <domain>")

            elif emistype == 'p':
                #mymap.scatter(lons_point,lats_point, s=2, c=tarvar, cmap=colormap, vmin=np.array(cbounds[0]), vmax=np.array(cbounds[-1]))
                mymap.scatter(lons_point,lats_point, s=markersize, c=tarvar, cmap=colormap, vmin=np.array(cbounds[0]), vmax=np.array(cbounds[-1]))
                mymap.readshapefile('./shapes/municipios', 'municipios',linewidth=0.5, color='k', antialiased=1)
                mymap.readshapefile('./shapes/espana', 'espana',linewidth=0.5, color='k', antialiased=1)
                mymap.readshapefile('./shapes/portugal', 'portugal',linewidth=0.5, color='k', antialiased=1)    
            else:
                raise Exception('check entry for <emistype>!')

            cbar = mymap.colorbar(shrink=0.8)
            cbar.set_label(species+' in '+VAR_units)
                        
            #figure finetuning
            #ax.set_title(str(species )+', month '+tarmonth+', day '+str(dayvec[dd])+', hour '+str(hourvec[tt])+' at '+str(height[hh])+' m')
            plt.title(str(species)+', month '+tarmonth+', day '+str(dayvec[dd])+', hour '+str(hourvec[tt])+' at '+str(height[hh])+' m')    

            if not os.path.exists(savepath+'/'+str(species)):
               os.makedirs(savepath+'/'+str(species))
                           
            savename=savepath+'/'+str(species)+'/'+str(domain)+'_'+str(species)+'_'+str(tarmonth)+'_'+str(dayvec[dd])+'_'+str(hourvec[tt])+'_'+str(height[hh])+'.png'
                        
            print(savename)
            fig1.savefig(savename, dpi=300)
            plt.close(fig1)

