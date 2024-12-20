# -*- coding: utf-8 -*-

#load python packages
import numpy as np
import os
#from netCDF4 import Dataset
import xarray as xr
from mpl_toolkits.basemap import Basemap
from utiles import crea_cmap
import pandas as pd
import os
import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import matplotlib.colors as colors
import Ngl

##set input variables
startdate=str(sys.argv[1]) # defines the start date and hour of the simulation
enddate=str(sys.argv[2]) # defines the end date and hour of the simulation
domain=str(sys.argv[3]) # defines the domain to be mapped, lowercase letters
species=str(sys.argv[4]) # defines first input argument, i.e. the chemical species

print(startdate)
print(enddate)
print(domain)
print(species)

#define the paths
homedir = os.getenv("HOME")
lustre = os.getenv("LUSTRE")
srcpath = lustre+'/OP/PRED/chimere2017r4'
#srcpath = lustre+'/OP/DATOS/RESULTADOS/chimere2017r4/h2v3e1p4/'+str(startdate)
rundir = homedir+'/OP/LANZAR/chimere2017r4'
scriptpath = homedir+'/OP/LANZAR/chimere2017r4/figures'
savepath = srcpath+'/figs' 
level = [15]
plimit = 0.5 #precipitation limit in mm
barbwidth = 0.5
leadtime = 25 #leadtime of the forecast in hours, 25 or 73
starthour = 03 #init of the forecasts
lat1 = 40.1 #start latitude of the SW to NE profile
lat2 = 47.9 #start latitutde of the NW to SE profile

## EXECUTE #############################################################
os.chdir(scriptpath)

#define pressure level for interpolation from model levels
plevs = np.array([1000., 990., 970., 950., 925., 900., 850., 800., 750., 700., 650., 600., 550., 500.]) #must be hPa

#define the colorbar
BLUE = '#6699cc'
GRAY = '#999999'
under = '#320032'
#over = '#fe00c8'
over = '#f40e29'

rgbs = [
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
 '#fa00fe', #violeta
 '#9600fe', #violeta
 '#960096'] #violeta
 
if species == "CO":
    cbounds = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400]
elif species == "O3":
    cbounds = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200]
elif species in ("PM10","PM10bio","PM10ant","PM25","PM25bio","PM25ant","NO2","SO2","pSALT","pDUST","pOCAR"):
    cbounds = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100]
    #cbounds = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200]
    ##cbounds = list(np.array(cbounds)*10)
elif species in ("pBCAR","pH2SO4"):
    cbounds = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
elif species in ("pNA","pH2SO4","pHCL","pWATER"):
    cbounds = [0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 47.5, 50]
elif species in ("swrd"):
    cbounds = [0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375, 400, 425, 450, 475, 500]
elif species in ('atte'):
    cbounds = [0.8, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.0]
else:
     print('ATTENTION: unknown entry for <species>')

# Creamos el custom colormap:
colormap = crea_cmap(cbounds, rgbs, under, over)
 
filename = srcpath+'/out.'+str(startdate)+'_'+str(enddate)+'_'+domain+'.nc'
#filename = srcpath+'/out.'+str(startdate)+'03_'+str(enddate)+'03_'+domain+'.nc'
print(filename)
nc = xr.open_dataset(filename)
lons = nc.lon.values
lats = nc.lat.values
u10 = nc.u10m.values
v10 = nc.v10m.values
precip = nc.topc.values
wspeed = (u10 ** 2 + v10 ** 2) ** 0.5
pres = nc.pres.values #load pressure, must be Pa 
pres = pres[:,0,:,:] #retain surface pressure
tarvar = nc.variables[species].values
tarvar_units = str(nc.variables[species].attrs.get('units'))
dates = nc.Times.values.tolist()
#acoeff = np.flipud(nc.a_vcoord.values)
acoeff = np.zeros(len(nc.a_vcoord.values)) #create fake hybrid a coefficients vector with P0 = 1000 as described at https://www.ncl.ucar.edu/Applications/vert_interp.shtml
#acoeff[-1] = 100000.
bcoeff = np.flipud(nc.b_vcoord.values)
nc.close()

#create reduced x and y coordinates for plotting wind vector
xred = range(0,lons.shape[0],2)
yred = range(0,lons.shape[1],2)
#set default variables for Basemap
minlat = np.min(lats)
maxlat = np.max(lats)
minlon = np.min(lons)
maxlon = np.max(lons)
lon_0=-9; lat_0=0

#define barbs and crosses
if domain in {'gal05r', 'gal0504r','gal3'}:
    barblength = 4
    msize = 3
elif domain in {'pib27','ib15r', 'ib16r', 'ib1914r', 'gal1511r', 'gal1511r2', 'gal16r', 'gal15r', 'gal15r2'}:
    barblength = 4
    msize = 3
elif domain in {'km36', 'km12'}:
    barblength = 3
    msize = 2
else:
    print("ATTENTION: Check entry for <domain>")

#get parallels and meridians
parallels = np.arange(np.ceil(minlat),np.floor(maxlat+1))
meridians = np.arange(np.ceil(minlon),np.floor(maxlon+1))

startlat1 = lats[np.where(lats[:,0] == lat1)[0],0]
startlat2 = lats[np.where(lats[:,0] == lat2)[0],0]
startlon = np.min(lons)
latlimit1 = np.max(lats) #limit at the NE edge of the domain
latlimit2 = np.min(lats) #limit at the SE edge of the domain
latres = np.round(np.nanmax(np.abs(np.gradient(lats))),2)
lonres = np.round(np.nanmax(np.abs(np.gradient(lons))),2)

#construct coordinate array for the first profile
stopper1 = startlat1
diag_ind = np.where((lats == startlat1) & (lons == startlon))
runind = 0
while stopper1 <= latlimit1:
    diag_ind_step = np.where((lats == startlat1+runind*latres) & (lons == startlon+runind*lonres))
    diag_ind = np.append(np.array(diag_ind),np.array(diag_ind_step),axis=1)
    stopper1 = stopper1 + latres
    runind = runind+1

lon_diag1 = np.diagonal(lons,offset=offset1,axis1=0,axis2=1)
lat_diag1 = np.diagonal(lats,offset=offset1,axis1=0,axis2=1)
Xvprof1,Yvprof1 = np.meshgrid(lon_diag1,plevs)

lon_diag2 = np.diagonal(np.fliplr(lons),offset=offset2,axis1=0,axis2=1)
lat_diag2 = np.diagonal(np.fliplr(lats),offset=offset2,axis1=0,axis2=1)
Xvprof2,Yvprof2 = np.meshgrid(lon_diag2,plevs)

#get startpoints for the profiles
lons_sorted = np.sort(np.unique(lons))
lats_sorted = np.sort(np.unique(lats))

if len(tarvar.shape) == 3: #for 2d variables
    print('INFO: '+species+' is a 2d variable.')
    tarvar_2d = tarvar
elif len(tarvar.shape) == 4: #for 3d variables
    print('INFO: '+species+' is a 3d variable, cut out model level '+str(level))
    tarvar_2d = np.squeeze(tarvar[:,level,:,:])
    tarvar_flip = np.flip(tarvar,axis=1) #for interpolation to pressure levels with ngl
    tarvar_pres = Ngl.vinth2p(tarvar_flip, acoeff, bcoeff, plevs, pres, intyp=1, p0=1000., ii=1, kxtrp=True) #p0 must be hPa
    tarvar_diag1 = tarvar_pres.diagonal(offset=offset1,axis1=2,axis2=3) #get the SW to NE diagonal
    
    #flip the lon lat matrix in the 4d array tarvar_pres left-right
    tarvar_pres_lr = np.zeros(tarvar_pres.shape)
    for tt in list(range(tarvar_pres_lr.shape[0])):
        for ll in list(range(tarvar_pres_lr.shape[1])):
            tarvar_pres_lr[tt,ll,:,:] = np.fliplr(tarvar_pres[tt,ll,:,:])            
    tarvar_diag2 = tarvar_pres_lr.diagonal(offset=offset2,axis1=2,axis2=3) #get the SW to NE diagonal #get the SW to NE diagonal #get the NW to SE diagonal
else:
    raise Exception('ERROR: invalid number of dimensions for <tarva>!')

for tt in np.arange(tarvar_2d.shape[0]):
    fig1 = plt.figure()
    ax = fig1.add_axes([0.1,0.1,0.8,0.8])
    if domain in {'gal05r', 'gal0504r','gal3'}:
        #mymap = Basemap(projection='tmerc',lon_0=lon_0,lat_0=lat_0, k_0=0.9996, rsphere=(6378137.00,6356752.314245179), llcrnrlon=minlon,llcrnrlat=minlat,urcrnrlon=maxlon,urcrnrlat=maxlat,resolution='c')
        mymap = Basemap(projection='cea',llcrnrlon=minlon,llcrnrlat=minlat,urcrnrlon=maxlon,urcrnrlat=maxlat,resolution='c',lat_ts=43)#lat_ts=43
        X, Y = mymap(lons,lats)        
        mymap.contourf(X,Y, np.squeeze(tarvar_2d[tt,:,:]), cbounds, cmap=colormap)
        mymap.readshapefile(scriptpath+'/shapes/municipios', 'municipios',linewidth=0.35, color='k', antialiased=1)          
        mymap.readshapefile(scriptpath+'/shapes/espana', 'espana',linewidth=0.35, color='k', antialiased=1)
        mymap.readshapefile(scriptpath+'/shapes/portugal', 'portugal',linewidth=0.35, color='k', antialiased=1)
    elif domain in {'pib27','ib15r', 'ib16r', 'km12', 'ib1914r', 'km36', 'gal1511r', 'gal1511r2', 'gal16r', 'gal15r', 'gal15r2'}:
        mymap = Basemap(projection='cea', resolution='l', llcrnrlon=np.min(lons),llcrnrlat=np.min(lats),urcrnrlon=np.max(lons),urcrnrlat=np.max(lats),lat_ts=44)
        X, Y = mymap(lons,lats) #activate if contourf is used
        mymap.contourf(X, Y, np.squeeze(tarvar_2d[tt,:,:]), cbounds, cmap=colormap)
        mymap.drawcoastlines()
    else:
        print("ATTENTION: Check entry for <domain>")
    #plot precipitation occurrence
    idp = np.where(precip[tt,] > plimit)
    plt.plot(X[idp], Y[idp], color='black', marker='x', linestyle='none', markersize=msize)
    #then plot the wind barbs
    X = X[xred,:]
    X = X[:,yred]
    Y = Y[xred,:]
    Y = Y[:,yred]
    u10step = u10[tt,xred,:]
    u10step = u10step[:,yred]
    v10step = v10[tt,xred,:]
    v10step = v10step[:,yred]
    mymap.barbs(X, Y, u10step, v10step, length=barblength, pivot='middle', barbcolor='grey', lw=barbwidth, barb_increments=dict(half=5, full=10, flag=46))
    
    cbar = plt.colorbar()
    cbar.set_label('ug/m3')
            
    ##figure finetuning
    ax.set_title(u'Concentracións de '+species+' no '+dates[tt])
    #mymap.drawparallels(parallels,labels=[True,False,False,True], color='None', size=4)
    #mymap.drawmeridians(meridians,labels=[True,False,False,True], color='None', size=4)
                      
    if not os.path.exists(savepath+'/'+str(domain)+'/'+str(species)):
       os.makedirs(savepath+'/'+str(domain)+'/'+str(species))
              
    savename=savepath+'/'+str(domain)+'/'+str(species)+'/'+str(domain)+'_'+str(species)+'_'+dates[tt]+'.png'
    #savename=str(domain)+'_'+str(species)+'_'+str(startdate)+'_'+str(int(tt)+3)+'.png'
    #savename=str(domain)+'_'+str(species)+'_'+dates[tt]+'.png'
            
    print(savename)
    fig1.savefig(savename, dpi=300)
    plt.close(fig1)
    
    #plot vertical profile 1 (from SW to NE)
    fig_vprof1 = plt.figure()
    plt.contourf(Xvprof1,Yvprof1,tarvar_diag1[tt,:,:],cmap=colormap,vmin=cbounds[0],vmax=cbounds[-1])
    plt.plot(lon_diag1,pres[tt,:,:].diagonal(offset=offset1,axis1=0,axis2=1)/100.,'-k') #plot first model level to mimic surface orography
    #figure finetuning
    plt.ylim(plevs[0],plevs[-1])
    #plt.gca().invert_yaxis()    
    cbar = plt.colorbar()
    cbar.set_label(tarvar_units)    
    xlabels1 = [str(lon_diag1[ii])+', '+str(lat_diag1[ii]) for ii in list(range(len(lon_diag1)))]
    plt.xticks(lon_diag1[0::4],xlabels1[0::4],rotation=45,fontsize=6.)
    plt.ylabel(u'Nivel de presión (hPa)')
    plt.xlabel('Coordenadas de suroeste a noreste')
    
    #save the figure
    savename = savepath+'/'+str(domain)+'/'+str(species)+'/vprof_diag1_'+str(domain)+'_'+str(species)+'_'+dates[tt]+'.png'
    fig_vprof1.savefig(savename)
    plt.close(fig_vprof1)
    
    #plot vertical profile 2 (from SW to NE)
    fig_vprof2 = plt.figure()
    plt.contourf(Xvprof2,Yvprof2,tarvar_diag2[tt,:,:],cmap=colormap,vmin=cbounds[0],vmax=cbounds[-1])
    plt.plot(lon_diag2,np.fliplr(pres[tt,:,:]).diagonal(offset=offset2,axis1=0,axis2=1)/100.,'-k') #plot first model level to mimic surface orography
    #figure finetuning
    plt.ylim(plevs[0],plevs[-1])
    plt.xlim(lon_diag2[0],lon_diag2[-1])
    #plt.gca().invert_yaxis()    
    cbar = plt.colorbar()
    cbar.set_label(tarvar_units)    
    xlabels2 = [str(lon_diag2[ii])+', '+str(lat_diag2[ii]) for ii in list(range(len(lon_diag2)))]
    plt.xticks(lon_diag2[0::4],xlabels2[0::4],rotation=45,fontsize=6.)
    plt.ylabel(u'Nivel de presión (hPa)')
    plt.xlabel('Coordenadas de sureste a noroeste')
    
    #save the figure
    savename = savepath+'/'+str(domain)+'/'+str(species)+'/vprof_diag2_'+str(domain)+'_'+str(species)+'_'+dates[tt]+'.png'
    fig_vprof2.savefig(savename)
    plt.close(fig_vprof2)

#write logfile
logfile=rundir+'/FLAG/map_'+str(domain)+'_'+str(startdate)+'_'+species+'.flag'
file = open(logfile,'w')
file.write('figures for '+domain+' and '+species+' have been generated')
file.close()

#plot map containing the lat lons of the two profiles
map_prof = plt.figure()
if domain in {'gal05r', 'gal0504r','gal3'}:
   mymap = Basemap(projection='cyl',llcrnrlon=minlon,llcrnrlat=minlat,urcrnrlon=maxlon,urcrnrlat=maxlat,resolution='c',lat_ts=43)#lat_ts=43
   X1, Y1 = mymap(lon_diag1,lat_diag1)
   X2, Y2 = mymap(lon_diag2,lat_diag2)
   plt.plot(X1,Y1)
   plt.plot(X2,Y2)
   mymap.readshapefile(scriptpath+'/shapes/municipios', 'municipios',linewidth=0.35, color='k', antialiased=1)          
   mymap.readshapefile(scriptpath+'/shapes/espana', 'espana',linewidth=0.35, color='k', antialiased=1)
   mymap.readshapefile(scriptpath+'/shapes/portugal', 'portugal',linewidth=0.35, color='k', antialiased=1)
elif domain in {'pib27','ib15r', 'ib16r', 'km12', 'ib1914r', 'km36', 'gal1511r', 'gal1511r2', 'gal16r', 'gal15r', 'gal15r2'}:
   mymap = Basemap(projection='cyl', resolution='l', llcrnrlon=np.min(lons),llcrnrlat=np.min(lats),urcrnrlon=np.max(lons),urcrnrlat=np.max(lats),lat_ts=44)
   X1, Y1 = mymap(lon_diag1,lat_diag1)
   X2, Y2 = mymap(lon_diag2,lat_diag2)
   plt.plot(X1,Y1)
   plt.plot(X2,Y2)
   mymap.drawcoastlines()
else:
   print("ATTENTION: Check entry for <domain>")

#save the figure
savename = savepath+'/'+str(domain)+'/vertical_profiles_'+str(domain)+'.png'
map_prof.savefig(savename)
plt.close(map_prof)

