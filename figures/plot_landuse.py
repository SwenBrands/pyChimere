# -*- coding: utf-8 -*-

#this scripts loads the CHIMERE output files and searches out the data at the grid-boxes nearest to those defined in <tarlat> and <tarlon>, the output is then save in nc format.

import xarray as xr
import numpy as np
import pandas as pd
import sys
import dask
import os
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap
from utiles import crea_cmap

homedir = os.getenv("HOME")
lustre = os.getenv("LUSTRE")
threshold = 2 #0.0000000001
domain = 'GAL3'
database = 'USGS' #USGS or GLOBCOVER

##EXECUTE ##############################################################
srcfile = homedir+'/OP/LANZAR/chimere2017r4/domains/'+domain+'/LANDUSE_'+database+'_'+domain+'.nc'
under = '#8B4513'
#over = '#fe00c8'
over = '#3200fe'

rgbs = [
 '#8B4513', #agricultural land
 '#7CFC00', #grassland
 '#FFE4E1', #barren land
 '#00BFFF', #deepskyblue, water inland
 '#2c2c2c', #grey, urban
 '#f37736', #orange, shrubs
 '#228B22', #green, needle
 '#808000', #green, broad
 '#3200fe'] #blue, ocean
nc = xr.open_dataset(srcfile)
lat = nc.lat.values
lon = nc.lon.values

agri = nc.Agricultural_land_crops.values
grass = nc.Grassland.values
barren = nc.Barren_land_bare_ground.values
water = nc.Inland_Water.values
urban = nc.Urban.values
urban_mask = np.where(urban>threshold)
urban[urban_mask] = np.array(2)
shrubs = nc.Shrubs.values
needleaf = nc.Needleaf_forest.values
broadleaf = nc.Broadleaf_forest.values
ocean = nc.Ocean.values
nc.close()
 
all_classes = np.stack((agri,grass,barren,water,urban,shrubs,needleaf,broadleaf,ocean))
dominant_class = all_classes.argmax(axis=0)+0.5

#PLOT MAP
cbounds = np.arange(np.max(dominant_class)+1)
cbounds = cbounds.tolist()
colormap = crea_cmap(cbounds, rgbs, under, over)
minlat = np.min(lat)
minlon = np.min(lon)
maxlat = np.max(lat)
maxlon = np.max(lon)
#get parallels and meridians
parallels = np.arange(np.ceil(minlat),np.floor(maxlat+1))
meridians = np.arange(np.ceil(minlon),np.floor(maxlon+1))

fig = plt.figure()
mymap = Basemap(projection='cea',llcrnrlon=minlon,llcrnrlat=minlat,urcrnrlon=maxlon,urcrnrlat=maxlat,resolution='i',lat_ts=43)#lat_ts=43
X, Y = mymap(lon,lat)
mymap.pcolormesh(X, Y, dominant_class, cmap=colormap, vmin=cbounds[0], vmax=cbounds[-1])
mymap.drawcoastlines()
#ax.set_title(u'Concentracións de '+species+' no '+dates[tt])
mymap.drawparallels(parallels,labels=[True,False,False,True], color='None', size=10)
mymap.drawmeridians(meridians,labels=[True,False,False,True], color='None', size=10)
plt.colorbar()
fig.savefig('landuse_'+database+'_'+domain+'.pdf', dpi=300)
plt.close(fig)

    
