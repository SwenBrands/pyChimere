# -*- coding: utf-8 -*-

#this scripts loads the CHIMERE output files and searches out the data at the grid-boxes nearest to those defined in <tarlat> and <tarlon>, the output is then save in nc format.

# import xarray as xr
# import numpy as np
# import pandas as pd
# import sys
# import dask
# import os
# import matplotlib
# matplotlib.use('Agg')
# from matplotlib import pyplot as plt
# from mpl_toolkits.basemap import Basemap
# from utiles import crea_cmap

import numpy as np
import os
import xarray as xr
from mpl_toolkits.basemap import Basemap
import pandas as pd
import os
import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import matplotlib.colors as colors
import Ngl
from utiles import crea_cmap

homedir = os.getenv("HOME")
lustre = os.getenv("LUSTRE")
threshold = 2 #0.0000000001
domain1 = 'GAL15R'
domain2 = 'GAL0504R'

##EXECUTE ##############################################################
dom1file = homedir+'/OP/LANZAR/chimere2017r4/domains/'+domain1+'/LANDUSE_USGS_'+domain1+'.nc'
dom2file = homedir+'/OP/LANZAR/chimere2017r4/domains/'+domain2+'/LANDUSE_USGS_'+domain2+'.nc'

nc = xr.open_dataset(dom1file)
lat1 = nc.lat.values
lon1 = nc.lon.values
#XX1, YY1 = np.meshgrid(lat1,lon1)
nc.close()

nc = xr.open_dataset(dom2file)
lat2 = nc.lat.values
lon2 = nc.lon.values
#XX2, YY2 = np.meshgrid(lat2,lon2)
nc.close()

data1 = np.zeros(lon1.shape)+1
data2 = np.zeros(lon2.shape)+2

#plot map
minlat = np.min(lat1)-5
minlon = np.min(lon1)-5
maxlat = np.max(lat1)+5
maxlon = np.max(lon1)+5
#get parallels and meridians
parallels = np.arange(np.ceil(minlat),np.floor(maxlat+1),5)
meridians = np.arange(np.ceil(minlon),np.floor(maxlon+1),5)

fig = plt.figure()
#cbounds = [0,0.5,1,1.5,2]
cbounds = [0,1,2]
mymap = Basemap(projection='cyl',llcrnrlon=-20,llcrnrlat=35,urcrnrlon=5,urcrnrlat=55,resolution='i',lat_ts=43)#lat_ts=43
mymap.drawlsmask(resolution='i',grid=1.25)
#XX1, YY1 = mymap(lon1,lat1)
#XX2, YY2 = mymap(lon2,lat2)
mymap.contourf(lon1, lat1, data1, cbounds, cmap='jet')
mymap.contourf(lon2, lat2, data2, cbounds, cmap='jet')
mymap.drawcoastlines(linewidth=0.5)
#ax.set_title(u'Concentracións de '+species+' no '+dates[tt])
mymap.drawparallels(parallels,labels=[True,False,False,True], color='grey', size=10)
mymap.drawmeridians(meridians,labels=[True,False,False,True], color='grey', size=10)
mymap.drawcountries()
#mymap.shadedrelief()
mymap.etopo()
#mymap.fillcontinents


#plt.colorbar()
fig.savefig('domains.pdf', dpi=300)
plt.close(fig)

    
