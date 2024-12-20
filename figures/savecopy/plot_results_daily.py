#load python packages
from netCDF4 import Dataset
import numpy as np
#import xarray as xr
#import cartopy.crs as crs
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap

srcpath = '/mnt/lustre/scratch/home/cesga/orballo/OP/PRED/chimere2016a'
savepath = '/mnt/lustre/scratch/home/cesga/orballo/OP/PRED/chimere2016a/figs' 
days = [20170413, 20170414, 20170415, 20170416, 20170417, 20170418, 20170419, 20170420, 20170421, 20170422, 20170423, 20170424]

for dd in np.arange(len(days)):
    filename = srcpath+'/out.'+str(days[dd])+'_00_km36.nc'

    ## EXECUTE #############################################################
    #dataset = xr.open_dataset(filename)

    nc = Dataset(filename, mode='r')
    lons = nc.variables['lon'][:]
    lats = nc.variables['lat'][:]
    pres = nc.variables['pres'][:]
    pm10 = nc.variables['PM10'][:]
    pm10_units = nc.variables['PM10'].units
    dates = nc.variables['Times']
    nc.close
    #v = np.linspace(-.1, 20, 200, endpoint=True)

    for tt in np.arange(pm10.shape[0]):
        for hh in np.arange(pm10.shape[1]):
            fig1 = plt.figure()
            map = Basemap(llcrnrlon=np.min(lons),llcrnrlat=np.min(lats),urcrnrlon=np.max(lons),urcrnrlat=np.max(lats))
            map.contourf(lons, lats, pm10[tt,hh,],vmin=0, vmax=1000)
            map.drawcoastlines()
            #plt.clim(0,1000)
            plt.colorbar()
            savename=savepath+'/lvl'+str(hh)+'/pm10_'+str(days[dd])+'_'+str(tt)
            fig1.savefig(savename, dpi=150)
            plt.close(fig1)

