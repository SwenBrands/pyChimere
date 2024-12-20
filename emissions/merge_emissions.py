# -*- coding: utf-8 -*-
'''reads csv files from the EMEP 01x01 dataset, compares the emissions therien with those obtained from the Galician
PRTR inventory (previously obtaiend with prtr2nc.py) and saves missing emissions in new csv files. This
is done for both the sector-species specific csv files, as well was for the national totals files. '''

#load python packages
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
from math import radians, cos, sin, asin, sqrt
import csv
import shutil

## set input paths #####################################################
lustre = os.getenv('LUSTRE')
home = os.getenv('HOME')
runpath = home+'/OP/LANZAR/chimere2017r4/emissions'
execfile(runpath+'/haversine.py')
#srcpath = lustre+'/SWEN/chimere2020r1/emiSURF2020/annual-EMEP01x01/preproc'
srcpath = lustre+'/SWEN/chimere2020r1/emisurf2020r3/annual-EMEP01x01/preproc'
prtrpath = lustre+'/OP/DATOS/CICC/chimere2017r4/EMISSIONS/prtr_galicia'
emisfolder = 'EMEP2019_merged' #folder containing the emission files
makeyear = 2021 #year the emission files were made
validyear = 2019 #year for which the emission files are valid
resolution = 0.1 #resolution of the EMEP grid in degrees

## set input parameters ################################################
variables = ['CO','NH3','NMVOC','NOx','PM2_5','PM10','PMcoarse','SOx']
#variables = ['NH3']
#activities = np.array(['B_Industry','A_PublicPower','B_Industry','A_PublicPower','B_Industry']) #set activity abbreviations as set in the filenames
#activities = np.array(['B_Industry','A_PublicPower']) #set activity abbreviations as set in the filenames
startind = 5 #index rank of the first data row in the input txt files, may change for different versions of the EMEP inventory
tarcountry = 'ES'

#mapping options
plotmap = 'yes'
latlims = [41.45,44.05] #must be coordinates of the emep grid!!
lonlims = [-9.45,-6.45] #must be coordinates of the coordinates!!
colormap = 'hot_r'
shapefile_path = home+'/OP/LANZAR/chimere2017r4/figures/shapes'

## EXECUTE #############################################################
print('Megagrams (Mg) /year are assumed for EMEP and kg/year for PRTR!')
os.chdir(runpath)
execfile('my_readcsv.py')

#load the pointwise PRTR emissions for <validyear>; the netCDF source file must have been previously generated with prtr2nc.py
prtrfile = prtrpath+'/prtr_emissions_'+str(validyear)+'.nc'
nc = xr.open_dataset(prtrfile)
lat_prtr = nc.lat.values
lon_prtr = nc.lon.values
sector_prtr = nc.sector_emep.values

#get all unique emep activities from the prtr nc file obtained with prtr2nc.py
activities = np.unique(nc.sector_emep.values)

#generate an empty lat-lon mesh to be used to fill the national totals file below
addmesh_lons = np.linspace(lonlims[0],lonlims[-1],np.round(np.subtract(lonlims[-1],lonlims[0]),2)/resolution+1)
addmesh_lats = np.linspace(latlims[0],latlims[-1],np.round(np.subtract(latlims[-1],latlims[0]),2)/resolution+1)
XX_red,YY_red = np.meshgrid(addmesh_lons,addmesh_lats)
XX_red = np.round(XX_red,2)
YY_red = np.round(YY_red,2)

#load the txt files
for vv in list(range(len(variables))):
    #generate empty mesh for each variable; used for the NT file below
    addmesh = np.zeros(XX_red.shape)
    #copy the national totals file for this variable
    ifile_tot_path_orig = srcpath+'/'+emisfolder+'/'+variables[vv]+'_'+str(makeyear)+'_GRID_'+str(validyear)+'/'+variables[vv]+'_NT_'+str(makeyear)+'_GRID_'+str(validyear)+'.txt'
    ifile_tot_path = srcpath+'/'+emisfolder+'/'+variables[vv]+'_'+str(makeyear)+'_GRID_'+str(validyear)+'/modified_'+variables[vv]+'_NT_'+str(makeyear)+'_GRID_'+str(validyear)+'.txt'
    shutil.copy(ifile_tot_path_orig,ifile_tot_path)
    #empty dictionares were the sector specific information is stored for later use in the national totals file
    diff_tont = {}
    lat_tont = {}
    lon_tont = {}
    #Then modify the sector files and subtract the difference from the national emissions file
    for aa in list(range(len(activities))):
        #load inventory file for this variable and activity
        print('Merge '+variables[vv]+' for '+activities[aa])
        print('A in activity file')
        ifile_path = srcpath+'/'+emisfolder+'/'+variables[vv]+'_'+str(makeyear)+'_GRID_'+str(validyear)+'/'+variables[vv]+'_'+activities[aa]+'_'+str(makeyear)+'_GRID_'+str(validyear)+'.txt'
        ifile = open(ifile_path, "rb")
        emisinv = my_readcsv(ifile_path,';')
        
        #get lats, lons, boundaries, units and values from the inventory file
        emisinv_lon = np.array([np.float(emisinv[ii][4]) for ii in list(range(startind,len(emisinv)))])
        emisinv_lon_elim = emisinv_lon+0.05
        emisinv_lon_wlim = emisinv_lon-0.05
        emisinv_lat = np.array([np.float(emisinv[ii][5]) for ii in list(range(startind,len(emisinv)))])
        emisinv_lat_nlim = emisinv_lat+0.05
        emisinv_lat_slim = emisinv_lat-0.05        
        emisinv_unit = [emisinv[ii][6] for ii in list(range(startind,len(emisinv)))]
        #emisinv_value = np.array([np.float(emisinv[ii][7]) for ii in list(range(startind,len(emisinv)))])
        emisinv_value = np.array([np.float(emisinv[ii][7]) for ii in list(range(startind,len(emisinv)))])/1000000 #convert from Mg to tons
        emisinv_country = np.array([emisinv[ii][0] for ii in list(range(startind,len(emisinv)))])
        emisinv_id = np.array(range(len(emisinv_value))) 

        #cut out a subdomain for the specific target region defined by <latlims> and <lonlims>
        getind = np.where((emisinv_lon >= lonlims[0]) & (emisinv_lon <= lonlims[1]) & (emisinv_lat >= latlims[0]) & (emisinv_lat <= latlims[1]) & (emisinv_country==tarcountry))[0]
        print(getind)
        emisinv_id_red = emisinv_id[getind]
        emisinv_lon_red = emisinv_lon[getind]
        emisinv_lon_elim_red = emisinv_lon_elim[getind]
        emisinv_lon_wlim_red = emisinv_lon_wlim[getind]
        emisinv_lat_nlim_red = emisinv_lat_nlim[getind]
        emisinv_lat_slim_red = emisinv_lat_slim[getind]
        emisinv_lat_red = emisinv_lat[getind]
        emisinv_value_red = emisinv_value[getind]
        #emisinv_unit_red = emisinv_unit[getind]
        #emisinv_country_red = emisinv_country[getind]
        
        #cut out prtr emissions for this variables
        emis_prtr = nc.emission.values[:,nc.species==variables[vv]]/1000 #all prtr point sources for this variable, devide by 1000 to obtain tons
        emis_prtr_sum = np.zeros(len(emisinv_id_red)) #prtr point sources summed up at each grid box of the reduced EMEP grid
        diff_val = list()
        diff_origval = list()
        diff_newval = list()
        diff_lat = list()
        diff_lon = list()
        diff_id = list()
        for pp in list(range(len(emisinv_id_red))):
            getind = np.where((lon_prtr >= emisinv_lon_wlim_red[pp]) & (lon_prtr < emisinv_lon_elim_red[pp]) \
             & (lat_prtr >= emisinv_lat_slim_red[pp]) & (lat_prtr < emisinv_lat_nlim_red[pp]) & (sector_prtr == activities[aa]))
            emis_prtr_sum[pp] = np.nansum(emis_prtr[getind])
            diff = emisinv_value_red[pp]-emis_prtr_sum[pp]
            if diff < 0:
                idind = int(np.where(emisinv_id == emisinv_id_red[pp])[0])
                emisinv[startind+idind][7] = str(emis_prtr_sum[pp])
                emisinv_value[idind] = emis_prtr_sum[pp]
                diff_val.append(diff)
                diff_origval.append(emisinv_value_red[pp])
                diff_newval.append(emis_prtr_sum[pp])
                diff_lat.append(emisinv_lat[idind])
                diff_lon.append(emisinv_lon[idind])
                diff_id.append(emisinv_id[idind])
                print('PRTR emissions exceed EMEP emissions at the gridbox centered at lon '+str(emisinv_lon_red[pp])+' and lat '+str(emisinv_lat_red[pp])+':')
                print('lons of PRTR point sources are:')
                print(lon_prtr[getind])
                print('lats of PRTR point sources are:')
                print(lat_prtr[getind])
                print('EMEP value is '+str(emisinv_value_red[pp]))
                print('PRTR sum is = '+str(emis_prtr_sum[pp]))
                print('Difference ='+str(emisinv_value_red[pp]-emis_prtr_sum[pp]))
                print('EMEP values are substituted by PRTR sum at this gridbox!')
                print('--------------------------------------------------')
        
        #then save the results in txt format
        outputfile = srcpath+'/'+emisfolder+'/'+variables[vv]+'_'+str(makeyear)+'_GRID_'+str(validyear)+'/modified_'+variables[vv]+'_'+activities[aa]+'_'+str(makeyear)+'_GRID_'+str(validyear)+'.txt'
        dataframe = pd.DataFrame(emisinv)
        dataframe.to_csv(outputfile,index=False,header=False,sep=';')
        del(outputfile,dataframe)
        #retain information for this activity sector for modifying the national totals file later
        diff_tont[activities[aa]] = diff_val
        lat_tont[activities[aa]] = diff_lat
        lon_tont[activities[aa]] = diff_lon
        
    ##return all values and coordinates of the dictonaries, np.sum is here used to stack the differences, not for summing !!
    #diff_all = np.array(np.sum(diff_tont.values()))*-1.0
    stacked_diff = np.array(np.sum(diff_tont.values()))
    #check if all values are negative
    if np.any(np.array(stacked_diff)>0):
        raise Exception('ERROR: positive values are present in stacked_diff!')
    diff_all = stacked_diff*-1.0 #switch the sign
    lat_all = np.round(np.array(np.sum(lat_tont.values())),2)
    lon_all = np.round(np.array(np.sum(lon_tont.values())),2)
    
    #assign the differences to the grid boxes of the empty mesh and sum them up
    for ll in range(len(diff_all)):
        coordind = np.where((XX_red == lon_all[ll]) & (YY_red == lat_all[ll]))
        addmesh[coordind] = addmesh[coordind]+diff_all[ll]

    #Then flatten the arrays and retain only the non-zero scalars
    addmesh_flat = addmesh.flatten()
    lon_flat = XX_red.flatten()
    lat_flat= YY_red.flatten()    
    retainind = np.where(addmesh_flat != 0)
    addmesh_flat = addmesh_flat[retainind]
    lon_flat = lon_flat[retainind]
    lat_flat = lat_flat[retainind]

    ##SWITCH TO NATIONAL TOTALS ###############################################################################################################################################
    #load national totals file, substract reduction obtained from the previous step and save as csv
    print('----------------------------------------------------------------')
    print('B in national totals file')
    nattotals = my_readcsv(ifile_tot_path,';')

    #get lats, lons, units and values from the national totals file
    nattotals_lon = np.array([np.float(nattotals[ii][4]) for ii in list(range(startind,len(nattotals)))])
    nattotals_lat = np.array([np.float(nattotals[ii][5]) for ii in list(range(startind,len(nattotals)))])
    nattotals_unit = [nattotals[ii][6] for ii in list(range(startind,len(nattotals)))]
    nattotals_value = np.array([np.float(nattotals[ii][7]) for ii in list(range(startind,len(nattotals)))])
    nattotals_country = np.array([nattotals[ii][0] for ii in list(range(startind,len(nattotals)))])
    
    #bis hierhin

    for lc in list(range(len(addmesh_flat))):
        coordind = np.where((nattotals_lat == lat_flat[lc]) & (nattotals_lon == lon_flat[lc]) & (nattotals_country==tarcountry))
        if np.size(coordind) == 0:
            raise Exception('ERROR: check entry for <coordind>')
        
        #search nearest neighbour cell in the EMEP grid
        distances = [haversine(lon_flat[lc],lat_flat[lc],nattotals_lon[ii],nattotals_lat[ii]) for ii in list(range(len(nattotals_lon)))]
        distances = np.array(distances)
        mindist = np.min(distances)
        mindist_id = int(np.where((distances == np.min(distances)) & (nattotals_country == tarcountry))[0])
        #adding prtr emissions to emep emissions at that grid-box
        origval = np.float(nattotals[startind+mindist_id][7])
        tarval = np.float(nattotals[startind+mindist_id][7])+addmesh_flat[lc]
        
        #print the output
        print('Variable ' + variables[vv])
        print('Grid box ' + str(lc)+':')
        print('Adding emission flux '+str(origval)+' + '+str(addmesh_flat[lc])+ ' = '+str(tarval))
        print('at lon '+str(nattotals_lon[mindist_id])+' and lat '+str(nattotals_lat[mindist_id])+' for '+nattotals_country[mindist_id])
        print('the grid box is located at lon '+str(lon_flat[lc])+' and lat '+str(lat_flat[lc]))
        nattotals[startind+mindist_id][7] = str(tarval)
        
        ##save the results for this stations
        #origval_vec_tot[lc] = origval
        #tarval_vec_tot[lc] = tarval
        #mindist_id_vec_tot[lc] = mindist_id

    #then save the results in txt format
    dataframe = pd.DataFrame(nattotals)
    dataframe.to_csv(ifile_tot_path,index=False,header=False,sep=';')
    print('----------------------------------------------------------------')
    print('----------------------------------------------------------------')

    #optionally map the target concentrations, also save the EMEP emission values and latlon coordinates over Galicia for this activity sector
    if plotmap == 'yes':
        fig1 = plt.figure()
        mymap = Basemap(projection='cyl', resolution='i', llcrnrlon=lonlims[0],llcrnrlat=latlims[0],urcrnrlon=lonlims[1],urcrnrlat=latlims[1])   
        #mymap.pcolormesh(X,Y, tarval, cmap=colormap, latlon=True, vmin=np.array(cbounds[0]), vmax=np.array(cbounds[-1]))
        mymap.pcolormesh(XX_red-resolution/2, YY_red-resolution/2, addmesh, cmap=colormap, latlon=True) #remove 0.05 degrees to center the gridbox at its corresponding coordinates
        mymap.readshapefile(shapefile_path+'/municipios', 'municipios',linewidth=0.5, color='k', antialiased=1)
        mymap.readshapefile(shapefile_path+'/espana', 'espana',linewidth=0.5, color='k', antialiased=1)
        mymap.readshapefile(shapefile_path+'/portugal', 'portugal',linewidth=0.5, color='k', antialiased=1)               
        #figure finetuning
        cbar = mymap.colorbar(shrink=0.8)
        cbar.ax.set_ylabel('Toneladas por ano', rotation=90)
        plt.title('Added to EMEP NT for '+str(validyear)+', '+variables[vv])
        savename = srcpath+'/figs/add2EMEP_NT_'+str(validyear)+'_'+variables[vv]+'.png'
        fig1.savefig(savename, dpi=300)
        plt.close(fig1)

nc.close()
