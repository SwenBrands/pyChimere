# -*- coding: utf-8 -*-

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
srcpath = lustre+'/SWEN/chimere2020r1/emiSURF2020/annual-EMEP01x01/preproc'
emisfolder = 'EMEP2018_reduced' #folder containing the emission files
makeyear = 2020 #year the emission files were made
validyear = 2018 #year for which the emission files are valid
resolution = 0.1 #resolution of the EMEP grid in degrees

## set input parameters ################################################
variables = ['CO','NH3','NMVOC','NOx','PM2_5','PM10','PMcoarse','SOx']

#facilities = ['Finsa','Meirama','Alcoa Coruna','As Pontes', 'Alcoa San-Ciprian']
#tarlat = [42.913333, 43.168406, 43.346918, 43.440292, 43.700505] #latitude of the location
#tarlon = [-8.500278, -8.409483, -8.438604, -7.864642, -7.466308] #longitude of the location
#activities = ['B_Industry','A_PublicPower','B_Industry','A_PublicPower','B_Industry'] #set activity abbreviations as set in the filenames

facilities = ['Meirama']
tarlat = [43.168406] #latitude of the location
tarlon = [-8.409483] #longitude of the location
activities = ['A_PublicPower'] #set activity abbreviations as set in the filenames

#multiply or set, multiply the existing emissions with the factors provided in <scalefactor> or simply set the emissions to this factor 
execmode = 'set'

##set scalefactor for each facility and emitted species
#scalefactor_CO = [100., 0.5, 0.5, 0.1, 0.2] #set estimated reduction factor for each location
#scalefactor_NH3 = [100., 0.5, 0.5, 0.1, 0.2] #set estimated reduction factor for each location
#scalefactor_NMVOC = [100., 0.5, 0.5, 0.1, 0.2] #set estimated reduction factor for each location
#scalefactor_NOx = [100., 0.5, 0.5, 0.1, 0.2] #set estimated reduction factor for each location
#scalefactor_PM2_5 = [100., 0.5, 0.5, 0.1, 0.2] #set estimated reduction factor for each location
#scalefactor_PM10 = [100., 0.5, 0.5, 0.1, 0.2] #set estimated reduction factor for each location
#scalefactor_PMcoarse = [100., 0.5, 0.5, 0.1, 0.2] #set estimated reduction factor for each location
#scalefactor_SOx = [100., 0.5, 0.5, 0.1, 0.2] #set estimated reduction factor for each location

scalefactor_CO = [3000.] #set estimated reduction factor for each location
scalefactor_NH3 = [0.5] #set estimated reduction factor for each location
scalefactor_NMVOC = [100.] #set estimated reduction factor for each location
scalefactor_NOx = [1000.] #set estimated reduction factor for each location
scalefactor_PM2_5 = [200.] #set estimated reduction factor for each location
scalefactor_PM10 = [200.] #set estimated reduction factor for each location
scalefactor_PMcoarse = [100.] #set estimated reduction factor for each location
scalefactor_SOx = [400.] #set estimated reduction factor for each location

startind = 5 #index rank of the first data row in the input txt files, may change for different versions of the EMEP inventory
tarcountry = 'ES'

#mapping options
plotmap = 'yes'
latlims = [41.5,44.] #[42.,44.]
lonlims = [-9.5,-6.5] #[-9.5,-7.]
colormap = 'hot_r'
shapefile_path = home+'/OP/LANZAR/chimere2017r4/figures/shapes'

## EXECUTE #############################################################
activities = np.array(activities)
tarlat = np.array(tarlat)
tarlon = np.array(tarlon)
facilities = np.array(facilities)
activities_unique = list(np.unique(activities))
scalefactor = np.stack((scalefactor_CO,scalefactor_NH3,scalefactor_NMVOC,scalefactor_NOx,scalefactor_PM2_5,scalefactor_PM10,scalefactor_PMcoarse,scalefactor_SOx))
#load the txt files
for vv in list(range(len(variables))):
    #copy the national totals file for this variable
    ifile_tot_path_orig = srcpath+'/'+emisfolder+'/'+variables[vv]+'_'+str(makeyear)+'_GRID_'+str(validyear)+'/'+variables[vv]+'_NT_'+str(makeyear)+'_GRID_'+str(validyear)+'.txt'
    ifile_tot_path = srcpath+'/'+emisfolder+'/'+variables[vv]+'_'+str(makeyear)+'_GRID_'+str(validyear)+'/modified_'+variables[vv]+'_NT_'+str(makeyear)+'_GRID_'+str(validyear)+'.txt'
    shutil.copy(ifile_tot_path_orig,ifile_tot_path)
    #Then modify the sector files and subtract the difference from the national emissions file
    for aa in list(range(len(activities_unique))):
        #load inventory file for this variable and activity
        print('Reduce '+variables[vv]+' for '+activities_unique[aa])
        print('A in activity file')
        ifile_path = srcpath+'/'+emisfolder+'/'+variables[vv]+'_'+str(makeyear)+'_GRID_'+str(validyear)+'/'+variables[vv]+'_'+activities_unique[aa]+'_'+str(makeyear)+'_GRID_'+str(validyear)+'.txt'
        ifile = open(ifile_path, "rb")
        reader = csv.reader(ifile, delimiter=';')
        #create a list containing the csv's rows 
        emisinv=[]
        for row in reader:
            #print row
            emisinv.append(row)
        ifile.close()
        del(row,reader,ifile,ifile_path)
        #get lats, lons, units and values from the inventory file
        emisinv_lon = np.array([np.float(emisinv[ii][4]) for ii in list(range(startind,len(emisinv)))])
        emisinv_lat = np.array([np.float(emisinv[ii][5]) for ii in list(range(startind,len(emisinv)))])
        emisinv_unit = [emisinv[ii][6] for ii in list(range(startind,len(emisinv)))]
        emisinv_value = np.array([np.float(emisinv[ii][7]) for ii in list(range(startind,len(emisinv)))])
        emisinv_country = np.array([emisinv[ii][0] for ii in list(range(startind,len(emisinv)))])
        
        #filter out which facilities pertain to the activity and filter relevant input variables
        facil_ind = np.where(activities == activities_unique[aa])[0]
        tarlat_step = tarlat[facil_ind]
        tarlon_step = tarlon[facil_ind]
        facilities_step = facilities[facil_ind]
        scalefactor_step = scalefactor[vv,facil_ind]
        
        origval_vec = np.zeros(len(tarlat_step))
        tardiff_vec = np.zeros(len(tarlat_step))
        tarval_vec = np.zeros(len(tarlat_step))
        mindist_id_vec = np.zeros(len(tarlat_step))
        
        for lc in list(range(len(tarlat_step))):
            #search nearest neighbour cell in the EMEP grid
            distances = [haversine(tarlon_step[lc],tarlat_step[lc],emisinv_lon[ii],emisinv_lat[ii]) for ii in list(range(len(emisinv_lon)))]
            distances = np.array(distances)
            mindist = np.min(distances)
            mindist_id = int(np.where((distances == np.min(distances)) & (emisinv_country == tarcountry))[0])
            #change emission flux for the identified nearest neighbour
            origval = np.float(emisinv[5+mindist_id][7])
            if execmode == 'multiply':
                tarval = np.float(emisinv[5+mindist_id][7])*scalefactor_step[lc]
            elif execmode == 'set':
                tarval = scalefactor_step[lc]
            else:
                raise Exception('ERROR: check entry for <execmode>!')
            tardiff = origval - tarval
        
            #print the output
            print('Location ' + str(lc)+':')
            print('Changing emission flux from '+str(origval)+' to '+str(tarval))
            print('at lon '+str(emisinv_lon[mindist_id])+' and lat '+str(emisinv_lat[mindist_id])+' for '+emisinv_country[mindist_id])
            print('the facility is located at lon '+str(tarlon_step[lc])+' and lat '+str(tarlat_step[lc]))
            emisinv[startind+mindist_id][7] = str(tarval)
            emisinv_value[mindist_id] = tarval
        
            #save the results for this station
            origval_vec[lc] = origval
            tardiff_vec[lc] = tardiff
            tarval_vec[lc] = tarval
            mindist_id_vec[lc] = mindist_id
        
        #then save the results in txt format
        outputfile = srcpath+'/'+emisfolder+'/'+variables[vv]+'_'+str(makeyear)+'_GRID_'+str(validyear)+'/modified_'+variables[vv]+'_'+activities[aa]+'_'+str(makeyear)+'_GRID_'+str(validyear)+'.txt'
        dataframe = pd.DataFrame(emisinv)
        dataframe.to_csv(outputfile,index=False,header=False,sep=';')
        del(outputfile,dataframe)

        ##SWITCH TO NATIONAL TOTALS ###############################################################################################################################################
        #load national totals file, substract reduction obtained from the previous step and save as csv
        print('----------------------------------------------------------------')
        print('B in national totals file')
        ifile_tot = open(ifile_tot_path, "rb")
        reader = csv.reader(ifile_tot, delimiter=';')
        #create a list containing the csv's rows 
        nattotals=[]
        for row in reader:
            #print row
            nattotals.append(row)
        ifile_tot.close()
        del(row,reader,ifile_tot)
        #get lats, lons, units and values from the national totals file
        nattotals_lon = np.array([np.float(nattotals[ii][4]) for ii in list(range(startind,len(nattotals)))])
        nattotals_lat = np.array([np.float(nattotals[ii][5]) for ii in list(range(startind,len(nattotals)))])
        nattotals_unit = [nattotals[ii][6] for ii in list(range(startind,len(nattotals)))]
        nattotals_value = np.array([np.float(nattotals[ii][7]) for ii in list(range(startind,len(nattotals)))])
        nattotals_country = np.array([nattotals[ii][0] for ii in list(range(startind,len(nattotals)))])

        origval_vec_tot = np.zeros(len(tarlat_step))
        tardiff_vec_tot = np.zeros(len(tarlat_step))
        tarval_vec_tot = np.zeros(len(tarlat_step))
        mindist_id_vec_tot = np.zeros(len(tarlat_step))
        
        for lc in list(range(len(tarlat_step))):
            #search nearest neighbour cell in the EMEP grid
            distances = [haversine(tarlon_step[lc],tarlat_step[lc],nattotals_lon[ii],nattotals_lat[ii]) for ii in list(range(len(nattotals_lon)))]
            distances = np.array(distances)
            mindist = np.min(distances)
            mindist_id = int(np.where((distances == np.min(distances)) & (nattotals_country == tarcountry))[0])
            #subtracting emission reduction at that grid-box
            origval = np.float(nattotals[startind+mindist_id][7])
            tarval = np.float(nattotals[startind+mindist_id][7])-tardiff_vec[lc]
        
            #print the output
            print('Location ' + str(lc)+':')
            print('Subtracting emission flux '+str(origval)+' - '+str(tardiff_vec[lc])+ ' = '+str(tarval))
            print('at lon '+str(nattotals_lon[mindist_id])+' and lat '+str(nattotals_lat[mindist_id])+' for '+nattotals_country[mindist_id])
            print('the facility is located at lon '+str(tarlon_step[lc])+' and lat '+str(tarlat_step[lc]))
            nattotals[startind+mindist_id][7] = str(tarval)
        
            #save the results for this stations
            origval_vec_tot[lc] = origval
            tarval_vec_tot[lc] = tarval
            mindist_id_vec_tot[lc] = mindist_id

        #then save the results in txt format
        dataframe = pd.DataFrame(nattotals)
        dataframe.to_csv(ifile_tot_path,index=False,header=False,sep=';')
        print('----------------------------------------------------------------')
        print('----------------------------------------------------------------')

        #optionally map the target concentrations
        if plotmap == 'yes':
            ##Then plot this grid (only once per variable and activity)
            getind = np.where((emisinv_lon >= lonlims[0]) & (emisinv_lon <= lonlims[1]) & (emisinv_lat >= latlims[0]) & (emisinv_lat <= latlims[1]) & (emisinv_country==tarcountry))[0]
            emisinv_lon_red = emisinv_lon[getind]
            emisinv_lat_red = emisinv_lat[getind]
            emisinv_value_red = emisinv_value[getind]
            XX,YY = np.meshgrid(np.unique(emisinv_lon_red),np.unique(emisinv_lat_red))
            #init regular emissions grid for this regions, look for matches and fill (values not found remain 0)
            ZZ = np.zeros(XX.shape)
            for ii in list(range(XX.shape[0])):
                for jj in list(range(XX.shape[1])):
                    tarind = np.where((emisinv_lon_red == XX[ii,jj]) & (emisinv_lat_red == YY[ii,jj]))[0]
                    if np.size(tarind) > 0:
                       ZZ[ii,jj] = emisinv_value_red[tarind]
                                
            fig1 = plt.figure()
            mymap = Basemap(projection='cyl', resolution='i', llcrnrlon=lonlims[0],llcrnrlat=latlims[0],urcrnrlon=lonlims[1],urcrnrlat=latlims[1])   
            #mymap.pcolormesh(X,Y, tarval, cmap=colormap, latlon=True, vmin=np.array(cbounds[0]), vmax=np.array(cbounds[-1]))
            mymap.pcolormesh(XX-resolution/2, YY-resolution/2, ZZ, cmap=colormap, latlon=True) #remove 0.05 degrees to center the gridbox at its corresponding coordinates
            mymap.readshapefile(shapefile_path+'/municipios', 'municipios',linewidth=0.5, color='k', antialiased=1)
            mymap.readshapefile(shapefile_path+'/espana', 'espana',linewidth=0.5, color='k', antialiased=1)
            mymap.readshapefile(shapefile_path+'/portugal', 'portugal',linewidth=0.5, color='k', antialiased=1)               
            ##plot locations of the emitting facilities
            plt.plot(tarlon,tarlat,'cX')
            plt.plot(emisinv_lon_red,emisinv_lat_red,'k.')
            #figure finetuning
            cbar = mymap.colorbar(shrink=0.8)
            plt.title('EMEP '+str(validyear)+', '+variables[vv]+', '+activities[aa])
            savename = srcpath+'/figs/EMEP_'+str(validyear)+'_'+variables[vv]+'_'+activities[aa]+'.png'
            fig1.savefig(savename, dpi=300)
            plt.close(fig1)
    
