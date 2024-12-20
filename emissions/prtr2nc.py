# -*- coding: utf-8 -*-
'''assigns latitude-longitude coordinates and other metadata to the prtr inventory for
Galicia and saves the dataset in netCDF format'''

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
srcpath = lustre+'/OP/DATOS/CICC/chimere2017r4/EMISSIONS/prtr_galicia'
emis_file = 'emisiones_galicia_2014_2019.csv' #folder containing the emission files
topo_file = 'prtr_calidadeAire_2021_v2.csv' #prtr_calidadeAire_2021_v2.csv, INSTALACIONES_PRTR_2020.csv folder containing the lat lon coordinates 

## set input parameters ################################################
species_prtr = ['Monoxido de carbono (CO)','Amoniaco (NH3)','Compuestos organicos volatiles distintos del metano (COVNM)','Oxidos de nitrogeno (NOx/NO2)','notused','Particulas (PM10)','notused','Oxidos de azufre (SOx/SO2)']
species_emep = ['CO','NH3','NMVOC','NOx','PM2_5','PM10','PMcoarse','SOx']
taryear = 'Validado 18' #as defined in <emis_file>, 'Validado 18' works
scale_pm25 = 0.4
scale_pmcoarse = 0.6

##headers of the input files, currently not used
variables_emissions = ['PPAL','Centro_id','Centro','Conteminante','Medio','Validado 14','Validado 15','Validado 16','Validado 17','Validado 18','Validado 19']
variables_topo = ['NATIONALID','Codigo','Nome','Concello','Latitude','Lonxitude','Altura','CNAE']
#variables_topo = ['NationalId','CodigoCentroCA','Nombre','ActividadEconomica','Id_Municipio','Id_Poblacion','Direccion','NumeroVia','CodPostal','Latitud','Longitud','Altitud','CNAECode','Id_Provincia','TITULAR']

## EXECUTE #############################################################
activity_emep = ['A_PublicPower','B_Industry','C_OtherStationaryComb','D_Fugitive','E_Solvents','F_RoadTransport','G_Shipping','H_Aviation','I_Offroad', \
                'J_Waste', 'K_AgriLivestock', 'L_AgriOther', 'M_Other']
activity_prtr_number = range(1,10)
activity_prtr_description = ['Sector de la energía','Produccion y transformacion de metales','Industria mineral', \
                            'Industria quimica','Gestion de residuos y aguas residuales', 'Fabricacion y transformación de papel y madera', \
                            'Ganaderia y acuicultura intensiva', 'Productos de origen animal y vegetal de la industria alimentaria y de las bebidas', \
                            'Otras actividades']
activity_emep_equival = ['A_PublicPower','B_Industry','B_Industry','B_Industry','J_Waste','B_Industry','K_AgriLivestock','B_Industry','M_Other']

execfile(runpath+'/my_readcsv.py')
emis_path = srcpath+'/'+emis_file
topo_path = srcpath+'/'+topo_file
output_path = srcpath+'/prtr_emissions_20'+taryear[-2:]+'.nc'

## load emission file and coordinates
emissions = my_readcsv(emis_path,';')
topo = my_readcsv(topo_path,';')

emissions_metadata = emissions[0]
topo_metadata = topo[0]

#get index of the column pertaining to 'PPAL' (the activity sector) in the emissions file and gather index of all point emissions in a list
ind_id = int(np.where(np.array(emissions_metadata)=='PPAL')[0])
sector_prtr = [emissions[ii][ind_id] for ii in range(1,len(emissions))]
sector_prtr = np.array(sector_prtr)
sector_prtr_short = np.array([sector_prtr[ii][0] for ii in range(len(sector_prtr))])
sector_prtr_short_unique = np.unique(sector_prtr_short)
#get equivalent emep sector for each prtr sector entry
sector_emep = np.array([None] * len(sector_prtr_short))
sector_prtr_long = np.array([None] * len(sector_prtr_short))
for se in range(len(sector_prtr_short_unique)):
    secind = np.where(sector_prtr_short == sector_prtr_short_unique[se])[0]
    sector_emep[secind] = activity_emep_equival[se]
    sector_prtr_long[secind] = activity_prtr_description[se]

#get index of the column pertaining to 'Conteminante' in the emissions file and gather index of all point emissions in a list
ind_id = int(np.where(np.array(emissions_metadata)=='Conteminante')[0])
species = [emissions[ii][ind_id] for ii in range(1,len(emissions))]
species = np.array(species)
#get index of the column pertaining to 'Centro_id' in the emissions file and gather index of all point emissions in a list
ind_id = int(np.where(np.array(emissions_metadata)=='Centro_id')[0])
id_emis = [emissions[ii][ind_id] for ii in range(1,len(emissions))]
id_emis = np.array(id_emis)
#get index of the column pertaining to the target year of the emission values defined in <taryear>
ind_id = int(np.where(np.array(emissions_metadata)==taryear)[0])
flux_emis = [emissions[ii][ind_id] for ii in range(1,len(emissions))]
flux_emis = np.array(flux_emis)
nanind = np.where(flux_emis == '')[0]
flux_emis[nanind] = np.nan

#get index of the column pertaining to 'CodigoCentroCA' in the topographic infos file and gather index of all topographic instances in a list
ind_id = int(np.where(np.array(topo_metadata)=='Codigo')[0])
id_topo = [topo[ii][ind_id] for ii in range(1,len(topo))]
id_topo = np.array(id_topo)
#get index of the column pertaining to 'Latitud' in the topographic infos file and gather index of all topographic instances in a list
ind_id = int(np.where(np.array(topo_metadata)=='Latitude')[0])
lat = [topo[ii][ind_id] for ii in range(1,len(topo))]
lat = np.array(lat)
#get index of the column pertaining to 'Longitud' in the topographic infos file and gather index of all topographic instances in a list
ind_id = int(np.where(np.array(topo_metadata)=='Lonxitude')[0])
lon = [topo[ii][ind_id] for ii in range(1,len(topo))]
lon = np.array(lon)
#get index of the column pertaining to 'Lonxitude' in the topographic infos file and gather index of all topographic instances in a list
ind_id = int(np.where(np.array(topo_metadata)=='Altura')[0])
altitude = [topo[ii][ind_id] for ii in range(1,len(topo))]
altitude = np.array(altitude)

#check whether the point emission and topo files have the same unique IDs, then assign lat and lon to each point emission using the ID
id_topo_unique = np.array(np.unique(id_topo))
id_emis_unique = np.array(np.unique(id_emis))
#check whether the topographic instances provided by Augustin are unique
if np.sum(id_topo_unique.astype('float')) == np.sum(id_topo.astype('float')):
    print('info: all topographic instances are unique, so we can proceed....')
else:
    raise Exception('ERROR: chech <id_topo> for repeated entries!')

lat_all = np.zeros(len(id_emis))
lon_all = np.zeros(len(id_emis))
altitude_all = np.zeros(len(id_emis))
common_ind = np.zeros(len(id_emis))
for ee in range(len(id_emis)):
    ind_common_id = np.where(id_topo == id_emis[ee])[0]
    #add exception for pointwise emissions without coordinates
    if ind_common_id.size == 0:
        print('warning: no coordinates were found for point emission with ID '+str(id_emis[ee])+'! Assinging NaN in this case....') 
        lat_all[ee] = np.nan
        lon_all[ee] = np.nan
        altitude_all[ee] = np.nan
        common_ind[ee] = np.nan
    else:
        lat_all[ee] = lat[ind_common_id]
        lon_all[ee] = lon[ind_common_id]
        altitude_all[ee] = altitude[ind_common_id]
        common_ind[ee] = ind_common_id

species_unique = np.unique(species)
emis_matrix = np.zeros((len(id_emis),len(species_emep)))
for vv in range(len(species_emep)):
    if species_emep[vv] == 'PM2_5':
        tarvar = 'Particulas (PM10)'
        scalefactor = scale_pm25
        print('info: '+species_emep[vv]+' does not exist in PRTR inventory, assuming '+tarvar+' values multiplied by '+str(scalefactor))
    elif species_emep[vv] == 'PMcoarse':
        tarvar = 'Particulas (PM10)'
        scalefactor = scale_pmcoarse
        print('info: '+species_emep[vv]+' does not exist in PRTR inventory, assuming '+tarvar+' values multiplied by '+str(scalefactor))
    else:
        tarvar = species_prtr[vv]
        scalefactor = 1.
        print('info: '+species_emep[vv]+' does exist in PRTR inventory...')
    varind = np.where(species == tarvar)[0]
    emis_matrix[varind,vv] = flux_emis.astype('float')[varind]*scalefactor

#save as xarray dataset
ds = xr.Dataset(data_vars=dict(emission=(["x", "y"], emis_matrix),lat=(["x"], lat_all),lon=(["x"], lon_all),altitude=(["x"], altitude_all), \
    sector_prtr=(["x"], sector_prtr), sector_prtr_long=(["x"], sector_prtr_long), sector_emep=(["x"], sector_emep)),coords=dict(species=(["y"], species_emep), id=(["x"], id_emis)), \
    attrs=dict(description="Pointwise emission data from the Galician PRTR report, provided by Agustin Diaz Alonso (Laboratorio de Medio Ambiente de Galicia) and re-organized in netCDF format by Swen Brands (MeteoGalicia)",units='kg/year',year=taryear,contact='swen.brands@gmail.com',scalefactor_pm25=str(scale_pm25),scalefactor_pmcoarse=str(scale_pmcoarse)))
ds.to_netcdf(output_path)
ds.close()

