# -*- coding: utf-8 -*-


"""this script drops variables and dimensions from the EMEP01x01 emission files generetated by emiSURF2020
that are not in use with chimere2017r4 and save the modified files with the same name, but in netcdf3 classic
format, which can be read by chimere2017r4. In the future, please test whether the model also runs with netcdf4
classic to save disk space!"""

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
domain = 'GAL15R'
outputformat = 'NETCDF3_CLASSIC' #output format as defined by <format> parameter of nc.to_netcdf, see http://xarray.pydata.org/en/stable/generated/xarray.Dataset.to_netcdf.html
dropvariables = 'no'

#set environmental variables
lustre = os.getenv("LUSTRE")
home = os.getenv("HOME")
mechanism = 'melchior'

#set path to the landuse file
rundir = home+'/OP/LANZAR/chimere2017r4/emissions'
emisdir_src = lustre+'/OP/DATOS/CICC/chimere2017r4/EMISSIONS/EMEP19_'+domain+'_pop1r1ag1sh0lsp0prtr1_globcover_'+mechanism+'/EMI_2019'
#emisdir_tar = lustre+'/OP/DATOS/CICC/chimere2017r4/EMISSIONS/myEMEP01x01_'+domain+'_pop1r1ag1sh0lsp0_globcover_'+mechanism+'/EMI_2018_v2017' # in this case, variables were dropped
emisdir_tar = lustre+'/OP/DATOS/CICC/chimere2017r4/EMISSIONS/EMEP19_'+domain+'_pop1r1ag1sh0lsp0prtr1_globcover_'+mechanism+'/EMI_2019_nc3' # in this case, variables were NOT dropped


##EXECUTE #############################################################
os.chdir(rundir)
print('INFO: the domain to be processed is '+domain+'...')
print('INFO: the chemical mechanism to be processed is '+mechanism+'...')
print('INFO: the format of the output netCDF files is '+outputformat+'...')
print('INFO: are variables dropped? '+dropvariables+'...')

#define species for the respective chemical mechanism, valid for emisurf2020r3
if mechanism == 'melchior':
    species = ['ALK4','ALK5','APINEN','ARO1','ARO2','BCAR_coa','BCAR_fin','C2H4','C2H5OH','C2H6', \
                'C3H6','C5H8','CH3CHO','CH3COE','CH3OH','CO','H2SO4_fin','HCHO','HONO','NC4H10', \
                'NH3', 'NO2', 'NO', 'OCAR_coa', 'OCAR_fin', 'OLE1', 'OLE2', 'OXYL', 'PPM_coa', 'PPM_fin', \
                'SO2', 'TMB', 'TOL']
elif mechanism == 'saprc':
    species = ['AACD','ACET','ACYE','ALK1','ALK2','ALK3','ALK4','ALK5', \
                'APINEN','ARO1','ARO2','BALD','BaP_fin','BbF_fin','BCAR_fin','BENZ','BkF_fin','C2H4', 'C5H8', \
                'CCHO', 'CH4', 'CO', 'CRES', 'GLY', 'H2SO4_fin', 'HCHO', 'HONO', \
                'IPRD', 'MACR', 'MEOH', 'NH3', 'NO2', 'N0', 'OCAR_fin', 'OLE1', 'OLE2', \
                'PPM_big', 'PPM_coa', 'PPM_fin', 'PRD2', 'RCHO', 'SO2', 'TMB', 'TOL']
else:
    raise Exception('ERROR: check entry for <mechanism>!')

tarmonths = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']

#load emission files for the old domain
for mm in list(range(len(tarmonths))):
    print('INFO: target month is '+tarmonths[mm]+' ######################################')
    for ii in list(range(len(species))):
        print('INFO: target species is '+species[ii])
        emisfile_src = emisdir_src+'/EMIS.'+domain+'.'+tarmonths[mm]+'.'+species[ii]+'.s.nc'        
        nc = xr.open_dataset(emisfile_src)
        
        #optionally drop variables
        if dropvariables == 'yes':
            print('INFO: SW_DATES, SW_ZONES, SW_SHIFT, TZ_PER_GRID and SW_SHIFT are dropped!')
            ncnew = nc.drop(['SW_DATES','SW_ZONES','SW_SHIFT','TZ_PER_GRID','SW_SHIFT'])
        elif dropvariables == 'no':
            print('INFO: no variables are dropped...')
            ncnew = nc
        else:
            raise Exception('ERROR: checke entry for <dropvariables>!')

        #then save the new file and close both source and new file
        emisfile_tar= emisdir_tar+'/EMIS.'+domain+'.'+tarmonths[mm]+'.'+species[ii]+'.s.nc'
        ncnew.to_netcdf(path=emisfile_tar,format=outputformat)
        nc.close()
        ncnew.close()

print('INFO: cut_emissions.py has been run successfully!')
