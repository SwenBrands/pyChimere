# -*- coding: utf-8 -*-

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
domain_old = 'GAL0504R'
domain_new = 'GAL0504S'

#set environmental variables
lustre = os.getenv("LUSTRE")
home = os.getenv("HOME")
mechanism = 'melchior'

#set path to the landuse file
rundir = home+'/OP/LANZAR/chimere2017r4/emissions'
lufile = home+'/OP/LANZAR/chimere2017r4/domains/'+domain_new+'/LANDUSE_GLOBCOVER_'+domain_new+'.nc'
emisdir_src = lustre+'/OP/DATOS/CICC/chimere2017r4/EMISSIONS/emep01x01_2017/'+mechanism+'/globcover/'+domain_old+'/pop1r1ag0sh0lsp1'
emisdir_tar = lustre+'/OP/DATOS/CICC/chimere2017r4/EMISSIONS/emep01x01_2017/'+mechanism+'/globcover/'+domain_new+'/pop1r1ag0sh0lsp1'

#species = ['APINEN','BaP_fin','BbF_fin','BCAR_coa','BCAR_fin','BkF_fin','C2H4','C2H5OH','C2H6', \
#            'C3H6','C5H8','CH3CHO','CH3COE','CH3OH','CH4','CO','H2SO4_fin','HCHO','HONO', 'NC4H10', \
#            'NH3', 'NO2', 'NO', 'OCAR_coa', 'OCAR_fin', 'OXYL', 'PPM_big', 'PPM_coa', 'PPM_fin', \
#            'SO2', 'TMB', 'TOL', 'TPPM_coa', 'TPPM_fin']

##EXECUTE #############################################################
os.chdir(rundir)
#define species for the respective chemical mechanism, ONLY VALID for Marta's files!
if mechanism == 'melchior':
	species = ['APINEN','BaP_fin','BbF_fin','BCAR_fin','BkF_fin','C2H4','C2H5OH','C2H6', \
				'C3H6','C5H8','CH3CHO','CH3COE','CH3OH','CH4','CO','H2SO4_fin','HCHO','HONO', 'NC4H10', \
				'NH3', 'NO2', 'NO', 'OCAR_fin', 'OXYL', 'PPM_big', 'PPM_coa', 'PPM_fin', \
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

#load landuse file of the new domain
nc = xr.open_dataset(lufile)
tarlat = nc.lat.values
tarlon = nc.lon.values
tarlat_vec = np.sort(np.unique(tarlat))
tarlon_vec = np.sort(np.unique(tarlon))
nc.close()

#load emission files for the old domain
for mm in list(range(len(tarmonths))):
	print('INFO: target month is '+tarmonths[mm])
	for ii in list(range(len(species))):
		print('INFO: target species is '+species[ii])
		emisfile_src = emisdir_src+'/EMIS.'+domain_old+'.'+tarmonths[mm]+'.'+species[ii]+'.s.nc'		
		nc = xr.open_dataset(emisfile_src)
		srclat = nc.lat.values
		srclon = nc.lon.values
		srclat_vec = np.sort(np.unique(srclat))
		srclon_vec = np.sort(np.unique(srclon))

		#find lat and lon indices and cut put
		latind = np.where((srclat_vec>=tarlat_vec.min()) & (srclat_vec<=tarlat_vec.max()))[0]
		lonind = np.where((srclon_vec>=tarlon_vec.min()) & (srclon_vec<=tarlon_vec.max()))[0]
		ncnew = nc.isel(south_north=latind,west_east=lonind)

		#rename domain
		ncnew.attrs['Domain'] = domain_new
		ncnew.attrs['Origin'] = 'cut out from '+domain_old

		#then save the new file and close both source and new file
		emisfile_tar= emisdir_tar+'/EMIS.'+domain_new+'.'+tarmonths[mm]+'.'+species[ii]+'.s.nc'
		ncnew.to_netcdf(emisfile_tar)
		nc.close()
		ncnew.close()

print('INFO: cut_emissions.py has been run successfully!')
