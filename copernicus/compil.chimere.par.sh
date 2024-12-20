# THIS FILE IS AUTO-Generated by ./chimere.sh from compil.chimere.par
# DO NOT CHANGE!
# Simulation No 1

export forecast='1'
export di=20180115
export hour=03
export nhours=24
export de=20180117
export endhour=03
export nested='no'
export dom='GAL3'
export nzdoms='6'
export nmdoms='4'
export nlevels='9'
export top_chimere_pressure='500'
export first_layer_pressure='997'
export meteo='WRF'
export metdom='d03'
export lab='gal3'
export clab='none'
export accur='low'
export nsconcs=24
export nsdepos=2
export sim='${di}${hour}_${dehour}_${lab}'
export emis='emep'
export surface_emissions='1'
export point_emissions='0'
export fire_emissions='0'
export iusebemis='1'
export iusedust='1'
export icuth='1'
export ifluxv='2'
export ifecan='1'
export iusebound='1'
export iuseini='1'
export glob_top_conc='0'
export mecachim='1'
export aero='1'
export nbins='10'
export seasalt='1'
export isaltp='0'
export carb='1'
export soatyp='2'
export trc='0'
export nequil='0'
export npeq='1'
export iadrydep='2'
export resusp='0'
export dtphys='10'
export dtchem='5'
export ngs='1'
export nsu='5'
export irs='1'
export iadv='2'
export iadvv='1'
export ideepconv='1'
export urbancorr='0'
export ilidar='0'
export imakechemprep='2'
export istopchemprep='0'
export imakemeteo='2'
export imakeaemis='1'
export imakefemis='0'
export imakerun='1'
export chimverb='5'
export simuldir='/mnt/lustre/scratch/newhome/xunta/rai/ba/OP/PRED/chimere2017r4'
export csimuldir='/mnt/lustre/scratch/newhome/xunta/rai/ba/OP/PRED/chimere2017r4'
export psimuldir='$simuldir'
export bigfilesdir='/mnt/lustre/scratch/newhome/xunta/rai/ba/OP/DATOS/CICC/chimere2017r4'
export iland_cover='3'
export landcover_dir='${bigfilesdir}/LANDCOVER'
export megan_data='${bigfilesdir}/MEGAN'
export bcdir='${bigfilesdir}'
export bcgas='MACC'
export bcgasdt='2'
export bcaer='MACC'
export bcaerdt='2'
export bcdust='MACC'
export bcdustdt='2'
export meteo_DIR='/mnt/lustre/scratch/newhome/xunta/rai/ba/SWEN/WRF_TEST/${di}/00'
export meteo_file='wrfout_d03_$(date -u -d "${di}" +%Y%m%d)_0000.nc'
export emissdir='${bigfilesdir}/EMISSIONS/galicia_2008'
export fire_emissdir='${bigfilesdir}/FIREMIS'
export datadir='${simuldir}/data_${dom}_${lab}'
export metdir='${meteo_DIR}/exdomout.${dom}'
export fnmeteo='${simuldir}/meteo.${sim}.nc'
export exdomout='${metdir}/exdomout.${sim}.nc'
export aemisdir='${datadir}'
export fnemisa='${aemisdir}/AEMISSIONS.${sim}.nc'
export fnemisf='${datadir}/FEMISSIONS.${sim}.nc'
export fnemisb='${datadir}/BEMISSIONS.${sim}.nc'
export fnemissalt='${datadir}/SEMISSIONS.${sim}.nc'
export fnemisd='${datadir}/DEMISSIONS.${sim}.nc'
export iniboundir='${simuldir}/INIBOUN.${nbins}'
export garbagedir='${chimere_root}/compilogs'
export clean='full'
