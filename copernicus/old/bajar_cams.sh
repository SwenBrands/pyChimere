#!/bin/bash

USERNAME=swen.brands
KEY=NjamQcPd

#load cdos and auxiliar function
module load gcc/5.3.0  openmpi/1.10.2 cdo/1.7.0
source validate_url.sh


#STORE=${LUSTRE}/OP/DATOS/CICC/chimere2016a/copernicus
#TARDATE=$(date +"%Y%m%d")
#TARYEAR=$(date +"%Y")
#TARMONTH=$(date +"%m")
TARHOUR=00

TARDATE=$(date +"%Y%m%d" --date="3 days ago")
TARYEAR=$(date +"%Y" --date="3 days ago")
TARMONTH=$(date +"%m" --date="3 days ago")

TARHOUR=00 #TAR = "target"
TYPE=fc # an or fc
#LEADTIME=(000 003 006 009 012 015 018 021 024 027 030 033 036 039 042 045 048 051 054 057 060 063 066 069 072 075)
LEADTIME=(000 003 006 009 012 015 018 021 024 027 030)
#VARIABLES=(c2h6 ch4 co hno3 no2 go3 pan so2 aermr04 aermr05 aermr06 lnsp)
VARIABLES=(aermr01 aermr02 aermr03 aermr04 aermr05 aermr06 aermr08 aermr10 aermr11 lnsp)
TARDIR=${LUSTRE}/OP/DATOS/CICC/chimere2016a/chimere_bigfiles_2016a/MACC/${TARDATE}
CDIR=${LUSTRE}/OP/DATOS/CICC/chimere2016a/chimere_bigfiles_2016a/MACC
RUNDIR=/home/cesga/orballo/OP/LANZAR/chimere2016a/copernicus

#create target directory for the specific date and copy the alpha and beta coefficients to this directory
mkdir ${TARDIR}
cp ${CDIR}/aux/coeffs.nc ${TARDIR}/

#download analysis
for VAR in ${VARIABLES[*]}
    do
    for LEAD in ${LEADTIME[*]}
        do
        TARFILE=z_cams_c_ecmf_${TARDATE}000000_prod_${TYPE}_ml_${LEAD}_${VAR}.nc
        URL=ftp://dissemination.ecmwf.int/DATA/CAMS_NREALTIME/${TARDATE}${TARHOUR}/${TARFILE}
        echo ${URL}        
        
        #retry to download every 30 seconds
        while [ ! -f "${TARDIR}/${TARFILE}" ]
        do        
        wget --user ${USERNAME} --password ${KEY} --tries=1 -P ${TARDIR} ${URL}
            if [ ! -f "${TARDIR}/${TARFILE}" ]
            then
               echo "${TARFILE} does not exist, wait 120 seconds and try again"
               sleep 120
            fi
        done        
        echo "${TARFILE} has been downloaded"
        
        ##cut out target region
        #cdo sellonlatbox,-20,10,30,53 ${TARDIR}/${TARFILE} ${TARDIR}/myregion.nc #myregion.nc is a temporary file
        #mv ${TARDIR}/myregion.nc ${TARDIR}/${TARFILE}
       
        #convert log. pressure to Pa
        if [ "$VAR" == "lnsp" ]; then
        mv ${TARDIR}/${TARFILE} ${TARDIR}/tempfile.nc
        cdo -b F64 exp ${TARDIR}/tempfile.nc ${TARDIR}/${TARFILE}
        fi
    done    
done

##postprocess cams output to format required by CHIMERE
#module load chimplot
#python postprocess_cams.py ${TARDATE} ${TARYEAR} ${TARMONTH} ${TARHOUR}

##set global varialbe 'refpres' and copy to final directory
#module load nco/4.5.5
#ncatted -a refpres,global,c,f,1.0 ${TARDIR}/gasmet_${TARYEAR}${TARMONTH}.nc
#cp ${TARDIR}/gasmet_${TARYEAR}${TARMONTH}.nc ${CDIR}/${TARYEAR}/gasmet_${TARYEAR}${TARMONTH}.nc
