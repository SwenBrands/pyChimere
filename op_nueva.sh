#!/bin/bash

ulimit -s unlimited

ARCH_VARIABLES=${1} #read environmental variables again
source ${ARCH_VARIABLES}
export WRFDIR
export PRED
export OP_STORE
export BCDIR
export CDIR
export CAMSDATE
export CAMSHOUR
export AUXDIR
export OP_CHIMERE

##postprocess cams output to the format required by CHIMERE

#create temporary directory where the results are saved in
mkdir ${AUXDIR}
#copy alpha and beta coefficients for model levels to AUXDIR
cp ${CDIR}/aux/coeffs.nc ${AUXDIR}

cd ${OP_CHIMERE}
echo "Info: start the postprocessing of the lateral boundary condition files from CAMS; postprocess_cams_regional.py is used for this."

#load all modules
module purge
module load meteogalicia/2021
module load intel/2021.3.0
module load impi/2021.3.0
module load python/2.7.18
module load pyngl/1.6.1-python-2.7.18
module load cdo/1.9.10

echo ${OP_CHIMERE}
echo "INFO: Checking whether the boundary data from CAMS have been downloaded to ${BCDIR}/${CAMSDATE}/${CAMSHOUR}..."
if [ -f ${BCDIR}/${CAMSDATE}/${CAMSHOUR}/download_complete_${CAMSDATE}_${CAMSHOUR}.flag ]; then
    echo "INFO: The CAMS data have been found and their postprocessing is now started..."
    python ${OP_CHIMERE}/copernicus/postprocess_cams_regional_137.py ${CAMSDATE} ${CAMSHOUR} > ./LOG/postprocess_cams_regional.log
    sleep 2
    ###join files, set global variable 'refpres' and copy to final directory
    #module purge
    #module load cdo/1.9.5
    #module load nco/4.7.7
    #sleep 5
    ##merge all variables (time-varying, climatological and coeffs) to one file
    #cdo -O merge ${AUXDIR}/C2H6.nc ${AUXDIR}/CH2O.nc ${AUXDIR}/CH4.nc ${AUXDIR}/CO.nc ${AUXDIR}/HNO3.nc ${AUXDIR}/ISOP.nc ${AUXDIR}/NO2.nc ${AUXDIR}/O3.nc ${AUXDIR}/PAN.nc ${AUXDIR}/SO2.nc ${AUXDIR}/SS1.nc ${AUXDIR}/SS2.nc ${AUXDIR}/SS3.nc ${AUXDIR}/DUST1.nc ${AUXDIR}/DUST2.nc ${AUXDIR}/DUST3.nc ${AUXDIR}/OM.nc ${AUXDIR}/BC.nc ${AUXDIR}/SO4.nc ${AUXDIR}/PSFC.nc ${AUXDIR}/BIGALK.nc ${AUXDIR}/BIGENE.nc ${AUXDIR}/CH3CHO.nc ${AUXDIR}/GLYOXAL.nc ${AUXDIR}/NH3.nc ${AUXDIR}/H2O2.nc ${AUXDIR}/TOLUENE.nc ${AUXDIR}/C10H16.nc ${AUXDIR}/TEMP.nc ${AUXDIR}/coeffs.nc ${AUXDIR}/gasmet_${TARYEAR}${TARMONTH}.nc
    #sleep 5
    #ncatted -a refpres,global,c,f,1.0 ${AUXDIR}/gasmet_${TARYEAR}${TARMONTH}.nc
    mkdir ${CDIR}/${TARYEAR}
    echo "change format of the boundary conditions file ${AUXDIR}/gasmet_${TARYEAR}${TARMONTH}.nc to 64-bit offest..."
    ncks -6 -O ${AUXDIR}/gasmet_${TARYEAR}${TARMONTH}.nc ${AUXDIR}/gasmet_${TARYEAR}${TARMONTH}.nc #converts boundary conditions file ($PRED/aux/gasmet*.nc) from netcdf4 to netcdf3 64 bit offset and overwrites; the ncks command brings the coordinates in the same order than those provided in MACC/clim
    mv ${AUXDIR}/gasmet_${TARYEAR}${TARMONTH}.nc ${CDIR}/${TARYEAR}/gasmet_${TARYEAR}${TARMONTH}.nc
    sleep 2
else
    echo "WARNING: The lateral boundary data from CAMS is missing in ${BCDIR}/${CAMSDATE}/${CAMSHOUR}."
    echo "CHIMERE will be run with climatological boundary conditions..."
fi

##launch CHIMERE
module purge
module load meteogalicia/2021
module load intel/2021.3.0
module load impi/2021.3.0
module load hdf5/1.10.7
module load netcdf/4.7.4
module load netcdf-fortran/4.5.3
module load pnetcdf/1.12.2
sleep 2
echo "launch CHIMERE2017r4 for ${TARDATE} ${TARHOUR}..."
./chimere.sh chimere.par f ${TARDATE} ${RUNHOURS} ${ENDDATE}
sleep 2

#postprocess the model output
module purge
module load meteogalicia/2021
module load intel/2021.3.0
module load impi/2021.3.0
module load python/2.7.18
module load pyngl/1.6.1-python-2.7.18
sleep 2

#step one of the postprocessing
for region in ${REG1} ${REG2}
do
    for species in O3 PM10 PM10bio PM10ant PM25 PM25bio PM25ant CO SO2 NO2 pBCAR #pSALT 
    do
    python -W ignore ${OP_CHIMERE}/figures/plot_results_new.py ${TARDATE}${INITHOUR} ${ENDDATE}${INITHOUR} ${region} ${species} > ${OP_CHIMERE}/LOG/mapping_${region}_${species}.log &
    sleep 1
    done
done

#check if step 1 of the postprocessing is complete
while [ ! -f ${OP_CHIMERE}/FLAG/map_${REG2}_${TARDATE}${INITHOUR}_pBCAR.flag ]
    do
    echo "INFO: the figures are still being generated, wait another 30 seconds..."
    sleep 30
done
echo "INFO: step 1 of the postprocessing has finished!"

#step 2 of the postprocessing
for region in ${REG1} ${REG2}
do
    for species in pOCAR pSALT pH2SO4 pDUST #pNA pH2SO4 pHCL pWATER
    #for species in pDUST
    do
    python -W ignore ${OP_CHIMERE}/figures/plot_results_new.py ${TARDATE}${INITHOUR} ${ENDDATE}${INITHOUR} ${region} ${species} > ${OP_CHIMERE}/LOG/mapping_${region}_${species}.log &
    sleep 1
    done
done

##check if step 2 of the postprocessing is complete
while [ ! -f ${OP_CHIMERE}/FLAG/map_${REG2}_${TARDATE}${INITHOUR}_pDUST.flag ]
    do
    echo "INFO: the figures are still being generated, wait another 30 seconds..."
    sleep 30
done
echo "INFO: step 2 of the postprocessing has finished!"

#Copy model output to "RESULTADOS"
mkdir ${RES}/${TARDATE}
cp -r ${PRED}/out.${TARDATE}${INITHOUR}_${ENDDATE}${INITHOUR}_${REG1}.nc  ${RES}/${TARDATE}
cp -r ${PRED}/end.${TARDATE}${INITHOUR}_${ENDDATE}${INITHOUR}_${REG1}.nc  ${RES}/${TARDATE}
cp -r ${PRED}/out.${TARDATE}${INITHOUR}_${ENDDATE}${INITHOUR}_${REG2}.nc  ${RES}/${TARDATE}
cp -r ${PRED}/end.${TARDATE}${INITHOUR}_${ENDDATE}${INITHOUR}_${REG2}.nc  ${RES}/${TARDATE}

##copy only the fine domain files in case a single nested run is defined in variables.sh
#if [ ${singlerun} == "no" ]
    #then
    #cp ${COARSERES}/${TARDATE}/out.${TARDATE}${INITHOUR}_${ENDDATE}${INITHOUR}_${REG1}.nc ${PRED}
    #cp ${COARSERES}/${TARDATE}/end.${TARDATE}${INITHOUR}_${ENDDATE}${INITHOUR}_${REG1}.nc ${PRED}    
#else
    #echo "INFO: This is a single nested run. Only the fine domain files are copied to the RESULTS folder."
#fi

cp -r ${PRED}/figs ${RES}/${TARDATE}
cp ${PRED}/chimere.par ${RES}/${TARDATE}
cp ${OP_CHIMERE}/LOG/CH_${TARDATE}_${TARHOUR}*.out ${RES}/${TARDATE}
sleep 2

##Clean the PRED folder
cd ${PRED}
rm -r data*
rm -r INIBOUN*
rm dep.${TARDATE}${INITHOUR}_${ENDDATE}${INITHOUR}*.nc
rm out.${TARDATE}${INITHOUR}_${ENDDATE}${INITHOUR}*.nc
rm chimer*
rm -r figs
rm -r exdomou*
rm end.*.nc

##Clean the CAMS folder
rm -r ${AUXDIR}

#write flag
echo "INFO: op_nueva.sh has been run completely on $(date)"
touch ${OP_CHIMERE}/FLAG/OP_${TARDATE}_end.flag
exit 0
