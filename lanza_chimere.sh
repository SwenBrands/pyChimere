#!/bin/bash -l

####################################################################################################
# lanza_chimere.sh:
####################################################################################################
#
# Description: This script launches the operative forecasts of the CHIMERE chemical transport model
#
####################################################################################################

#--------------------------------------------------------------------------------------------------
# Check if the script is called correctly
#--------------------------------------------------------------------------------------------------
ARCH_VARIABLES=${1}

NOME_TAREFA="$(basename "$0")"
MENSAXE="Uso: ${NOME_TAREFA} [arquivo de variables]"

if [ ! -f "${ARCH_VARIABLES}" ]
then 
    echo "Non se atopa o arquivo ${ARCH_VARIABLES}"
    echo ${MENSAXE}
    exit 1
fi

#--------------------------------------------------------------------------------------------------
# Load the variables
#--------------------------------------------------------------------------------------------------
#source /etc/profile.d/zz-cesga.sh
source ${ARCH_VARIABLES}

#--------------------------------------------------------------------------------------------------
# Check if the WRF simulation is complete
#--------------------------------------------------------------------------------------------------
while [ ! -f ${WRFDIR}/${TARDATE}/${TARHOUR}/fin_operativa_${TARDATE}${TARHOUR} ]
    do
    echo "INFO: CHIMERE is waiting for the meteorological input data from WRF, wait another 60 seconds and try again..."
    sleep 60
done
echo "INFO: The flag for the WRF run has been found! Proceed to see whether the boundary data from CAMS have been downloaded..."
sleep 1

#--------------------------------------------------------------------------------------------------
# Delete namelist (called chimere.par) from the previous execution and replace with namelist for boundary conditions from CAMS. If these are absent, copy namlist for using climatology instead.
#--------------------------------------------------------------------------------------------------
echo "INFO: Checking whether the boundary data from CAMS have been downloaded to ${BCDIR}/${CAMSDATE}/${CAMSHOUR}..."
rm ${OP_CHIMERE}/chimere.par
sleep 1
if [ -f ${BCDIR}/${CAMSDATE}/${CAMSHOUR}/download_complete_${CAMSDATE}_${CAMSHOUR}.flag ]; then
    echo "INFO: The CAMS data have been found. Copying chimere.par.cams file for execution with CAMS data..."
    cp ${OP_CHIMERE}/namelists/chimere.par.cams ${OP_CHIMERE}/chimere.par
    echo "INFO: CHIMERE run for ${STARTDATE}${TARHOUR} will be run with boundary conditions from CAMS from ${CAMSDATE}${CAMSHOUR}..."
else
    echo "WARNING: The flag for the download of the CAMS data cannot be found! Most probably, there are files missing in ${BCDIR}/${CAMSDATE}/${CAMSHOUR} that could not be downloaded. Check also if the path to the files exists..."
    echo "Copying chimere.par.clim file for execution with CAMS data..."
    cp ${OP_CHIMERE}/namelists/chimere.par.clim ${OP_CHIMERE}/chimere.par
    echo "INFO: CHIMERE run for ${STARTDATE}${TARHOUR} will be exectued with boundary conditions from climatololgical files..."
fi

##--------------------------------------------------------------------------------------------------
## Check if the download of the CAMS data (boundary conditions) is complete
##--------------------------------------------------------------------------------------------------
#while [ ! -f ${BCDIR}/${TARDATE}/${TARHOUR}/download_complete_${TARDATE}_${TARHOUR}.flag ]
    #do
    #echo "INFO: CHIMERE is waiting for the boundary data from CAMS, wait another 60 seconds and try again..."
    #sleep 60
#done
#echo "INFO: The CAMS data have been completely downloaded, their postprocessing is now initiated..."
#sleep 1

## delete all postprocessibng flags from previous simulations
echo "INFO: Delete all postprocessing flags from previous simulations..."
cd ${OP_CHIMERE}/FLAG/
rm map*.flag
#rm *flag

##copy the initial conditions file from the previous day
if [ ${contrun} == "yes" ]
    then
    echo "INFO: This is a contintuation run; copy the initial conditions file from the previous day..."
    cp ${RES}/${PREVDATE}/end.${PREVDATE}*.nc ${PRED}/
else
    echo "INFO: This is NOT a continuation run; the init files are not copied..."
fi

#copy the coarse resolution run in case a single nested run is defined in chimere.par
if [ ${singlerun} == "yes" ]
    then
    echo "INFO: This is a single nested run. Copy the files of the coarse resolution domain" ${REG1}
    cp ${COARSERES}/${TARDATE}/end.${TARDATE}${INITHOUR}_${ENDDATE}${INITHOUR}_${REG1}.nc ${PRED}
    cp ${COARSERES}/${TARDATE}/out.${TARDATE}${INITHOUR}_${ENDDATE}${INITHOUR}_${REG1}.nc ${PRED}
else
    echo "INFO: This is NOT a single nested run. Therefore there is no need to copy any coarse run files."
fi

##--------------------------------------------------------------------------------------------------
## Send to queue
##--------------------------------------------------------------------------------------------------
echo "INFO: CHIMERE2017r4 operative model is sent to queue for " ${TARDATE}" "${TARHOUR} 

#     --nodes=${NumNodos} \
cd ${OP_CHIMERE} #folder where op_nueva.sh is located
QSUB="sbatch \
     --time=${exectime} \
     --job-name=chimere_testrun \
     --export=ALL \
     --begin=now \
     --output=${OP_CHIMERE_LOG}/CH_${TARDATE}_${TARHOUR}_%j.out \
     --ntasks-per-node=${NumProcesosNodo} \
     --ntasks=${NumMPI} \
     --cpus-per-task=${NumOMP} \
     --mem=${Memoria} \
     --mail-type=ALL \
     --mail-user=${DESTINATARIO} \
     ./op_nueva.sh variables_run.sh" #o variables_00_137.sh
echo ${QSUB}
${QSUB}


sleep 3
    
exit 0
