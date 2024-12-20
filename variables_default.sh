#!/bin/bash -l

#this script defines the default variables for the CHIMERE sensitivity experiments

TARDATE=19000101 #is overwritten by runme.sh, DO NOT CHANGE
ENDDATE=19000102 # dito
PREVDATE=18991231 # dito
CAMSDATE=19000101
TARHOUR=00
CAMSHOUR=00
INITHOUR=03 #init hour for all runs
RUNHOURS=24 #run duration in hours for all runs, 24 for one day
exectime=00:53:00 #requested run time (send to queue), 00:47:00 is sufficient for summer runs, max. time consumption on 20190730, 00:43:00 is sufficient for the rest of the summer days
singlerun=no #is this a single nested run?
contrun=yes

#additional variables
TARYEAR="${TARDATE:0:4}"
TARMONTH="${TARDATE:4:2}"

#domains
REG1='gal15r'
REG2='gal3'

# Variables para a execucion
OP=${HOME}/OP                                #added by Swen
OP_CHIMERE=${HOME}/OP/LANZAR/chimere2017r4
OP_STORE=${STORE_METEO}/OP
OP_CHIMERE_LOG=${OP_CHIMERE}/LOG

#rutas de los ficheros de salida
PRED=${LUSTRE}/OP/PRED/chimere2017r4
RES=${STORE2}/OP/DATOS/RESULTADOS/chimere2017r4/test
COARSERES=${STORE2}/DATOS/RESULTADOS/chimere2017r4/exp2 #directory to the coarse resolution file, only needed in case a single nested run is executed
GRAF=${OP_STORE}/DATOS/RESULTADOS/chimere2017r4 #currently not used in the chimere2017r4 folder
PROD=${OP_STORE}/DATOS/PRODUCTOS/chimere2017r4 #currently not used in the chimere2017r4 folder
SCRIPTS=${OP_CHIMERE}/SCRIPTS
AUXDIR=${PRED}/aux

#variables needed by postprocess_cams.py
WRFDIR=${STORE2}/OP/DATOS/RESULTADOS/WRF_ARW_DET
#WRFDIR=${OP_STORE}/DATOS/RESULTADOS/WRF_ARW_DET
BCDIR=${STORE2}/OP/DATOS/CICC/CAMS
#BCDIR=${OP_STORE}/DATOS/CICC/CAMS
CDIR=${OP_STORE}/DATOS/CICC/chimere2017r4/MACC

#Variables para mandar a la cola
Particion="short" #Particion="thin-shared" # partition name
NumNodos=1                                      # Número de nodos que se van a usar
NumProcesosNodo=24                               # Número de procesos por nodo
NumMPI=$((${NumNodos}*${NumProcesosNodo}))       # Número total de procesos (MPI)
NumOMP=1                                         # Número de CPUs por proceso (OMP)
Memoria=128GB
Logfile=${OP_CHIMERE}/LOG/slurm_${TARDATE}_${TARHOUR}.txt

#Variables para la reservea
MODELO="CH"
REMITENTE="saraiba@ft2.cesga.es"                # Remitente del envio de comprobacion
DESTINATARIO="swen.brands@gmail.com"  # Destinatario
LISTA_ARID="/tmp/lista_CH.arid"
NOME_TAREFA=CH2017
