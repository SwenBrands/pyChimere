#!/bin/bash
export LANG=en_US
export LC_NUMERIC=C
export LC_ALL=C
ulimit -s unlimited

#source /etc/profile.d/zz-cesga.sh
source ${ARCH_VARIABLES}

#---------------------------------------------------------------------------------------
#  Main script running the CHIMERE model
#
#  http://www.lmd.polytechnique.fr/chimere
#  Questions: chimere@lmd.polytechnique.fr
#
#  Usage: ./chimere.sh [param_file [task [start_date [nhours [end_date [run_number] ] ] ] ] ]
#     task is one of : (c)ompilation, (s)equential, (p)arallel, (f)lags
#     run_number : if specified overrides values given in the 'runs' row of chimere.par
#
#---------------------------------------------------------------------------------------

echo
echo -e "\033[1;47m       _________________________________      \033[0m"
echo -e "\033[1;47m                                              \033[0m"
echo -e "\033[1;47m       CHIMERE chemistry-transport model      \033[0m"
echo -e "\033[1;47m                  2017                        \033[0m"
echo -e "\033[1;47m       _________________________________      \033[0m"
echo -e "\033[1;47m                                              \033[0m"

#---------------------------------------------------------------------------------------
# To make system settings visible by the scripts

# chimere root directory
export chimere_root=`pwd`

source src/mychimere.sh || { echo '   => Try running ./config.sh mychimere.sh.<my_configuration>'  ;  exit 1 ; }
rm -f src/Makefile.hdr
ln -s ../makefiles.hdr/${my_hdr} src/Makefile.hdr

# MPI configuration
export mpirun=${my_mpirun}    # path to your mpirun command
export mpiframe=${my_mpiframe}   # "lam" or "openmpi"
export hostfile=${my_hostfile}   # only for OpenMPI

# Required utilities. You may have to define manually their full path
export MAKE=${my_make}
export AWK=${my_awk}
export NCDUMP=${my_ncdump}   # should be a netcdf4 compatible ncdump
export compmode=${my_mode}   # only for OpenMPI

# model version
version=`basename ${chimere_root}`

# MAKE,AWK and NCDUMP are set by the calling script. We check it again
${MAKE} --version 2>/dev/null >/dev/null || \
    { echo "You need gmake to run CHIMERE. Bye ..."; exit 1; }
${AWK} --version 2>/dev/null >/dev/null || \
    { echo "You need gawk to run CHIMERE. Bye ..."; exit 1; }
which ${NCDUMP} 2>/dev/null >/dev/null || \
    { echo "You need ncdump to run CHIMERE. Bye ..."; exit 1; }

export MAKE
export AWK
export NCDUMP

#---------------------------------------------------------------------------------------
# chimere command line arguments

# default values
task=c # f for running, c for compiling
chimparams="chimere.par"

# Check if (default) or (modified by user) chimere.par file
if [ $# -gt 0 ] ; then
   chimparams=$1
fi

runlist=$(gawk '$1=="runs"{s="";for (i=3;i<=NF;i++) if (substr($i,1,1)=="#") break; else s=s" "$i; print s}' ${chimparams})

# Check flags for compilation tasks: (s)equential, (p)arallel, (c)ompilation
if [ $# -gt 1 ] ; then
   task=$2
   if [ ${task:0:1} != "s" -a ${task:0:1} != "p" -a ${task:0:1} != "c"  -a ${task:0:1} != "f" ] ; then
      echo "   Unknown task ${task} for compilation."
      echo "   Valid tasks are: (s)equential, (p)arallel, (c)ompilation, (f)params"
      exit 1
   fi
fi

export task

# Simulation dates as arguments to the script: override those in chimere.par
di=$3
nhours=$4
de=$5

# which column of chimere.par should we run? This overrides the "runs" line of chimere.par
if [ $# -ge 6 ] ; then
   runlist=$6
fi

# Get the 1st column to run to set params for the 1st time
for r in ${runlist} ; do
   ru=${r}
   break
done

#---------------------------------------------------------------------------------------
# First read of the namelist just to have main paths of the simulation before compilation
. ${chimere_root}/scripts/define_params.sh ${chimparams} ${ru} ${di} ${nhours} ${de}

# Define TMP directory on time for all simulations
export tmplab=$(date +%Y%m%d_%H-%M)
export chimere_tmp=${simuldir}/tmp${di}-${lab}_${tmplab}

echo "   TMP directory: "${chimere_tmp}
mkdir -p ${garbagedir}
mkdir -p ${chimere_tmp}
# Dir to store executables
export exedir=${chimere_root}/exe_${compmode}

# Compilation (if necessary or requested by the user)

if [ ${task:0:1} == "c" ] ; then
   . ${chimere_root}/scripts/chimere-compil.sh
   ${chimere_root}/scripts/cleanup.sh
   exit 0
fi


# Copy of the exe dir without compilation
if [[ ${task:0:1} == "s" || ${task:0:1} == "p" || ${task:0:1} == "f" ]] ; then
   . ${chimere_root}/scripts/chimere-copyexe.sh
fi

# Make chemistry input files

. ${chimere_root}/scripts/make-chemistry.sh

#---------------------------------------------------------------------------------------

for runs in ${runlist} ; do

   echo
   echo -e "\033[1;47m o   Simulation directory: ${simuldir} \033[0m"

   echo "   Verbose mode level "${chimverb}

   # Get simulation parameters for the current run, redefine dates, and define some common derived params
   cd ${chimere_root}
   . ${chimere_root}/scripts/define_params.sh ${chimparams} ${runs} ${di}${hour} ${nhours} ${de}${endhour}

   # Check flags
   . ${chimere_root}/scripts/chimere-flags-define.sh ${task:0:1}

   # There is no online version implemented.
   online=0

   # Print screen of all flags status
   if [[ "${task:0:1}" = "s" || "${task:0:1}" = "p" || "${task:0:1}" = "f" ]] ; then
      . ${chimere_root}/scripts/chimere-flags-synthesis.sh
   fi

   # Make domains

   . ${chimere_root}/scripts/chimere-domain.sh
   [ $? == 0 ] || { echo "Abnormal termination of chimere-domain.sh"; exit 1; }

   # export to following steps
   export arome_root AWK bins boundaer boundgas chimere_root chemprepfic cmd coarsefile_info \
      day DIRO firedir flatreduc GRIBEX hour ifires iusedust imakeaemis imakefemis iusebemis \
      readaemis readfemis ipointsources lsm_suffix MAKE meteo_file_2d meteo_file_3d meteo_pfile \
      month NCDUMP nhours oro_suffix pemissdir reactive runlist runs scratch \
      tmplab tmp_meteo_pfile version year mpirun dustdir chimparams metdom lunom online

   # step 1 : sequential
   ${chimere_root}/scripts/chimere-step1.sh
   [ $? == 0 ] || { echo "Abnormal termination of step1.sh"; exit 1; }

   ## step 2 : parallel
   #lamparams="${my_lamparams}" #modified by Aurelio
   #ompiparams="${my_ompiparams}" #modified by Aurelio

   #if [ ${mpiframe} == "lam" ]; then
     #mpiparams="${lamparams}"
   #elif [ ${mpiframe} == "openmpi" ]; then
     #mpiparams="${ompiparams}"
   #else
     #echo "Unknown MPI frame $mpiframe. Bye"
     #exit 1
   #fi
   #export mpiparams

   # Run parallel part (Chimere core)
   if [ $imakerun == 1 ] || [ $iusebemis == 1 ]  || [ $iusedust == 1 ] ; then
      if [ "${task:0:1}" != "s" ] ; then
         ${chimere_root}/scripts/chimere-step2.sh
         [ $? == 0 ] || { echo "Abnormal termination of step2.sh"; exit 1; }
      fi
   fi

   cd ${chimere_root}

done

# write flag # added by Swen
date > ${chimere_root}/FLAG/end_chimere_${di}_${hour}.flag

# clean-up : sequential
${chimere_root}/scripts/cleanup.sh
