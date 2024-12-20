#!/bin/bash

# CHIMERE installation script
# Run this script from the directory where you want Chimere to be installed
# Usage: ./configure.sh [mychimere-label.sh]

#================================================================
# function to choose from a list of installed library versions
# usage: makeChoice <search_expr> <library_name> <out_var_name> [strict (1/0)]
makeChoice(){
   search_expr="$1"
   lib_name="$2"
   out_var_name=$3
   strict=$4
   choicelst=`locate ${search_expr}`
   if [ ${#choicelst} -eq 0 ] ; then
      if [ ${strict} = "1" ] ; then
         echo "You need to install a ${lib_name}"'! Bye.'
         exit 1;
      else
         echo 'Warning! '"Doing without ${lib_name}. Make sure you really don't need it"'!'
         eval $out_var_name=""
         return
      fi
   fi
   i=0
   echo "Available ${lib_name} versions:"
   for f in ${choicelst}; do
      if [ "${strict}" = "2" ] ; then
         [ -x ${f} ] || continue
         choice_array[$i]=${f}
      else
         choice_array[$i]=`dirname $(dirname ${f})`
      fi
      echo "[${i}] ${choice_array[$i]}"
      (( i++ ))
   done
   n=0
   if [ $i -gt 1 ]; then
      echo -n "Choose your version, please: "
      read n
      [ ${#choice_array[$n]} = "0" ] && { echo 'Wrong choice!' ; exit 1 ; }
   fi
   eval $out_var_name=${choice_array[$n]}
}
#================================================================
echo
echo "-------------------------------"
echo "   CHIMERE Configuration utility"
echo "-------------------------------"
echo

mychimere=$1

if [ -z $1 ] ; then
   echo "Format : $0 <mychimere.sh.version>"
   echo
   exit 1
elif [ $1 = "mychimere.sh" ]; then
   echo "Please, specify a name other than mychimere.sh"
   echo "It should be intuitive, corresponding to your particular machine and software configuration,"
   echo "   e.g., $0 mychimere.sh.cluster1.gfortran"
   echo
   echo "After you are done configuring your specified file,"
   echo "   mychimere.sh will be a symbolic link referring to this file"
   echo
   exit 1
fi
#-----------------------
mychimex=mychimere.sh.example
[ -s ${mychimex} ] || { echo ${mychimex}' not found! You can download from http://www.lmd.polytechnique.fr/chimere' ; exit 1; }
rm -f ${mychimere} || { echo "No right to remove ${mychimere} ?" ; exit 1; }
cp ${mychimex} ${mychimere} || { echo "No right to copy ${mychimex} => ${mychimere} ?" ; exit 1; }
echo
echo "Configuring ${mychimere}"
echo

# 1. Try to determine whether the system is 32 or 64 bit
syst=`arch 2>/dev/null` ||\
    { echo "arch command is not available. Is your system 64-bit? (y/[n])" ;
      read stat ;
      [[ "${stat}" = "y" || "${stat}" = "Y" ]] && syst="-64"; }
suf=""
echo $syst | grep "64" >/dev/null && suf="-64"

# 2. Choose Fortran compiler
makeChoice "bin/g95 bin/gfortran bin/ifort" "Fortran compiler" REALFC 2

# 3. Netcdf library
echo "----------------------------"
makeChoice libnetcdf.a "NetCDF C library" netcdfdirc 1
stat=$(nm ${netcdfdirc}/lib/libnetcdf.a  | grep -q nf90_open && echo "yes")
if [ "${stat}" = "yes" ] ; then
   netcdfdir=${netcdfdirc}
else
   echo "----------------------------"
   echo "Fortran functions are not found in ${netcdfdirc}/lib. Trying to locate in a separate lib"
   makeChoice libnetcdff.a "NetCDF Fortran library" netcdfdir 1
fi   

# 3.1. HDF5 (if use netcdf4)
if [ -e ${netcdfdirc}/bin/nc-config ] && ${netcdfdirc}/bin/nc-config --has-hdf5 ; then
   echo "----------------------------"
   makeChoice libhdf5.a "HDF5 library" hdfdir
else
   echo "Warning: HDF5 does not seem to be linked to ${netcdfdirc}"
fi

# 3.2. grib_api
echo "----------------------------"
makeChoice libgrib_api.a GRIB_API gribapidir 0
#[ -z ${gribexdir} ] || gribexdir=${gribexdir}/lib


# 4. MPI
#MF77=`which mpif77 2>/dev/null` || { echo 'You need to have MPI installed!'; exit 1; }
echo "----------------------------"
makeChoice "--regexp bin/mpif77$" "MPI library" MPIDIR 1
MPIBIN=${MPIDIR}/bin
# check if chosen MPi is LAM
[ -e ${MPIBIN}/lamboot ] && { sufmpi="lam" ; mpiframe="lam" ; }
mpif77=${MPIBIN}/mpif77
mpif90=
# check if mpif90 is available
[ -e ${MPIBIN}/mpif90 ] && mpif90=${MPIBIN}/mpif90
mpirun=${MPIBIN}/mpirun

# 4.1. IF several MPI versions are installed => need right dynamic libraries in LD_LIBRARY_PATH

#  4.1.1. remove current MPI lib from the beginning if already there
#         suppose only two lib directories
#~ lib_path=$LD_LIBRARY_PATH
#~ for i in 1 2 ; do
   #~ lib_path=$(sed 's,^'${MPIDIR}'/lib[/openmpi]*:,,' <<< ${lib_path})
#~ done

#  4.1.2. add current MPI lib to the beginning
#~ lib_path=${MPIDIR}/lib:${MPIDIR}/lib/openmpi:${lib_path}

# 5. pnetcdf
echo "----------------------------"
makeChoice "--regexp pnetcdf.inc" "Pnetcdf Library" pnetcdfdir


# 6. Resume
echo "----------------------------"
echo " ---- Your chosen options: ----"
echo "   Fortran Compiler       : ${REALFC}"
echo "   NetCDF Fortran Library : ${netcdfdir}"
echo "   NetCDF C Library       : ${netcdfdirc}"
echo "   HDF5 Library           : ${hdfdir}"
echo "   GRIB_API Library       : ${gribapidir}"
echo "   MPI Compilers          : ${mpif77}   ${mpif90}"
echo "   Pnetcdf Library        : ${pnetcdfdir}"
#~ echo "   Look for MPI libs : ${lib_path}"
echo -n "   Your system is "
if [ "${suf}" = "-64" ]; then
   echo "64 bit"
else
   echo "32 bit"
fi
if [ "${sufmpi}" = "lam" ]; then
   echo "   Your MPI version is LAM"
else
   if [ -e ${MPIBIN}/opal_wrapper ] ; then
      echo "   Your MPI version is Open MPI"
      sufmpi="ompi"
      mpiframe="openmpi"
   else
      echo 'Your MPI version is neither OpenMPI nor LAM!'
      echo '   Attention! You might have to manually adjust your Makefile.hdr'
      echo '   for your MPI realization!'
   fi
fi
echo -n "   Is it OK? ([y]/n) : "
read stat
[[ "${stat}" = "" || "${stat}" = "y" || "${stat}" = "Y" ]] || { echo 'Canceled! Come again!'; exit 1 ; }

# 6. Change PATHs in mychimere.sh

# netcdf & hdf5
sed -i 's,\(^export my_netcdfdir=\).*,\1'${netcdfdir}',' ${mychimere}
sed -i 's,\(^export my_netcdfdirc=\).*,\1'${netcdfdirc}',' ${mychimere}
sed -i 's,\(^export my_hdfdir=\).*,\1'${hdfdir}',' ${mychimere}

# gribex
sed -i 's,\(^export my_gribapi=\).*,\1'${gribapidir}',' ${mychimere}

# Fortran
for v in my_g95 my_gfortran my_pgf90 my_ifort ; do
   sed -i 's,\(^export '${v}'=\).*,\1'${REALFC}',' ${mychimere}
done

# MPI
sed -i 's,\(^export my_mpif90=\).*,\1'${mpif90}',' ${mychimere}
sed -i 's,\(^export my_mpif77=\).*,\1'${mpif77}',' ${mychimere}
sed -i 's,\(^export my_mpirun=\).*,\1'${mpirun}',' ${mychimere}
sed -i 's,\(^export my_mpiframe=\).*,\1'${mpiframe}',' ${mychimere}
#~ sed -i 's,\(^export LD_LIBRARY_PATH=\).*,\1'${lib_path}',' ${mychimere}
sed -i 's,\(^export my_mpilib=\).*,\1'${MPIDIR}/lib',' ${mychimere}

# my_hdr
sed -i 's,\(^export my_hdr=\).*,\1Makefile.hdr.'$(basename ${REALFC})${suf}-${sufmpi}',' ${mychimere}

# pnetcdfdir
sed -i 's,\(^export my_pnetcdflib=\).*,\1'${pnetcdfdir}/lib',' ${mychimere}
sed -i 's,\(^export my_pnetcdfinc=\).*,\1'${pnetcdfdir}/include',' ${mychimere}

# 8. If the <myfile> was not mychimere.sh => make symbolic link mychimere.sh -> <myfile>
[ ${mychimere} = "mychimere.sh" ] || ln -sf ../${mychimere} src/mychimere.sh

echo "New ${mychimere} is written"

# 8. make symbolic link makefiles.hdr/${my_hdr} src/Makefile.hdr
ln -sf makefiles.hdr/${my_hdr} src/Makefile.hdr

