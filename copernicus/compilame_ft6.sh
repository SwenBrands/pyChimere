#!/bin/bash -l

module purge
module load intel/2016
module load impi/5.1
#module load openmpi/1.10.2
module load szip/2.1
module load hdf5/1.8.16
module load netcdf-fortran/4.4.3
module load pnetcdf/1.7.0

./chimere.sh chimere.par c

