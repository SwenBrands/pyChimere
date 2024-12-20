#!/bin/bash -l

#final, most rapid configuration for ft7
module purge
#module load cesga/2018
module load meteogalicia/2018
module load intel/2018.5.274
module load openmpi/2.1.1
module load netcdf-fortran/4.4.4
module load pnetcdf/1.7.0
#optionally include NCOs to work with emiSURF2020; this changes netcdf and hdf versions
module load nco/4.7.7

## test with cesga/2018
#module load impi/2018.3.222
#module load impi/5.1.3.223
#module load impi/2018.4.274

## tests with meteogalicia/2018
#module load intel/2016.4.258
#module load impi/2018.4.274
#module load openmpi/2.1.1
#module load netcdf-fortran/4.4.4
#module load pnetcdf/1.7.0

./chimere.sh chimere.par c

