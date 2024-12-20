#!/bin/bash -l

#startdate=(20170601 20170602 20170603 20170604 20170605 20170606 20170607 20170608 20170609 20170610 20170611 20170612 20170613 20170614 20170615 20170616 20170617 20170618 20170619 20170620 20170621 20170622 20170623 20170624 20170625) #array containing the start dates of the forecast runs
startdate=(20170601) #array containing the start dates of the forecast runs
#variables=(PM10 PM25 NO2 O3 CO SO2)
variables=(PM10 NO2 O3)
#--------------------------------------------------------------------------------------------------
# Cargaomos os modulos
#--------------------------------------------------------------------------------------------------
module load gcc/5.3.0
module load python/2.7.11

#--------------------------------------------------------------------------------------------------
# Send to queue, loop through all days in <startdate> (defined variables_chimere00.sh)
#--------------------------------------------------------------------------------------------------
for tarvar in ${variables[*]}
    do
    for tardate in ${startdate[*]}
    do
        python plot_results.py ${tardate} af12 ${tarvar} af12
    done
done
exit 0
