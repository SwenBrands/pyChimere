#!/bin/bash -l
echo "Renames the <modified...txt> files created by merge_emissions.py to match the name of the original files. Prior to this, savecopies are made in the <origfiles> directory"

species=(CO NH3 NMVOC NOx PM10 PM2_5 PMcoarse SOx)
activities=(A_PublicPower B_Industry J_Waste K_AgriLivestock M_Other NT)
srcpath=${LUSTRE}/SWEN/chimere2020r1/emisurf2020r3/annual-EMEP01x01/preproc/EMEP2019_merged
prodyear=2021 #year the EMEP files were produced
valyear=2019 #year the EMEP files are valid for

## EXECUTE #############################################################
cd ${srcpath}

#for i in `seq 0 4`;
COUNT1=`expr ${#species[@]} - $#1`
echo "COUNT1 is ${COUNT1}"
COUNT2=`expr ${#activities[@]} - $#1`
echo "COUNT2 is ${COUNT2}"

for i in `seq 0 ${COUNT1}`;
    do
    mkdir ${srcpath}/${species[i]}_${prodyear}_GRID_${valyear}/origfiles
    cd ${srcpath}/${species[i]}_${prodyear}_GRID_${valyear}
    
    for j in `seq 0 ${COUNT2}`;
        do
        origfile=${srcpath}/${species[i]}_${prodyear}_GRID_${valyear}/${species[i]}_${activities[j]}_${prodyear}_GRID_${valyear}.txt
        savecopy=${srcpath}/${species[i]}_${prodyear}_GRID_${valyear}/origfiles/${species[i]}_${activities[j]}_${prodyear}_GRID_${valyear}.txt
        modfile=${srcpath}/${species[i]}_${prodyear}_GRID_${valyear}/modified_${species[i]}_${activities[j]}_${prodyear}_GRID_${valyear}.txt
        ##make a savecopy in origfiles folder
        cp ${origfile} ${savecopy}
        ##then overwrite original files with modified files
        mv ${modfile} ${origfile}
    done
done
echo "rename_modified_emepfiles.sh has been executed successfully"
