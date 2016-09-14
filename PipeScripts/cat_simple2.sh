#!/bin/bash 
#$ -l sps=1
#$ -o /sps/lsst/data/rcecile/TJP_noBAO/output_cat_simple
#$ -M rcecile@in2p3.fr
#$ -q long
#$ -m eas
#! /usr/local/bin/bash -l

# for noBAO case, just change the dir & dirg paths & nCase

function do_rdlss {
    echo "==========================================================================="
    echo "==========================================================================="
    # R = randomize in cell, z = z-axis of the cube is los axis = redshift
    command="${code1}rdlss -C ${dir}simu$2_$1_r.fits -O ${dir}cat_G2$3_$1.fits -i 1 -R -M -z -g ${filecut} -a 1.04720 "
    echo "I launch " $command
    $command
}


function do_sel_simple {
    echo "==========================================================================="

    command="${code1}getsf -F $2 -O $1 -o ${dirg}SelFunc_G2$4_full -z $3"
    echo "I launch " $command
    rm -f tmptmp
    $command

}

#dir="/sps/lsst/data/rcecile/TJP_BAO/"
#dirg="/sps/lsst/data/rcecile/TJP_BAO_grids/"
#nCase=1

dir="/sps/lsst/data/rcecile/TJP_noBAO/"
dirg="/sps/lsst/data/rcecile/TJP_noBAO_grids/"
nCase=4

code1="/sps/lsst/dev/rcecile/BAOProgs/DirectSimPerso/exe/"
#filecut="/sps/lsst/PhotozBAO/ricol/SIMU50deg/newgrid_absmag/goldencut.txt"
filecut="/sps/lsst/PhotozBAO/ricol/SIMU50deg/grid_absmag_atm/goldencut.txt"

Nslice=70
 
j=13
while [ ${j} -lt  ${nCase} ]
do
    if [ ${nCase} -le  1 ]
    then
	sj=""
	sj2=""
    else
	sj=$j
	sj2=_$j
     fi

    i0=38
    while [ ${i0} -lt  ${Nslice} ]
    do
	simu="do_rdlss Slice$i0 ${sj} ${sj2}  "
	$simu ${sj} ${sj2}
	((i0++ ))
    done
    ((j++ ))
done

#########################################################################################################################
 
j=2
while [ ${j} -lt  ${nCase} ]
do
    if [ ${nCase} -le  1 ]
    then
	sj=""
    else
	sj=_$j
    fi

    i0=0
    i_max=$((${Nslice} -1 ))
    obs_list5=${dir}cat_G2${sj}_Slice${i0}.fits
    full_list=${dir}/cat_G2${sj}_Slice${i0}_ZONLY.fits
   
    i1=$((${i0} + 1))
    
    while [ ${i1} -le  ${i_max} ]
    do
	obs=${dir}cat_G2${sj}_Slice${i1}.fits    
	full=${dir}/cat_G2${sj}_Slice${i1}_ZONLY.fits
	
	obs_list5=${obs_list5},$obs
	full_list=${full_list},$full
	
	(( i1++ ))
    done
    
    do_sel_simple $obs_list5 $full_list ZS,ZS,zs ${sj}
    ((j++ ))
done
