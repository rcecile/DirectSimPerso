#!/bin/bash 
#$ -l sps=1
#$ -o /sps/lsst/data/rcecile/TJP_BAO/output_self_err
#$ -M rcecile@in2p3.fr
#$ -q long
#$ -m eas
#! /usr/local/bin/bash -l

# for noBAO case, just change the dir & dirg paths & nCase


function do_sel {
    echo "==========================================================================="
    echo "===========================================================================%"

    command="${code1}getsf -F $2 -O $4 -o ${dirg}SelFunc_G2$5_errPpodds -z $3"
    echo "I launch " $command
    rm -f tmptmp
    $command

    command="${code1}getsf -F $2 -O $1 -o ${dirg}SelFunc_G2$5_errP -z $3"
    echo "I launch " $command
    rm -f tmptmp
    $command
}


function do_sel_gauss {
    echo "==========================================================================="

    command="${code1}getsf -F $2 -O $1 -o ${dirg}SelFunc_G2$4_Gauss0.03 -z $3"
    echo "I launch " $command
    rm -f tmptmp
    $command
}

dir="/sps/lsst/data/rcecile/TJP_BAO/"
dirg="/sps/lsst/data/rcecile/TJP_BAO_grids/"
nCase=1

#dir="/sps/lsst/data/rcecile/TJP_noBAO/"
#dirg="/sps/lsst/data/rcecile/TJP_noBAO_grids/"
#nCase=10

code1="/sps/lsst/dev/rcecile/BAOProgs/DirectSim/exe/"

Nslice=70
 
j=0
while [ ${j} -lt  ${nCase} ]
do
    if [ ${nCase} -le  1 ]
    then
	sj=""
    else
	sj=_$j
    fi

    i_min=0
    i_max=$((${Nslice} -1 ))
    full_list=${dir}/cat_G2${sj}_Slice${i_min}_ZONLY.fits
    i1=$((${i_min} + 1))
    while [ ${i1} -le  ${i_max} ]
    do
	full=${dir}/cat_G2${sj}_Slice${i1}_ZONLY.fits
	full_list=${full_list},$full
	(( i1++ ))
    done
    
    
    i_min=3
    i_max=$((${Nslice} -1 ))
    obs_list=${dir}cat_G2${sj}_AllSlice_errP_cataLT10.fits
    obs_list2=${dir}cat_G2${sj}_AllSlice_errPpodds_cataLT10.fits
    obs_list4=${dir}cat_G2${sj}_AllSlice_err0.03_cataLT10.fits

    i1=$((${i_min}))
    while [ ${i1} -le  ${i_max} ]
    do
	obs=${dir}cat_G2${sj}_Slice${i1}_errP.fits
	obs_list=${obs_list},$obs
	
	obs=${dir}cat_G2${sj}_Slice${i1}_errPpodds.fits
	obs_list2=${obs_list2},$obs
	
	obs=${dir}cat_G2${sj}_Slice${i1}_err0.03.fits
	obs_list4=${obs_list4},$obs
	
	(( i1++ ))
    done
    
    do_sel_gauss $obs_list4 $full_list ZG,ZS,zs ${sj}
#    do_sel $obs_list $full_list ZP,ZS,zs $obs_list2 ${sj}
    ((j++ ))
done
