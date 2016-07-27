#!/bin/bash 
#$ -l sps=1
#$ -o /sps/lsst/data/rcecile/TJP_BAO/output_self_err
#$ -M rcecile@in2p3.fr
#$ -q long
#$ -m eas
#! /usr/local/bin/bash -l


function do_sel {
    echo "==========================================================================="
    echo "===========================================================================%"

    command="${code1}getsf -F $2 -R -O $4 -o ${dirg}SelFunc_zOrd_errPpodds -z $3"
    echo "I launch " $command
    rm -f tmptmp
    $command

  command="${code1}getsf -F $2 -R -O $1 -o ${dirg}SelFunc_zOrd_errP -z $3"
   echo "I launch " $command
    rm -f tmptmp
    $command
}


function do_sel_gauss {
    echo "==========================================================================="

    command="${code1}getsf -F $2 -R -O $1 -o ${dirg}SelFunc_zOrd_Gauss0.03 -z $3"
    echo "I launch " $command
    rm -f tmptmp
    $command
}


code1="/sps/lsst/dev/rcecile/BAOProgs/DirectSim/exe/"
dir="/sps/lsst/data/rcecile/TJP_BAO/"
dirg="/sps/lsst/data/rcecile/TJP_BAO_grids/"

Nslice=70
Nslicez=100
 
i_min=0
i_max=$((${Nslice} -1 ))
obs_list5=${dir}cat_zOrd_Slice${i_min}.fits
full_list=${dir}/cat_zOrd_Slice${i_min}_ZONLY.fits
i1=$((${i_min} + 1))
while [ ${i1} -le  ${i_max} ]
do
    full=${dir}/cat_zOrd_Slice${i1}_ZONLY.fits
    full_list=${full_list},$full
    (( i1++ ))
done


i_min=0
i_max=$((${Nslicez} -1 ))
obs_list=${dir}cat_zOrd_AllSlice_errP_cataLT10.fits
obs_list2=${dir}cat_zOrd_AllSlice_errPpodds_cataLT10.fits
obs_list4=${dir}cat_zOrd_AllSlice_err0.03_cataLT10.fits

i1=$((${i_min}))
while [ ${i1} -le  ${i_max} ]
do
    obs=${dir}cat_zOrd_Slice${i1}_errP.fits
    obs_list=${obs_list},$obs

    obs=${dir}cat_zOrd_Slice${i1}_errPpodds.fits
    obs_list2=${obs_list2},$obs
	
    obs=${dir}cat_zOrd_Slice${i1}_err0.03.fits
    obs_list4=${obs_list4},$obs

    (( i1++ ))
done

do_sel $obs_list $full_list ZP,ZS,zs $obs_list2
do_sel_gauss $obs_list4 $full_list ZG,ZS,zs 
