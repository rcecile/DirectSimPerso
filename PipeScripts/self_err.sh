#!/bin/bash 
#$ -l sps=1
#$ -o /sps/lsst/data/rcecile/Planck_BAO2/output_self_err_Zucca
#$ -M rcecile@in2p3.fr
#$ -q long
#$ -m eas
#! /usr/local/bin/bash -l

# for noBAO case, just change the dir & dirg paths & nCase


function do_sel {
    echo "==========================================================================="
    echo "===========================================================================%"
    command="${code1}getsf -O $4 -o ${dirg}SelFunc_lfZuccaAllFalse_errPBDT8 -H ${dirg}SelFunc_lfZuccaAllFalse -z $2"
    rm -f tmptmp
    $command

    command="${code1}getsf -O $5 -o ${dirg}SelFunc_lfZuccaAllFalse_errPBDT9 -H ${dirg}SelFunc_lfZuccaAllFalse -z $2"
    rm -f tmptmp
    $command

    command="${code1}getsf -O $1 -o ${dirg}SelFunc_lfZuccaAllFalse_errP -H ${dirg}SelFunc_lfZuccaAllFalse -z $2"
    rm -f tmptmp
    $command
}


function do_sel_gauss {
    echo "==========================================================================="

    command="${code1}getsf -O $1 -o ${dirg}SelFunc_lfZuccaAllFalse_err0.03 -H ${dirg}SelFunc_lfZuccaAllFalse -z $2"
    echo "I launch " $command
    rm -f tmptmp
    $command
}

dir="/sps/lsst/data/rcecile/Planck_BAO2/"
dirg="/sps/lsst/data/rcecile/Planck_BAO2_grids/"

code1="/sps/lsst/users/rcecile/BAOProgs/DirectSimPerso/exe/"

Nslice=100
 
i_min=3
i_max=$((${Nslice} -1 ))
obs_list=${dir}cat_lfZuccaAllFalse_zOrd_AllSlice_errP_cataLT10.fits
obs_list2=${dir}cat_lfZuccaAllFalse_zOrd_AllSlice_errPpodds_cataLT10.fits
obs_list4=${dir}cat_lfZuccaAllFalse_zOrd_AllSlice_err0.03_cataLT10.fits
obs_list5=${dir}cat_lfZuccaAllFalse_zOrd_AllSlice_errPBDT8_cataLT10.fits
obs_list6=${dir}cat_lfZuccaAllFalse_zOrd_AllSlice_errPBDT9_cataLT10.fits

i1=$((${i_min}))
while [ ${i1} -le  ${i_max} ]
do
    obs=${dir}cat_lfZuccaAllFalse_zOrd_Slice${i1}_errP.fits
    obs_list=${obs_list},$obs
    
    obs=${dir}cat_lfZuccallFalse_zOrd_Slice${i1}_errPpodds.fits
    obs_list2=${obs_list2},$obs
    
    obs=${dir}cat_lfZuccaAllFalse_zOrd_Slice${i1}_err0.03.fits
    obs_list4=${obs_list4},$obs
    
    obs=${dir}cat_lfZuccaAllFalse_zOrd_Slice${i1}_errPBDT8.fits
    obs_list5=${obs_list5},$obs
    
    obs=${dir}cat_lfZuccaAllFalse_zOrd_Slice${i1}_errPBDT9.fits
    obs_list6=${obs_list6},$obs
    
    (( i1++ ))
done

do_sel_gauss $obs_list4 ZG,ZS,zs 
#do_sel $obs_list ZP,ZS,zs $obs_list2 $obs_list5 $obs_list6
