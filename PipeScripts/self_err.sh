#!/bin/bash 
#$ -l sps=1
#$ -o /sps/lsst/data/rcecile/Planck_BAO/output_self_err
#$ -M rcecile@in2p3.fr
#$ -q long
#$ -m eas
#! /usr/local/bin/bash -l

# for noBAO case, just change the dir & dirg paths & nCase


function do_sel {
    echo "==========================================================================="
    echo "===========================================================================%"
    command="${code1}getsf -O $4 -o ${dirg}SelFunc_gold$6_errPBDT8_bin -H ${dirg}SelFunc_gold_$6 -z $2"
    rm -f tmptmp
#     $command

    command="${code1}getsf -O $5 -o ${dirg}SelFunc_gold$6_errPBDT9_bin -H ${dirg}SelFunc_gold_$6 -z $2"
    rm -f tmptmp
    $command

    command="${code1}getsf -O $1 -o ${dirg}SelFunc_gold$6_errP_bin -H ${dirg}SelFunc_gold_$6 -z $2"
    rm -f tmptmp
#    $command
}


function do_sel_gauss {
    echo "==========================================================================="

    command="${code1}getsf -O $1 -o ${dirg}SelFunc_gold$3_ -H ${dirg}SelFunc_gold_$3 -z $2"
    echo "I launch " $command
    rm -f tmptmp
    $command
}

#dir="/sps/lsst/data/rcecile/Planck_BAO/"
#dirg="/sps/lsst/data/rcecile/Planck_BAO_grids/"
#nCase=1

dir="/sps/lsst/data/rcecile/Planck_noBAO/"
dirg="/sps/lsst/data/rcecile/Planck_noBAO_grids/"
nCase=10

code1="/sps/lsst/dev/rcecile/BAOProgs/DirectSimPerso/exe/"

Nslice=100
 
j=9
while [ ${j} -lt  ${nCase} ]
do
    if [ ${nCase} -le  1 ]
    then
	sj=""
    else
	sj=_$j
    fi
   
    i_min=3
    i_max=$((${Nslice} -1 ))
    obs_list=${dir}cat_gold${sj}_AllSlice_errP_bin_cataLT10.fits
    obs_list2=${dir}cat_gold${sj}_AllSlice_errPpodds_cataLT10.fits
    obs_list4=${dir}cat_gold${sj}_AllSlice_err0.03_cataLT10.fits
    obs_list5=${dir}cat_gold${sj}_AllSlice_errPBDT8_bin_cataLT10.fits
    obs_list6=${dir}cat_gold${sj}_AllSlice_errPBDT9_bin_cataLT10.fits

    i1=$((${i_min}))
    while [ ${i1} -le  ${i_max} ]
    do
	obs=${dir}cat_zOrd${sj}_Slice${i1}_errP_bin.fits
	obs_list=${obs_list},$obs
	
	obs=${dir}cat_zOrd${sj}_Slice${i1}_errPpodds.fits
	obs_list2=${obs_list2},$obs
	
	obs=${dir}cat_zOrd${sj}_Slice${i1}_err0.03.fits
	obs_list4=${obs_list4},$obs
	
	obs=${dir}cat_zOrd${sj}_Slice${i1}_errPBDT8_bin.fits
	obs_list5=${obs_list5},$obs
	
	obs=${dir}cat_zOrd${sj}_Slice${i1}_errPBDT9_bin.fits
	obs_list6=${obs_list6},$obs
	
	(( i1++ ))
    done
    
#    do_sel_gauss $obs_list4 ZG,ZS,zs ${sj}
    do_sel $obs_list ZP,ZS,zs $obs_list2 $obs_list5 $obs_list6 ${sj} 
    ((j++ ))
done
