#!/bin/bash 
#$ -l sps=1
#$ -o /sps/lsst/data/rcecile/TJP_BAO/output_cat_simple
#$ -M rcecile@in2p3.fr
#$ -q long
#$ -m eas
#! /usr/local/bin/bash -l

# for noBAO case, just change the dir & dirg paths

function do_rdlss {
    echo "==========================================================================="
    echo "==========================================================================="
    command="${code1}rdlss -C ${dir}simu_$1_r.fits -O ${dir}cat_$1.fits -i 1 -R -M -g ${filecut} -a 1.04720 "
    echo "I launch " $command
    $command
}


function do_sel_simple {
    echo "==========================================================================="

    command="${code1}getsf -F $2 -R -O $1 -o ${dirg}SelFunc_ -z $3"
    echo "I launch " $command
    rm -f tmptmp
    $command

}

dir="/sps/lsst/data/rcecile/TJP_BAO/"
dirg="/sps/lsst/data/rcecile/TJP_BAO_grids/"

code1="/sps/lsst/dev/rcecile/BAOProgs/DirectSimPerso/exe/"
filecut="/sps/lsst/PhotozBAO/ricol/SIMU50deg/newgrid_absmag/goldencut.txt"

Nslice=70
 
i0=0
while [ ${i0} -lt  ${Nslice} ]
do
    simu="do_rdlss Slice$i0"
#    echo "Simu lancee : " $simu
#    $simu
    (( i0++ ))
done
 
i_min=0
i_max=$((${Nslice} -1 ))
obs_list5=${dir}cat_Slice${i_min}.fits
full_list=${dir}/cat_Slice${i_min}_ZONLY.fits
i1=$((${i_min} + 1))

while [ ${i1} -le  ${i_max} ]
do
    obs=${dir}cat_Slice${i1}.fits
    obs_list5=${obs_list5},$obs

    full=${dir}/cat_Slice${i1}_ZONLY.fits
    full_list=${full_list},$full

    (( i1++ ))
done

do_sel_simple $obs_list5 $full_list ZS,ZS,zs


