#!/bin/bash 
#$ -l sps=1
#$ -o /sps/lsst/data/rcecile/Planck_BAO2_SP/output_cat_All
#$ -M rcecile@in2p3.fr
#$ -q long
#$ -m eas
#! /usr/local/bin/bash -l


function do_rdlss {
    echo "==========================================================================="
    echo "==========================================================================="
    # R = randomize in cell; G = goldencut; a = opening angle cut ( 1.04720 = 60 deg)
    command="${code1}rdlss -C ${dir}simu$2_$1_r.fits -O ${dir2}cat${2}_$1.fits -i 1 -R -G -a 1.04720 "
    echo $command
    $command
}


function do_sel_simple {
    echo "==========================================================================="

    command="${code1}getsf -F $2 -O $1 -o ${dirg}SelFunc${4} -z $3"
    echo "I launch " $command
    rm -f tmptmp
    $command

}
code1="/sps/lsst/dev/rcecile/BAOProgs/DirectSimPerso/exe/"
Nslice=41

dir="/sps/lsst/data/rcecile/Planck_BAO2_SP/"
dir2="/sps/lsst/data/rcecile/Planck_BAO2_SP/"
dirg="/sps/lsst/data/rcecile/Planck_BAO2_SP/"
nCase=10

j=0
while [ ${j} -lt  ${nCase} ]
do  
    i0=40
    while [ ${i0} -lt  ${Nslice} ]
    do
	simu="do_rdlss Slice$i0 $j"
	$simu
	((i0++ ))
    done
	
    ((j++ ))
done
