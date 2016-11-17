#!/bin/bash 
#$ -l sps=1
#$ -o /sps/lsst/data/rcecile/Planck_BAO/output_cat_err
#$ -M rcecile@in2p3.fr
#$ -q long
#$ -m eas
#! /usr/local/bin/bash -l

# for noBAO case, just change the dir & dirg paths & nCase

function do_err {
     command="${code1}addGausszerr -C ${dir}cat_zOrd$3_$1.fits -E $2 -c ZS -d 0.09 -O ${dir}cat_zOrd$3_$1_err$2.fits"
     $command

     command="${code1}addGausszerr -C ${dir}cat_zOrd$3_$1.fits -c ZS -p $pdf  -d 0.09 -O ${dir}cat_zOrd$3_$1_errP.fits"
     $command

     command="${code1}addGausszerr -C ${dir}cat_zOrd$3_$1.fits -c ZS -p $podds -d 0.09 -O ${dir}cat_zOrd$3_$1_errPpodds.fits"
     $command

     command="${code1}addGausszerr -C ${dir}cat_zOrd$3_$1.fits -c ZS -p $bdt -d 0.09 -O ${dir}cat_zOrd$3_$1_errPBDT.fits"
     $command
}

dir="/sps/lsst/data/rcecile/Planck_BAO/"
dirg="/sps/lsst/data/rcecile/Planck_BAO_grids/"
nCase=1

#dir="/sps/lsst/data/rcecile/Planck_noBAO/"
#dirg="/sps/lsst/data/rcecile/Planck_noBAO_grids/"
#nCase=10

code1="/sps/lsst/dev/rcecile/BAOProgs/DirectSimPerso/exe/"
#pdf="/sps/lsst/PhotozBAO/ricol/SIMU50deg/newgrid_absmag/pdfz"
#podds="/sps/lsst/PhotozBAO/ricol/SIMU50deg/newgrid_absmag/pdfz_podds"
pdf="/sps/lsst/PhotozBAO/ricol/SIMU50deg/grid_absmag_atm/pdfz"
podds="/sps/lsst/PhotozBAO/ricol/SIMU50deg/grid_absmag_atm/pdfz_podds"
bdt="/sps/lsst/PhotozBAO/ricol/SIMU50deg/grid_absmag_atm/pdfz_bdt"

Nslicez=100
 
 
j=0
while [ ${j} -lt  ${nCase} ]
do
    i0=3
    while [ ${i0} -lt  ${Nslicez} ]
    do
	if [ ${nCase} -le  1 ]
	then
	    simu="do_err Slice$i0 0.03"	
	else
	    simu="do_err Slice$i0 0.03 _$j"	
	fi

	echo "Simu lancee : " $simu
	$simu
	
	(( i0++ ))
    done
        ((j++ ))

done
