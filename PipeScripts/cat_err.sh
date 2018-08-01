#!/bin/bash 
#$ -l sps=1
#$ -o /sps/lsst/data/rcecile/Planck_BAO2/output_cat_err_lfZuccaAllFalse
#$ -M rcecile@in2p3.fr
#$ -q long
#$ -m eas
#! /usr/local/bin/bash -l

function do_err {
     command="${code1}addGausszerr -C ${dir}cat_lfZuccaAllFalse_zOrd$3_$1.fits -E $2 -c ZS -d 0.15 -O ${dir}cat_lfZuccaAllFalse_zOrd$3_$1_err$2.fits"
     $command

     command="${code1}addGausszerr -C ${dir}cat_lfZuccaAllFalse_zOrd$3_$1.fits -c ZS -p $pdf  -d 0.15 -O ${dir}cat_lfZuccaAllFalse_zOrd$3_$1_errP.fits"
     $command

     command="${code1}addGausszerr -C ${dir}cat_lfZuccaAllFalse_zOrd$3_$1.fits -c ZS -p ${bdt9} -d 0.15 -O ${dir}cat_lfZuccaAllFalse_zOrd$3_$1_errPBDT9.fits"
     $command

     command="${code1}addGausszerr -C ${dir}cat_lfZuccaAllFalse_zOrd$3_$1.fits -c ZS -p ${bdt8} -d 0.15 -O ${dir}cat_lfZuccaAllFalse_zOrd$3_$1_errPBDT8.fits"
     $command
}


dir="/sps/lsst/data/rcecile/Planck_BAO2/"
dirg="/sps/lsst/data/rcecile/Planck_BAO2_grids/"

code1="/sps/lsst/users/rcecile/BAOProgs/DirectSimPerso/exe/"

pdf="/sps/lsst/PhotozBAO/ricol/SIMU50deg/grid_absmag_atm/pdfz"
#pdf="/sps/lsst/data/rcecile/PhotozBAO/FastPz/pdfz"
bdt8="/sps/lsst/PhotozBAO/ricol/SIMU50deg/grid_absmag_atm/pdfz_bdt_0p8"
bdt9="/sps/lsst/PhotozBAO/ricol/SIMU50deg/grid_absmag_atm/pdfz_bdt_0p9"


Nslicez=100
 
i0=0
while [ ${i0} -lt  ${Nslicez} ]
do
    simu="do_err Slice$i0 0.03"	
    
    echo "Simu lancee : " $simu
    $simu
    
    (( i0++ ))
done
