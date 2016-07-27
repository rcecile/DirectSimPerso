#!/bin/bash 
#$ -l sps=1
#$ -o /sps/lsst/data/rcecile/TJP_BAO/output_cat_err
#$ -M rcecile@in2p3.fr
#$ -q long
#$ -m eas
#! /usr/local/bin/bash -l


function do_err {
     echo "==========================================================================="
     echo "==========================================================================="
     command="${code1}addGausszerr -C ${dir}cat_zord_$1.fits -E $2 -r -c ZS -d 0.09 -O ${dir}cat_zOrd_$1_err$2.fits"
     echo "I launch " $command
     $command
 }

function do_perr { 
    echo "==========================================================================="
    echo "==========================================================================="
    command="${code1}addGausszerr -C ${dir}cat_zOrd_$1.fits -c ZS -r -p $pdf  -d 0.09 -O ${dir}cat_zOrd_$1_errP.fits"
    echo "I launch " $command
    $command
    command="${code1}addGausszerr -C ${dir}cat_zOrd_$1.fits -c ZS -r -p $podds -d 0.09 -O ${dir}cat_zOrd_$1_errPpodds.fits"
    echo "I launch " $command
    $command
}

dir="/sps/lsst/data/rcecile/TJP_BAO/"
dirg="/sps/lsst/data/rcecile/TJP_BAO_grids/"

code1="/sps/lsst/dev/rcecile/BAOProgs/DirectSimPerso/exe/"
pdf="/sps/lsst/PhotozBAO/ricol/SIMU50deg/newgrid_absmag/pdfz"
podds="/sps/lsst/PhotozBAO/ricol/SIMU50deg/newgrid_absmag/pdfz_podds"

Nslicez=100
 
i0=0
while [ ${i0} -lt  ${Nslicez} ]
do
    simu="do_err Slice$i0 0.03"
    echo "Simu lancee : " $simu
    $simu
    
    simu="do_perr Slice$i0"
    echo "Simu lancee : " $simu
    $simu
    
    (( i0++ ))
done
 