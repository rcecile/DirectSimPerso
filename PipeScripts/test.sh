#!/bin/bash 
#$ -l sps=1
#$ -o /sps/lsst/data/rcecile/TJP_BAO/output_self_err
#$ -M rcecile@in2p3.fr
#$ -q long
#$ -m eas
#! /usr/local/bin/bash -l

# for noBAO case, just change the dir & dirg paths & nCase


function do_sel {
    command="${code1}getsf -O $1 -o ${dirg}SelFunc_Test -F  $1 -z $2"
    rm -f tmptmp
    $command
    command="${code1}getsf -O $1 -o ${dirg}SelFunc_TestBDT -F  $1 -z $2 -b 0."
    rm -f tmptmp
    $command
    command="${code1}getsf -O $1 -o ${dirg}SelFunc_TestPODDS -F  $1 -z $2 -p 0.5"
    rm -f tmptmp
    $command
}
function do_sel0 {
   command="${code1}getsf -O $1 -o ${dirg}SelFunc_TestOld -F $3 -z $2 "
    rm -f tmptmp
    $command
}
dirg="/sps/lsst/data/rcecile/TJP_BAO_grids/"

code1="/sps/lsst/dev/rcecile/BAOProgs/DirectSimPerso/exe/"
file="/sps/lsst/data/rcecile/TJP_BAO/file0.fits"
file0="/sps/lsst/data/rcecile/TJP_BAO/cat_G2_Slice60_errPpodds.fits"
file00="/sps/lsst/data/rcecile/TJP_BAO/cat_G2_Slice60_ZONLY.fits"

#do_sel $file ZP,ZS,ZS 

do_sel0 $file0 ZP,ZS,zs $file00
