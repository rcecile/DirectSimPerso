#!/bin/bash 
#$ -l sps=1
#$ -o /sps/lsst/data/rcecile/Planck_BAO/output_cube_no
#$ -M rcecile@in2p3.fr
#$ -q long
#$ -m eas
#! /usr/local/bin/bash -l

function do_big {
     # replace -u ... by -U 15 for the last Planck 2015 results
     # here with the TJP convention
#     command="${code0}cmvginit3d -o ${dir}simu -u 0.7,0.3,0.05,0.7,-1.,0.,9.06581172724113E-05 -8 0.8154,-8. -G 2 -s 1000,0.001,0.5 -x $1,$3 -y $1,$3 -z $2,$3 -D $4 -O 2,2"
    command="${code0}cmvginit3d -o ${dir}simuRef -U 15 -G 1 -s 1000,0.001,0.5 -x $1,$3 -y $1,$3 -z $2,$3 -D $4 -O 2,2"
    $command   
}

function do_big_no {
     # add -a for random
     # add -L to not simulate oscillations
#     command="${code0}cmvginit3d -o ${dirn}simu$5 -u 0.7,0.3,0.05,0.7,-1.,0.,9.06581172724113E-05 -8 0.8154,-8. -G 2 -s 1000,0.001,0.5 -x $1,$3 -y $1,$3 -z $2,$3 -D $4 -O 2,2 -L -a"
    command="${code0}cmvginit3d -o ${dirn}simu$5 -U 15 -G 1 -s 1000,0.001,0.5 -x $1,$3 -y $1,$3 -z $2,$3 -D $4 -O 0,2  -a"
    $command
}

dir="/sps/lsst/data/rcecile/Planck_BAO/"
dirn="/sps/lsst/data/rcecile/Planck_BAO2_SP/"
code0="/sps/lsst/users/rcecile/BAOProgs/SimLSS/exe/"

cell=8
Nz_tot=700
Nx_tot=1250
h0=3200

simu="do_big  $Nx_tot $Nz_tot $cell $h0"
echo "Simu lancee : " $simu
$simu

nCase=0
j=0
while [ ${j} -lt  ${nCase} ]
do
    simu="do_big_no  $Nx_tot $Nz_tot $cell $h0 $j"
    echo "Simu lancee : " $simu
    $simu
    ((j++ ))
done
