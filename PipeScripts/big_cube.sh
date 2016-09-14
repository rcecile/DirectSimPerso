#!/bin/bash 
#$ -l sps=1
#$ -o /sps/lsst/data/rcecile/TJP_BAO/output_cube_no
#$ -M rcecile@in2p3.fr
#$ -q huge
#$ -m eas
#! /usr/local/bin/bash -l

 function do_big {
     # replace -u ... by -U 15 for the last Planck 2015 results
     # here with th TJP convention
     command="${code0}cmvginit3d -o ${dir}simu2200 -u 0.7,0.3,0.05,0.7,-1.,0.,9.06581172724113E-05 -8 0.8154,-8. -G 2 -s 1000,0.001,0.5 -x $1,$3 -y $1,$3 -z $2,$3 -D $4 -O 2,2"
     echo "I launch " $command
     $command

}
 function do_big_no {
     # add -a for random
     # add -L to not simulate oscillations
     command="${code0}cmvginit3d -o ${dirn}simu2200$5 -u 0.7,0.3,0.05,0.7,-1.,0.,9.06581172724113E-05 -8 0.8154,-8. -G 2 -s 1000,0.001,0.5 -x $1,$3 -y $1,$3 -z $2,$3 -D $4 -O 2,2 -L -a"

     echo "I launch " $command
     $command
 }

dir="/sps/lsst/data/rcecile/TJP_BAO/"
dirn="/sps/lsst/data/rcecile/TJP_noBAO/"
code0="/sps/lsst/dev/rcecile/BAOProgs/SimLSS/exe/"

cell=8
Nz_tot=700
Nx_tot=1600
Nx_tot=2200
h0=3200

simu="do_big  $Nx_tot $Nz_tot $cell $h0"
echo "Simu lancee : " $simu
$simu

nCase=10
j=0
while [ ${j} -lt  ${nCase} ]
do
    simu="do_big_no  $Nx_tot $Nz_tot $cell $h0 $j"
    $simu
    ((j++ ))
done
