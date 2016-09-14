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
     command="${code0}cmvginit3d -o ${dir}simuTJP_G1 -u 0.7,0.3,0.05,0.7,-1.,0.,9.06581172724113E-05 -8 0.8154,-8. -G 1 -s 1000,0.001,0.5 -x $1,$3 -y $1,$3 -z $2,$3 -D $4 -O 2,2"
 #    echo "I launch " $command
#    $command
    command="${code0}cmvginit3d -o ${dir}simuTJP_G2 -u 0.7,0.3,0.05,0.7,-1.,0.,9.06581172724113E-05 -8 0.8154,-8. -G 2 -s 1000,0.001,0.5 -x $1,$3 -y $1,$3 -z $2,$3 -D $4 -O 2,2"
     echo "I launch " $command
    $command

     command="${code0}cmvginit3d -o ${dir}simuU15_G1 -U 15. -G 1 -s 1000,0.001,0.5 -x $1,$3 -y $1,$3 -z $2,$3 -D $4 -O 2,2"
#     echo "I launch " $command
#    $command
 
     command="${code0}cmvginit3d -o ${dir}simuU15_G2 -U 15. -G 2 -s 1000,0.001,0.5 -x $1,$3 -y $1,$3 -z $2,$3 -D $4 -O 2,2"
 #    echo "I launch " $command
 #   $command
 
     command="${code0}cmvginit3d -o ${dir}simuU15 -U 15. -s 1000,0.001,0.5 -x $1,$3 -y $1,$3 -z $2,$3 -D $4 -O 2,2"
#     echo "I launch " $command
#    $command
 }
 function do_big_no {
     # add -a for random
     # add -L to not simulate oscillations
     command="${code0}cmvginit3d -o ${dirn}simu$5 -u 0.7,0.3,0.05,0.7,-1.,0.,9.06581172724113E-05 -8 0.8154,-8. -G 1 -s 1000,0.001,0.5 -x $1,$3 -y $1,$3 -z $2,$3 -D $4 -O 2,2 -L -a"

     echo "I launch " $command
     $command
 }

dir="/sps/lsst/data/rcecile/TJP_BAO/"
dirn="/sps/lsst/data/rcecile/TJP_noBAO/"
code0="/sps/lsst/dev/rcecile/BAOProgs/SimLSS/exe/"

cell=8
Nz_tot=750
Nx_tot=1350
h0=3200

simu="do_big  $Nx_tot $Nz_tot $cell $h0"
echo "Simu lancee : " $simu
$simu

for i in  0 1 2 3 4 5 6 7 8 9
do
    simu="do_big_no  $Nx_tot $Nz_tot $cell $h0 $i"
#    echo "Simu lancee : " $simu
#    $simu
done
