#!/bin/bash 
#$ -l sps=1
#$ -o /sps/lsst/data/rcecile/Planck_BAO2_SP/output_cube_no
#$ -M rcecile@in2p3.fr
#$ -q long
#$ -m eas
#! /usr/local/bin/bash -l

 function do_big {
     # replace -u ... by -U 15 for the last Planck 2015 results
     # here with the TJP convention
#     command="${code0}cmvginit3d -o ${dir}simu -u 0.7,0.3,0.05,0.7,-1.,0.,9.06581172724113E-05 -8 0.8154,-8. -G 2 -s 1000,0.001,0.5 -x $1,$3 -y $1,$3 -z $2,$3 -D $4 -O 2,2"
     command="${code0}cmvginit3d -o ${dir}simu -U 15 -G 1 -s 1000,0.001,0.5 -x $1,$3 -y $1,$3 -z $2,$3 -D $4 -O 2,2"
     $command

}

dir="/sps/lsst/data/rcecile/Planck_BAO/"
dirn="/sps/lsst/data/rcecile/Planck_BAO2_SP/"
code0="/sps/lsst/dev/rcecile/BAOProgs/SimLSS/exe/"

cell=8
Nz_tot=700
Nx_tot=1250
h0=3200

simu="do_big  $Nx_tot $Nz_tot $cell $h0"
echo "Simu lancee : " $simu
