#!/bin/bash 
#$ -l sps=1
#$ -o /sps/lsst/data/rcecile/Planck_BAO2_PS/output_ps_err
#$ -M rcecile@in2p3.fr
#$ -m eas
#$ -q long
#! /usr/local/bin/bash -l


#########################################################################################################
 function do_ps_err {
     echo "==========================================================================="
     command="${code1}computepsfromgrids -C ${dirg}grids_superSN5_${nx_i}_z${z_i}_err0.03_0${1}_${2}_5cubes.fits -N ${norm} -O ${dirps}PS_superSN5_${nx_i}_z${z_i}_err0.03_0${1}_${2}"
      $command
 }

#########################################################################################################

code0="/sps/lsst/users/rcecile/BAOProgs/SimLSS/exe/"     
code1="/sps/lsst/users/rcecile/BAOProgs/DirectSimPerso/exe/"     
err2=0.03

cell=(8 8 8)
z=('0.5' '0.9' '1.5')
nx=(120 225 320)
nz=(125 125 175)
zmean=('0.51' '0.93' '1.57')

# coeff Cecile sqrt(fit rapport spectres th / spectres avec Norm=1)
NormNg=(1.30 1.18 1.08)


#########################################################################################################
dirg="/sps/lsst/data/rcecile/Planck_BAO2_grids/"
dirps="/sps/lsst/data/rcecile/Planck_BAO2_PS/"


for j in 0 2 3 4
do
    err=${err[j]}
    for i in 1
    do
	nx_i=$((${nx[i]}))
	nz_i=$((${nz[i]}))
	z_i=${z[i]}
	zmean_i=${zmean[i]}
	cell_i=$((${cell[i]}))
	norm=${NormNg[i]}
	for k in 0 1 2 3 4 5 6 7 8 9
	do
	    do_ps_err $j $k
	done
    done
done
