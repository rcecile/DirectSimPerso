#!/bin/bash 
#$ -l sps=1
#$ -o /sps/lsst/data/rcecile/Planck_BAO2_PS/output_ps_thin
#$ -M rcecile@in2p3.fr
#$ -m eas
#$ -q long
#! /usr/local/bin/bash -l

# just for plots, to compare with

 function do_ps_base {
     echo "==========================================================================="
     command="${code1}computepsfromgrids -C ${dirg}gridSN_lfZuccaAllFalse_${nx_i}_z${z_i}_1cubes.fits          -N ${norm} -O ${dirps}PS_SNthin_${nx_i}_z${z_i}"
     $command
  }

code0="/sps/lsst/users/rcecile/BAOProgs/SimLSS/exe/"     
code1="/sps/lsst/users/rcecile/BAOProgs/DirectSimPerso/exe/"     
err2=0.03

cell=(8 8 8)
z=('0.5' '0.9' '1.3')
nx=(60 100 150)
nz=(125 125 125)

# coeff Cecile sqrt(fit rapport spectres th / spectres avec Norm=1)
NormNg=(1.30 1.18 1.08)


#########################################################################################################
dirg="/sps/lsst/data/rcecile/Planck_BAO2_grids/"
dirps="/sps/lsst/data/rcecile/Planck_BAO2_PS/"


for i in 0 1 2
do
    nx_i=$((${nx[i]}))
    nz_i=$((${nz[i]}))
    z_i=${z[i]}
    cell_i=$((${cell[i]}))
    norm=${NormNg[i]}
    do_ps_base 
     
done
