#!/bin/bash 
#$ -l sps=1
#$ -o /sps/lsst/data/rcecile/Planck_BAO2_PS/output_ps_noErr
#$ -M rcecile@in2p3.fr
#$ -m eas
#$ -q long
#! /usr/local/bin/bash -l

# just for plots, to compare with

 function do_ps_base {
     echo "==========================================================================="
     command="${code1}computepsfromgrids -C ${dirg}grids_testSN_${nx_i}_z${z_i}_5cubes.fits          -N ${norm} -O ${dirps}PS_testSN_${nx_i}_z${z_i}"
     $command
  }

#########################################################################################################
 function do_ps_err {
     echo "==========================================================================="
     command="${code1}computepsfromgrids -C ${dirg}grids_testSN_${nx_i}_z${z_i}_err$1_5cubes.fits -N ${norm} -O ${dirps}PS_testSN_${nx_i}_z${z_i}_err${1}"
     $command
 }

 function do_ps_perr {
     command="${code1}computepsfromgrids -C ${dirg}grids_testSN_${nx_i}_z${z_i}_errP_5cubes.fits     -N ${norm} -O ${dirps}PS_testSN_${nx_i}_z${z_i}_errP"
     $command
     command="${code1}computepsfromgrids -C ${dirg}grids_testSN_${nx_i}_z${z_i}_errPBDT9_5cubes.fits -N ${norm} -O ${dirps}PS_testSN_${nx_i}_z${z_i}_errPBDT9"
     $command
     command="${code1}computepsfromgrids -C ${dirg}grids_testSN_${nx_i}_z${z_i}_errPBDT8_5cubes.fits -N ${norm} -O ${dirps}PS_testSN_${nx_i}_z${z_i}_errPBDT8"
     $command
 }
#########################################################################################################

code0="/sps/lsst/users/rcecile/BAOProgs/SimLSS/exe/"     
code1="/sps/lsst/users/rcecile/BAOProgs/DirectSimPerso/exe/"     
err2=0.03

cell=(8 8 8)
z=('0.5' '0.9' '1.5')
#nx=(120 120 120)
#nz=(125 125 125)
nx=(80 80 80)
nz=(80 80 80)
zmean=('0.51' '0.93' '1.57')

# coeff Cecile sqrt(fit rapport spectres th / spectres avec Norm=1)
NormNg=(1.30 1.18 1.08)


#########################################################################################################
dirg="/sps/lsst/data/rcecile/Planck_BAO2_grids/"
dirps="/sps/lsst/data/rcecile/Planck_BAO2_PS/"

err2=0.03

for i in 0 1 2
do
    nx_i=$((${nx[i]}))
    nz_i=$((${nz[i]}))
    z_i=${z[i]}
    zmean_i=${zmean[i]}
    cell_i=$((${cell[i]}))
    norm=${NormNg[i]}
    do_ps_base 
    do_ps_err $err2
    do_ps_perr 
    
done
