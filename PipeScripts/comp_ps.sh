#!/bin/bash 
#$ -l sps=1
#$ -o /sps/lsst/data/rcecile/Planck_BAO2_PS/output_ps_Zucca
#$ -M rcecile@in2p3.fr
#$ -m eas
#$ -q long
#! /usr/local/bin/bash -l

# just for plots, to compare with
 function do_big {
     echo "==========================================================================="
     command="${code0}cmvginit3d -o ${dirps}simu_ps_${nx_i}_z${zmean_i} -U 15 -G 1 -s 1000,0.001,0.5 -x ${nx_i},${cell_i} -y ${nx_i},${cell_i} -z ${nz_i},${cell_i} -Z ${zmean_i} -O 2,2"
     $command
 }

 function do_ps_base {
     echo "==========================================================================="
     command="${code1}computepsfromgrids -C ${dirg}grids${suff}_${nx_i}_z${z_i}_5cubes.fits          -N ${norm} -O ${dirps}PS${suff}_${nx_i}_z${z_i}"
     $command
  }

#########################################################################################################
 function do_ps_err {
     echo "==========================================================================="
     command="${code1}computepsfromgrids -C ${dirg}grids${suff}_${nx_i}_z${z_i}_err$1_5cubes.fits -N ${norm} -O ${dirps}PS${suff}_${nx_i}_z${z_i}_err${1}"
     $command
 }

 function do_ps_perr {
     command="${code1}computepsfromgrids -C ${dirg}grids${suff}_${nx_i}_z${z_i}_errP_5cubes.fits     -N ${norm} -O ${dirps}PS${suff}_${nx_i}_z${z_i}_errP"
     $command
     command="${code1}computepsfromgrids -C ${dirg}grids${suff}_${nx_i}_z${z_i}_errPBDT9_5cubes.fits -N ${norm} -O ${dirps}PS${suff}_${nx_i}_z${z_i}_errPBDT9"
     $command
     command="${code1}computepsfromgrids -C ${dirg}grids${suff}_${nx_i}_z${z_i}_errPBDT8_5cubes.fits -N ${norm} -O ${dirps}PS${suff}_${nx_i}_z${z_i}_errPBDT8"
     $command
 }
#########################################################################################################

code0="/sps/lsst/users/rcecile/BAOProgs/SimLSS/exe/"     
code1="/sps/lsst/users/rcecile/BAOProgs/DirectSimPerso/exe/"     
err2=0.03

cell=(8 8 8 8)
z=('0.5' '0.9' '1.5' '1.3')
nx=(120 225 320 300)
nz=(125 125 175 125)
zmean=('0.51' '0.93' '1.57' '1.36')

# coeff Cecile sqrt(fit rapport spectres th / spectres avec Norm=1)
NormNg=(1.30 1.18 1.08 1.10)


#########################################################################################################
dirg="/sps/lsst/data/rcecile/Planck_BAO2_grids/"
dirps="/sps/lsst/data/rcecile/Planck_BAO2_PS/"

suff="_lfZuccaAllFalse"

err2=0.03

for i in 3
do
    nx_i=$((${nx[i]}))
    nz_i=$((${nz[i]}))
    z_i=${z[i]}
    zmean_i=${zmean[i]}
    cell_i=$((${cell[i]}))
    norm=${NormNg[i]}

    do_big 
    do_ps_base 

    do_ps_err $err2
    do_ps_perr 
    
done
