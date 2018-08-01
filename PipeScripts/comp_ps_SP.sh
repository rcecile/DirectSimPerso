#!/bin/bash 
#$ -l sps=1
#$ -o /sps/lsst/data/rcecile/Planck_BAO2_PS/output_ps
#$ -M rcecile@in2p3.fr
#$ -m eas
#$ -q long
#! /usr/local/bin/bash -l

 function do_big {
     echo "==========================================================================="
     command="${code0}cmvginit3d -o ${dirps}simu_ps_${nx_i}_z${zmean_i} -U 15 -G 1 -s 1000,0.001,0.5 -x ${nx_i},${cell_i} -y ${nx_i},${cell_i} -z ${nz_i},${cell_i} -Z ${zmean_i} -O 2,2"
     $command
 }


 function do_ps_base {
     echo "==========================================================================="
     command="${code1}computepsfromarray -C ${dirg}grids_nCR$1_${nx_i}_z${z_i}_5cubes.fits -S ${dirps}simu_ps_${nx_i}_z${zmean_i}_r.fits -x -N ${norm} -O ${dirps}PS_nCR$1_${nx_i}_z${z_i}"
     $command
  }

#########################################################################################################

code0="/sps/lsst/dev/rcecile/BAOProgs/SimLSS/exe/"     
code1="/sps/lsst/dev/rcecile/BAOProgs/DirectSimPerso/exe/"     

cell=(8 8 8)
z=('0.5' '0.9' '1.5')
nx=(120 225 320)
nz=(125 125 175)
zmean=('0.51' '0.93' '1.57')

#coef fournis par Reza, prends en compte la mise à 0 des pixels négatifs avant conversion en galaxies
#NormNg=(1.13 1.08 1.034)
#NormNg=(1.0 1.0 1.0)
# coeff Cecile sqrt(fit rapport spectres th / spectres avec Norm=1)
NormNg=(1.30 1.18 1.08)


#########################################################################################################
dirg="/sps/lsst/data/rcecile/Planck_BAO2_SP/"
dirps="/sps/lsst/data/rcecile/Planck_BAO2_SP/"

nCase=10
j=3
while [ ${j} -lt  ${nCase} ]
do  
    for i in 0 1 2
    do
	nx_i=$((${nx[i]}))
	nz_i=$((${nz[i]}))
	z_i=${z[i]}
	zmean_i=${zmean[i]}
	cell_i=$((${cell[i]}))
	norm=${NormNg[i]}
#	do_big
	do_ps_base $j
		
    done
    ((j++ ))
done

