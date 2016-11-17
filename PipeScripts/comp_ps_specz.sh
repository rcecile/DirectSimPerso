#!/bin/bash 
#$ -l sps=1
#$ -o /sps/lsst/data/rcecile/TJP_BAO_PS/output_ps
#$ -M rcecile@in2p3.fr
#$ -m eas
#$ -q long
#! /usr/local/bin/bash -l

#########################################################################################################

 function do_ps_base {
     echo "==========================================================================="
     command="${code1}computepsfromarray -C ${dirg}grid_gauss${sj}_${nx_i}_z${z_i}.fits -S ${dirps}simu_ps_${nx_i}_z${z_i}_r.fits  -O ${dirps}PS_gauss${sj}_${nx_i}_z${z_i}"
     $command
  }

 function do_ps_perr {
     command="${code1}computepsfromarray -C ${dirg}grid_gauss${sj}_${nx_i}_z${z_i}_errP.fits -S ${dirps}simu_ps_${nx_i}_z${z_i}_r.fits  -O ${dirps}PS_gauss${sj}_${nx_i}_z${z_i}_errP"
     $command
     command="${code1}computepsfromarray -C ${dirg}grid_gauss${sj}_${nx_i}_z${z_i}_errPpodds.fits -S ${dirps}simu_ps_${nx_i}_z${z_i}_r.fits  -O ${dirps}PS_gauss${sj}_${nx_i}_z${z_i}_errPpodds"
     $command
 }

 function do_ps_perr_k {
     command="${code1}computepsfromarray -C ${dirg}grid_gauss${sj}_${nx_i}_z${z_i}_errPpodds.fits -S ${dirps}simu_ps_${nx_i}_z${z_i}_r.fits -m $1  -O ${dirps}PS_gauss${sj}_${nx_i}_z${z_i}_k${1}_errPpodds"
     $command
 }
#########################################################################################################
nCase=10

code0="/sps/lsst/dev/rcecile/BAOProgs/SimLSS/exe/"     
code1="/sps/lsst/dev/rcecile/BAOProgs/DirectSimPerso/exe/"     
err2=0.03

z=('0.9' '1.3' '1.8' '0.5')
nz=(125 75 65 140)
nx=(640 900 1024 350)
nk=(0.12 0.12 0.12 -1)
cell=(8 8 8 8)


#########################################################################################################
dirg="/sps/lsst/data/rcecile/TJP_BAO_grids/"
dirps="/sps/lsst/data/rcecile/TJP_BAO_PS/"
sj=""
for i in 3
do
    nx_i=$((${nx[i]}))
    nz_i=$((${nz[i]}))
    z_i=${z[i]}
    k_i=${nk[i]}
    cell_i=$((${cell[i]}))
#    do_ps_base
    if [ ${k_i} -lt 0 ]
    then 
	do_ps_perr 
    else
	do_ps_perr_k  ${k_i}
    fi
done
