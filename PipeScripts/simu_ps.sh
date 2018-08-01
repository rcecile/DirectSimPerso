#!/bin/bash 
#$ -l sps=1
#$ -o /sps/lsst/data/rcecile/TJP_BAO_PS/output_ps
#$ -M rcecile@in2p3.fr
#$ -m eas
#$ -q long
#! /usr/local/bin/bash -l


 function do_ps_err {
     echo "==========================================================================="
     command="${code1}computepsfromarray -C ${dirg}grid_G2${sj}_${nx_i}_z${z_i}_G$1_nocata.fits -S ${dirps}simu_ps_${nx_i}_z${z_i}_r.fits -O ${dirps}PS_G2${sj}_${nx_i}_z${z_i}_err${1}_nocata"
     $command
  }

 function do_ps_perr {
     command="${code1}computepsfromarray -C ${dirg}grid_G2${sj}_${nx_i}_z${z_i}_errP_nocata.fits -S ${dirps}simu_ps_${nx_i}_z${z_i}_r.fits  -O ${dirps}PS_G2${sj}_${nx_i}_z${z_i}_errP_nocata"
     $command
     command="${code1}computepsfromarray -C ${dirg}grid_G2${sj}_${nx_i}_z${z_i}_errPpodds_nocata.fits -S ${dirps}simu_ps_${nx_i}_z${z_i}_r.fits  -O ${dirps}PS_G2${sj}_${nx_i}_z${z_i}_errPpodds_nocata"
     $command
 }

 function do_ps_perr_k {
     command="${code1}computepsfromarray -C ${dirg}grid_G2${sj}_${nx_i}_z${z_i}_errPpodds_nocata.fits -S ${dirps}simu_ps_${nx_i}_z${z_i}_r.fits -m $1  -O ${dirps}PS_G2${sj}_${nx_i}_z${z_i}_k${1}_errPpodds_nocata"
     $command
 }
#########################################################################################################
nCase=10

code0="/sps/lsst/users/rcecile/BAOProgs/SimLSS/exe/"     
code1="/sps/lsst/users/rcecile/BAOProgs/DirectSimPerso/exe/"     
err2=0.03

cell=8
	
z=('0.7' '1.4' '0.9' '1.3' '1.8')
D=(2507 4176 3055 3975 4875)
nz=(200 200 125 75 65)
nx=(450 875 640 900 1024)
nk=(-1 -1 -1 -1 0.03)

z=('0.7' '1.4' '0.9' '1.3' '1.8')
nz=(200 200 125 75 75)
nx=(450 875 640 900 500)
nk=(-1 -1 -1 -1 -1)
cell=(8 8 8 8 16)


#########################################################################################################
dirg="/sps/lsst/data/rcecile/TJP_BAO_grids/"
dirps="/sps/lsst/data/rcecile/TJP_BAO_PS/"
sj=""
for i in 4 2 3 
do
    nx_i=$((${nx[i]}))
    nz_i=$((${nz[i]}))
    z_i=${z[i]}
    k_i=$((${nk[i]}))
    cell_i=$((${cell[i]}))
    do_ps_err $err2
    if [ ${k_i} -lt 0 ]
    then 
	do_ps_perr 
    else
	do_ps_perr_k  ${k_i}
    fi
done
