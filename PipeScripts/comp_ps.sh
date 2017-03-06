#!/bin/bash 
#$ -l sps=1
#$ -o /sps/lsst/data/rcecile/Planck_BAO_PS/output_ps
#$ -M rcecile@in2p3.fr
#$ -m eas
#$ -q long
#! /usr/local/bin/bash -l

 function do_big {
     echo "==========================================================================="
     command="${code0}cmvginit3d -o ${dirps}simu_ps_${nx_i}_z${zmean_i} -U 15 -G 1 -s 1000,0.001,0.5 -x ${nx_i},${cell_i} -y ${nx_i},${cell_i} -z ${nz_i},${cell_i} -Z ${zmean_i} -O 2,2"
     $command
 }

 function do_big_no {
     echo "==========================================================================="
     command="${code0}cmvginit3d -o ${dirps}simu_ps_${nx_i}_z${zmean_i} -U 15 -L -G 1 -s 1000,0.001,0.5 -x ${nx_i},${cell_i} -y ${nx_i},${cell_i} -z ${nz_i},${cell_i} -Z ${zmean_i} -O 2,2"
     $command
 }

 function do_ps_base {
     echo "==========================================================================="
     command="${code1}computepsfromarray -C ${dirg}grids${sj}_${nx_i}_z${z_i}_5cubes.fits -S ${dirps}simu_ps_${nx_i}_z${zmean_i}_r.fits  -O ${dirps}PS_newErr${sj}_${nx_i}_z${z_i}"
      $command
  }

#########################################################################################################
 function do_ps_err {
     echo "==========================================================================="
     command="${code1}computepsfromarray -C ${dirg}grids${sj}_${nx_i}_z${z_i}_err$1_5cubes.fits -S ${dirps}simu_ps_${nx_i}_z${zmean_i}_r.fits -O ${dirps}PS_newErr${sj}_${nx_i}_z${z_i}_err${1}"
     $command
  }

 function do_ps_perr {
     command="${code1}computepsfromarray -C ${dirg}grids${sj}_${nx_i}_z${z_i}_errP_5cubes.fits -S ${dirps}simu_ps_${nx_i}_z${zmean_i}_r.fits      -O ${dirps}PS_newErr${sj}_${nx_i}_z${z_i}_errP"
     $command
     command="${code1}computepsfromarray -C ${dirg}grids${sj}_${nx_i}_z${z_i}_errPBDT9_5cubes.fits -S ${dirps}simu_ps_${nx_i}_z${zmean_i}_r.fits  -O ${dirps}PS_newErr${sj}_${nx_i}_z${z_i}_errPBDT9"
     $command
     command="${code1}computepsfromarray -C ${dirg}grids${sj}_${nx_i}_z${z_i}_errPBDT8_5cubes.fits -S ${dirps}simu_ps_${nx_i}_z${zmean_i}_r.fits  -O ${dirps}PS_newErr${sj}_${nx_i}_z${z_i}_errPBDT8"
     $command
 }
#########################################################################################################
nCase=10

code0="/sps/lsst/dev/rcecile/BAOProgs/SimLSS/exe/"     
code1="/sps/lsst/dev/rcecile/BAOProgs/DirectSimPerso/exe/"     
err2=0.03
# erreur photoZ
err_tot=(0.023 0.021 0.029 0.0198 0.020 0.025 0.0173 0.0191 0.0228)

#nx=(160 300 225)
#nz=(140 125 150)
#zmean :    0.52496901    0.96134880   1.5432835
#zmedian :  0.52200496    0.95816419   1.5342830


cell=(8 8 8)
z=('0.5' '0.9' '1.5')
nx=(120 225 320)
nz=(125 125 175)
#zmean :    0.51408011 0.93618589  1.5776486
#zmedian : 0.51176616  0.93265427  1.5653133
zmean=('0.51' '0.93' '1.57')


#########################################################################################################
for k in 0 1 2 3 4 5 6 7 8
do
    dirg="/sps/lsst/data/rcecile/Planck_BAO_grids/"
    dirps="/sps/lsst/data/rcecile/Planck_BAO_PS/"
    sj=""
    err2=${err_tot[k]}
    
    for i in 0 1 2
    do
	nx_i=$((${nx[i]}))
	nz_i=$((${nz[i]}))
	z_i=${z[i]}
	zmean_i=${zmean[i]}
	cell_i=$((${cell[i]}))
#    do_big 
#    do_ps_base 
	do_ps_err $err2
#    do_ps_perr 
	
    done
    
    dirg="/sps/lsst/data/rcecile/Planck_noBAO_grids/"
    dirps="/sps/lsst/data/rcecile/Planck_noBAO_PS/"
    for i in 0 1 2
    do
	nx_i=$((${nx[i]}))
	nz_i=$((${nz[i]}))
	z_i=${z[i]}
	zmean_i=${zmean[i]}
	cell_i=$((${cell[i]}))
#    do_big_no
    done
    
    j=100
    while [ ${j} -lt ${nCase} ]
    do
	sj=_$j
	for i in 0 1 2
	do
	    nx_i=$((${nx[i]}))
	    nz_i=$((${nz[i]}))
	    z_i=${z[i]}
	    zmean_i=${zmean[i]}
	    cell_i=$((${cell[i]}))
#	do_ps_base
	    do_ps_err $err2
#	do_ps_perr 
	done
	((j++ ))
    done
    
done