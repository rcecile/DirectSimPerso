#!/bin/bash 
#$ -l sps=1
#$ -o /sps/lsst/data/rcecile/TJP_BAO_PS/output_ps
#$ -M rcecile@in2p3.fr
#$ -m eas
#$ -q huge
#! /usr/local/bin/bash -l

 function do_big {
     echo "==========================================================================="
     command="${code0}cmvginit3d -o ${dirps}simu_ps_${nx_i}_z${z_i} -u 0.7,0.3,0.05,0.7,-1.,0.,9.06581172724113E-05 -8 0.8154,-8. -G 2 -s 1000,0.001,0.5 -x ${nx_i},${cell_i} -y ${nx_i},${cell_i} -z ${nz_i},${cell_i} -Z ${z_i} -O 2,2"
    $command
 }

 function do_big_no {
     echo "==========================================================================="
     command="${code0}cmvginit3d -o ${dirps}simu_ps_${nx_i}_z${z_i} -u 0.7,0.3,0.05,0.7,-1.,0.,9.06581172724113E-05 -8 0.8154,-8. -L -G 2 -s 1000,0.001,0.5 -x ${nx_i},${cell_i} -y ${nx_i},${cell_i} -z ${nz_i},${cell_i} -Z ${z_i} -O 2,2"
     $command
 }

 function do_ps_base {
     echo "==========================================================================="
     command="${code1}computepsfromarray -C ${dirg}grid_G2${sj}_${nx_i}_z${z_i}.fits -S ${dirps}simu_ps_${nx_i}_z${z_i}_r.fits  -O ${dirps}PS_G2${sj}_${nx_i}_z${z_i}"
     $command
  }

#########################################################################################################
 function do_ps_err {
     echo "==========================================================================="
     command="${code1}computepsfromarray -C ${dirg}grid_G2${sj}_${nx_i}_z${z_i}_G$1.fits -S ${dirps}simu_ps_${nx_i}_z${z_i}_r.fits -O ${dirps}PS_G2${sj}_${nx_i}_z${z_i}_err${1}"
     $command
  }

 function do_ps_perr {
     command="${code1}computepsfromarray -C ${dirg}grid_G2${sj}_${nx_i}_z${z_i}_errP.fits -S ${dirps}simu_ps_${nx_i}_z${z_i}_r.fits  -O ${dirps}PS_G2${sj}_${nx_i}_z${z_i}_errP"
     $command
     command="${code1}computepsfromarray -C ${dirg}grid_G2${sj}_${nx_i}_z${z_i}_errPpodds.fits -S ${dirps}simu_ps_${nx_i}_z${z_i}_r.fits  -O ${dirps}PS_G2${sj}_${nx_i}_z${z_i}_errPpodds"
     $command
 }

 function do_ps_perr_k {
     command="${code1}computepsfromarray -C ${dirg}grid_G2${sj}_${nx_i}_z${z_i}_errPpodds.fits -S ${dirps}simu_ps_${nx_i}_z${z_i}_r.fits -m $1  -O ${dirps}PS_G2${sj}_${nx_i}_z${z_i}_k${1}_errPpodds"
     $command
 }
#########################################################################################################
nCase=10

code0="/sps/lsst/dev/rcecile/BAOProgs/SimLSS/exe/"     
code1="/sps/lsst/dev/rcecile/BAOProgs/DirectSimPerso/exe/"     
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
for i in 4
do
    nx_i=$((${nx[i]}))
    nz_i=$((${nz[i]}))
    z_i=${z[i]}
    k_i=$((${nk[i]}))
    cell_i=$((${cell[i]}))
#    do_big 
    do_ps_base 
    do_ps_err $err2
    if [ ${k_i} -lt 0 ]
    then 
	do_ps_perr 
    else
	do_ps_perr_k  ${k_i}
    fi
done

j=0
dirg="/sps/lsst/data/rcecile/TJP_noBAO_grids/"
dirps="/sps/lsst/data/rcecile/TJP_noBAO_PS/"
while [ ${j} -lt ${nCase} ]
do
    sj=_$j
    for i in 4
    do
	nx_i=$((${nx[i]}))
	nz_i=$((${nz[i]}))
	z_i=${z[i]}
	k_i=$((${nk[i]}))
	cell_i=$((${cell[i]}))
	
	if [ ${j} -eq 0 ]
	then 
	    echo do_big_no
#	    do_big_no
	fi

	do_ps_base
	do_ps_err $err2
	if [ ${k_i} -lt 0 ]
	then 
	    do_ps_perr 
	else
	    do_ps_perr_k ${k_i} 
	fi
	
    done
    ((j++ ))
done
