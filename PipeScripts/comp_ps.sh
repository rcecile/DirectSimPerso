#!/bin/bash 
#$ -l sps=1
#$ -o /sps/lsst/data/rcecile/PaperBAO_cell8_PS/output_ps_cell8_G
#$ -M rcecile@in2p3.fr
#$ -m eas
#$ -q huge
#! /usr/local/bin/bash -l

 function do_big {
     echo "==========================================================================="
     command="${code0}cmvginit3d -o ${dirp}simu_ps_${1}_z$4 -U 15 -G 2 -s 1000,0.001,0.5 -x $1,$3 -y $1,$3 -z $2,$3 -Z $4 -O 2,2"
     echo "I launch " $command
     $command
 }

 function do_ps_base {
     echo "==========================================================================="
     command="${code1}computepsfromarray -C ${dir}grid_G2sfmean_${1}_z$2.fits -S ${dirp}simu_ps_${1}_z$2_r.fits  -O ${dirp}PS_G2sfmean_${1}_z$2"
     echo "I launch " $command
     $command
  }

#########################################################################################################
 function do_ps_err {
     echo "==========================================================================="
     command="${code1}computepsfromarray -C ${dir}grid_G2sfmean_${1}_z${2}_G$3.fits -S ${dirp}simu_ps_${1}_z$2_r.fits -O ${dirp}PS_G2sfmean_${1}_z${2}_err${3}"
     echo "I launch " $command
     $command
  }

 function do_ps_perr {
     command="${code1}computepsfromarray -C ${dir}grid_G2sfmean_${1}_z${2}_errP.fits -S ${dirp}simu_ps_${1}_z$2_r.fits  -O ${dirp}PS_G2sfmean$4_${1}_z${2}_errP"
     echo "I launch " $command
     $command
     command="${code1}computepsfromarray -C ${dir}grid_G2sfmean_${1}_z${2}_errPpodds.fits -S ${dirp}simu_ps_${1}_z$2_r.fits  -O ${dirp}PS_G2sfmean$4_${1}_z${2}_errPpodds"
     echo "I launch " $command
     $command
 }
#########################################################################################################

 function do_big_w {
     echo "==========================================================================="    
     command="${code0}cmvginit3d -o ${dirp}simu_ps$8_${1}_z$4_w${6}${7}  -P -t 2 -u 0.6935,0.3065,0.0483,0.679,${6},${7},9.06581172724113E-05 -8 0.8154,-8. -G 2 -s 1000,0.001,0.5 -x $1,$3 -y $1,$3 -z $2,$3 -Z $4 -O 0,2"
     echo "I launch " $command
     $command
 }

 function do_ps_w {
     echo "==========================================================================="
     command="${code1}computepsfromarray -C ${dir}grid_G2sfmean_${1}_z$2.fits -S ${dirp}simu_ps_${1}_z$2_w${3}${4}_r.fits -w $3,$4 -O ${dirp}PS_G2sfmean_${1}_z$2_w$3$4 "
     echo "I launch " $command
     $command
  }

dir="/sps/lsst/data/rcecile/PaperBAO_cell8_grids/"
dirp="/sps/lsst/data/rcecile/PaperBAO_cell8_PS/"
code0="/sps/lsst/dev/rcecile/BAOProgs/SimLSS/exe/"     
code1="/sps/lsst/dev/rcecile/BAOProgs/DirectSim/exe/"     
err2=0.02

cell=(8 8 8)

z=('0.5' '1.0' '1.5' '0.65' '1.3' )
D=(1945 3400 4485 2428 4085)
nz=(300 200 200 200 20)
nx=(300 675 960 432 864)


#########################################################################################################

for i in 0 1 2 3 4
do
 nx_i=$((${nx[i]}))
 nz_i=$((${nz[i]}))
 z_i=${z[i]}
 cell_i=${cell[i]}

# do_big ${nx_i} ${nz_i} ${cell_i} ${z_i} 

# do_ps_base ${nx_i} ${z_i}
 do_ps_err ${nx_i} ${z_i} $err2
# do_ps_perr ${nx_i} ${z_i} 
 
 for w0 in -2 -1.5 -1.0 -0.5
 do
     for wa in -1. -0.5 0. 0.5 1.
     do
	 echo  ${w0} ${wa}
# 	 do_big_w ${nx_i} ${nz_i} ${cell_i} ${z_i} $dirp ${w0} ${wa}
#	 do_ps_w ${nx_i} ${z_i} ${w0} ${wa}
     done
 done
done
