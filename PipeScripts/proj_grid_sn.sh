#!/bin/bash 
#$ -l sps=1
#$ -o /sps/lsst/data/rcecile/Planck_BAO2_grids/output_sn
#$ -M rcecile@in2p3.fr
#$ -m eas
#$ -q long
#! /usr/local/bin/bash -l


 function do_big_base {
     echo "==========================================================================="
     command="${code1}grid_data -C $4 -s ${dirsf}SelFunc${suff}_speczmean  -P $1,$1,$3,$2,$cell,$cell  -a 1.04720 -S -A 0.,0. -z ZS,ZS -N ${norm} -O ${dirg}gridSN${suff}_$1_z$2_mini"
     echo $command
     $command
 }

dir="/sps/lsst/data/rcecile/Planck_BAO2/"
dirg="/sps/lsst/data/rcecile/Planck_BAO2_grids/"

dirsf="/sps/lsst/data/rcecile/Planck_BAO2_grids/"
code1="/sps/lsst/users/rcecile/BAOProgs/DirectSimPerso/exe/"     

suff="_lfZuccaAllFalse"

z0=(0.5 0.9 1.5 1.3)
# Nx_tot_err=(120 225 320 300)
Nx_tot_err=(60 100 150 300)
Nz_tot=(125 125 175 125)

cell_tot=(8 8 8 8 8)

# coeff Cecile sqrt(fit rapport spectres th / spectres avec Norm=1)
NormNg=(1.30 1.18 1.08 1.10)

i_mid=(16 29 49 43)

########## without error

#i_min=(5 16 29 26)
#i_max=(37 65 99 67)

#i_min=(5 16 29 26)
#i_max=(37 65 99 61)
i_min=(5 16 29 26)
i_max=(37 65 99 96)

    ################ grid choice ###########
for i in 3
    ########################################
do
    cell=$((${cell_tot[i]}))
    norm=${NormNg[i]}

    i0=$((${i_mid[i]}))
    obs_list=${dir}cat${suff}_zOrd_Slice${i0}.fits
    i1=$(($i0 -1 ))
    while [ ${i1} -ge  ${i_min[i]} ]
    do
	obs=${dir}cat${suff}_zOrd_Slice${i1}.fits
#	obs_list=${obs_list},$obs
	(( i1-- ))
    done
    
    i1=$(($i0 +1 ))
    while [ ${i1} -le  ${i_max[i]} ]
    do
	obs=${dir}cat${suff}_zOrd_Slice${i1}.fits
#	obs_list=${obs_list},$obs
	(( i1++ ))
    done
    
    simu="do_big_base ${Nx_tot_err[i]}  ${z0[i]} ${Nz_tot[i]} $obs_list "
    $simu
done
 
