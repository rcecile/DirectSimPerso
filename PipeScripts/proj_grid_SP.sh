#!/bin/bash 
#$ -l sps=1
#$ -o /sps/lsst/data/rcecile/Planck_BAO2_SP/output_grid_5cubes_SNSP
#$ -M rcecile@in2p3.fr
#$ -m eas
#$ -q long
#! /usr/local/bin/bash -l


 function do_big_base {
     echo "==========================================================================="
     command="${code1}grid_data -C $4 -s ${dirsf}SelFunc_speczmean  -P $1,$1,$3,$2,$cell,$cell  -a 1.04720 -S -A 0.,0. -A 40.,0. -A 40.,90. -A 40.,180. -A 40.,270. -z ZS,ZS -N ${norm} -O ${dirg}grids_nCR$5_$1_z$2"
     $command
 }


dir="/sps/lsst/data/rcecile/Planck_BAO2_SP/"
dirg="/sps/lsst/data/rcecile/Planck_BAO2_SP/"

dirsf="/sps/lsst/data/rcecile/Planck_BAO2_grids/"
code1="/sps/lsst/dev/rcecile/BAOProgs/DirectSimPerso/exe/"     

z0=(0.5 0.9 1.5)
Nx_tot_err=(120 225 320)
Nz_tot=(125 125 175)

cell_tot=(8 8 8 8)

#coef fournis par Reza, prends en compte la mise à 0 des pixels négatifs avant conversion en galaxies
#NormNg=(1.13 1.08 1.034)
# coeff Cecile sqrt(fit rapport spectres th / spectres avec Norm=1)
NormNg=(1.30 1.18 1.08)

nCase=10

i_mid=(16 29 49)

########## without error

i_min=(5 16 26)
i_max=(39 61 92)

j=0
while [ ${j} -lt  ${nCase} ]
do  
    ################ grid choice ###########
    for i in 0 1 2
    ########################################
    do
	cell=$((${cell_tot[i]}))
	norm=${NormNg[i]}
	
	i0=$((${i_mid[i]}))
	obs_list=${dir}cat${j}_Slice${i0}.fits
	i1=$(($i0 -1 ))
	while [ ${i1} -ge  ${i_min[i]} ]
	do
	    obs=${dir}cat${j}_Slice${i1}.fits
	    obs_list=${obs_list},$obs
	    (( i1-- ))
	done
	
	i1=$(($i0 +1 ))
	while [ ${i1} -le  ${i_max[i]} ]
	do
	    obs=${dir}cat${j}_Slice${i1}.fits
	    obs_list=${obs_list},$obs
	    (( i1++ ))
	done
	
	simu="do_big_base ${Nx_tot_err[i]}  ${z0[i]} ${Nz_tot[i]} $obs_list $j"
	$simu
    done
    ((j++ ))
done
