#!/bin/bash 
#$ -l sps=1
#$ -o /sps/lsst/data/rcecile/Planck_BAO_grids/output_grid_5cubes_spectro
#$ -M rcecile@in2p3.fr
#$ -m eas
#$ -q long
#! /usr/local/bin/bash -l

# for noBAO case, just change the dir & dirg paths & nCase


 function do_big_err {
     echo "==========================================================================="
     command="${code1}grid_data -C $4 -s ${dirsf}SelFunc_gaussmean -R -P $1,$1,$3,$2,$cell,$cell -A 0.,0. -A 40.,0. -A 40.,90. -A 40.,180. -A 40.,270. -z ZS,ZS -e $5 -O ${dirg}grids$6_$1_z$2_err$5"
     $command
}


dir="/sps/lsst/data/rcecile/Planck_BAO/"
dirg="/sps/lsst/data/rcecile/Planck_BAO_grids/"
 
dirn="/sps/lsst/data/rcecile/Planck_noBAO/"
dirgn="/sps/lsst/data/rcecile/Planck_noBAO_grids/"
nCase=10

dirsf="/sps/lsst/data/rcecile/Planck_BAO_grids/"
code1="/sps/lsst/dev/rcecile/BAOProgs/DirectSimPerso/exe/"     

z0=(0.5 0.9 1.5)
Nx_tot_err=(120 225 320)
Nz_tot=(125 125 175)

cell_tot=(8 8 8 8)

i_mid=(16 29 49)

########## computed in flight

# erreur photoZ
#err_tot=(0.023 0.021 0.029)

# erreur photoZ BDT90
#err_tot=(0.0198 0.020 0.025)

# erreur photoZ BDT80
err_tot=(0.0173 0.0191 0.0228)

i_min=(3 8 19)
i_max=(48 80 99)

################ grid choice ###########
for i in 0 1 2
########################################
do
    err=${err_tot[i]}
    cell=$((${cell_tot[i]}))
    
    i0=$((${i_mid[i]}))
    obs_list=${dir}cat_zOrd_Slice${i0}.fits
    i1=$(($i0 -1 ))
    while [ ${i1} -ge  ${i_min[i]} ]
    do
	obs=${dir}cat_zOrd_Slice${i1}.fits
	obs_list=${obs_list},$obs
	(( i1-- ))
    done
    
    i1=$(($i0 +1 ))
    while [ ${i1} -le  ${i_max[i]} ]
    do
	obs=${dir}cat_zOrd_Slice${i1}.fits
	obs_list=${obs_list},$obs
	(( i1++ ))
    done
    
    simu="do_big_err ${Nx_tot_err[i]}  ${z0[i]} ${Nz_tot[i]} $obs_list $err "
    $simu
done

j=0
while [ ${j} -lt  ${nCase} ]
do
    sj=_$j
    ################ grid choice ###########
    for i in 0 1 2
    ########################################
    do
	err=${err_tot[i]}
	cell=$((${cell_tot[i]}))

	i0=$((${i_mid[i]}))
	obs_list=${dirn}cat_zOrd${sj}_Slice${i0}.fits
	i1=$(($i0 -1 ))
	while [ ${i1} -ge  ${i_min[i]} ]
	do
	    obs=${dirn}cat_zOrd${sj}_Slice${i1}.fits
	    obs_list=${obs_list},$obs
	    (( i1-- ))
	done
    
	i1=$(($i0 +1 ))
	while [ ${i1} -le  ${i_max[i]} ]
	do
	    obs=${dirn}cat_zOrd${sj}_Slice${i1}.fits
	    obs_list=${obs_list},$obs
	    (( i1++ ))
	done
	
	simu="do_big_err ${Nx_tot_err[i]}  ${z0[i]} ${Nz_tot[i]} $obs_list $err ${sj}"
	$simu
    done
    ((j++ ))

done
