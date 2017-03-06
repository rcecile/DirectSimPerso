#!/bin/bash 
#$ -l sps=1
#$ -o /sps/lsst/data/rcecile/Planck_BAO_grids/output_grid_sigma
#$ -M rcecile@in2p3.fr
#$ -m eas
#$ -q long
#! /usr/local/bin/bash -l

# for noBAO case, just change the dir & dirg paths & nCase



 function do_big_err {
     echo "==========================================================================="
     command="${code1}grid_data -C $4 -s ${dirsf}SelFunc_gaussmean -e $5 -R -P $1,$1,$3,$2,$cell,$cell -z ZS,ZS -O ${dirg}grid_cube$6_$1_z$2_G$5.fits"
     $command
}

dir="/sps/lsst/data/rcecile/Planck_BAO/"
dirg="/sps/lsst/data/rcecile/Planck_BAO_grids/"
nCase=1
 
dirsf="/sps/lsst/data/rcecile/Planck_BAO_grids/"
code1="/sps/lsst/dev/rcecile/BAOProgs/DirectSimPerso/exe/"     

z0=(0.5 0.9 1.5)
Nx_tot_err=(160 300 225)
Nz_tot=(140 125 150)
cell_tot=(8 8 8 8)

i_mid=(16 29 49)

########## with error 0.03

#erreur mediane photoZ
err_tot=(0.030515025 0.028397875 0.037388461)

#erreur mediane photoZ BDT90
err90_tot=(0.025927096 0.027347027 0.032890236)

#erreur mediane photoZ BDT80
err80_tot=(0.023046306 0.026009080 0.030297966)

i_min=(3 8 20)
i_max=(58 94 99)

j=0
while [ ${j} -lt  ${nCase} ]
do
    if [ ${nCase} -le  1 ]
    then
	sj=""
    else
	sj=_$j
    fi
    ################ grid choice ###########
    for i in 0 1 2
    ########################################
    do
	err=${err_tot[i]}
	err90=${err90_tot[i]}
	err80=${err80_tot[i]}
	cell=$((${cell_tot[i]}))
	i0=$((${i_mid[i]}))
	obs_list=${dir}cat_zOrd${sj}_Slice${i0}.fits
	i1=$(($i0 -1 ))
	while [ ${i1} -ge  ${i_min[i]} ]
	do
	    obs=${dir}cat_zOrd${sj}_Slice${i1}.fits
	    obs_list=${obs_list},$obs
	    (( i1-- ))
	done
	
	i1=$(($i0 +1 ))
	while [ ${i1} -le  ${i_max[i]} ]
	do
	    obs=${dir}cat_zOrd${sj}_Slice${i1}.fits
	    obs_list=${obs_list},$obs
	    (( i1++ ))
	done
	simu="do_big_err ${Nx_tot_err[i]}  ${z0[i]} ${Nz_tot[i]} $obs_list $err ${sj}"
	$simu
	simu="do_big_err ${Nx_tot_err[i]}  ${z0[i]} ${Nz_tot[i]} $obs_list ${err90} ${sj}"
	$simu
	simu="do_big_err ${Nx_tot_err[i]}  ${z0[i]} ${Nz_tot[i]} $obs_list ${err80} ${sj}"
	$simu
    done
    ((j++ ))
    
done
