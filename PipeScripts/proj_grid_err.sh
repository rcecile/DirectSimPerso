#!/bin/bash 
#$ -l sps=1
#$ -o /sps/lsst/data/rcecile/Planck_BAO2_grids/output_grid_err
#$ -M rcecile@in2p3.fr
#$ -m eas
#$ -q long
#! /usr/local/bin/bash -l


 function do_big_err {
     echo "==========================================================================="
     command="${code1}grid_data -C $4 -s ${dirsf}SelFunc_gaussmean  -R -P $1,$1,$3,$2,$cell,$cell  -E /sps/lsst/data/rcecile/Planck_BAO2_grids/Sigma_Gauss0.00.txt -a 1.04720 -S -A 0.,0. -A 40.,0. -A 40.,90. -A 40.,180. -A 40.,270. -z ZG,ZS -N ${norm} -O ${dirg}grids_superSN5_$1_z$2_err$5_00_$6"
     $command
     command="${code1}grid_data -C $4 -s ${dirsf}SelFunc_gaussmean  -R -P $1,$1,$3,$2,$cell,$cell  -E /sps/lsst/data/rcecile/Planck_BAO2_grids/Sigma_Gauss0.02.txt -a 1.04720 -S -A 0.,0. -A 40.,0. -A 40.,90. -A 40.,180. -A 40.,270. -z ZG,ZS -N ${norm} -O ${dirg}grids_superSN5_$1_z$2_err$5_02_$6"
     $command
     command="${code1}grid_data -C $4 -s ${dirsf}SelFunc_gaussmean  -R -P $1,$1,$3,$2,$cell,$cell  -E /sps/lsst/data/rcecile/Planck_BAO2_grids/Sigma_Gauss0.03.txt -a 1.04720 -S -A 0.,0. -A 40.,0. -A 40.,90. -A 40.,180. -A 40.,270. -z ZG,ZS -N ${norm} -O ${dirg}grids_superSN5_$1_z$2_err$5_03_$6"
     $command
     command="${code1}grid_data -C $4 -s ${dirsf}SelFunc_gaussmean  -R -P $1,$1,$3,$2,$cell,$cell  -E /sps/lsst/data/rcecile/Planck_BAO2_grids/Sigma_Gauss0.04.txt -a 1.04720 -S -A 0.,0. -A 40.,0. -A 40.,90. -A 40.,180. -A 40.,270. -z ZG,ZS -N ${norm} -O ${dirg}grids_superSN5_$1_z$2_err$5_04_$6"
     $command
}


dir="/sps/lsst/data/rcecile/Planck_BAO2/"
dirg="/sps/lsst/data/rcecile/Planck_BAO2_grids/"

dirsf="/sps/lsst/data/rcecile/Planck_BAO2_grids/"
code1="/sps/lsst/dev/rcecile/BAOProgs/DirectSimPerso/exe/"     

z0=(0.5 0.9 1.5)
Nx_tot_err=(120 225 320)
Nz_tot=(125 125 175)

cell_tot=(8 8 8 8)

# coeff Cecile sqrt(fit rapport spectres th / spectres avec Norm=1)
NormNg=(1.30 1.18 1.08)

i_mid=(16 29 49)


########## with error 0.03 or photoZ
err=0.03

i_min=(3 8 19)
i_max=(48 80 99)

    ################ grid choice ###########
for i in 1
    ########################################
do
    cell=$((${cell_tot[i]}))
    norm=${NormNg[i]}

    i0=$((${i_mid[i]}))
    obs_list=${dir}cat_AllzOrd_Slice${i0}_err${err}.fits
    i1=$(($i0 -1 ))
    while [ ${i1} -ge  ${i_min[i]} ]
    do
	obs=${dir}cat_AllzOrd_Slice${i1}_err${err}.fits
	obs_list=${obs_list},$obs
	(( i1-- ))
    done
    
    i1=$(($i0 +1 ))
    while [ ${i1} -le  ${i_max[i]} ]
    do
  	obs=${dir}cat_AllzOrd_Slice${i1}_err${err}.fits
  	obs_list=${obs_list},$obs
 	(( i1++ ))
    done
    obs_list=${obs_list},${dir}cat_AllSlice_err${err}_cataLT10.fits
    
    for isim in 0 1 2 3 4 5 6 7 8 9
    do
	simu="do_big_err ${Nx_tot_err[i]}  ${z0[i]} ${Nz_tot[i]} $obs_list $err $isim"
	$simu
    done

   done
