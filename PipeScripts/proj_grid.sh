#!/bin/bash 
#$ -l sps=1
#$ -o /sps/lsst/data/rcecile/Planck_BAO2_grids/output_grid_5cubes_Zucca
#$ -M rcecile@in2p3.fr
#$ -m eas
#$ -q long
#! /usr/local/bin/bash -l


 function do_big_base {
     echo "==========================================================================="
     command="${code1}grid_data -C $4 -s ${dirsf}SelFunc${suff}_speczmean  -P $1,$1,$3,$2,$cell,$cell  -a 1.04720 -S -A 0.,0. -A 40.,0. -A 40.,90. -A 40.,180. -A 40.,270. -z ZS,ZS -N ${norm} -O ${dirg}grids${suff}_$1_z$2"
     $command
 }

 function do_big_err {
     echo "==========================================================================="
     command="${code1}grid_data -C $4 -s ${dirsf}SelFunc${suff}_gaussmean  -R -P $1,$1,$3,$2,$cell,$cell  -E /sps/lsst/data/rcecile/Planck_BAO2_grids/Sigma_Gauss0.03.txt -a 1.04720 -S -A 0.,0. -A 40.,0. -A 40.,90. -A 40.,180. -A 40.,270. -z ZG,ZS -N ${norm} -O ${dirg}grids${suff}_$1_z$2_err$5"
     $command
}

 function do_big_perr {
     echo "==========================================================================="
     command="${code1}grid_data -C $4 -s ${dirsf}SelFunc${suff}_photzmean -R -P $1,$1,$3,$2,$cell,$cell   -E /sps/lsst/data/rcecile/Planck_BAO2_grids/Sigma_photoZ.txt -a 1.04720 -S -A 0.,0. -A 40.,0. -A 40.,90. -A 40.,180. -A 40.,270. -z ZP,ZS -N ${norm} -O ${dirg}grids${suff}_$1_z$2_errP"
     $command

     command="${code1}grid_data -C $5 -s ${dirsf}SelFunc${suff}_pbdt9mean  -R -P $1,$1,$3,$2,$cell,$cell  -E /sps/lsst/data/rcecile/Planck_BAO2_grids/Sigma_photoZ_BDT90.txt -a 1.04720 -S -A 0.,0. -A 40.,0. -A 40.,90. -A 40.,180. -A 40.,270. -z ZP,ZS -N ${norm} -O ${dirg}grids${suff}_$1_z$2_errPBDT9"
     $command

     command="${code1}grid_data -C $6 -s ${dirsf}SelFunc${suff}_pbdt8mean  -R -P $1,$1,$3,$2,$cell,$cell  -E /sps/lsst/data/rcecile/Planck_BAO2_grids/Sigma_photoZ_BDT80.txt -a 1.04720 -S -A 0.,0. -A 40.,0. -A 40.,90. -A 40.,180. -A 40.,270. -z ZP,ZS -N ${norm} -O ${dirg}grids${suff}_$1_z$2_errPBDT8"
    $command

}

dir="/sps/lsst/data/rcecile/Planck_BAO2/"
dirg="/sps/lsst/data/rcecile/Planck_BAO2_grids/"

dirsf="/sps/lsst/data/rcecile/Planck_BAO2_grids/"
code1="/sps/lsst/users/rcecile/BAOProgs/DirectSimPerso/exe/"     

suff="_lfZuccaAllFalse"

z0=(0.5 0.9 1.5 1.3)
Nx_tot_err=(120 225 320 300)
Nz_tot=(125 125 175 125)

cell_tot=(8 8 8 8 8)

# coeff Cecile sqrt(fit rapport spectres th / spectres avec Norm=1)
NormNg=(1.30 1.18 1.08 1.10)

i_mid=(16 29 49 43)

########## without error

i_min=(5 16 29 26)
i_max=(37 65 99 67)

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
	obs_list=${obs_list},$obs
	(( i1-- ))
    done
    
    i1=$(($i0 +1 ))
    while [ ${i1} -le  ${i_max[i]} ]
    do
	obs=${dir}cat${suff}_zOrd_Slice${i1}.fits
	obs_list=${obs_list},$obs
	(( i1++ ))
    done
    
    simu="do_big_base ${Nx_tot_err[i]}  ${z0[i]} ${Nz_tot[i]} $obs_list "
    $simu
done
 
########## with error 0.03 or photoZ
err=0.03

i_min=(3 8 19 17)
i_max=(48 80 99 82)

    ################ grid choice ###########
for i in 3
    ########################################
do
    cell=$((${cell_tot[i]}))
    norm=${NormNg[i]}

    i0=$((${i_mid[i]}))
    obs_list=${dir}cat${suff}_zOrd_Slice${i0}_err${err}.fits
    obs_list1=${dir}cat${suff}_zOrd_Slice${i0}_errP.fits
    obs_list2=${dir}cat${suff}_zOrd_Slice${i0}_errPBDT9.fits
    obs_list3=${dir}cat${suff}_zOrd_Slice${i0}_errPBDT8.fits
    i1=$(($i0 -1 ))
    while [ ${i1} -ge  ${i_min[i]} ]
    do
	obs=${dir}cat${suff}_zOrd_Slice${i1}_err${err}.fits
	obs_list=${obs_list},$obs
	obs=${dir}cat${suff}_zOrd_Slice${i1}_errP.fits
	obs_list1=${obs_list1},$obs
	obs=${dir}cat${suff}_zOrd_Slice${i1}_errPBDT9.fits
	obs_list2=${obs_list2},$obs
	obs=${dir}cat${suff}_zOrd_Slice${i1}_errPBDT8.fits
	obs_list3=${obs_list3},$obs
	(( i1-- ))
    done
    
    i1=$(($i0 +1 ))
    while [ ${i1} -le  ${i_max[i]} ]
    do
  	obs=${dir}cat${suff}_zOrd_Slice${i1}_err${err}.fits
  	obs_list=${obs_list},$obs
  	obs=${dir}cat${suff}_zOrd_Slice${i1}_errP.fits
  	obs_list1=${obs_list1},$obs
  	obs=${dir}cat${suff}_zOrd_Slice${i1}_errPBDT9.fits
  	obs_list2=${obs_list2},$obs
  	obs=${dir}cat${suff}_zOrd_Slice${i1}_errPBDT8.fits
  	obs_list3=${obs_list3},$obs
 	(( i1++ ))
    done
    obs_list=${obs_list},${dir}cat${suff}_zOrd_AllSlice_err${err}_cataLT10.fits
    obs_list1=${obs_list1},${dir}cat${suff}_zOrd_AllSlice_errP_cataLT10.fits
    obs_list2=${obs_list2},${dir}cat${suff}_zOrd_AllSlice_errPBDT9_cataLT10.fits
    obs_list3=${obs_list3},${dir}cat${suff}_zOrd_AllSlice_errPBDT8_cataLT10.fits
    
    simu="do_big_err ${Nx_tot_err[i]}  ${z0[i]} ${Nz_tot[i]} $obs_list $err"
    $simu
    
    simu="do_big_perr ${Nx_tot_err[i]} ${z0[i]} ${Nz_tot[i]} $obs_list1 $obs_list2 $obs_list3"
    $simu
done
