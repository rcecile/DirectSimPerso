#!/bin/bash 
#$ -l sps=1
#$ -o /sps/lsst/data/rcecile/TJP_BAO_grids/output_grid
#$ -M rcecile@in2p3.fr
#$ -m eas
#$ -q huge
#! /usr/local/bin/bash -l

 function do_big_base {
     echo "==========================================================================="
     echo "==========================================================================="
     command="${code1}grid_data -C $4 -s ${dirg}SelFunc_speczmean -P $1,$1,$3,$2,$cell -z zs,zs -r -O ${dirg}grid_$1_z$2.fits"
     echo "I launch " $command
     $command
 }

 function do_big_err {
     echo "==========================================================================="
     echo "==========================================================================="
     command="${code1}grid_data -C $4 -s ${dirg}SelFunc_gaussmean -R -P $1,$1,$3,$2,$cell -z zG,zs -r -O ${dirg}grid_$1_z$2_G$5.fits"
     echo "I launch " $command
     $command
}

 function do_big_perr {
     echo "==========================================================================="
     echo "==========================================================================="
     command="${code1}grid_data -C $4 -s ${dirg}SelFunc_photzmean -R -P $1,$1,$3,$2,$cell -z zP,zs -r -O ${dirg}grid_$1_z$2_errP.fits"
     echo "I launch " $command
     $command

     command="${code1}grid_data -C $5 -s ${dirg}SelFunc_poddsmean -R -P $1,$1,$3,$2,$cell -z zP,zs -r -O ${dirg}grid_$1_z$2_errPodds.fits"
     echo "I launch " $command
     $command
}
 
code1="/sps/lsst/dev/rcecile/BAOProgs/DirectSimPerso/exe/"     
dir="/sps/lsst/data/rcecile/TJP_BAO/"
dirg="/sps/lsst/data/rcecile/TJP_BAO_grids/"

cell=8

z0=(0.7 1.4)
Nx_tot_err=(450 900)
Nz_tot=(200 200)

i_mid=(23 47)


########## without error
i_min=(14 33)
i_max=(34 63)

for i in 0 
do
    i0=$((${i_mid[i]}))
    obs_list=${dir}cat_zOrd_Slice${i0}.fits
    echo $obs
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
 
    simu="do_big_base ${Nx_tot_err[i]}  ${z0[i]} ${Nz_tot[i]} $obs_list"
#    echo "Simu lancee : " $simu
#    $simu
done

########## with error 0.03
err=0.03

i_min=(10 27)
i_max=(40 72)

for i in 0 1
do
    i0=$((${i_mid[i]}))
    obs_list=${dir}cat_zOrd_Slice${i0}_err${err}.fits
    obs_list1=${dir}cat_zOrd_Slice${i0}_errP.fits
    obs_list2=${dir}cat_zOrd_Slice${i0}_errPodds.fits
    echo $obs
    i1=$(($i0 -1 ))
    while [ ${i1} -ge  ${i_min[i]} ]
    do
	obs=${dir}cat_zOrd_Slice${i1}_err${err}.fits
	obs_list=${obs_list},$obs
	obs=${dir}cat_zOrd_Slice${i1}_errP.fits
	obs_list1=${obs_list1},$obs
	obs=${dir}cat_zOrd_Slice${i1}_errPodds.fits
	obs_list2=${obs_list2},$obs
	(( i1-- ))
    done
    
    i1=$(($i0 +1 ))
    while [ ${i1} -le  ${i_max[i]} ]
    do
	obs=${dir}cat_zOrd_Slice${i1}_err${err}.fits
	obs_list=${obs_list},$obs
	obs=${dir}cat_zOrd_Slice${i1}_errP.fits
	obs_list1=${obs_list1},$obs
	obs=${dir}cat_zOrd_Slice${i1}_errPodds.fits
	obs_list2=${obs_list2},$obs
	(( i1++ ))
    done
    obs_list=${obs_list},${dir}cat_zOrd_AllSlice_err${err}.fits
    obs_list1=${obs_list1},${dir}cat_zOrd_AllSlice_errP.fits
    obs_list2=${obs_list},${dir}cat_zOrd_AllSlice_errPodds.fits

    simu="do_big_err ${Nx_tot_err[i]}  ${z0[i]} ${Nz_tot[i]} $obs_list $err"
#    echo "Simu lancee : " $simu
#    $simu

    simu="do_big_perr ${Nx_tot_err[i]} ${z0[i]} ${Nz_tot[i]} $obs_list1 $obs_list2"
#    echo "Simu lancee : " $simu
#    $simu

done
