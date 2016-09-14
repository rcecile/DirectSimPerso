#!/bin/bash 
#$ -l sps=1
#$ -o /sps/lsst/data/rcecile/TJP_BAO_grids/output_grid
#$ -M rcecile@in2p3.fr
#$ -m eas
#$ -q huge
#! /usr/local/bin/bash -l

# for noBAO case, just change the dir & dirg paths & nCase


 function do_big_base {
     echo "==========================================================================="
     echo "==========================================================================="
     command="${code1}grid_data -C $4 -s ${dirg}SelFunc__specz -P $1,$1,$3,$2,$cell -z ZS,ZS -O ${dirg}grid_G1$5_$1_z$2.fits"

     echo "I launch " $command
     $command
 }

 function do_big_err {
     echo "==========================================================================="
     echo "==========================================================================="
     command="${code1}grid_data -C $4 -s ${dirg}SelFunc_gaussmean -R -P $1,$1,$3,$2,$cell -z ZG,ZS -O ${dirg}grid_G2$6_$1_z$2_G$5.fits"
     echo "I launch " $command
     $command
}

 function do_big_perr {
     echo "==========================================================================="
     echo "==========================================================================="
     command="${code1}grid_data -C $4 -s ${dirg}SelFunc_photzmean -R -P $1,$1,$3,$2,$cell -z ZP,ZS -O ${dirg}grid_G2$6_$1_z$2_errP.fits"
     echo "I launch " $command
     $command

     command="${code1}grid_data -C $5 -s ${dirg}SelFunc_poddsmean -R -P $1,$1,$3,$2,$cell -z ZP,ZS -O ${dirg}grid_G2$6_$1_z$2_errPodds.fits"
     echo "I launch " $command
     $command
}
 
dir="/sps/lsst/data/rcecile/TJP_BAO/"
dirg="/sps/lsst/data/rcecile/TJP_BAO_grids/"
nCase=1
 
#dir="/sps/lsst/data/rcecile/TJP_noBAO/"
#dirg="/sps/lsst/data/rcecile/TJP_noBAO_grids/"
#nCase=10

code1="/sps/lsst/dev/rcecile/BAOProgs/DirectSimPerso/exe/"     
dirgref="/sps/lsst/data/rcecile/PaperBAO_cell8_grids/"

cell=8

z0=(0.7 1.4 0.65 1.3)
Nx_tot_err=(450 900 432 864)
Nz_tot=(200 200 200 200)

i_mid=(23 47 22 43)


########## without error
i_min=(14 33 12 31)
i_max=(34 63 32 59)
i_min=(14 33 0 31)
i_max=(34 63 40 59)

j=0
while [ ${j} -lt  ${nCase} ]
do
    if [ ${nCase} -le  1 ]
    then
	sj=""
    else
	sj=_$j
    fi
    for i in 2 
    do

	i0=$((${i_mid[i]}))
	obs_list=${dir}cat_Slice${i0}.fits
	i1=$(($i0 -1 ))
	while [ ${i1} -ge  ${i_min[i]} ]
	do
	    obs=${dir}cat_Slice${i1}.fits
	    obs_list=${obs_list},$obs
	    (( i1-- ))
	done
    
	i1=$(($i0 +1 ))
	while [ ${i1} -le  ${i_max[i]} ]
	do
	    obs=${dir}cat_Slice${i1}.fits
	    obs_list=${obs_list},$obs
	    (( i1++ ))
	done
	
	simu="do_big_base ${Nx_tot_err[i]}  ${z0[i]} ${Nz_tot[i]} $obs_list ${sj}"
	$simu
    done
    ((j++ ))

done

########## with error 0.03
err=0.03

i_min=(10 27 8 25)
i_max=(40 72 38 67)

j=0
while [ ${j} -lt  ${nCase} ]
do
    if [ ${nCase} -le  1 ]
    then
	sj=""
    else
	sj=_$j
    fi
    for i in 2
    do
	i0=$((${i_mid[i]}))
	obs_list=${dir}cat_zOrd${sj}_Slice${i0}_err${err}.fits
	obs_list1=${dir}cat_zOrd${sj}_Slice${i0}_errP.fits
	obs_list2=${dir}cat_zOrd${sj}_Slice${i0}_errPodds.fits
	i1=$(($i0 -1 ))
	while [ ${i1} -ge  ${i_min[i]} ]
	do
	    obs=${dir}cat_zOrd${sj}_Slice${i1}_err${err}.fits
	    obs_list=${obs_list},$obs
	    obs=${dir}cat_zOrd${sj}_Slice${i1}_errP.fits
	    obs_list1=${obs_list1},$obs
	    obs=${dir}cat_zOrd${sj}_Slice${i1}_errPodds.fits
	    obs_list2=${obs_list2},$obs
	    (( i1-- ))
	done
	
	i1=$(($i0 +1 ))
	while [ ${i1} -le  ${i_max[i]} ]
	do
	    obs=${dir}cat_zOrd${sj}_Slice${i1}_err${err}.fits
	    obs_list=${obs_list},$obs
	    obs=${dir}cat_zOrd${sj}_Slice${i1}_errP.fits
	    obs_list1=${obs_list1},$obs
	    obs=${dir}cat_zOrd${sj}_Slice${i1}_errPodds.fits
	    obs_list2=${obs_list2},$obs
	    (( i1++ ))
	done
	obs_list=${obs_list},${dir}cat_zOrd${sj}_AllSlice_err${err}.fits
	obs_list1=${obs_list1},${dir}cat_zOrd${sj}_AllSlice_errP.fits
	obs_list2=${obs_list},${dir}cat_zOrd${sj}_AllSlice_errPodds.fits
	
	simu="do_big_err ${Nx_tot_err[i]}  ${z0[i]} ${Nz_tot[i]} $obs_list $err ${sj}"
#    $simu
	
	simu="do_big_perr ${Nx_tot_err[i]} ${z0[i]} ${Nz_tot[i]} $obs_list1 $obs_list2 ${sj}"
#    $simu
    done
    ((j++ ))
	
done
