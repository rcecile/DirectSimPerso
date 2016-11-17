#!/bin/bash 
#$ -l sps=1
#$ -o /sps/lsst/data/rcecile/TJP_BAO_grids/output_grid
#$ -M rcecile@in2p3.fr
#$ -m eas
#$ -q long
#! /usr/local/bin/bash -l

# for noBAO case, just change the dir & dirg paths & nCase


 function do_big_base {
     echo "==========================================================================="
     echo "==========================================================================="
#     command="${code1}grid_data -C $4 -s ${dirsf}SelFunc_speczmean -P $1,$1,$3,$2,$cell,$cell -z ZS,ZS -O ${dirg}grid_specz$5_$1_z$2.fits"
     command="${code1}grid_data -C $4 -s ${dirsf}SelFunc_gaussmean -P $1,$1,$3,$2,$cell,$cell -z ZS,ZS -O ${dirg}grid_gauss$5_$1_z$2.fits"

     echo "I launch " $command
     $command
 }

 function do_big_perr {
     echo "==========================================================================="
     echo "==========================================================================="
#     command="${code1}grid_data -C $4 -s ${dirsf}SelFunc_photzmean -R -P $1,$1,$3,$2,$cell,$cell -z ZP,ZS -O ${dirg}grid_G2$6_$1_z$2_errP.fits"
     command="${code1}grid_data -C $4 -s ${dirsf}SelFunc_gaussmean -R -P $1,$1,$3,$2,$cell,$cell -z ZP,ZS -O ${dirg}grid_gauss$6_$1_z$2_errP.fits"
     $command

#     command="${code1}grid_data -C $5 -s ${dirsf}SelFunc_poddsmean -R -P $1,$1,$3,$2,$cell,$cell -z ZP,ZS -O ${dirg}grid_G2$6_$1_z$2_errPpodds.fits"
     command="${code1}grid_data -C $5 -s ${dirsf}SelFunc_gaussmean -R -P $1,$1,$3,$2,$cell,$cell -z ZP,ZS -O ${dirg}grid_gauss$6_$1_z$2_errPpodds.fits"
     $command
}
 
dir="/sps/lsst/data/rcecile/TJP_BAO/"
dirg="/sps/lsst/data/rcecile/TJP_BAO_grids/"
nCase=1
 
#dir="/sps/lsst/data/rcecile/TJP_noBAO/"
#dirg="/sps/lsst/data/rcecile/TJP_noBAO_grids/"
#nCase=10

dirsf="/sps/lsst/data/rcecile/TJP_BAO_grids/"
code1="/sps/lsst/dev/rcecile/BAOProgs/DirectSimPerso/exe/"     
dirgref="/sps/lsst/data/rcecile/PaperBAO_cell8_grids/"

z0=(0.5 0.9 1.3 1.8 1.8)
Nx_tot_err=(350 640 900 1024 500)
Nz_tot=(140 125 75 65 75)
cell_tot=(8 8 8 8 16)

i_mid=(18 33 44 55 55)


########## without error
i_min=(11 26 40 52 48)
i_max=(26 40 49 60 64)

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
    for i in 0
    ########################################
    do
	cell=$((${cell_tot[i]}))

	i0=$((${i_mid[i]}))
	obs_list=${dir}cat_G2${sj}_Slice${i0}.fits
	i1=$(($i0 -1 ))
	while [ ${i1} -ge  ${i_min[i]} ]
	do
	    obs=${dir}cat_G2${sj}_Slice${i1}.fits
	    obs_list=${obs_list},$obs
	    (( i1-- ))
	done
    
	i1=$(($i0 +1 ))
	while [ ${i1} -le  ${i_max[i]} ]
	do
	    obs=${dir}cat_G2${sj}_Slice${i1}.fits
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

i_min=(6 21 35 47 42)
i_max=(31 45 54 64 69)

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
    for i in 0 
    ########################################
    do
	cell=$((${cell_tot[i]}))
	i0=$((${i_mid[i]}))
	obs_list=${dir}cat_G2${sj}_Slice${i0}_err${err}.fits
	obs_list1=${dir}cat_G2${sj}_Slice${i0}_errP.fits
	obs_list2=${dir}cat_G2${sj}_Slice${i0}_errPpodds.fits
	i1=$(($i0 -1 ))
	while [ ${i1} -ge  ${i_min[i]} ]
	do
	    obs=${dir}cat_G2${sj}_Slice${i1}_err${err}.fits
	    obs_list=${obs_list},$obs
	    obs=${dir}cat_G2${sj}_Slice${i1}_errP.fits
	    obs_list1=${obs_list1},$obs
	    obs=${dir}cat_G2${sj}_Slice${i1}_errPpodds.fits
	    obs_list2=${obs_list2},$obs
	    (( i1-- ))
	done
	
	i1=$(($i0 +1 ))
	while [ ${i1} -le  ${i_max[i]} ]
	do
	    obs=${dir}cat_G2${sj}_Slice${i1}_err${err}.fits
	    obs_list=${obs_list},$obs
	    obs=${dir}cat_G2${sj}_Slice${i1}_errP.fits
	    obs_list1=${obs_list1},$obs
	    obs=${dir}cat_G2${sj}_Slice${i1}_errPpodds.fits
	    obs_list2=${obs_list2},$obs
	    (( i1++ ))
	done
	obs_list=${obs_list},${dir}cat_G2${sj}_AllSlice_err${err}_cataLT10.fits
	obs_list1=${obs_list1},${dir}cat_G2${sj}_AllSlice_errP_cataLT10.fits
	obs_list2=${obs_list2},${dir}cat_G2${sj}_AllSlice_errPpodds_cataLT10.fits
	
	simu="do_big_err ${Nx_tot_err[i]}  ${z0[i]} ${Nz_tot[i]} $obs_list $err ${sj}"
#	$simu
	
	simu="do_big_perr ${Nx_tot_err[i]} ${z0[i]} ${Nz_tot[i]} $obs_list1 $obs_list2 ${sj}"
	$simu
    done
    ((j++ ))
    
done
