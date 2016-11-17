#!/bin/bash 
#$ -q long
#$ -l sps=1
#$ -cwd
#$ -o output -j y
#$ -M rcecile@in2p3.fr
#$ -m eas

 function do_fit {
     echo "==========================================================================="
     command="${code1}fitkbao -P ${dirp}PS_G2_2fitErr_$1_z$2$4_wngal.txt   -U ${dirg}grid_G2_$1_z$2.fits -R ${dirp}simu_ps_$1_z$2_ntpk.txt -k ${maxk} -O ${dirp}fit_G2_Err_k${maxk}_$1_z$2$4" 
     echo "I launch" $command
     $command
     
       rm -f *.ppf
 }

 function do_fit_simple {
     echo "==========================================================================="
     command="${code1}fitkbao -P ${dirp}PS_G2_2fitErr_$1_z$2$4_wngal.txt   -U ${dirg}grid_G2_$1_z$2.fits -R ${dirp}simu_ps_$1_z$2_ntpk.txt  -O ${dirp}fit_G2_Err_$1_z$2$4" 
     echo "I launch" $command
     $command
     
      rm -f *.ppf
 }
#-U : cosmological parameters
#-P : power spectrum to fit
#-R : theoretical case compared to the previous spectrum

code1="/sps/lsst/dev/rcecile/BAOProgs/DirectSimPerso/exe/"
dirg="/sps/lsst/data/rcecile/TJP_BAO_grids/"
dirp="/sps/lsst/data/rcecile/TJP_BAO_PS/"

Nx_list=(350 640 900 1024 500)
z_list=("0.5" "0.9" "1.3" "1.8" "1.8")
grid_err_list=("" "_G0.03" "_errPpodds" "_errPBDT")
ps_err_list=("" "_err0.03" "_errPpodds" "_errPBDT")

for i in 0 1 2 3
do
    Nx="${Nx_list[$i]}"
    z="${z_list[$i]}"
    for j in 3
    do
	grid_err="${grid_err_list[$j]}"
	ps_err="${ps_err_list[$j]}"
	do_fit_simple $Nx $z $grid_err $ps_err
    done

done

maxk=0.12
for i in 0 1 2 3
do
    Nx="${Nx_list[$i]}"
    z="${z_list[$i]}"
    for j in 3
    do
	grid_err="${grid_err_list[$j]}"
	ps_err="${ps_err_list[$j]}"
	do_fit $Nx $z $grid_err $ps_err
    done

done

