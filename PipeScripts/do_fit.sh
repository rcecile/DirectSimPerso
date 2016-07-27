#!/bin/bash 
#$ -q huge
#$ -l sps=1
#$ -cwd
#$ -o output -j y
#$ -M rcecile@in2p3.fr
#$ -m eas

 function do_fit {
     echo "==========================================================================="
     command="${code1}fitkbao -P ${dirp}PS_G2sfmean2fitErr_$1_z$2$4_wngal.txt   -U ${dirg}grid_G2sfmean_$1_z$2.fits -R ${dirp}simu_ps_$1_z$2_ntpk.txt -O ${dirp}fit_G2sfmeanErr_$1_z$2$4" 
     echo "I launch" $command
     $command
     
     rm -f *.ppf
 }
#-U : cosmological parameters
#-P : power spectrum to fit
#-R : theoretical case compared to the previous spectrum

code1="/sps/lsst/dev/rcecile/BAOProgs/DirectSim/exe/"
#code1="/sps/lsst/dev/mmoneuse/BAOProgs/DirectSim/exe/"
dirg="/sps/lsst/data/rcecile/PaperBAO_cell8_grids/"
dirp="/sps/lsst/data/rcecile/PaperBAO_cell8_PS/"

Nx_list=(300 675 960 432 864)
z_list=("0.5" "1.0" "1.5" "0.65" "1.3")
grid_err_list=("" "_G0.03" "_errPpodds")
ps_err_list=("" "_err0.03" "_errPpodds")

for i in 0 1 2 3 4
do
    Nx="${Nx_list[$i]}"
    z="${z_list[$i]}"
    for j in 0 1 2
    do
	grid_err="${grid_err_list[$j]}"
	ps_err="${ps_err_list[$j]}"
	do_fit $Nx $z $grid_err $ps_err
    done

done
