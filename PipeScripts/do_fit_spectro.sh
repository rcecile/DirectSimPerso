#!/bin/bash 
#$ -q long
#$ -l sps=1
#$ -cwd
#$ -o output -j y
#$ -M rcecile@in2p3.fr
#$ -m eas

 function do_fit {
     echo "==========================================================================="
     command="${code1}fitkbao -P ${dirp}PS_2fit_spec_$1_z$2$4_wngal.txt -U ${dirg}grids_$1_z$2_5cubes.fits -R ${dirp}simu_ps_$1_z${zmean}_ntpk.txt -k ${maxk} -O ${dirp}fit_spec_k${maxk}_$1_z$2$4" 
     echo "I launch" $command
     $command
     
       rm -f *.ppf
 }

 function do_fit_simple {
     echo "==========================================================================="
     command="${code1}fitkbao -P ${dirp}PS_2fit_spec_$1_z$2$4_wngal.txt -U ${dirg}grids_$1_z$2_5cubes.fits -R ${dirp}simu_ps_$1_z${zmean}_ntpk.txt -O ${dirp}fit_spec_$1_z$2$4" 
     echo "I launch" $command
     $command
     
      rm -f *.ppf
 }
#-U : cosmological parameters
#-P : power spectrum to fit 
#-R : theoretical case compared to the previous spectrum

code1="/sps/lsst/dev/rcecile/BAOProgs/DirectSimPerso/exe/"
dirg="/sps/lsst/data/rcecile/Planck_BAO_grids/"
dirp="/sps/lsst/data/rcecile/Planck_BAO_PS/"

#Nx_list=(160 300 225)
Nx_list=(120 225 320)
z_list=("0.5" "0.9" "1.5")
zmean_list=('0.51' '0.93' '1.57')
grid_err_list=("" "_err0.03" "_errP" "_errPBDT9" "_errPBDT8")
ps_err_list=("" "_err0.03" "_errP" "_errPBDT9" "_errPBDT8")

for i in 0 1 2
do
    Nx="${Nx_list[$i]}"
    z="${z_list[$i]}"
    zmean="${zmean_list[$i]}"
    for j in 0 
    do
	grid_err="${grid_err_list[$j]}"
	ps_err="${ps_err_list[$j]}"
	do_fit_simple $Nx $z $grid_err $ps_err
    done

done

maxk_list=(0.20 0.18 0.16 0.14 0.12 0.10 0.08)
for k in 0 1 2 3 4 5 6
do
    maxk="${maxk_list[$k]}"
    for i in 0 1 2
    do
	Nx="${Nx_list[$i]}"
	z="${z_list[$i]}"
	zmean="${zmean_list[$i]}"
	
	for j in 0
	do
	    grid_err="${grid_err_list[$j]}"
	    ps_err="${ps_err_list[$j]}"
	    do_fit $Nx $z $grid_err $ps_err
	done	
    done
done

