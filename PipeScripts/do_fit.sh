#!/bin/bash 
#$ -q long
#$ -l sps=1
#$ -cwd
#$ -o output -j y
#$ -M rcecile@in2p3.fr
#$ -m eas

function do_fit {
     echo "==========================================================================="
    command="${code1}fitkbao -P ${dirp}PS_fromgrids_$1_z$2$3_wngal.txt -O ${dirp}fit_fromgrids_k${maxk}_$1_z$2$4" 
     $command
     
       rm -f *.ppf
}

code1="/sps/lsst/dev/rcecile/BAOProgs/DirectSimPerso/exe/"
dirg="/sps/lsst/data/rcecile/Planck_BAO2_grids/"
dirp="/sps/lsst/data/rcecile/Planck_BAO2_PS/"

Nx_list=(120 225 320)
z_list=("0.5" "0.9" "1.5")
zmean_list=('0.51' '0.93' '1.57')
grid_err_list=("" "_err0.03" "_errP" "_errPBDT9" "_errPBDT8")
ps_err_list=("" "_err0.03" "_errP" "_errPBDT9" "_errPBDT8")


maxk_list=(0.20 0.16 0.14 0.10)
for k in 0 #1 2 3
do
    maxk="${maxk_list[$k]}"
    for i in 0 #1 2
    do
	Nx="${Nx_list[$i]}"
	z="${z_list[$i]}"
	zmean="${zmean_list[$i]}"
	
	for j in 0 #1 2 3 4
	do
	    grid_err="${grid_err_list[$j]}"
	    ps_err="${ps_err_list[$j]}"
	    do_fit $Nx $z $grid_err $ps_err
	done	
    done
done

