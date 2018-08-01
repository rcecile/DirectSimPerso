#!/bin/bash 
#$ -q long
#$ -l sps=1
#$ -cwd
#$ -o output -j y
#$ -M rcecile@in2p3.fr
#$ -m eas

#fichier -P : 4 col [1]:k [2]:Pkwosc [3]:Shot-Noise [4]: Error
 function do_fit {
     echo "==========================================================================="
     command="${code1}fitkbaobaseline -P ${dir}PS${suff}_$1_z$2$3_wngal.txt -k ${mink}${maxk} -O ${dirp}fit${suff}_k${mink}_${maxk}_$1_z$2$3" 
     echo "I launch" $command
     $command
     rm -f *.ppf
 }

suff="_lfZuccaAllFalse"
 
code1="/sps/lsst/users/rcecile/BAOProgs/DirectSimPerso/exe/"
dir="/sps/lsst/data/rcecile/Planck_BAO2_PS/"    
dirp="/sps/lsst/data/rcecile/Planck_BAO2_PS/"    

Nx_list=(120 225 320 300)
z_list=("0.5" "0.9" "1.5" "1.3")
zmean_list=('0.51' '0.93' '1.57' '1.35')
grid_err_list=("" "_err0.03" "_errP" "_errPBDT9" "_errPBDT8")
ps_err_list=("" "_err0.03" "_errP" "_errPBDT9" "_errPBDT8")


maxk_list=(0.20 0.15 0.10 0.12)
mink_list=(0.02 0.03 0.04)
for kmax in  0 1 2 3
do
    maxk="${maxk_list[$kmax]}"
    for kmin in 0 1 2
    do
	mink="${mink_list[$kmin]}"
	for i in 3  
	do
	    Nx="${Nx_list[$i]}"
	    z="${z_list[$i]}"
	    zmean="${zmean_list[$i]}"
	    
            for j in 0 1 2 3 4
	    do
		grid_err="${grid_err_list[$j]}"
		ps_err="${ps_err_list[$j]}"
		run="${run_list[$m]}"
		do_fit $Nx $z $ps_err
	    done
	done
    done
done
