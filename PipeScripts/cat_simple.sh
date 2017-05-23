bg
#!/bin/bash 
#$ -l sps=1
#$ -o /sps/lsst/data/rcecile/Planck_BAO2/output_cat_All
#$ -M rcecile@in2p3.fr
#$ -q long
#$ -m eas
#! /usr/local/bin/bash -l


function do_rdlss {
    echo "==========================================================================="
    echo "==========================================================================="
    # R = randomize in cell; G = goldencut; a = opening angle cut ( 1.04720 = 60 deg)
    command="${code1}rdlss -C ${dir}simu_$1_r.fits -O ${dir2}cat${suff}_$1.fits -i 1 -R -G -a 1.04720 "
    echo $command
    $command
}


function do_sel_simple {
    echo "==========================================================================="

    command="${code1}getsf -F $2 -O $1 -o ${dirg}SelFunc${suff}_$4 -z $3"
    echo "I launch " $command
    rm -f tmptmp
    $command

}
suff="_All"
code1="/sps/lsst/dev/rcecile/BAOProgs/DirectSimPerso/exe/"
Nslice=70

dir="/sps/lsst/data/rcecile/Planck_BAO/"
dir2="/sps/lsst/data/rcecile/Planck_BAO2/"
dirg="/sps/lsst/data/rcecile/Planck_BAO2_grids/"
nCase=1

i0=0
while [ ${i0} -lt  ${Nslice} ]
do
    simu="do_rdlss Slice$i0 "
    $simu
    ((i0++ ))
    
done

#########################################################################################################################
 
i0=0
i_max=$((${Nslice} -1 ))
obs_list5=${dir2}cat${suff}_Slice${i0}.fits
full_list=${dir2}/cat${suff}_Slice${i0}_ZONLY.fits

i1=$((${i0} + 1))

while [ ${i1} -le  ${i_max} ]
do
    obs=${dir2}cat${suff}_Slice${i1}.fits    
    full=${dir2}/cat${suff}_Slice${i1}_ZONLY.fits
    
    obs_list5=${obs_list5},$obs
    full_list=${full_list},$full
    
    (( i1++ ))
done

do_sel_simple $obs_list5 $full_list ZS,ZS,zs 
