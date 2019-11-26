#!/bin/bash


# make all the directories

run_path='/scratch/wally/FAC/FBM/DBC/cdessim2/default/kgilbert/'
curr_path='/scratch/wally/FAC/FBM/DBC/cdessim2/default/kgilbert/Genes_Length100/'


${curr_path}SetupAndSimulate.sh $curr_path $run_path

# simulations are done, input files are reformatted and put into folders, now run each set with partitionfinder and IQtree


######### DO ALL THE PARTITIONING AND IQTREE analyses

pop_size=10000

basedirname=SlimTestJan21_


### s = -0.001

dirname=s0p001_h0p5_rep
for i in `seq $1 $2`;
    do
        fulldirname=${basedirname}${dirname}$i
        sbatch "$fulldirname/Run.slurm"        
done

dirname=s0p001_h0p25_rep
for i in `seq $1 $2`;
    do
        fulldirname=${basedirname}${dirname}$i
        sbatch "$fulldirname/Run.slurm"        
done

dirname=s0p001_h0p0_rep
for i in `seq $1 $2`;
    do
        fulldirname=${basedirname}${dirname}$i
        sbatch "$fulldirname/Run.slurm"        
done



### s = -0.01

dirname=s0p01_h0p5_rep
for i in `seq $1 $2`;
    do
        fulldirname=${basedirname}${dirname}$i
        sbatch "$fulldirname/Run.slurm"        
done

dirname=s0p01_h0p25_rep
for i in `seq $1 $2`;
    do
        fulldirname=${basedirname}${dirname}$i
        sbatch "$fulldirname/Run.slurm"        
done

dirname=s0p01_h0p0_rep
for i in `seq $1 $2`;
    do
        fulldirname=${basedirname}${dirname}$i
        sbatch "$fulldirname/Run.slurm"        
done



### s = -0.1

dirname=s0p1_h0p5_rep
for i in `seq $1 $2`;
    do
        fulldirname=${basedirname}${dirname}$i
        sbatch "$fulldirname/Run.slurm"        
done

dirname=s0p1_h0p25_rep
for i in `seq $1 $2`;
    do
        fulldirname=${basedirname}${dirname}$i
        sbatch "$fulldirname/Run.slurm"        
done

dirname=s0p1_h0p0_rep
for i in `seq $1 $2`;
    do
        fulldirname=${basedirname}${dirname}$i
        sbatch "$fulldirname/Run.slurm"        
done
