#!/bin/bash


# make all the directories
current_path='/scratch/axiom/FAC/FBM/DBC/cdessim2/default/kgilbert/Genes_Length100/'
##__________________________________##
##__________________ /\  ___________##
##______MODIFY ME___ ||  ___________##
##__________________ ||  ___________##
##__________________________________##


run_path='/scratch/wally/FAC/FBM/DBC/cdessim2/default/kgilbert/'
ax_run_path='/scratch/axiom/FAC/FBM/DBC/cdessim2/default/kgilbert/'

${current_path}SetupAndSimulate.sh $current_path $run_path $ax_run_path

# simulations are done, input files are reformatted and put into folders, 
#     now run each set with partitionfinder and IQtree


######### DO ALL THE PARTITIONING AND IQTREE analyses

input_folder=$current_path'sampledMSAs/'
files=${input_folder}*.phy

for f in $files;
    do
    # sample size of genes
    sample_size=$(echo $f | tr "_" "\n" | tail -n 1 | tr "." "\n" | head -n 1)
    #raw numeric sample size for directory label
    raw_samp=$(echo $sample_size | tr "g" "\n" | head -n 1)
    # species:
    species=$(echo $f | tr "_" "\n" | tail -n 2 | head -n 1)
    
    sbatch ${current_path}'part_N'$raw_samp'_'$species'_r-a/Run.slurm'
	sbatch ${current_path}'part_N'$raw_samp'_'$species'_r-b/Run.slurm'
	sbatch ${current_path}'part_N'$raw_samp'_'$species'_r-c/Run.slurm'
	sbatch ${current_path}'part_N'$raw_samp'_'$species'_g-a/Run.slurm'
	sbatch ${current_path}'part_N'$raw_samp'_'$species'_g-b/Run.slurm'
	sbatch ${current_path}'part_N'$raw_samp'_'$species'_g-c/Run.slurm'
	sbatch ${current_path}'part_N'$raw_samp'_'$species'_noPart/Run.slurm'

done

