#!/bin/bash


curr_path=$1
run_path=$2
ax_run_path=$3


mkdir ${curr_path}'sampledMSAs'



# do the sims
${run_path}ALF_standalone/bin/alfsim alf-params_lopho.drw
${run_path}ALF_standalone/bin/alfsim alf-params_myria.drw
echo "alf simulation complete"

# make a concatenated MSA of output genes
${curr_path}ConcatenateGenes_MSAs.sh
echo "MSA concatenation complete"

# convert the fasta format to phylip for partitioning (and treebuilding)
source /scratch/wally/FAC/FBM/DBC/cdessim2/default/kgilbert/miniconda3/etc/profile.d/conda.sh
conda activate my_python3
python ${run_path}convert_fa_to_phylipSequential.py
conda deactivate
echo "python conversion to Phylip complete"

# reformat the sequential phylip format to have the correct whitespace separating species from bases
${curr_path}ReformatPhylip.sh $curr_path
echo "reformatting done to Phylip files"


# make all the directories and put the files in the right place for running partitionfinder and iqtree
# MAKE ALSO THE SLURM FILES!

input_folder=$curr_path'sampledMSAs/'
files=${input_folder}*.phy


for f in $files;
    do

    # sample size of genes
    sample_size=$(echo $f | tr "_" "\n" | tail -n 1 | tr "." "\n" | head -n 1)
    #raw numeric sample size for directory label
    raw_samp=$(echo $sample_size | tr "g" "\n" | head -n 1)

    # species name:
    species=$(echo $f | tr "_" "\n" | tail -n 2 | head -n 1)
 
	mkdir ${curr_path}'part_N'$raw_samp'_'$species'_g-a'
	cp $input_folder'MSA_'$species'_'$raw_samp'genes.phy' ${curr_path}'part_N'$raw_samp'_'$species'_g-a/temp_align.phy' 	
	
	mkdir ${curr_path}'part_N'$raw_samp'_'$species'_g-b'
	cp $input_folder'MSA_'$species'_'$raw_samp'genes.phy' ${curr_path}'part_N'$raw_samp'_'$species'_g-b/temp_align.phy' 

	mkdir ${curr_path}'part_N'$raw_samp'_'$species'_g-c'
	cp $input_folder'MSA_'$species'_'$raw_samp'genes.phy' ${curr_path}'part_N'$raw_samp'_'$species'_g-c/temp_align.phy' 

	mkdir ${curr_path}'part_N'$raw_samp'_'$species'_r-a'
	cp $input_folder'MSA_'$species'_'$raw_samp'genes.phy' ${curr_path}'part_N'$raw_samp'_'$species'_r-a/temp_align.phy' 

	mkdir ${curr_path}'part_N'$raw_samp'_'$species'_r-b'
	cp $input_folder'MSA_'$species'_'$raw_samp'genes.phy' ${curr_path}'part_N'$raw_samp'_'$species'_r-b/temp_align.phy' 

	mkdir ${curr_path}'part_N'$raw_samp'_'$species'_r-c'
	cp $input_folder'MSA_'$species'_'$raw_samp'genes.phy' ${curr_path}'part_N'$raw_samp'_'$species'_r-c/temp_align.phy' 

	mkdir ${curr_path}'part_N'$raw_samp'_'$species'_noPart'
	cp $input_folder'MSA_'$species'_'$raw_samp'genes.phy' ${curr_path}'part_N'$raw_samp'_'$species'_noPart/temp_align.phy' 

done


# make config files for partitionfinder
${curr_path}MakeConfigs.sh $curr_path $run_path $ax_run_path
echo "config and slurm files done"

