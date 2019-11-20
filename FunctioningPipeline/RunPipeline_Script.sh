#!/bin/bash

run_path='/scratch/wally/FAC/FBM/DBC/cdessim2/default/kgilbert/'
curr_path='/scratch/wally/FAC/FBM/DBC/cdessim2/default/kgilbert/Genes_Length100/'

mkdir ${curr_path}'sampledMSAs'

mkdir ${curr_path}'part_N10'
mkdir ${curr_path}'part_N20'
mkdir ${curr_path}'part_N40'
mkdir ${curr_path}'part_N80'


# do the sims
${path}ALF_standalone/bin/alfsim alf-params_lopho.drw
${path}ALF_standalone/bin/alfsim alf-params_myria.drw

# make a concatenated MSA of output genes
${curr_path}ConcatenateGenes_MSAs.sh

# convert the fasta format to phylip for partitioning (and treebuilding)
source /scratch/wally/FAC/FBM/DBC/cdessim2/default/kgilbert/miniconda3/etc/profile.d/conda.sh
conda activate my_python3
python ${path}convert_fa_to_phylipSequential.py
conda deactivate

# reformat the sequential phylip format to have the correct whitespace separating species from bases
${curr_path}ReformatPhylip.sh $path

# partition and treebuild
conda activate my_python2
${curr_path}PartitionMSAs.sh $curr_path $run_path
conda deactivate
