#!/bin/bash

run_path='/scratch/wally/FAC/FBM/DBC/cdessim2/default/kgilbert/'
curr_path='/scratch/wally/FAC/FBM/DBC/cdessim2/default/kgilbert/Genes_Length100/'

mkdir ${curr_path}'sampledMSAs'

mkdir ${curr_path}'part_N10'
mkdir ${curr_path}'part_N20'
mkdir ${curr_path}'part_N40'
mkdir ${curr_path}'part_N80'


# do the sims
#${run_path}ALF_standalone/bin/alfsim alf-params_lopho.drw
#${run_path}ALF_standalone/bin/alfsim alf-params_myria.drw
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
${curr_path}ReformatPhylip.sh $path
echo "reformatting done to phylip files"

# partition and treebuild
conda activate my_python2
${curr_path}PartitionMSAs.sh $curr_path $run_path
conda deactivate
echo "full analysis complete"
