#!/bin/bash


# do the sims

alfsim alf-params_lopho.drw

alfsim alf-params_myria.drw

# make a concatenated MSA of output genes
sh ./ConcatenateGenes_MSAs.sh

# convert the fasta format to phylip for partitioning (and treebuilding)

python3.7 convert_fa_to_phylipSequential.py

# reformat the sequential phylip format to have the correct whitespace separating species from bases

sh ./ReformatPhylip.sh

# partition and treebuild

sh ./PartitionMSAs.sh