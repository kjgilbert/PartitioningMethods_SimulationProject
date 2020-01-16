#!/bin/bash


current_path='/scratch/axiom/FAC/FBM/DBC/cdessim2/default/kgilbert/Rep1_genes100/'
##__________________________________##
##__________________ /\  ___________##
##______MODIFY ME___ ||  ___________##
##__________________ ||  ___________##
##__________________________________##

run_path='/scratch/axiom/FAC/FBM/DBC/cdessim2/default/kgilbert/'
outdir='AnalysisResults/'
mkdir $outdir


# get all inferred tree outputs into a good format for comparison and summarizing


module load Bioinformatics/Software/vital-it
module load R/latest

Rscript ExtractOutputs_InDirectories_Partitioning.R $current_path






##__________________________________##
# get the real trees into a format for comparison to the inferred trees


# 62 lopho species
raw_real_lopho="results_lopho/sim_match_lophotroch/RealTree.nwk"
# 25 myria species
raw_real_myria="results_myria/sim_match_myriapod/RealTree.nwk"


# lophotrochozoa
for i in `seq 1 62`;
	do
	if (( $i < 10 ))
	then
	    num='0'$i
	    new_ID='S00'$i'_00001'
	else
		num=$i
	    new_ID='S0'$i'_00001'
	fi
	curr_ID='SE0'$num
	sed -i -- "s/$curr_ID/$new_ID/g" $raw_real_lopho
done
cp $raw_real_lopho $outdir'RealTree_lopho.nwk'


# myriapoda

for i in `seq 1 25`;
	do
	if (( $i < 10 ))
	then
	    num='0'$i
	    new_ID='S00'$i'_00001'
	else
		num=$i
	    new_ID='S0'$i'_00001'
	fi
	curr_ID='SE0'$num
	sed -i -- "s/$curr_ID/$new_ID/g" $raw_real_myria
done
cp $raw_real_myria $outdir'RealTree_myria.nwk'




##__________________________________##
# use python to get RF distance and Euclidean distance between trees

lophoFiles=${outdir}maxtree_lopho_*.nwk
myriaFiles=${outdir}maxtree_myria_*.nwk

## CHECK THAT THEY ARE READ RIGHT INTO THE PYTHON SCRIPT - quotes and final directory not being doubled

echo "filename rf_dist eucl_dist" > RF_Eucl_results.txt

source /scratch/wally/FAC/FBM/DBC/cdessim2/default/kgilbert/miniconda3/etc/profile.d/conda.sh
conda activate my_python3

# lopho
for f in $lophoFiles;
    do
    python ${run_path}OutputAnalysis_RF-Euc.py "${outdir}RealTree_lopho.nwk" "${outdir}${f}" >> RF_Eucl_results.txt
done

# myria
for f in $myriaFiles;
    do
    python ${run_path}OutputAnalysis_RF-Euc.py "${outdir}RealTree_myria.nwk" "${outdir}${f}" >> RF_Eucl_results.txt
done

conda deactivate
