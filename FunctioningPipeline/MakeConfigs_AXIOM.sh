#!/bin/bash

current_path=$1
run_path=$2
ax_run_path=$3
# actually it's only used for iqtree and partitionfinder, so keep name same and modify to correct path:
run_path=$ax_run_path

# put the cfg file in each corresponding directory

aic_string=$'## ALIGNMENT FILE ##
alignment = temp_align.phy;

## BRANCHLENGTHS: linked | unlinked ##
branchlengths = linked;

## MODELS OF EVOLUTION: all | allx | mrbayes | beast | gamma | gammai | <list> ##
models = all;    # allx uses ML instead of empirical estimation of amino acid frequencies

# MODEL SELECCTION: AIC | AICc | BIC #
model_selection = aic;

## DATA BLOCKS: see manual for how to define ##
[data_blocks]
'

aicc_string=$'## ALIGNMENT FILE ##
alignment = temp_align.phy;

## BRANCHLENGTHS: linked | unlinked ##
branchlengths = linked;

## MODELS OF EVOLUTION: all | allx | mrbayes | beast | gamma | gammai | <list> ##
models = all;    # allx uses ML instead of empirical estimation of amino acid frequencies

# MODEL SELECCTION: AIC | AICc | BIC #
model_selection = aicc;

## DATA BLOCKS: see manual for how to define ##
[data_blocks]
'

bic_string=$'## ALIGNMENT FILE ##
alignment = temp_align.phy;

## BRANCHLENGTHS: linked | unlinked ##
branchlengths = linked;

## MODELS OF EVOLUTION: all | allx | mrbayes | beast | gamma | gammai | <list> ##
models = all;    # allx uses ML instead of empirical estimation of amino acid frequencies

# MODEL SELECCTION: AIC | AICc | BIC #
model_selection = bic;

## DATA BLOCKS: see manual for how to define ##
[data_blocks]
'

rcluster_string=$'
## SCHEMES, search: all | user | greedy | rcluster | rclusterf | kmeans ##
[schemes]
search = rcluster;  # only works with raxml command line specification
'

greedy_string=$'
## SCHEMES, search: all | user | greedy | rcluster | rclusterf | kmeans ##
[schemes]
search = greedy;
'


slurm_string_1=$'#!/bin/bash
#SBATCH --mail-user=kimberly.gilbert@unil.ch
#SBATCH --mail-type=fail
#SBATCH --job-name="'
slurm_string_2=$'"
#SBATCH --workdir=.
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=2-00:00
#SBATCH --partition=ax-long
#SBATCH --account=cdessim2_default


source /scratch/axiom/FAC/FBM/DBC/cdessim2/default/kgilbert/miniconda3/etc/profile.d/conda.sh
conda activate /scratch/axiom/FAC/FBM/DBC/cdessim2/default/kgilbert/miniconda3/envs
'

slurm_string_noPart=$'"
#SBATCH --workdir=.
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1-00:00
#SBATCH --partition=ax-long
#SBATCH --account=cdessim2_default
'

#________go through each gene sample size set________#

input_folder=$current_path'sampledMSAs/'
files=${input_folder}*.phy


for f in $files;
    do

    # sample size of genes
    sample_size=$(echo $f | tr "_" "\n" | tail -n 1 | tr "." "\n" | head -n 1)
    #raw numeric sample size for directory label
    raw_samp=$(echo $sample_size | tr "g" "\n" | head -n 1)

    # make the partition file based on the gene sites stored in config_genes__.txt
    # format the gene sites
    # species:
    species=$(echo $f | tr "_" "\n" | tail -n 2 | head -n 1)
    gene_sites_file=$input_folder'config_'$species'_'$sample_size'.txt'
    gene_string=''
    gene_end='0'
    prev_gene_end='0'
    it='0'
    while read line
        do
        it=$(expr $it + 1)
        gene_start=$(expr $gene_end + 1)
        gene_end=$(expr $line - 1 + $prev_gene_end)
        prev_gene_end=$(expr $line - 1 + $prev_gene_end)
        gene_string=$gene_string'gene'$it' = '$gene_start'-'$gene_end$';\n'
    done <$gene_sites_file

    # put all the strings together into the final .cfg file
    echo "$aic_string" "$gene_string" "$rcluster_string" > ${curr_path}'part_N'$raw_samp'_'$species'_r-a/partition_finder.cfg'
    echo "$aic_string" "$gene_string" "$greedy_string" > ${curr_path}'part_N'$raw_samp'_'$species'_g-a/partition_finder.cfg'
    echo "$bic_string" "$gene_string" "$rcluster_string" > ${curr_path}'part_N'$raw_samp'_'$species'_r-b/partition_finder.cfg'
    echo "$bic_string" "$gene_string" "$greedy_string" > ${curr_path}'part_N'$raw_samp'_'$species'_g-b/partition_finder.cfg'
    echo "$aicc_string" "$gene_string" "$rcluster_string" > ${curr_path}'part_N'$raw_samp'_'$species'_r-c/partition_finder.cfg'
    echo "$aicc_string" "$gene_string" "$greedy_string" > ${curr_path}'part_N'$raw_samp'_'$species'_g-c/partition_finder.cfg'

	# make slurm run files for each as well

	echo "$slurm_string_1"$raw_samp"_"$species'_r-a'"$slurm_string_2" > ${curr_path}'part_N'$raw_samp'_'$species'_r-a/Run.slurm'
	echo "$slurm_string_1"$raw_samp"_"$species'_g-a'"$slurm_string_2" > ${curr_path}'part_N'$raw_samp'_'$species'_g-a/Run.slurm'
	echo "$slurm_string_1"$raw_samp"_"$species'_r-b'"$slurm_string_2" > ${curr_path}'part_N'$raw_samp'_'$species'_r-b/Run.slurm'
	echo "$slurm_string_1"$raw_samp"_"$species'_g-b'"$slurm_string_2" > ${curr_path}'part_N'$raw_samp'_'$species'_g-b/Run.slurm'
	echo "$slurm_string_1"$raw_samp"_"$species'_r-c'"$slurm_string_2" > ${curr_path}'part_N'$raw_samp'_'$species'_r-c/Run.slurm'
	echo "$slurm_string_1"$raw_samp"_"$species'_g-c'"$slurm_string_2" > ${curr_path}'part_N'$raw_samp'_'$species'_g-c/Run.slurm'
	echo "$slurm_string_1"$raw_samp"_"$species'_noPart'"$slurm_string_noPart" > ${curr_path}'part_N'$raw_samp'_'$species'_noPart/Run.slurm'

	# make output file in each directory

	touch ${curr_path}'part_N'$raw_samp'_'$species'_r-a/partitioning_treebuilding_output.txt'
	touch ${curr_path}'part_N'$raw_samp'_'$species'_g-a/partitioning_treebuilding_output.txt'
	touch ${curr_path}'part_N'$raw_samp'_'$species'_r-b/partitioning_treebuilding_output.txt'
	touch ${curr_path}'part_N'$raw_samp'_'$species'_g-b/partitioning_treebuilding_output.txt'
	touch ${curr_path}'part_N'$raw_samp'_'$species'_r-c/partitioning_treebuilding_output.txt'
	touch ${curr_path}'part_N'$raw_samp'_'$species'_g-c/partitioning_treebuilding_output.txt'
	touch ${curr_path}'part_N'$raw_samp'_'$species'_noPart/partitioning_treebuilding_output.txt'


	# raxml AIC
	echo "python ${run_path}partitionfinder-2.1.1/PartitionFinderProtein.py ${curr_path}part_N${raw_samp}_${species}_r-a/ --raxml -p -1
conda deactivate
# copy nexus best scheme to a partition file
awk '/nexus/,/end/' ${curr_path}part_N${raw_samp}_${species}_r-a/analysis/best_scheme.txt > ${curr_path}part_N${raw_samp}_${species}_r-a/best_scheme_nexus.txt
# save information criteria output from partition finder
echo '${species}_${raw_samp}_raxml-aic' >> ${curr_path}part_N${raw_samp}_${species}_r-a/partitioning_treebuilding_output.txt
echo \$(grep search ${curr_path}part_N${raw_samp}_${species}_r-a/analysis/best_scheme.txt) >> ${curr_path}part_N${raw_samp}_${species}_r-a/partitioning_treebuilding_output.txt
echo \$(grep Scheme ${curr_path}part_N${raw_samp}_${species}_r-a/analysis/best_scheme.txt) >> ${curr_path}part_N${raw_samp}_${species}_r-a/partitioning_treebuilding_output.txt

#________run iqtree with partition file - raxml aic________#
${run_path}iqtree-1.6.12-Linux/bin/iqtree -s ${curr_path}part_N${raw_samp}_${species}_r-a/temp_align.phy -seed 123456789 -nt 8 -spp ${curr_path}part_N${raw_samp}_${species}_r-a/best_scheme_nexus.txt -bb 1000
echo \$(grep -A 12 MAXIMUM ${curr_path}part_N${raw_samp}_${species}_r-a/best_scheme_nexus.txt.iqtree) >> ${curr_path}part_N${raw_samp}_${species}_r-a/partitioning_treebuilding_output.txt
echo \$(grep -A 2 newick ${curr_path}part_N${raw_samp}_${species}_r-a/best_scheme_nexus.txt.iqtree | head -n 3 | tail -n 1) >> ${curr_path}part_N${raw_samp}_${species}_r-a/partitioning_treebuilding_output.txt
echo \$(grep -A 5 CONSENSUS ${curr_path}part_N${raw_samp}_${species}_r-a/best_scheme_nexus.txt.iqtree) >> ${curr_path}part_N${raw_samp}_${species}_r-a/partitioning_treebuilding_output.txt
echo \$(grep -A 2 newick ${curr_path}part_N${raw_samp}_${species}_r-a/best_scheme_nexus.txt.iqtree | tail -n 1) >> ${curr_path}part_N${raw_samp}_${species}_r-a/partitioning_treebuilding_output.txt" >> ${curr_path}'part_N'$raw_samp'_'$species'_r-a/Run.slurm'


	# raxml BIC
	echo "python ${run_path}partitionfinder-2.1.1/PartitionFinderProtein.py ${curr_path}part_N${raw_samp}_${species}_r-b/ --raxml -p -1
conda deactivate
# copy nexus best scheme to a partition file
awk '/nexus/,/end/' ${curr_path}part_N${raw_samp}_${species}_r-b/analysis/best_scheme.txt > ${curr_path}part_N${raw_samp}_${species}_r-b/best_scheme_nexus.txt
# save information criteria output from partition finder
echo '${species}_${raw_samp}_raxml-bic' >> ${curr_path}part_N${raw_samp}_${species}_r-b/partitioning_treebuilding_output.txt
echo \$(grep search ${curr_path}part_N${raw_samp}_${species}_r-b/analysis/best_scheme.txt) >> ${curr_path}part_N${raw_samp}_${species}_r-b/partitioning_treebuilding_output.txt
echo \$(grep Scheme ${curr_path}part_N${raw_samp}_${species}_r-b/analysis/best_scheme.txt) >> ${curr_path}part_N${raw_samp}_${species}_r-b/partitioning_treebuilding_output.txt

#________run iqtree with partition file - raxml aic________#
${run_path}iqtree-1.6.12-Linux/bin/iqtree -s ${curr_path}part_N${raw_samp}_${species}_r-b/temp_align.phy -seed 123456789 -nt 8 -spp ${curr_path}part_N${raw_samp}_${species}_r-b/best_scheme_nexus.txt -bb 1000
echo \$(grep -A 12 MAXIMUM ${curr_path}part_N${raw_samp}_${species}_r-b/best_scheme_nexus.txt.iqtree) >> ${curr_path}part_N${raw_samp}_${species}_r-b/partitioning_treebuilding_output.txt
echo \$(grep -A 2 newick ${curr_path}part_N${raw_samp}_${species}_r-b/best_scheme_nexus.txt.iqtree | head -n 3 | tail -n 1) >> ${curr_path}part_N${raw_samp}_${species}_r-b/partitioning_treebuilding_output.txt
echo \$(grep -A 5 CONSENSUS ${curr_path}part_N${raw_samp}_${species}_r-b/best_scheme_nexus.txt.iqtree) >> ${curr_path}part_N${raw_samp}_${species}_r-b/partitioning_treebuilding_output.txt
echo \$(grep -A 2 newick ${curr_path}part_N${raw_samp}_${species}_r-b/best_scheme_nexus.txt.iqtree | tail -n 1) >> ${curr_path}part_N${raw_samp}_${species}_r-b/partitioning_treebuilding_output.txt" >> ${curr_path}'part_N'$raw_samp'_'$species'_r-b/Run.slurm'

	# raxml AICc
	echo "python ${run_path}partitionfinder-2.1.1/PartitionFinderProtein.py ${curr_path}part_N${raw_samp}_${species}_r-c/ --raxml -p -1
conda deactivate
# copy nexus best scheme to a partition file
awk '/nexus/,/end/' ${curr_path}part_N${raw_samp}_${species}_r-c/analysis/best_scheme.txt > ${curr_path}part_N${raw_samp}_${species}_r-c/best_scheme_nexus.txt
# save information criteria output from partition finder
echo '${species}_${raw_samp}_raxml-aicc' >> ${curr_path}part_N${raw_samp}_${species}_r-c/partitioning_treebuilding_output.txt
echo \$(grep search ${curr_path}part_N${raw_samp}_${species}_r-c/analysis/best_scheme.txt) >> ${curr_path}part_N${raw_samp}_${species}_r-c/partitioning_treebuilding_output.txt
echo \$(grep Scheme ${curr_path}part_N${raw_samp}_${species}_r-c/analysis/best_scheme.txt) >> ${curr_path}part_N${raw_samp}_${species}_r-c/partitioning_treebuilding_output.txt

#________run iqtree with partition file - raxml aic________#
${run_path}iqtree-1.6.12-Linux/bin/iqtree -s ${curr_path}part_N${raw_samp}_${species}_r-c/temp_align.phy -seed 123456789 -nt 8 -spp ${curr_path}part_N${raw_samp}_${species}_r-c/best_scheme_nexus.txt -bb 1000
echo \$(grep -A 12 MAXIMUM ${curr_path}part_N${raw_samp}_${species}_r-c/best_scheme_nexus.txt.iqtree) >> ${curr_path}part_N${raw_samp}_${species}_r-c/partitioning_treebuilding_output.txt
echo \$(grep -A 2 newick ${curr_path}part_N${raw_samp}_${species}_r-c/best_scheme_nexus.txt.iqtree | head -n 3 | tail -n 1) >> ${curr_path}part_N${raw_samp}_${species}_r-c/partitioning_treebuilding_output.txt
echo \$(grep -A 5 CONSENSUS ${curr_path}part_N${raw_samp}_${species}_r-c/best_scheme_nexus.txt.iqtree) >> ${curr_path}part_N${raw_samp}_${species}_r-c/partitioning_treebuilding_output.txt
echo \$(grep -A 2 newick ${curr_path}part_N${raw_samp}_${species}_r-c/best_scheme_nexus.txt.iqtree | tail -n 1) >> ${curr_path}part_N${raw_samp}_${species}_r-c/partitioning_treebuilding_output.txt" >> ${curr_path}'part_N'$raw_samp'_'$species'_r-c/Run.slurm'

	# greedy AIC
	echo "python ${run_path}partitionfinder-2.1.1/PartitionFinderProtein.py ${curr_path}part_N${raw_samp}_${species}_g-a/ -p -1
conda deactivate
# copy nexus best scheme to a partition file
awk '/nexus/,/end/' ${curr_path}part_N${raw_samp}_${species}_g-a/analysis/best_scheme.txt > ${curr_path}part_N${raw_samp}_${species}_g-a/best_scheme_nexus.txt
# save information criteria output from partition finder
echo '${species}_${raw_samp}_greedy-aic' >> ${curr_path}part_N${raw_samp}_${species}_g-a/partitioning_treebuilding_output.txt
echo \$(grep search ${curr_path}part_N${raw_samp}_${species}_g-a/analysis/best_scheme.txt) >> ${curr_path}part_N${raw_samp}_${species}_g-a/partitioning_treebuilding_output.txt
echo \$(grep Scheme ${curr_path}part_N${raw_samp}_${species}_g-a/analysis/best_scheme.txt) >> ${curr_path}part_N${raw_samp}_${species}_g-a/partitioning_treebuilding_output.txt

#________run iqtree with partition file - raxml aic________#
${run_path}iqtree-1.6.12-Linux/bin/iqtree -s ${curr_path}part_N${raw_samp}_${species}_g-a/temp_align.phy -seed 123456789 -nt 8 -spp ${curr_path}part_N${raw_samp}_${species}_g-a/best_scheme_nexus.txt -bb 1000
echo \$(grep -A 12 MAXIMUM ${curr_path}part_N${raw_samp}_${species}_g-a/best_scheme_nexus.txt.iqtree) >> ${curr_path}part_N${raw_samp}_${species}_g-a/partitioning_treebuilding_output.txt
echo \$(grep -A 2 newick ${curr_path}part_N${raw_samp}_${species}_g-a/best_scheme_nexus.txt.iqtree | head -n 3 | tail -n 1) >> ${curr_path}part_N${raw_samp}_${species}_g-a/partitioning_treebuilding_output.txt
echo \$(grep -A 5 CONSENSUS ${curr_path}part_N${raw_samp}_${species}_g-a/best_scheme_nexus.txt.iqtree) >> ${curr_path}part_N${raw_samp}_${species}_g-a/partitioning_treebuilding_output.txt
echo \$(grep -A 2 newick ${curr_path}part_N${raw_samp}_${species}_g-a/best_scheme_nexus.txt.iqtree | tail -n 1) >> ${curr_path}part_N${raw_samp}_${species}_g-a/partitioning_treebuilding_output.txt" >> ${curr_path}'part_N'$raw_samp'_'$species'_g-a/Run.slurm'

	# greedy BIC
	echo "python ${run_path}partitionfinder-2.1.1/PartitionFinderProtein.py ${curr_path}part_N${raw_samp}_${species}_g-b/ -p -1
conda deactivate
# copy nexus best scheme to a partition file
awk '/nexus/,/end/' ${curr_path}part_N${raw_samp}_${species}_g-b/analysis/best_scheme.txt > ${curr_path}part_N${raw_samp}_${species}_g-b/best_scheme_nexus.txt
# save information criteria output from partition finder
echo '${species}_${raw_samp}_greedy-bic' >> ${curr_path}part_N${raw_samp}_${species}_g-b/partitioning_treebuilding_output.txt
echo \$(grep search ${curr_path}part_N${raw_samp}_${species}_g-b/analysis/best_scheme.txt) >> ${curr_path}part_N${raw_samp}_${species}_g-b/partitioning_treebuilding_output.txt
echo \$(grep Scheme ${curr_path}part_N${raw_samp}_${species}_g-b/analysis/best_scheme.txt) >> ${curr_path}part_N${raw_samp}_${species}_g-b/partitioning_treebuilding_output.txt

#________run iqtree with partition file - raxml aic________#
${run_path}iqtree-1.6.12-Linux/bin/iqtree -s ${curr_path}part_N${raw_samp}_${species}_g-b/temp_align.phy -seed 123456789 -nt 8 -spp ${curr_path}part_N${raw_samp}_${species}_g-b/best_scheme_nexus.txt -bb 1000
echo \$(grep -A 12 MAXIMUM ${curr_path}part_N${raw_samp}_${species}_g-b/best_scheme_nexus.txt.iqtree) >> ${curr_path}part_N${raw_samp}_${species}_g-b/partitioning_treebuilding_output.txt
echo \$(grep -A 2 newick ${curr_path}part_N${raw_samp}_${species}_g-b/best_scheme_nexus.txt.iqtree | head -n 3 | tail -n 1) >> ${curr_path}part_N${raw_samp}_${species}_g-b/partitioning_treebuilding_output.txt
echo \$(grep -A 5 CONSENSUS ${curr_path}part_N${raw_samp}_${species}_g-b/best_scheme_nexus.txt.iqtree) >> ${curr_path}part_N${raw_samp}_${species}_g-b/partitioning_treebuilding_output.txt
echo \$(grep -A 2 newick ${curr_path}part_N${raw_samp}_${species}_g-b/best_scheme_nexus.txt.iqtree | tail -n 1) >> ${curr_path}part_N${raw_samp}_${species}_g-b/partitioning_treebuilding_output.txt" >> ${curr_path}'part_N'$raw_samp'_'$species'_g-b/Run.slurm'

	# greedy AICc
	echo "python ${run_path}partitionfinder-2.1.1/PartitionFinderProtein.py ${curr_path}part_N${raw_samp}_${species}_g-c/ -p -1
conda deactivate
# copy nexus best scheme to a partition file
awk '/nexus/,/end/' ${curr_path}part_N${raw_samp}_${species}_g-c/analysis/best_scheme.txt > ${curr_path}part_N${raw_samp}_${species}_g-c/best_scheme_nexus.txt
# save information criteria output from partition finder
echo '${species}_${raw_samp}_greedy-aicc' >> ${curr_path}part_N${raw_samp}_${species}_g-c/partitioning_treebuilding_output.txt
echo \$(grep search ${curr_path}part_N${raw_samp}_${species}_g-c/analysis/best_scheme.txt) >> ${curr_path}part_N${raw_samp}_${species}_g-c/partitioning_treebuilding_output.txt
echo \$(grep Scheme ${curr_path}part_N${raw_samp}_${species}_g-c/analysis/best_scheme.txt) >> ${curr_path}part_N${raw_samp}_${species}_g-c/partitioning_treebuilding_output.txt

#________run iqtree with partition file - raxml aic________#
${run_path}iqtree-1.6.12-Linux/bin/iqtree -s ${curr_path}part_N${raw_samp}_${species}_g-c/temp_align.phy -seed 123456789 -nt 8 -spp ${curr_path}part_N${raw_samp}_${species}_g-c/best_scheme_nexus.txt -bb 1000
echo \$(grep -A 12 MAXIMUM ${curr_path}part_N${raw_samp}_${species}_g-c/best_scheme_nexus.txt.iqtree) >> ${curr_path}part_N${raw_samp}_${species}_g-c/partitioning_treebuilding_output.txt
echo \$(grep -A 2 newick ${curr_path}part_N${raw_samp}_${species}_g-c/best_scheme_nexus.txt.iqtree | head -n 3 | tail -n 1) >> ${curr_path}part_N${raw_samp}_${species}_g-c/partitioning_treebuilding_output.txt
echo \$(grep -A 5 CONSENSUS ${curr_path}part_N${raw_samp}_${species}_g-c/best_scheme_nexus.txt.iqtree) >> ${curr_path}part_N${raw_samp}_${species}_g-c/partitioning_treebuilding_output.txt
echo \$(grep -A 2 newick ${curr_path}part_N${raw_samp}_${species}_g-c/best_scheme_nexus.txt.iqtree | tail -n 1) >> ${curr_path}part_N${raw_samp}_${species}_g-c/partitioning_treebuilding_output.txt" >> ${curr_path}'part_N'$raw_samp'_'$species'_g-c/Run.slurm'

	# run IQtree withOUT partitioning	
	echo "
${run_path}iqtree-1.6.12-Linux/bin/iqtree -s ${curr_path}part_N${raw_samp}_${species}_noPart/temp_align.phy -seed 123456789 -nt 8 -mrate G -bb 1000
echo \$(grep -A 12 MAXIMUM ${curr_path}part_N${raw_samp}_${species}_noPart/temp_align.phy.iqtree) >> ${curr_path}part_N${raw_samp}_${species}_noPart/partitioning_treebuilding_output.txt
echo \$(grep -A 2 newick ${curr_path}part_N${raw_samp}_${species}_noPart/temp_align.phy.iqtree | head -n 3 | tail -n 1) >> ${curr_path}part_N${raw_samp}_${species}_noPart/partitioning_treebuilding_output.txt
echo \$(grep -A 5 CONSENSUS ${curr_path}part_N${raw_samp}_${species}_noPart/temp_align.phy.iqtree) >> ${curr_path}part_N${raw_samp}_${species}_noPart/partitioning_treebuilding_output.txt
echo \$(grep -A 2 newick ${curr_path}part_N${raw_samp}_${species}_noPart/temp_align.phy.iqtree | tail -n 1) >> ${curr_path}part_N${raw_samp}_${species}_noPart/partitioning_treebuilding_output.txt" >> ${curr_path}'part_N'$raw_samp'_'$species'_noPart/Run.slurm'
	
done
