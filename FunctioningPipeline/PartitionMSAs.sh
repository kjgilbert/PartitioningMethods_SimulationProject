#!/bin/bash

input_folder='sampledMSAs/'
files=${input_folder}*.phy
###files=`ls ${input_folder}*.phy`

#mkdir part_N20
#mkdir part_N40
#mkdir part_N100
#mkdir part_Nall

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



# do all N10

for f in files;
	do
	cp $f part_N10/temp_align.phy
	# make the partition file based on the gene sites stored in config_genes__.txt
		# format the gene sites
		# species:
		species=$(echo $f | tr "_" "\n" | tail -n 2 | head -n 1)
		# sample size of genes
		sample_size=$(echo $f | tr "_" "\n" | tail -n 1 | tr "." "\n" | head -n 1)
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
		echo "$aic_string" "$gene_string" "$greedy_string" > partition_finder_N10_raxml_aic.cfg	
		echo "$aic_string" "$gene_string" "$rcluster_string" > partition_finder_N10_greedy_aic.cfg	
		echo "$bic_string" "$gene_string" "$greedy_string" > partition_finder_N10_raxml_bic.cfg	
		echo "$bic_string" "$gene_string" "$rcluster_string" > partition_finder_N10_greedy_bic.cfg	
		echo "$aicc_string" "$gene_string" "$greedy_string" > partition_finder_N10_raxml_aicc.cfg	
		echo "$aicc_string" "$gene_string" "$rcluster_string" > partition_finder_N10_greedy_aicc.cfg	

	# run all combinations of analysis on that file
	cp partition_finder_N10_raxml_aic.cfg part_N10/partition_finder.cfg	
	python2.7 ../../partition_finder_2/partitionfinder-2.1.1/PartitionFinderProtein.py part_N10/ --raxml -p 4

	cp partition_finder_N10_greedy_aic.cfg part_N10/partition_finder.cfg	
	python2.7 ../../partition_finder_2/partitionfinder-2.1.1/PartitionFinderProtein.py part_N10/ -p 4

	cp partition_finder_N10_raxml_bic.cfg part_N10/partition_finder.cfg	
	python2.7 ../../partition_finder_2/partitionfinder-2.1.1/PartitionFinderProtein.py part_N10/ --raxml -p 4

	cp partition_finder_N10_greedy_bic.cfg part_N10/partition_finder.cfg	
	python2.7 ../../partition_finder_2/partitionfinder-2.1.1/PartitionFinderProtein.py part_N10/ -p 4

	cp partition_finder_N10_raxml_aicc.cfg part_N10/partition_finder.cfg	
	python2.7 ../../partition_finder_2/partitionfinder-2.1.1/PartitionFinderProtein.py part_N10/ --raxml -p 4

	cp partition_finder_N10_greedy_aicc.cfg part_N10/partition_finder.cfg	
	python2.7 ../../partition_finder_2/partitionfinder-2.1.1/PartitionFinderProtein.py part_N10/ -p 4

done


# do all N20


for f in files;
	do
	cp $f part_N20/temp_align.phy
	# make the partition file based on the gene sites stored in config_genes__.txt
		# format the gene sites
		# species:
		species=$(echo $f | tr "_" "\n" | tail -n 2 | head -n 1)
		# sample size of genes
		sample_size=$(echo $f | tr "_" "\n" | tail -n 1 | tr "." "\n" | head -n 1)
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
		echo "$aic_string" "$gene_string" "$greedy_string" > partition_finder_N20_raxml_aic.cfg	
		echo "$aic_string" "$gene_string" "$rcluster_string" > partition_finder_N20_greedy_aic.cfg	
		echo "$bic_string" "$gene_string" "$greedy_string" > partition_finder_N20_raxml_bic.cfg	
		echo "$bic_string" "$gene_string" "$rcluster_string" > partition_finder_N20_greedy_bic.cfg	
		echo "$aicc_string" "$gene_string" "$greedy_string" > partition_finder_N20_raxml_aicc.cfg	
		echo "$aicc_string" "$gene_string" "$rcluster_string" > partition_finder_N20_greedy_aicc.cfg	

	# run all combinations of analysis on that file
	cp partition_finder_N20_raxml_aic.cfg part_N20/partition_finder.cfg	
	python2.7 ../../partition_finder_2/partitionfinder-2.1.1/PartitionFinderProtein.py part_N20/ --raxml -p 4

	cp partition_finder_N20_greedy_aic.cfg part_N20/partition_finder.cfg	
	python2.7 ../../partition_finder_2/partitionfinder-2.1.1/PartitionFinderProtein.py part_N20/ -p 4

	cp partition_finder_N20_raxml_bic.cfg part_N20/partition_finder.cfg	
	python2.7 ../../partition_finder_2/partitionfinder-2.1.1/PartitionFinderProtein.py part_N20/ --raxml -p 4

	cp partition_finder_N20_greedy_bic.cfg part_N20/partition_finder.cfg	
	python2.7 ../../partition_finder_2/partitionfinder-2.1.1/PartitionFinderProtein.py part_N20/ -p 4

	cp partition_finder_N20_raxml_aicc.cfg part_N20/partition_finder.cfg	
	python2.7 ../../partition_finder_2/partitionfinder-2.1.1/PartitionFinderProtein.py part_N20/ --raxml -p 4

	cp partition_finder_N20_greedy_aicc.cfg part_N20/partition_finder.cfg	
	python2.7 ../../partition_finder_2/partitionfinder-2.1.1/PartitionFinderProtein.py part_N20/ -p 4

done


# do all N40


for f in files;
	do
	cp $f part_N40/temp_align.phy
	# make the partition file based on the gene sites stored in config_genes__.txt
		# format the gene sites
		# species:
		species=$(echo $f | tr "_" "\n" | tail -n 2 | head -n 1)
		# sample size of genes
		sample_size=$(echo $f | tr "_" "\n" | tail -n 1 | tr "." "\n" | head -n 1)
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
		echo "$aic_string" "$gene_string" "$greedy_string" > partition_finder_N40_raxml_aic.cfg	
		echo "$aic_string" "$gene_string" "$rcluster_string" > partition_finder_N40_greedy_aic.cfg	
		echo "$bic_string" "$gene_string" "$greedy_string" > partition_finder_N40_raxml_bic.cfg	
		echo "$bic_string" "$gene_string" "$rcluster_string" > partition_finder_N40_greedy_bic.cfg	
		echo "$aicc_string" "$gene_string" "$greedy_string" > partition_finder_N40_raxml_aicc.cfg	
		echo "$aicc_string" "$gene_string" "$rcluster_string" > partition_finder_N40_greedy_aicc.cfg	

	# run all combinations of analysis on that file
	cp partition_finder_N40_raxml_aic.cfg part_N40/partition_finder.cfg	
	python2.7 ../../partition_finder_2/partitionfinder-2.1.1/PartitionFinderProtein.py part_N40/ --raxml -p 4

	cp partition_finder_N40_greedy_aic.cfg part_N40/partition_finder.cfg	
	python2.7 ../../partition_finder_2/partitionfinder-2.1.1/PartitionFinderProtein.py part_N40/ -p 4

	cp partition_finder_N40_raxml_bic.cfg part_N40/partition_finder.cfg	
	python2.7 ../../partition_finder_2/partitionfinder-2.1.1/PartitionFinderProtein.py part_N40/ --raxml -p 4

	cp partition_finder_N40_greedy_bic.cfg part_N40/partition_finder.cfg	
	python2.7 ../../partition_finder_2/partitionfinder-2.1.1/PartitionFinderProtein.py part_N40/ -p 4

	cp partition_finder_N40_raxml_aicc.cfg part_N40/partition_finder.cfg	
	python2.7 ../../partition_finder_2/partitionfinder-2.1.1/PartitionFinderProtein.py part_N40/ --raxml -p 4

	cp partition_finder_N40_greedy_aicc.cfg part_N40/partition_finder.cfg	
	python2.7 ../../partition_finder_2/partitionfinder-2.1.1/PartitionFinderProtein.py part_N40/ -p 4

done


# do all N 80

rm partitioning_treebuilding_output.txt
touch partitioning_treebuilding_output.txt

#________go through each gene sample set________#


dir='N80'

#________go through each species file for a given gene sample________#

for f in files;
	do
	cp $f 'part_'$dir'/temp_align.phy'
	# make the partition file based on the gene sites stored in config_genes__.txt
		# format the gene sites
		# species:
		species=$(echo $f | tr "_" "\n" | tail -n 2 | head -n 1)
		# sample size of genes
		sample_size=$(echo $f | tr "_" "\n" | tail -n 1 | tr "." "\n" | head -n 1)
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
		echo "$aic_string" "$gene_string" "$greedy_string" > 'partition_finder_'$dir'_raxml_aic.cfg'	
		echo "$aic_string" "$gene_string" "$rcluster_string" > 'partition_finder_'$dir'_greedy_aic.cfg'	
		echo "$bic_string" "$gene_string" "$greedy_string" > 'partition_finder_'$dir'_raxml_bic.cfg'
		echo "$bic_string" "$gene_string" "$rcluster_string" > 'partition_finder_'$dir'_greedy_bic.cfg'	
		echo "$aicc_string" "$gene_string" "$greedy_string" > 'partition_finder_'$dir'_raxml_aicc.cfg'
		echo "$aicc_string" "$gene_string" "$rcluster_string" > 'partition_finder_'$dir'_greedy_aicc.cfg'	

	# run all combinations of partitioning analysis on that file and the associated iqtree analysis
	# label in the output I'm saving which file this is
	echo $f >> partitioning_treebuilding_output.txt
	
#________run partition finder - raxml aic________#
	# clean out partitioning results before doing new partitioning
	rm -r 'part_'$dir'/analysis'
	rm  'part_'$dir'/log.txt'
	cp 'partition_finder_'$dir'_raxml_aic.cfg' 'part_'$dir'/partition_finder.cfg'
	python2.7 ../../partition_finder_2/partitionfinder-2.1.1/PartitionFinderProtein.py 'part_'$dir'/' --raxml -p 4
	# copy nexus best scheme to a partition file
	awk '/nexus/,/end/' 'part_'$dir'/analysis/best_scheme.txt' > 'part_'$dir'/best_scheme_nexus.txt'
	# save information criteria output from partition finder
	echo $(grep search 'part_'$dir'/analysis/best_scheme.txt') >> partitioning_treebuilding_output.txt
	echo $(grep Scheme 'part_'$dir'/analysis/best_scheme.txt') >> partitioning_treebuilding_output.txt
	
#________run iqtree with partition file - raxml aic________#
	#IQTREE_PartitionFinder
	iqtree -s 'part_'$dir'/temp_align.phy' -seed 123456789 -nt 4 -spp 'part_'$dir'/best_scheme_nexus.txt' -bb 1000
	echo $(grep -A 12 MAXIMUM 'part_'$dir'/best_scheme_nexus.txt.iqtree') >> partitioning_treebuilding_output.txt
	echo $(grep -A 2 newick 'part_'$dir'/best_scheme_nexus.txt.iqtree' | head -n 3 | tail -n 1) >> partitioning_treebuilding_output.txt
	echo $(grep -A 5 CONSENSUS 'part_'$dir'/best_scheme_nexus.txt.iqtree') >> partitioning_treebuilding_output.txt
	echo $(grep -A 2 newick 'part_'$dir'/best_scheme_nexus.txt.iqtree' | tail -n 1) >> partitioning_treebuilding_output.txt
	
#________run partition finder - raxml bic________#
	# clean out partitioning results before doing new partitioning
	rm -r 'part_'$dir'/analysis'
	rm  'part_'$dir'/log.txt'
	cp 'partition_finder_'$dir'_raxml_bic.cfg' 'part_'$dir'/partition_finder.cfg'
	python2.7 ../../partition_finder_2/partitionfinder-2.1.1/PartitionFinderProtein.py 'part_'$dir'/' --raxml -p 4
	# copy nexus best scheme to a partition file
	awk '/nexus/,/end/' 'part_'$dir'/analysis/best_scheme.txt' > 'part_'$dir'/best_scheme_nexus.txt'
	# save information criteria output from partition finder
	echo $(grep search 'part_'$dir'/analysis/best_scheme.txt') >> partitioning_treebuilding_output.txt
	echo $(grep Scheme 'part_'$dir'/analysis/best_scheme.txt') >> partitioning_treebuilding_output.txt

#________run iqtree with partition file - raxml bic________#
	#IQTREE_PartitionFinder
	iqtree -s 'part_'$dir'/temp_align.phy' -seed 123456789 -nt 4 -spp 'part_'$dir'/best_scheme_nexus.txt' -bb 1000
	echo $(grep -A 12 MAXIMUM 'part_'$dir'/best_scheme_nexus.txt.iqtree') >> partitioning_treebuilding_output.txt
	echo $(grep -A 2 newick 'part_'$dir'/best_scheme_nexus.txt.iqtree' | head -n 3 | tail -n 1) >> partitioning_treebuilding_output.txt
	echo $(grep -A 5 CONSENSUS 'part_'$dir'/best_scheme_nexus.txt.iqtree') >> partitioning_treebuilding_output.txt
	echo $(grep -A 2 newick 'part_'$dir'/best_scheme_nexus.txt.iqtree' | tail -n 1) >> partitioning_treebuilding_output.txt

#________run partition finder - raxml aicc________#
	# clean out partitioning results before doing new partitioning
	rm -r 'part_'$dir'/analysis'
	rm  'part_'$dir'/log.txt'
	cp 'partition_finder_'$dir'_raxml_aicc.cfg' 'part_'$dir'/partition_finder.cfg'
	python2.7 ../../partition_finder_2/partitionfinder-2.1.1/PartitionFinderProtein.py 'part_'$dir'/' --raxml -p 4
	# copy nexus best scheme to a partition file
	awk '/nexus/,/end/' 'part_'$dir'/analysis/best_scheme.txt' > 'part_'$dir'/best_scheme_nexus.txt'
	# save information criteria output from partition finder
	echo $(grep search 'part_'$dir'/analysis/best_scheme.txt') >> partitioning_treebuilding_output.txt
	echo $(grep Scheme 'part_'$dir'/analysis/best_scheme.txt') >> partitioning_treebuilding_output.txt

#________run iqtree with partition file - raxml aicc________#
	#IQTREE_PartitionFinder
	iqtree -s 'part_'$dir'/temp_align.phy' -seed 123456789 -nt 4 -spp 'part_'$dir'/best_scheme_nexus.txt' -bb 1000
	echo $(grep -A 12 MAXIMUM 'part_'$dir'/best_scheme_nexus.txt.iqtree') >> partitioning_treebuilding_output.txt
	echo $(grep -A 2 newick 'part_'$dir'/best_scheme_nexus.txt.iqtree' | head -n 3 | tail -n 1) >> partitioning_treebuilding_output.txt
	echo $(grep -A 5 CONSENSUS 'part_'$dir'/best_scheme_nexus.txt.iqtree') >> partitioning_treebuilding_output.txt
	echo $(grep -A 2 newick 'part_'$dir'/best_scheme_nexus.txt.iqtree' | tail -n 1) >> partitioning_treebuilding_output.txt
	

	
#________run partition finder - greedy aic________#
	# clean out partitioning results before doing new partitioning
	rm -r 'part_'$dir'/analysis'
	rm  'part_'$dir'/log.txt'
	cp 'partition_finder_'$dir'_greedy_aic.cfg' 'part_'$dir'/partition_finder.cfg'
	python2.7 ../../partition_finder_2/partitionfinder-2.1.1/PartitionFinderProtein.py 'part_'$dir'/' -p 4
	# copy nexus best scheme to a partition file
	awk '/nexus/,/end/' 'part_'$dir'/analysis/best_scheme.txt' > 'part_'$dir'/best_scheme_nexus.txt'
	# save information criteria output from partition finder
	echo $(grep search 'part_'$dir'/analysis/best_scheme.txt') >> partitioning_treebuilding_output.txt
	echo $(grep Scheme 'part_'$dir'/analysis/best_scheme.txt') >> partitioning_treebuilding_output.txt
	
#________run iqtree with partition file - greedy aic________#
	#IQTREE_PartitionFinder
	iqtree -s 'part_'$dir'/temp_align.phy' -seed 123456789 -nt 4 -spp 'part_'$dir'/best_scheme_nexus.txt' -bb 1000
	echo $(grep -A 12 MAXIMUM 'part_'$dir'/best_scheme_nexus.txt.iqtree') >> partitioning_treebuilding_output.txt
	echo $(grep -A 2 newick 'part_'$dir'/best_scheme_nexus.txt.iqtree' | head -n 3 | tail -n 1) >> partitioning_treebuilding_output.txt
	echo $(grep -A 5 CONSENSUS 'part_'$dir'/best_scheme_nexus.txt.iqtree') >> partitioning_treebuilding_output.txt
	echo $(grep -A 2 newick 'part_'$dir'/best_scheme_nexus.txt.iqtree' | tail -n 1) >> partitioning_treebuilding_output.txt
	
#________run partition finder - greedy bic________#
	# clean out partitioning results before doing new partitioning
	rm -r 'part_'$dir'/analysis'
	rm  'part_'$dir'/log.txt'
	cp 'partition_finder_'$dir'_greedy_bic.cfg' 'part_'$dir'/partition_finder.cfg'
	python2.7 ../../partition_finder_2/partitionfinder-2.1.1/PartitionFinderProtein.py 'part_'$dir'/' -p 4
	# copy nexus best scheme to a partition file
	awk '/nexus/,/end/' 'part_'$dir'/analysis/best_scheme.txt' > 'part_'$dir'/best_scheme_nexus.txt'
	# save information criteria output from partition finder
	echo $(grep search 'part_'$dir'/analysis/best_scheme.txt') >> partitioning_treebuilding_output.txt
	echo $(grep Scheme 'part_'$dir'/analysis/best_scheme.txt') >> partitioning_treebuilding_output.txt

#________run iqtree with partition file - greedy bic________#
	#IQTREE_PartitionFinder
	iqtree -s 'part_'$dir'/temp_align.phy' -seed 123456789 -nt 4 -spp 'part_'$dir'/best_scheme_nexus.txt' -bb 1000
	echo $(grep -A 12 MAXIMUM 'part_'$dir'/best_scheme_nexus.txt.iqtree') >> partitioning_treebuilding_output.txt
	echo $(grep -A 2 newick 'part_'$dir'/best_scheme_nexus.txt.iqtree' | head -n 3 | tail -n 1) >> partitioning_treebuilding_output.txt
	echo $(grep -A 5 CONSENSUS 'part_'$dir'/best_scheme_nexus.txt.iqtree') >> partitioning_treebuilding_output.txt
	echo $(grep -A 2 newick 'part_'$dir'/best_scheme_nexus.txt.iqtree' | tail -n 1) >> partitioning_treebuilding_output.txt

#________run partition finder - greedy aicc________#
	# clean out partitioning results before doing new partitioning
	rm -r 'part_'$dir'/analysis'
	rm  'part_'$dir'/log.txt'
	cp 'partition_finder_'$dir'_greedy_aicc.cfg' 'part_'$dir'/partition_finder.cfg'
	python2.7 ../../partition_finder_2/partitionfinder-2.1.1/PartitionFinderProtein.py 'part_'$dir'/' -p 4
	# copy nexus best scheme to a partition file
	awk '/nexus/,/end/' 'part_'$dir'/analysis/best_scheme.txt' > 'part_'$dir'/best_scheme_nexus.txt'
	# save information criteria output from partition finder
	echo $(grep search 'part_'$dir'/analysis/best_scheme.txt') >> partitioning_treebuilding_output.txt
	echo $(grep Scheme 'part_'$dir'/analysis/best_scheme.txt') >> partitioning_treebuilding_output.txt

#________run iqtree with partition file - greedy aicc________#
	#IQTREE_PartitionFinder
	iqtree -s 'part_'$dir'/temp_align.phy' -seed 123456789 -nt 4 -spp 'part_'$dir'/best_scheme_nexus.txt' -bb 1000
	echo $(grep -A 12 MAXIMUM 'part_'$dir'/best_scheme_nexus.txt.iqtree') >> partitioning_treebuilding_output.txt
	echo $(grep -A 2 newick 'part_'$dir'/best_scheme_nexus.txt.iqtree' | head -n 3 | tail -n 1) >> partitioning_treebuilding_output.txt
	echo $(grep -A 5 CONSENSUS 'part_'$dir'/best_scheme_nexus.txt.iqtree') >> partitioning_treebuilding_output.txt
	echo $(grep -A 2 newick 'part_'$dir'/best_scheme_nexus.txt.iqtree' | tail -n 1) >> partitioning_treebuilding_output.txt
	




#________iqtree without partition file________#
	#IQTREE_Default, can be run once on the alignment without all the partition combinations:
	iqtree -s 'part_'$dir'/temp_align.phy' -seed 123456789 -nt 4 -mrate G -bb 1000
	echo $(grep -A 12 MAXIMUM 'part_'$dir'/temp_align.phy.iqtree') >> partitioning_treebuilding_output.txt
	echo $(grep -A 2 newick 'part_'$dir'/temp_align.phy.iqtree' | head -n 3 | tail -n 1) >> partitioning_treebuilding_output.txt
	echo $(grep -A 5 CONSENSUS 'part_'$dir'/temp_align.phy.iqtree') >> partitioning_treebuilding_output.txt
	echo $(grep -A 2 newick 'part_'$dir'/temp_align.phy.iqtree' | tail -n 1) >> partitioning_treebuilding_output.txt
	

	

done # done all combos on that file






## #RAxML_Default
## raxmlHPC-PTHREADS -T 4 -f a -k -m PROTGAMMAAUTO ­­auto­prot=bic -s ${CONCAT_ALIGNMENT} -p 15826 -x 15826 -w ${DIR} -n ${TREENAME} -# 100

## #RAxML_PartitionFinder
## raxmlHPC-PTHREADS -T 4 -f a -m PROTGAMMAWAG -s ${CONCAT_ALIGNMENT} -p 15826 -x 15826 -q ${PARTITION_FILE} -w ${DIR} -n ${TREENAME} -# 100
