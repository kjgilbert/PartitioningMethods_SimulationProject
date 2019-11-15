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

# aic rcluster
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
#			echo $gene_start
#			echo $gene_end
#			echo $prev_gene_end
		done <$gene_sites_file
		
		# put all the strings together into the final .cfg file
		echo "$bic_string" "$gene_string" "$greedy_string" > test.txt

	cp partition_finder_N10_raxml_aic.cfg part_N10/partition_finder.cfg	
	python2.7 ../../partition_finder_2/partitionfinder-2.1.1/PartitionFinderProtein.py part_N10/ --raxml -p 4
done


# do all N20

# aic rcluster
for f in files;
	do
	cp $f part_N20/temp_align.phy
	cp partition_finder_N20_raxml_aic.cfg part_N20/partition_finder.cfg	
	python2.7 ../../partition_finder_2/partitionfinder-2.1.1/PartitionFinderProtein.py part_N20/ --raxml -p 4
done

# bic rcluster
for f in files;
	do
	cp $f part_N20/temp_align.phy
	cp partition_finder_N20_raxml_bic.cfg part_N20/partition_finder.cfg	
	python2.7 ../../partition_finder_2/partitionfinder-2.1.1/PartitionFinderProtein.py part_N20/ --raxml -p 4
done

# aicc rcluster
for f in files;
	do
	cp $f part_N20/temp_align.phy
	cp partition_finder_N20_raxml_aicc.cfg part_N20/partition_finder.cfg	
	python2.7 ../../partition_finder_2/partitionfinder-2.1.1/PartitionFinderProtein.py part_N20/ --raxml -p 4

!iqtree -s ../partition_finder_2/PartitionAlfSims/test1_24oct/test_MSA1.phy -spp ../partition_finder_2/PartitionAlfSims/test1_partitionFile.nex

	cp part_N20/best_scheme.txt ../
	cp part_N20/
done






# aic greedy
for f in files;
	do
	cp $f part_N20/temp_align.phy
	cp partition_finder_N20_greedy_aic.cfg part_N20/partition_finder.cfg	
	python2.7 ../../partition_finder_2/partitionfinder-2.1.1/PartitionFinderProtein.py part_N20/ -p 4
done

# bic greedy
for f in files;
	do
	cp $f part_N20/temp_align.phy
	cp partition_finder_N20_greedy_bic.cfg part_N20/partition_finder.cfg	
	python2.7 ../../partition_finder_2/partitionfinder-2.1.1/PartitionFinderProtein.py part_N20/ -p 4
done

# aicc greedy
for f in files;
	do
	cp $f part_N20/temp_align.phy
	cp partition_finder_N20_greedy_aicc.cfg part_N20/partition_finder.cfg	
	python2.7 ../../partition_finder_2/partitionfinder-2.1.1/PartitionFinderProtein.py part_N20/ -p 4
done



# go through each file in sampled_MSAs/ and rename into the correpsonding partitioning
#     folder as temp_align.phy so that the cfg file runs on that input

partitionfinder-2.1.1/PartitionFinderProtein.py PartitionAlfSims/test3_13nov/ --raxml

# then save the desired output from that run into new folder with naming convention







#IQTREE_Default
iqtree -s ${CONCAT_ALIGNMENT} -seed 12345 -nt 4 -mrate G -bb 1000

#IQTREE_PartitionFinder
iqtree -s ${CONCAT_ALIGNMENT} -seed 12345 -nt 4 -spp ${PARTITION_FILE} -bb 1000



#RAxML_Default
raxmlHPC-PTHREADS -T 4 -f a -k -m PROTGAMMAAUTO ­­auto­prot=bic -s ${CONCAT_ALIGNMENT} -p 15826 -x 15826 -w ${DIR} -n ${TREENAME} -# 100

#RAxML_PartitionFinder
raxmlHPC-PTHREADS -T 4 -f a -m PROTGAMMAWAG -s ${CONCAT_ALIGNMENT} -p 15826 -x 15826 -q ${PARTITION_FILE} -w ${DIR} -n ${TREENAME} -# 100
