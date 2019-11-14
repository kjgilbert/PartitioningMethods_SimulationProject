#!/bin/bash

input_folder='sampledMSAs/'
files=${input_folder}*.phy
###files=`ls ${input_folder}*.phy`

#mkdir part_N20
#mkdir part_N40
#mkdir part_N100
#mkdir part_Nall

# put the cfg file in each corresponding directory

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
