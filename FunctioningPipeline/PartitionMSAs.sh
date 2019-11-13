#!/bin/bash

input_folder='sampledMSAs'

#mkdir part_N20
#mkdir part_N40
#mkdir part_N100
#mkdir part_Nall

# put the cfg file in each corresponding directory

# do all N20

# aic rcluster

for f in glob.glob(os.path.join(input_folder, '*.phy')):
    do
	cp f part_N20/temp_align.phy
	cp partition_finder_N20_raxml_aic.cfg part_N20/partition_finder_cfg
	
	python2.7 PartitionFinderProtein.py ${OUTDIR} --raxml -p 4

# bic rcluster
# aicc rcluster
# aic greedy
# bic greedy
# aicc greedy


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
