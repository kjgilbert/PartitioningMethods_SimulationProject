#!/bin/bash

# raxmlHPC-PTHREADS -T 4 -f a -m PROTGAMMAWAG -s ${CONCAT_ALIGNMENT} -p 15826 -x 15826 -q ${PARTITION_FILE} -w ${DIR} -n ${TREENAME} -# 100
# raxmlHPC-PTHREADS -T 4 -f a -m PROTGAMMAWAG -s part_N10_lopho_r-a/temp_align.phy -p 15826 -x 15826 -q part_N10_lopho_r-a/best_scheme_nexus.txt -w part_N10_lopho_r-a/ -n raxml_tree.out -# 100

#current_replicate='New_Rep1_100-1000/'

current_replicate=$1
full_path='/scratch/axiom/FAC/FBM/DBC/cdessim2/default/kgilbert/'


#num_threads=12
num_threads=$2

part_folders=part_N*

for part_folder in ${part_folders};
	do
	if echo $part_folder | grep -q noPart ; then
    	../standard-RAxML/raxmlHPC-PTHREADS -T 4 -f a -k -m PROTGAMMAAUTO ­­auto­prot=bic -s ${full_path}${current_replicate}${part_folder}/temp_align.phy -p 15826 -x 15826 -w ${full_path}${current_replicate}${part_folder}/ -n raxml_tree.out -# 100
	else
    	../standard-RAxML/raxmlHPC-PTHREADS -T ${num_threads} -f a -m PROTGAMMAWAG -s ${full_path}${current_replicate}${part_folder}/temp_align.phy -p 15826 -x 15826 -q ${full_path}${current_replicate}${part_folder}/best_scheme_nexus.txt -w ${full_path}${current_replicate}${part_folder}/ -n raxml_tree.out -# 100
	fi
done

