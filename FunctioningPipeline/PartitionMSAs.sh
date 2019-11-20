#!/bin/bash

current_path=$1
input_folder=$current_path'sampledMSAs/'
files=${input_folder}*.phy


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









#________go through each gene sample size set________#

rm -f partitioning_treebuilding_output.txt
touch partitioning_treebuilding_output.txt


for f in $files;
    do

    # sample size of genes
    sample_size=$(echo $f | tr "_" "\n" | tail -n 1 | tr "." "\n" | head -n 1)
    #raw numeric sample size for directory label
    raw_samp=$(echo $sample_size | tr "g" "\n" | head -n 1)

    cp $f 'part_N'$raw_samp'/temp_align.phy'
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
    echo "$aic_string" "$gene_string" "$rcluster_string" > 'partition_finder_N'$raw_samp'_raxml_aic.cfg'
    echo "$aic_string" "$gene_string" "$greedy_string" > 'partition_finder_N'$raw_samp'_greedy_aic.cfg'
    echo "$bic_string" "$gene_string" "$rcluster_string" > 'partition_finder_N'$raw_samp'_raxml_bic.cfg'
    echo "$bic_string" "$gene_string" "$greedy_string" > 'partition_finder_N'$raw_samp'_greedy_bic.cfg'
    echo "$aicc_string" "$gene_string" "$rcluster_string" > 'partition_finder_N'$raw_samp'_raxml_aicc.cfg'
    echo "$aicc_string" "$gene_string" "$greedy_string" > 'partition_finder_N'$raw_samp'_greedy_aicc.cfg'


    # run all combinations of partitioning analysis on that file and the associated iqtree analysis
    # label in the output I'm saving which file this is
    echo $f >> partitioning_treebuilding_output.txt

#________run partition finder - raxml aic________#
    # clean out partitioning results before doing new partitioning
    rm -rf 'part_N'$raw_samp'/analysis'
    rm -f 'part_N'$raw_samp'/log.txt'
    cp 'partition_finder_N'$raw_samp'_raxml_aic.cfg' 'part_N'$raw_samp'/partition_finder.cfg'
    python $current_path'partitionfinder-2.1.1/PartitionFinderProtein.py' 'part_N'$raw_samp'/' --raxml -p 4
    # copy nexus best scheme to a partition file
    awk '/nexus/,/end/' 'part_N'$raw_samp'/analysis/best_scheme.txt' > 'part_N'$raw_samp'/best_scheme_nexus.txt'
    # save information criteria output from partition finder
    echo $(grep search 'part_N'$raw_samp'/analysis/best_scheme.txt') >> partitioning_treebuilding_output.txt
    echo $(grep Scheme 'part_N'$raw_samp'/analysis/best_scheme.txt') >> partitioning_treebuilding_output.txt
#________run iqtree with partition file - raxml aic________#
    rm -f 'part_N'$raw_samp'/best_scheme.txt.*'
    #IQTREE_PartitionFinder
    $current_path'iqtree-1.6.12-Linux/bin/iqtree' -s 'part_N'$raw_samp'/temp_align.phy' -seed 123456789 -nt 4 -spp 'part_N'$raw_samp'/best_scheme_nexus.txt' -bb 1000
    echo $(grep -A 12 MAXIMUM 'part_N'$raw_samp'/best_scheme_nexus.txt.iqtree') >> partitioning_treebuilding_output.txt
    echo $(grep -A 2 newick 'part_N'$raw_samp'/best_scheme_nexus.txt.iqtree' | head -n 3 | tail -n 1) >> partitioning_treebuilding_output.txt
    echo $(grep -A 5 CONSENSUS 'part_N'$raw_samp'/best_scheme_nexus.txt.iqtree') >> partitioning_treebuilding_output.txt
    echo $(grep -A 2 newick 'part_N'$raw_samp'/best_scheme_nexus.txt.iqtree' | tail -n 1) >> partitioning_treebuilding_output.txt


#________run partition finder - raxml aicc________#
    # clean out partitioning results before doing new partitioning
    rm -rf 'part_N'$raw_samp'/analysis'
    rm -f 'part_N'$raw_samp'/log.txt'
    cp 'partition_finder_N'$raw_samp'_raxml_aicc.cfg' 'part_N'$raw_samp'/partition_finder.cfg'
    python $current_path'partitionfinder-2.1.1/PartitionFinderProtein.py' 'part_N'$raw_samp'/' --raxml -p 4
    # copy nexus best scheme to a partition file
    awk '/nexus/,/end/' 'part_N'$raw_samp'/analysis/best_scheme.txt' > 'part_N'$raw_samp'/best_scheme_nexus.txt'
    # save information criteria output from partition finder
    echo $(grep search 'part_N'$raw_samp'/analysis/best_scheme.txt') >> partitioning_treebuilding_output.txt
#________run iqtree with partition file - raxml aicc________#
    rm -f 'part_N'$raw_samp'/best_scheme.txt.*'
    #IQTREE_PartitionFinder
    $current_path'iqtree-1.6.12-Linux/bin/iqtree' -s 'part_N'$raw_samp'/temp_align.phy' -seed 123456789 -nt 4 -spp 'part_N'$raw_samp'/best_scheme_nexus.txt' -bb 1000
    echo $(grep -A 12 MAXIMUM 'part_N'$raw_samp'/best_scheme_nexus.txt.iqtree') >> partitioning_treebuilding_output.txt
    echo $(grep -A 2 newick 'part_N'$raw_samp'/best_scheme_nexus.txt.iqtree' | head -n 3 | tail -n 1) >> partitioning_treebuilding_output.txt
    echo $(grep -A 5 CONSENSUS 'part_N'$raw_samp'/best_scheme_nexus.txt.iqtree') >> partitioning_treebuilding_output.txt
    echo $(grep -A 2 newick 'part_N'$raw_samp'/best_scheme_nexus.txt.iqtree' | tail -n 1) >> partitioning_treebuilding_output.txt


#________run partition finder - raxml bic________#
    # clean out partitioning results before doing new partitioning
    rm -rf 'part_N'$raw_samp'/analysis'
    rm -f 'part_N'$raw_samp'/log.txt'
    cp 'partition_finder_N'$raw_samp'_raxml_bic.cfg' 'part_N'$raw_samp'/partition_finder.cfg'
    python $current_path'partitionfinder-2.1.1/PartitionFinderProtein.py' 'part_N'$raw_samp'/' --raxml -p 4
    # copy nexus best scheme to a partition file
    awk '/nexus/,/end/' 'part_N'$raw_samp'/analysis/best_scheme.txt' > 'part_N'$raw_samp'/best_scheme_nexus.txt'
    # save information criteria output from partition finder
    echo $(grep search 'part_N'$raw_samp'/analysis/best_scheme.txt') >> partitioning_treebuilding_output.txt
#________run iqtree with partition file - raxml aicc________#
    rm -f 'part_N'$raw_samp'/best_scheme.txt.*'
    #IQTREE_PartitionFinder
    $current_path'iqtree-1.6.12-Linux/bin/iqtree' -s 'part_N'$raw_samp'/temp_align.phy' -seed 123456789 -nt 4 -spp 'part_N'$raw_samp'/best_scheme_nexus.txt' -bb 1000
    echo $(grep -A 12 MAXIMUM 'part_N'$raw_samp'/best_scheme_nexus.txt.iqtree') >> partitioning_treebuilding_output.txt
    echo $(grep -A 2 newick 'part_N'$raw_samp'/best_scheme_nexus.txt.iqtree' | head -n 3 | tail -n 1) >> partitioning_treebuilding_output.txt
    echo $(grep -A 5 CONSENSUS 'part_N'$raw_samp'/best_scheme_nexus.txt.iqtree') >> partitioning_treebuilding_output.txt
    echo $(grep -A 2 newick 'part_N'$raw_samp'/best_scheme_nexus.txt.iqtree' | tail -n 1) >> partitioning_treebuilding_output.txt



#________run partition finder - greedy aic________#
    # clean out partitioning results before doing new partitioning
    rm -rf 'part_N'$raw_samp'/analysis'
    rm -f 'part_N'$raw_samp'/log.txt'
    cp 'partition_finder_N'$raw_samp'_greedy_aic.cfg' 'part_N'$raw_samp'/partition_finder.cfg'
    python $current_path'partitionfinder-2.1.1/PartitionFinderProtein.py' 'part_N'$raw_samp'/' -p 4
    # copy nexus best scheme to a partition file
    awk '/nexus/,/end/' 'part_N'$raw_samp'/analysis/best_scheme.txt' > 'part_N'$raw_samp'/best_scheme_nexus.txt'
    # save information criteria output from partition finder
    echo $(grep search 'part_N'$raw_samp'/analysis/best_scheme.txt') >> partitioning_treebuilding_output.txt
    echo $(grep Scheme 'part_N'$raw_samp'/analysis/best_scheme.txt') >> partitioning_treebuilding_output.txt
#________run iqtree with partition file - greedy aic________#
    rm -f 'part_N'$raw_samp'/best_scheme.txt.*'
    #IQTREE_PartitionFinder
    $current_path'iqtree-1.6.12-Linux/bin/iqtree' -s 'part_N'$raw_samp'/temp_align.phy' -seed 123456789 -nt 4 -spp 'part_N'$raw_samp'/best_scheme_nexus.txt' -bb 1000
    echo $(grep -A 12 MAXIMUM 'part_N'$raw_samp'/best_scheme_nexus.txt.iqtree') >> partitioning_treebuilding_output.txt
    echo $(grep -A 2 newick 'part_N'$raw_samp'/best_scheme_nexus.txt.iqtree' | head -n 3 | tail -n 1) >> partitioning_treebuilding_output.txt
    echo $(grep -A 5 CONSENSUS 'part_N'$raw_samp'/best_scheme_nexus.txt.iqtree') >> partitioning_treebuilding_output.txt
    echo $(grep -A 2 newick 'part_N'$raw_samp'/best_scheme_nexus.txt.iqtree' | tail -n 1) >> partitioning_treebuilding_output.txt


#________run partition finder - greedy aicc________#
    # clean out partitioning results before doing new partitioning
    rm -rf 'part_N'$raw_samp'/analysis'
    rm -f 'part_N'$raw_samp'/log.txt'
    cp 'partition_finder_N'$raw_samp'_greedy_aicc.cfg' 'part_N'$raw_samp'/partition_finder.cfg'
    python $current_path'partitionfinder-2.1.1/PartitionFinderProtein.py' 'part_N'$raw_samp'/' -p 4
    # copy nexus best scheme to a partition file
    awk '/nexus/,/end/' 'part_N'$raw_samp'/analysis/best_scheme.txt' > 'part_N'$raw_samp'/best_scheme_nexus.txt'
    # save information criteria output from partition finder
    echo $(grep search 'part_N'$raw_samp'/analysis/best_scheme.txt') >> partitioning_treebuilding_output.txt
#________run iqtree with partition file - greedy aicc________#
    rm -f 'part_N'$raw_samp'/best_scheme.txt.*'
    #IQTREE_PartitionFinder
    $current_path'iqtree-1.6.12-Linux/bin/iqtree' -s 'part_N'$raw_samp'/temp_align.phy' -seed 123456789 -nt 4 -spp 'part_N'$raw_samp'/best_scheme_nexus.txt' -bb 1000
    echo $(grep -A 12 MAXIMUM 'part_N'$raw_samp'/best_scheme_nexus.txt.iqtree') >> partitioning_treebuilding_output.txt
    echo $(grep -A 2 newick 'part_N'$raw_samp'/best_scheme_nexus.txt.iqtree' | head -n 3 | tail -n 1) >> partitioning_treebuilding_output.txt
    echo $(grep -A 5 CONSENSUS 'part_N'$raw_samp'/best_scheme_nexus.txt.iqtree') >> partitioning_treebuilding_output.txt
    echo $(grep -A 2 newick 'part_N'$raw_samp'/best_scheme_nexus.txt.iqtree' | tail -n 1) >> partitioning_treebuilding_output.txt


#________run partition finder - greedy bic________#
    # clean out partitioning results before doing new partitioning
    rm -rf 'part_N'$raw_samp'/analysis'
    rm -f 'part_N'$raw_samp'/log.txt'
    cp 'partition_finder_N'$raw_samp'_greedy_bic.cfg' 'part_N'$raw_samp'/partition_finder.cfg'
    python $current_path'partitionfinder-2.1.1/PartitionFinderProtein.py' 'part_N'$raw_samp'/' -p 4
    # copy nexus best scheme to a partition file
    awk '/nexus/,/end/' 'part_N'$raw_samp'/analysis/best_scheme.txt' > 'part_N'$raw_samp'/best_scheme_nexus.txt'
    # save information criteria output from partition finder
    echo $(grep search 'part_N'$raw_samp'/analysis/best_scheme.txt') >> partitioning_treebuilding_output.txt
#________run iqtree with partition file - greedy aicc________#
    rm -f 'part_N'$raw_samp'/best_scheme.txt.*'
    #IQTREE_PartitionFinder
    $current_path'iqtree-1.6.12-Linux/bin/iqtree' -s 'part_N'$raw_samp'/temp_align.phy' -seed 123456789 -nt 4 -spp 'part_N'$raw_samp'/best_scheme_nexus.txt' -bb 1000
    echo $(grep -A 12 MAXIMUM 'part_N'$raw_samp'/best_scheme_nexus.txt.iqtree') >> partitioning_treebuilding_output.txt
    echo $(grep -A 2 newick 'part_N'$raw_samp'/best_scheme_nexus.txt.iqtree' | head -n 3 | tail -n 1) >> partitioning_treebuilding_output.txt
    echo $(grep -A 5 CONSENSUS 'part_N'$raw_samp'/best_scheme_nexus.txt.iqtree') >> partitioning_treebuilding_output.txt
    echo $(grep -A 2 newick 'part_N'$raw_samp'/best_scheme_nexus.txt.iqtree' | tail -n 1) >> partitioning_treebuilding_output.txt






#________iqtree without partition file________#
    rm -f 'part_N'$raw_samp'/best_scheme.txt.*'
    #IQTREE_Default, can be run once on the alignment without all the partition combinations:
    $current_path'iqtree-1.6.12-Linux/bin/iqtree' -s 'part_N'$raw_samp'/temp_align.phy' -seed 123456789 -nt 4 -mrate G -bb 1000
    echo $(grep -A 12 MAXIMUM 'part_'$dir'/temp_align.phy.iqtree') >> partitioning_treebuilding_output.txt
    echo $(grep -A 2 newick 'part_'$dir'/temp_align.phy.iqtree' | head -n 3 | tail -n 1) >> partitioning_treebuilding_output.txt
    echo $(grep -A 5 CONSENSUS 'part_'$dir'/temp_align.phy.iqtree') >> partitioning_treebuilding_output.txt
    echo $(grep -A 2 newick 'part_'$dir'/temp_align.phy.iqtree' | tail -n 1) >> partitioning_treebuilding_output.txt




done # done all combos on that file