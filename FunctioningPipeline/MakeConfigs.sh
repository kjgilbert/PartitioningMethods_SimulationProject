current_path=$1
run_path=$2
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


done