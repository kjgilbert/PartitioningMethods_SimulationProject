#!/bin/bash

outfile='partitioning_treebuilding_output.txt'


while read line
	do
	
	if grep -c search $outfile; then
        search_algo=$(head -n $line | tail -n 1 | tr " : " "\n" | tail -n 1)
    else
        echo not found
    fi
done <$outfile




while read line
	do
	if grep -c 'search' $line; then
        search_algo=$(head -n $line | tail -n 1 | tr " : " "\n" | tail -n 1)
    else
        echo 'not found'
    fi
done <partitioning_treebuilding_output.txt




    STR=''
    rm -f 'sampledMSAs/config_lopho_'$samp'genes.txt'
    touch 'sampledMSAs/config_lopho_'$samp'genes.txt'
    for i in `seq 1 $samp`;
        do
        file=${lopho_dir}'MSA_'$i'_aa.fa '
        STR=$STR$file
        head -n 2 $file | tail -n 1 | wc -m >> 'sampledMSAs/config_lopho_'$samp'genes.txt'
    done
    paste $STR | awk -v f=1 -v t=$samp '{ if (NR%2==0) { for(i=f; i<=t;i++) {printf("%s%s",$i,(i==t)?"\n":OFS="")}} else {print $1}}' > 'sampledMSAs/MSA_lopho_'$samp'genes.fa'



grep -c Scheme




while read line
	do
	
	if grep -q search $outfile; then
        search_algo=$(echo $f | tr "_" "\n" | tail -n 2 | head -n 1)
    else
        echo not found
    fi
	
	it=$(expr $it + 1)
	gene_start=$(expr $gene_end + 1)
	gene_end=$(expr $line - 1 + $prev_gene_end)
	prev_gene_end=$(expr $line - 1 + $prev_gene_end)
	gene_string=$gene_string'gene'$it' = '$gene_start'-'$gene_end$';\n'
done <$outfile
