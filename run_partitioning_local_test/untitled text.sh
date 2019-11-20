#!/bin/bash


lopho_dir='results_lopho/sim_match_lophotroch/MSA/'
myria_dir='results_myria/sim_match_myriapod/MSA/'


## concatenate genes

sample_sizes=('10' '20' '40' '80')

for samp in ${sample_sizes[*]};
	do	
	
	# lophotrochozoa
	STR=''
	rm -f 'sampledMSAs/config_lopho_'$samp'genes.txt'
	touch 'sampledMSAs/config_lopho_'$samp'genes.txt'
	for i in `seq 1 $samp`;
		do
		file=${lopho_dir}'MSA_'$i'_aa.fa '
		STR=$STR$file
		head -n 2 $file | tail -n 1 | wc -m >> 'sampledMSAs/config_lopho_'$samp'genes.txt'
	done
	paste $STR | awk -v f=1 -v t=10 '{ if (NR%2==0) { for(i=f; i<=t;i++) {printf("%s%s",$i,(i==t)?"\n":OFS)}} else {print $1}}' > 'sampledMSAs/MSA_lopho_'$samp'genes.fa'
	
	# myriapoda
	STR=''
	rm -f 'sampledMSAs/config_myria_'$samp'genes.txt'
	touch 'sampledMSAs/config_myria_'$samp'genes.txt'
	for i in `seq 1 $samp`;
		do
		file=${myria_dir}'MSA_'$i'_aa.fa '
		STR=$STR$file
		head -n 2 $file | tail -n 1 | wc -m >> 'sampledMSAs/config_myria_'$samp'genes.txt'
	done
	paste $STR | awk -v f=1 -v t=10 '{ if (NR%2==0) { for(i=f; i<=t;i++) {printf("%s%s",$i,(i==t)?"\n":OFS)}} else {print $1}}' > 'sampledMSAs/MSA_myria_'$samp'genes.fa'

done