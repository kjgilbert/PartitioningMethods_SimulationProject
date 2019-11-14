#!/bin/bash

### take each MSA file and sample the first 20, 40, 100, or all character columns
### rename and put into a new directory for the next step (partitioning and treebuilding)

inDirLopho="results_lopho/sim_match_lophotroch/MSA/"
inDirMyria="results_myria/sim_match_myriapod/MSA/"
outDir="sampledMSAs/"
msaName1="MSA_"
msaName2="_aa.fa"

## mkdir $outDir

for i in `seq 1 10`;
    do
    	
        inFileNameLo=${inDirLopho}${msaName1}$i${msaName2}
        inFileNameMy=${inDirMyria}${msaName1}$i${msaName2}
		
		# do the sampling
		
		# N=20
		outLo="_aa_lopho_N20.fa"
		outMy="_aa_myria_N20.fa"
        outFileNameLo=${outDir}${msaName1}$i${outLo}
        outFileNameMy=${outDir}${msaName1}$i${outMy}
		cut -c-20 $inFileNameLo > $outFileNameLo
		cut -c-20 $inFileNameMy > $outFileNameMy

		# N=40
		outLo="_aa_lopho_N40.fa"
		outMy="_aa_myria_N40.fa"
        outFileNameLo=${outDir}${msaName1}$i${outLo}
        outFileNameMy=${outDir}${msaName1}$i${outMy}
		cut -c-40 $inFileNameLo > $outFileNameLo
		cut -c-40 $inFileNameMy > $outFileNameMy

		# N=100
		outLo="_aa_lopho_N100.fa"
		outMy="_aa_myria_N100.fa"
        outFileNameLo=${outDir}${msaName1}$i${outLo}
        outFileNameMy=${outDir}${msaName1}$i${outMy}
		cut -c-100 $inFileNameLo > $outFileNameLo
		cut -c-100 $inFileNameMy > $outFileNameMy

		# N = all
		outLo="_aa_lopho_Nall.fa"
		outMy="_aa_myria_Nall.fa"
        outFileNameLo=${outDir}${msaName1}$i${outLo}
        outFileNameMy=${outDir}${msaName1}$i${outMy}
		cat $inFileNameLo > $outFileNameLo
		cat $inFileNameMy > $outFileNameMy

done
