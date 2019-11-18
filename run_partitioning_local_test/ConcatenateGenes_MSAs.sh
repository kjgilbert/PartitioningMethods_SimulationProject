#!/bin/bash


lopho_dir='results_lopho/sim_match_lophotroch/MSA/'
myria_dir='results_myria/sim_match_myriapod/MSA/'


## concatenate 10 genes

STR=''
##STR2=''
rm sampledMSAs/config_lopho_10genes.txt; touch sampledMSAs/config_lopho_10genes.txt
for i in `seq 1 10`;
	do
	file=${lopho_dir}'MSA_'$i'_aa.fa '
    STR=$STR$file
##    STR2=$STR2'$'$i' '
    head -n 2 $file | tail -n 1 | wc -m >> sampledMSAs/config_lopho_10genes.txt
done
##echo $STR 
##echo $STR2 
paste $STR | awk '{if (NR%2==0) {print $1 $2 $3 $4 $5 $6 $7 $8 $9 $10} else {print $1}}' > sampledMSAs/MSA_lopho_10genes.fa
##paste $STR | awk '{if (NR%2==0) {for (i=1; i<=10; i++) print $i} else {print $1}}' > sampledMSAs/MSA_lopho_10genes.fa
##paste -d '\0' $STR | sed 's/>[A-Z]*//' > sampledMSAs/MSA_lopho_10genes.fa
 
STR=''
rm sampledMSAs/config_myria_10genes.txt; touch sampledMSAs/config_myria_10genes.txt
for i in `seq 1 10`;
	do
	file=${myria_dir}'MSA_'$i'_aa.fa '
    STR=$STR$file
    head -n 2 $file | tail -n 1 | wc -m >> sampledMSAs/config_myria_10genes.txt
done
paste $STR | awk '{if (NR%2==0) {print $1 $2 $3 $4 $5 $6 $7 $8 $9 $10} else {print $1}}' > sampledMSAs/MSA_myria_10genes.fa
 
 
 
 
## concatenate 20 genes

STR=''
rm sampledMSAs/config_lopho_20genes.txt; touch sampledMSAs/config_lopho_20genes.txt
for i in `seq 1 20`;
	do
	file=${lopho_dir}'MSA_'$i'_aa.fa '
    STR=$STR$file
    head -n 2 $file | tail -n 1 | wc -m >> sampledMSAs/config_lopho_20genes.txt
done
paste $STR | awk '{if (NR%2==0) {print $1 $2 $3 $4 $5 $6 $7 $8 $9 $10 $1 $2 $3 $4 $5 $6 $7 $8 $9 $10 $1 $2 $3 $4 $5 $6 $7 $8 $9 $10 $11 $12 $13 $14 $15 $16 $17 $18 $19 $20} else {print $1}}' > sampledMSAs/MSA_lopho_20genes.fa
 
 
STR=''
rm sampledMSAs/config_myria_20genes.txt; touch sampledMSAs/config_myria_20genes.txt
for i in `seq 1 20`;
	do
	file=${myria_dir}'MSA_'$i'_aa.fa '
    STR=$STR$file
    head -n 2 $file | tail -n 1 | wc -m >> sampledMSAs/config_myria_20genes.txt
done
paste $STR | awk '{if (NR%2==0) {print $1 $2 $3 $4 $5 $6 $7 $8 $9 $10 $1 $2 $3 $4 $5 $6 $7 $8 $9 $10 $1 $2 $3 $4 $5 $6 $7 $8 $9 $10 $11 $12 $13 $14 $15 $16 $17 $18 $19 $20} else {print $1}}' > sampledMSAs/MSA_myria_20genes.fa
 




## concatenate 40 genes


STR=''
rm sampledMSAs/config_lopho_40genes.txt; touch sampledMSAs/config_lopho_40genes.txt
for i in `seq 1 40`;
	do
	file=${lopho_dir}'MSA_'$i'_aa.fa '
    STR=$STR$file
    head -n 2 $file | tail -n 1 | wc -m >> sampledMSAs/config_lopho_40genes.txt
done
paste $STR | awk '{if (NR%2==0) {print $1 $2 $3 $4 $5 $6 $7 $8 $9 $10 $1 $2 $3 $4 $5 $6 $7 $8 $9 $10 $1 $2 $3 $4 $5 $6 $7 $8 $9 $10 $11 $12 $13 $14 $15 $16 $17 $18 $19 $20 $21 $22 $23 $24 $25 $26 $27 $28 $29 $30 $31 $32 $33 $34 $35 $36 $37 $38 $39 $40} else {print $1}}' > sampledMSAs/MSA_lopho_40genes.fa
 
 
STR=''
rm sampledMSAs/config_myria_40genes.txt; touch sampledMSAs/config_myria_40genes.txt
for i in `seq 1 40`;
	do
	file=${myria_dir}'MSA_'$i'_aa.fa '
    STR=$STR$file
    head -n 2 $file | tail -n 1 | wc -m >> sampledMSAs/config_myria_40genes.txt
done
paste $STR | awk '{if (NR%2==0) {print $1 $2 $3 $4 $5 $6 $7 $8 $9 $10 $1 $2 $3 $4 $5 $6 $7 $8 $9 $10 $1 $2 $3 $4 $5 $6 $7 $8 $9 $10 $11 $12 $13 $14 $15 $16 $17 $18 $19 $20 $21 $22 $23 $24 $25 $26 $27 $28 $29 $30 $31 $32 $33 $34 $35 $36 $37 $38 $39 $40} else {print $1}}' > sampledMSAs/MSA_myria_40genes.fa


 
 
## concatenate 80 genes


STR=''
rm sampledMSAs/config_lopho_80genes.txt; touch sampledMSAs/config_lopho_80genes.txt
for i in `seq 1 80`;
	do
	file=${lopho_dir}'MSA_'$i'_aa.fa '
    STR=$STR$file
    head -n 2 $file | tail -n 1 | wc -m >> sampledMSAs/config_lopho_80genes.txt
done
paste $STR | awk '{if (NR%2==0) {print $1 $2 $3 $4 $5 $6 $7 $8 $9 $10 $1 $2 $3 $4 $5 $6 $7 $8 $9 $10 $1 $2 $3 $4 $5 $6 $7 $8 $9 $10 $11 $12 $13 $14 $15 $16 $17 $18 $19 $20 $21 $22 $23 $24 $25 $26 $27 $28 $29 $30 $31 $32 $33 $34 $35 $36 $37 $38 $39 $40 $41 $42 $43 $44 $45 $46 $47 $48 $49 $50 $51 $52 $53 $54 $55 $56 $57 $58 $59 $60 $61 $62 $63 $64 $65 $66 $67 $68 $69 $70 $71 $72 $73 $74 $75 $76 $77 $78 $79 $80} else {print $1}}' > sampledMSAs/MSA_lopho_80genes.fa
 
 
STR=''
rm sampledMSAs/config_myria_80genes.txt; touch sampledMSAs/config_myria_80genes.txt
for i in `seq 1 80`;
	do
	file=${myria_dir}'MSA_'$i'_aa.fa '
    STR=$STR$file
    head -n 2 $file | tail -n 1 | wc -m >> sampledMSAs/config_myria_80genes.txt
done
paste $STR | awk '{if (NR%2==0) {print $1 $2 $3 $4 $5 $6 $7 $8 $9 $10 $1 $2 $3 $4 $5 $6 $7 $8 $9 $10 $1 $2 $3 $4 $5 $6 $7 $8 $9 $10 $11 $12 $13 $14 $15 $16 $17 $18 $19 $20 $21 $22 $23 $24 $25 $26 $27 $28 $29 $30 $31 $32 $33 $34 $35 $36 $37 $38 $39 $40 $41 $42 $43 $44 $45 $46 $47 $48 $49 $50 $51 $52 $53 $54 $55 $56 $57 $58 $59 $60 $61 $62 $63 $64 $65 $66 $67 $68 $69 $70 $71 $72 $73 $74 $75 $76 $77 $78 $79 $80} else {print $1}}' > sampledMSAs/MSA_myria_80genes.fa
