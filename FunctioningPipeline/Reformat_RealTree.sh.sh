#!/bin/bash


# 62 lopho species
raw_real_lopho="results_lopho/sim_match_lophotroch/RealTree.nwk"
# 25 myria species
raw_real_myria="results_myria/sim_match_myriapod/RealTree.nwk"

# lophotrochozoa
for i in `seq 1 62`;
	do
	if (( $i < 10 ))
	then
	    num='0'$i
	    new_ID='S0000'$i'_00001'
	else
		num=$i
	    new_ID='S000'$i'_00001'
	fi
	curr_ID='SE0'$num
	sed -i 's/$curr_ID/$new_ID/g' $raw_real_lopho
done

# myriapoda

for i in `seq 1 25`;
	do
	if (( $i < 10 ))
	then
	    num='0'$i
	    new_ID='S0000'$i'_00001'
	else
		num=$i
	    new_ID='S000'$i'_00001'
	fi
	curr_ID='SE0'$num
	sed -i 's/$curr_ID/$new_ID/g' $raw_real_myria
done