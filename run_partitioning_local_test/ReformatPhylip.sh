#!/bin/bash

### reformat the sequential phylip format to have a space between species name and base pairs, 
###  so that it works in partitionfinder 2

inDir="sampledMSAs/"

for i in ${inDir}*.phy;
    do
    sed "s/.\{10\}/& /" $i > temp_phylip.txt
    cp temp_phylip.txt $i
done
