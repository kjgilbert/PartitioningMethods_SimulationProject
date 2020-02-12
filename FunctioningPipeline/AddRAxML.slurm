#!/bin/bash
#SBATCH --mail-user=kimberly.gilbert@unil.ch
#SBATCH --mail-type=fail
#SBATCH --job-name="add raxml"
#SBATCH --workdir=.
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=10G
#SBATCH --time=10-00:00
#SBATCH --partition=ax-long
#SBATCH --account=cdessim2_default



./Add_raxML_Analysis.sh 'New_Rep1_100-1000/' 12