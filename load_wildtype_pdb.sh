#!/bin/bash
#SBATCH --partition phyq
#SBATCH --job-name=ec
#SBATCH --output=ec_clean_2.txt
# Script to download files from RCSB http file download services.
# Use the -h switch to get help on usage.
#SBATCH --nodes=1 --cpus-per-task=30



cd  /Users/robsonlab/Teetly/get_data_pyrosetta
source /Users/robsonlab/Teetly/.teetlyrc
conda activate Slusky
conda run -n Slusky python load_wildtype_ec.py

## OUTFILES: ec_clean.txt and non_enz_clean.txt