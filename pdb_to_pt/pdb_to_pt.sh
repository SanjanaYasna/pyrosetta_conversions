#!/bin/bash
#SBATCH --partition phyq
#SBATCH --job-name=non_enz
#SBATCH --output=non_enz_pt.txt
# Script to download files from RCSB http file download services.
# Use the -h switch to get help on usage.
#SBATCH --nodes=1 --cpus-per-task=30

cd  /Users/robsonlab/Teetly/get_data_pyrosetta/pdb_to_pt
source /Users/robsonlab/Teetly/.teetlyrc
conda activate gcn_heal
conda run -n gcn_heal python non_enz.py
