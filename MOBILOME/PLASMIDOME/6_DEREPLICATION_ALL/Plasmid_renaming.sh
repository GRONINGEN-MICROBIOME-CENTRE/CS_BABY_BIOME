#!/bin/bash
#SBATCH --job-name=Plasmid_renaming
#SBATCH --output=Plasmid_renaming.out
#SBATCH --mem=40gb
#SBATCH --time=00:19:00
#SBATCH --cpus-per-task=16
#SBATCH --export=NONE
#SBATCH --open-mode=truncate
#SBATCH --partition=himem

plasmids=$1 #path to FASTA file with all the plasmids

# Clean environment, load modules
module purge; ml R; module list

# Execute the R script
Rscript Rename_plasmid_sequences.R $plasmids

# Set permissions
chmod 440 rep_seqs_all_plasmids_renamed.fna Renamed_plasmids.txt
