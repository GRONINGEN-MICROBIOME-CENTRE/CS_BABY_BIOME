#!/bin/bash
#SBATCH --job-name=combine_circular_plasmids
#SBATCH --output=combine_circular_plasmids.out
#SBATCH --mem=20gb
#SBATCH --time=00:59:00
#SBATCH --cpus-per-task=24
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

IMG_DB=$1 #path to FASTA file all the plasmids from the IMG/PR database
IMG_DB_data=$2 #path to TSV file with the metadata of all the plasmids from the IMG/PR database
circular_plasmids=$3 #path to FASTA file with all the circular plasmids identified by metaPlasmidSPades and SCAPP
circular_plasmids_geNomad=$4 #path to FASTA file with all the circular plasmids identified by geNomad

# Clean environment, load modules
module purge; ml SeqKit; module list

# Select only circular plasmids (with DTR) from human origin from the IMG DB
head -n 1 $IMG_DB_data && grep -h "Human" $IMG_DB_data | grep "Direct terminal repeat" > IMGPR_human_circular_plasmid_data.tsv
cut -f1 IMGPR_human_circular_plasmid_data.tsv > IMGPR_human_circular_plasmid_list.txt
awk -F '|' '/^>/ {sub(/\|.*/, "", $1); print $1; next} {print}' $IMG_DB > IMG_DB_formatted.fna
seqkit grep -f IMGPR_human_circular_plasmid_list.txt IMG_DB_formatted.fna > IMGPR_human_circular_plasmids.fna

# Combine circular plasmids from the different analyses
cat IMGPR_human_circular_plasmids.fna $circular_plasmids $circular_plasmids_geNomad > all_plasmid_complete_sequences.fna

