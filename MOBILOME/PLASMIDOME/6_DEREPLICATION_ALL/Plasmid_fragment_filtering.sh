#!/bin/bash
#SBATCH --job-name=Plasmid_fragment_filtering
#SBATCH --output=Plasmid_fragment_filtering.out
#SBATCH --mem=10gb
#SBATCH --time=00:19:00
#SBATCH --cpus-per-task=16
#SBATCH --export=NONE
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

plasmid_clusters=$1 #path to TSV file with all the plasmid clusters
all_plasmids=$2 #path to FASTA file with all combined plasmids (before dereplication)

# Clean environment, load modules
module purge; ml R; module list

# Select clusters of circular plasmids with plasmid fragments
cat rep_seqs_all_plasmids_renamed_clusters.tsv | cut -f2 | grep "IMG\|Circular" | grep fragment > clustering_fragments.txt

# Execute the Python script to filter out plasmid fragments
Rscript Plasmid_fragment_filtering.R clustering_fragments.txt

# Clean environment, load modules
module purge; ml SeqKit; module list

# Get the FASTA and TXT files with the final set of plasmids 
seqkit grep -n -v -f Plasmid_fragments_to_remove.txt $all_plasmids > plasmids.fna
grep ">" plasmids.fna | cut -f2 -d ">"  > plasmids.txt

# Set permissions 
chmod 440 plasmids.fna plasmids.txt

