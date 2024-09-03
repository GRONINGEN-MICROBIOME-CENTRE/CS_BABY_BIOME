#!/bin/bash
#SBATCH --job-name=Combine_all_plasmids
#SBATCH --output=Combine_all_plasmids.out
#SBATCH --mem=2gb
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=4
#SBATCH --export=NONE
#SBATCH --open-mode=truncate
#SBATCH --partition=himem

circular_plasmids=$1 #path to FASTA file with all the circular plasmids
circular_plasmids_clusters=$2 #path to TSV file with clusters after dereplication
plasmid_fragments=$3 #path to FASTA file with all the plasmi fragments
plasmid_fragments_clusters=$4 #path to TSV file with clusters after dereplication

# Clean environment, load modules
module purge; ml SeqKit; module list

# Extract the representative sequences of the circular and plasmid fragments after dereplication (separately)
cat $circular_plasmids_clusters | cut -f1 > rep_seqs_circular_plasmids.txt
seqkit grep -f rep_seqs_circular_plasmids.txt $circular_plasmids > rep_seqs_circular_plasmids.fna
cat $plasmid_fragments_clusters | cut -f1 > rep_seqs_plasmid_fragments.txt
seqkit grep -f rep_seqs_plasmid_fragments.txt $plasmid_fragments > rep_seqs_plasmid_fragments.fna

# Combine the circular and linear plasmids (dereplicated separately)
cat rep_seqs_circular_plasmids.fna rep_seqs_plasmid_fragments.fna > rep_seqs_all_plasmids.fna
