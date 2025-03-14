#!/bin/bash
#SBATCH --job-name=CRISPR_matching_blastn
#SBATCH --output=CRISPR_matching_blastn.out
#SBATCH --mem=80gb
#SBATCH --time=3-0
#SBATCH --cpus-per-task=32
#SBATCH --export=NONE
#SBATCH --open-mode=truncate
#SBATCH --partition=himem

CRISPR_spacers=$1 # path to FASTA file with CRISPR spacers
plasmids=$2 #path to FASTA file with PTU representative genomes

mkdir -p CRISPR_DB

# Load BLAST
module purge; ml BLAST+/2.13.0-gompi-2022a SeqKit; ml list

#First, create a blast+ database:
makeblastdb \
	-in $CRISPR_spacers \
	-dbtype nucl \
	-out CRISPR_DB/crispr_spacers_db

#Next, perform blastn comparisons (plasmids vs crispr spacers)
blastn \
	-query $plasmids \
	-db CRISPR_DB/crispr_spacers_db \
	-outfmt "6 qseqid sseqid pident length mismatch qstart qend sstart send" \
	-max_target_seqs 1000 \
	-word_size 8 \
	-dust no \
	-out PTU_crispr_results.tsv \
	-num_threads ${SLURM_CPUS_PER_TASK}

# Estimate length of CRISPR spacers 
seqkit fx2tab --length $CRISPR_spacers | cut -f 1,4 > CRISPR_length.tsv

