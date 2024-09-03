#!/bin/bash
#SBATCH --job-name=Align_plasmids_BLAST
#SBATCH --output=Align_plasmids_BLAST.out
#SBATCH --mem=100gb
#SBATCH --time=04:59:00
#SBATCH --cpus-per-task=24
#SBATCH --open-mode=truncate
#SBATCH --partition=himem

plasmids=$1 #path to FASTA file with all the plasmids
plasmids_filename="$(basename "${plasmids}")" #extract filename
plasmids_filename="${plasmids_filename%.*}" #extract filename without the extension

# Load BLAST
module purge; ml BLAST+/2.13.0-gompi-2022a; ml list

#First, create a blast+ database:
makeblastdb \
	-in $plasmids \
	-dbtype nucl \
	-out $plasmids_filename

#Next, use megablast from blast+ package to perform all-vs-all blastn of sequences:
blastn \
	-query $plasmids \
	-db $plasmids_filename \
	-outfmt '6 std qlen slen' \
	-out ${plasmids_filename}.tsv \
	-num_threads ${SLURM_CPUS_PER_TASK}

