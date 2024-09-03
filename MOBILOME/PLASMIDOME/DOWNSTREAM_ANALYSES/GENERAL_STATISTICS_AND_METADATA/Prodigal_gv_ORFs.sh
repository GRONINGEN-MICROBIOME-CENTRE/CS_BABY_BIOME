#!/bin/bash
#SBATCH --job-name=Prodigal_ORFs
#SBATCH --output=Prodigal_ORFs.out
#SBATCH --mem=8gb
#SBATCH --time=00:29:00
#SBATCH --cpus-per-task=16
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

seqs=$1 #path to FASTA file with sequences
seqs_name="${vOTU_seqs%.*}" #extract filename without the extension

# Clean environment, load modules and activate conda environment
module purge; ml prodigal; module list

# Run Prodigal-gv
prodigal \
	-i $seqs \
	-d ${seqs_name}_genes.fna \
	-a ${seqs_name}_proteins.faa \
	-o ${seqs_name}_output.csv
