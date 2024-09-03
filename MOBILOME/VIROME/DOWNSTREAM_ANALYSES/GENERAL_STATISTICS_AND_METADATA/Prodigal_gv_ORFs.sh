#!/bin/bash
#SBATCH --job-name=Prodigal-gv_ORFs
#SBATCH --output=Prodigal-gv_ORFs.out
#SBATCH --mem=8gb
#SBATCH --time=00:29:00
#SBATCH --cpus-per-task=16
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

vOTU_seqs=$1 #path to FASTA file with viral sequences
vOTU_seqs_name="${vOTU_seqs%.*}" #extract filename without the extension

# Clean environment, load modules and activate conda environment
module purge; ml prodigal-gv; module list

# Run Prodigal-gv
prodigal-gv \
	-i $vOTU_seqs \
	-d ${vOTU_seqs_name}_genes.fna \
	-a ${vOTU_seqs_name}_proteins.faa \
	-o ${vOTU_seqs_name}_output.csv
