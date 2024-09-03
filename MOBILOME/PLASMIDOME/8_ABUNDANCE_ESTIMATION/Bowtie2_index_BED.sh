#!/bin/bash
#SBATCH --job-name=Bowtie2_index
#SBATCH --output=Bowtie2_index.out
#SBATCH --mem=8gb
#SBATCH --time=00:29:00
#SBATCH --cpus-per-task=4
#SBATCH --export=NONE
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

plasmids=$1 # FASTA file with all the plasmids

# Build bowtie2 index
module purge; ml Bowtie2; module list

bowtie2-build \
    $plasmids \
    plasmids_DB \
    --large-index \
    --threads ${SLURM_CPUS_PER_TASK}

# Clean environment and load modules
module purge; ml bioawk/1.0-GCC-11.2.0; module list

# Generate BED file
bioawk -c fastx '{print $name"\t0\t"length($seq)}' $plasmids > plasmids_DB.bed

chmod 440 plasmids_DB.bed
