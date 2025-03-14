#!/bin/bash
#SBATCH --job-name=HMMER_annotation_PHROGs
#SBATCH --output=HMMER_annotation_PHROGs.out
#SBATCH --mem=40gb
#SBATCH --time=19:59:00
#SBATCH --cpus-per-task=32
#SBATCH --open-mode=truncate
#SBATCH --partition=himem

HMM_DB=$1 #directory with HMM profiles of PHROGs
viral_proteins=$2 #FASTA file with predicted viral proteins

mkdir -p HMMER_RESULTS_PHROGs OUTPUT_files

# Clean environment and load modules 
module purge; ml HMMER; module list

# Analyze the viral protein sequences using our HMM DB
hmmsearch \
    -o HMMER_RESULTS_PHROGs/HMMER_result \
    --tblout HMMER_RESULTS_PHROGs/HMMER_table \
    -E 1e-5 \
    --cpu ${SLURM_CPUS_PER_TASK} \
    ${HMM_DB}/PHROGS_new.hmm \
    $viral_proteins
