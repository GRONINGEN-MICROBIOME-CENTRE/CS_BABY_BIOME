#!/bin/bash
#SBATCH --job-name=AMR_annotation
#SBATCH --output=AMR_annotation.out
#SBATCH --mem=50gb
#SBATCH --time=02:59:00
#SBATCH --cpus-per-task=24
#SBATCH --open-mode=truncate
#SBATCH --partition=himem

contig_file=$1 #path to FASTA file with contigs

echo -e '\n---- RUNNING RGI on CARD Database ----'
module purge; ml Anaconda3; module list
source activate /scratch/hb-llnext/conda_envs/RGI_env; conda list 

rgi load \
    --card_json /scratch/hb-tifn/DBs/CARD_2023_06/card.json \
    --local

rgi main \
    --input_sequence ${contig_file} \
    --output_file RGI_CARD_annotation \
    --local \
    --clean \
    -n ${SLURM_CPUS_PER_TASK}
