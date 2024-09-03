#!/bin/bash
#SBATCH --job-name=AMR_annotation_virus
#SBATCH --output=AMR_annotation_virus.out
#SBATCH --mem=80gb
#SBATCH --time=07:59:00
#SBATCH --cpus-per-task=24
#SBATCH --open-mode=truncate
#SBATCH --partition=himem

protein_file=$1 #path to FASTA file with proteins predicted by prodigal-gv from viral contigs

echo -e '\n---- RUNNING RGI on CARD Database ----'
module purge; ml Anaconda3; module list
source activate /scratch/hb-llnext/conda_envs/RGI_env; conda list 

rgi load \
    --card_json /scratch/hb-tifn/DBs/CARD_2023_06/card.json \
    --local

rgi main \
    --input_sequence ${protein_file} \
    --output_file RGI_CARD \
    -t protein \
    --local \
    --clean \
    -n ${SLURM_CPUS_PER_TASK}
