#!/bin/bash
#SBATCH --job-name=AMR_annotation
#SBATCH --output=AMR_annotation_%A_%a.out
#SBATCH --mem=50gb
#SBATCH --time=04:59:00
#SBATCH --cpus-per-task=24
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

der_bins_dir=$1 #directory with the dereplicated bins/MAGs
output_dir=$2 #output directory
sample_list=${der_bins_dir}/$3 #file with the list of all bins
SAMPLE_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${sample_list})

mkdir -p ${output_dir}/AMRFinderPlus ${output_dir}/RGI_CARD

echo -e '\n-------------------- WORKING WITH '${SAMPLE_ID}' SAMPLE --------------------'


echo -e '\n---- RUNNING AMR Finder Plus ----'

# Clean environment, load modules and activate conda environment
module purge; ml Anaconda3; module list
source activate /home2/p304845/Conda_envs/AMRFinderPlus_env/; conda list 

# Run AMR Finder (without plus)
amrfinder \
    -n ${der_bins_dir}/${SAMPLE_ID}.fa \
    --output ${output_dir}/AMRFinderPlus/${SAMPLE_ID}_AMRFinderPlus_annotation.csv 

conda deactivate

echo -e '\n---- RUNNING RGI on CARD Database ----'
module purge; ml Anaconda3; module list
source activate /scratch/hb-llnext/conda_envs/RGI_env; conda list 

rgi load \
    --card_json /scratch/hb-tifn/DBs/CARD_2023_06/card.json \
    --local

rgi main \
    --input_sequence ${der_bins_dir}/${SAMPLE_ID}.fa \
    --output_file ${output_dir}/RGI_CARD/${SAMPLE_ID}_RGI_CARD_annotation \
    --local \
    --clean \
    --include_loose \
    -n ${SLURM_CPUS_PER_TASK}

# no --low-quality argument here as we have high-quality genomes
