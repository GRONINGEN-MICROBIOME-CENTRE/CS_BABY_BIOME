#!/bin/bash
#SBATCH --job-name=Bin_abundance_estimation
#SBATCH --output=Bin_abundance_estimation_%A_%a.out
#SBATCH --mem=60gb
#SBATCH --time=01:59:00
#SBATCH --cpus-per-task=24
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --partition=regular

clean_FASTQs=$1 #directory with clean FASTQ files
der_bins_dir=$2 #directory with the dereplicated bins/MAGs
output_dir=$3 #output directory
sample_list=${clean_FASTQs}/$4 #file with the list of all samples in the directory

SAMPLE_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${sample_list})

# Create output directory if it doesn't exist
mkdir -p $output_dir

echo -e '\n-------------------- WORKING WITH '${SAMPLE_ID}' SAMPLE --------------------'


echo -e '\n-------------------- EESTIMATION BIN abundances --------------------'

# Clean environment, load modules and activate conda environment
module purge; ml Anaconda3; module list
source activate /scratch/hb-llnext/conda_envs/CoverM_env; conda list 

coverm genome \
    -c ${clean_FASTQs}/${SAMPLE_ID}_kneaddata_paired_1.fastq.gz ${clean_FASTQs}/${SAMPLE_ID}_kneaddata_paired_2.fastq.gz \
    -f ${der_bins_dir}/*fa \
    -m rpkm \
    -o $output_dir/${SAMPLE_ID}_CoverM_bin_abundance.tsv \
    -t ${SLURM_CPUS_PER_TASK} \
    -v

