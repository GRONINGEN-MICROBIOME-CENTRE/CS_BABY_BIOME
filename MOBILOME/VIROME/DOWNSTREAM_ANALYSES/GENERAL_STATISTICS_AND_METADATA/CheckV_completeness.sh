#!/bin/bash
#SBATCH --job-name=CheckV_completeness
#SBATCH --output=CheckV_completeness_%A_%a.out
#SBATCH --mem=24gb
#SBATCH --time=00:59:00
#SBATCH --cpus-per-task=16
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --partition=regular

contig_dir=$1 #path to directory with the split FASTA files with predicted viral contigs
file_list=${contig_dir}/$2 #file with the list of all split FASTA filenames in the directory
contig_file=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${file_list})
contig_file_name=$(basename "$contig_file" | cut -d . -f1)

echo '-------------------- WORKING WITH '${contig_file}' FILE --------------------'

# Clean environment, load modules and activate conda environment
module purge; ml CheckV/1.0.1-foss-2021b-DIAMOND-2.1.8; ml list

mkdir -p CheckV_results CheckV_results/${contig_file_name}

# Run CheckV
checkv \
    end_to_end \
	  ${contig_dir}/${contig_file} \
	  CheckV_results/${contig_file_name} \
	  -t ${SLURM_CPUS_PER_TASK} \
	  -d /scratch/hb-llnext/databases/checkv-db-v1.5

