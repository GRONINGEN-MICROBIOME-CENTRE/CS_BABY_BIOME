#!/bin/bash
#SBATCH --job-name=iPHOP_host_assignment
#SBATCH --output=iPHOP_host_assignment_%A_%a.out
#SBATCH --mem=50gb
#SBATCH --time=04:59:00
#SBATCH --cpus-per-task=32
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

contig_dir=$1 #path to directory with the split FASTA files with predicted viral contigs
file_list=${contig_dir}/$2 #file with the list of all split FASTA filenames in the directory
contig_file=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${file_list})
contig_file_name=$(basename "$contig_file" | cut -d . -f1)

echo '-------------------- WORKING WITH '${contig_file}' FILE --------------------'

# Clean environment, load modules and activate conda environment
module purge; ml Anaconda3; module list
source activate /scratch/hb-llnext/conda_envs/iphop133_env; conda list

mkdir -p iPHOP_RESULTS iPHOP_RESULTS/${contig_file} iPHOP_RESULTS_ENRICHED iPHOP_RESULTS_ENRICHED/${contig_file}

iphop predict \
	--fa_file ${contig_dir}/${contig_file} \
	--db_dir /scratch/hb-llnext/databases/iphop133_db/Aug_2023_pub_rw \
	--out_dir iPHOP_RESULTS/${contig_file} \
	--num_threads ${SLURM_CPUS_PER_TASK}

iphop predict \
	--fa_file ${contig_dir}/${contig_file} \
	--db_dir /scratch/hb-llnext/databases/iphop133_db/Aug_2023_pub_rw_CS_Baby_Biome \
	--out_dir iPHOP_RESULTS_ENRICHED/${contig_file} \
	--num_threads ${SLURM_CPUS_PER_TASK}
	 
conda deactivate
