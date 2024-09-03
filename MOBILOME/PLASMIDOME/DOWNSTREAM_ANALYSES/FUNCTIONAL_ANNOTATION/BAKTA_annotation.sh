#!/bin/bash
#SBATCH --job-name=BAKTA_ORF_annotation
#SBATCH --output=BAKTA_ORF_annotation.out
#SBATCH --mem=80gb
#SBATCH --time=12:59:00
#SBATCH --cpus-per-task=32
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

plasmid_genomes=$1 #path to FASTA file with vOTU sequences 
output_dir=$2 #output directory

# Clean environment, load modules and activate conda environment
module purge; ml Anaconda3; module list
source activate /scratch/hb-llnext/conda_envs/bakta_env/; conda list 

# Run Prodigal-gv (paralellized)
bakta $plasmid_genomes \
	-d /scratch/hb-llnext/databases/bakta_db/db \
	-o $output_dir \
	-p Plasmids \
	--keep-contig-headers \
	-v \
	-t ${SLURM_CPUS_PER_TASK}

conda deactivate
