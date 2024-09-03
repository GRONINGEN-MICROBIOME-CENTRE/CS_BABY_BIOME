#!/bin/bash
#SBATCH --job-name=Viral_detection_CS_Baby_Biome
#SBATCH --output=Viral_detection_CS_%A_%a.out
#SBATCH --mem=16gb
#SBATCH --time=16:59:00
#SBATCH --cpus-per-task=16
#SBATCH --export=NONE
#SBATCH --open-mode=truncate
#SBATCH --partition=himem

sample_dir=$1 #directory with the metaSPADEs scaffolds
sample_list=${sample_dir}$2 #file with the list of all samples in the directory
SAMPLE_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${sample_list} | cut -d "_" -f1)

echo '-------------------- WORKING WITH '${SAMPLE_ID}' SAMPLE --------------------'

echo '---- RUNNING VirSorter2 WITH '${SAMPLE_ID}' SAMPLE ----'

mkdir -p ${sample_dir}${SAMPLE_ID}/virome_discovery

# Clean environment, load modules and activate conda environment
module purge; ml Anaconda3; module list
source activate /scratch/hb-llnext/conda_envs/VirSorter2_env; conda list

# Run VirSorter2
virsorter run \
	-w ${sample_dir}${SAMPLE_ID}/virome_discovery/VirSorter2 \
	-i ${sample_dir}${SAMPLE_ID}_metaspades_scaffolds.fa \
	--min-length 10000 \
 	--keep-original-seq \
	--include-groups "dsDNAphage,NCLDV,ssDNA,lavidaviridae" \
 	--db-dir /scratch/hb-llnext/conda_envs/VirSorter2_env/db \
  	-j ${SLURM_CPUS_PER_TASK} \
	all

rm -r ${sample_dir}${SAMPLE_ID}/virome_discovery/VirSorter2/iter-0
rm -r ${sample_dir}${SAMPLE_ID}/virome_discovery/VirSorter2/log
rm ${sample_dir}${SAMPLE_ID}/virome_discovery/VirSorter2/config.yaml

echo '---- RUNNING DeepVirFinder WITH '${SAMPLE_ID}' SAMPLE ----'

mkdir -p ${sample_dir}${SAMPLE_ID}/virome_discovery

conda deactivate

# Clean environment, load modules and activate conda environment
module purge; ml Anaconda3; module list
source activate /scratch/hb-llnext/conda_envs/DeepVirFinder_env; conda list

# Run DeepVirFinder
python3 /scratch/hb-llnext/conda_envs/DeepVirFinder_env/DeepVirFinder/dvf.py \
	-i ${sample_dir}${SAMPLE_ID}_metaspades_scaffolds.fa \
	-o ${sample_dir}${SAMPLE_ID}/virome_discovery/DeepVirFinder \
	-l 10000

conda deactivate
