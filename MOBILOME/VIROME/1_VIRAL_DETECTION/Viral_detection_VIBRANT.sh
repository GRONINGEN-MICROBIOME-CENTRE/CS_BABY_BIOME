#!/bin/bash
#SBATCH --job-name=Viral_detection_VIBRANT 
#SBATCH --output=Viral_detection_VIBRANT_%A_%a.out
#SBATCH --mem=16gb
#SBATCH --time=01:59:00
#SBATCH --cpus-per-task=16
#SBATCH --export=NONE
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

sample_dir=$1 #directory with the metaSPAdes contigs/scaffolds
sample_list=${sample_dir}/$2 #file with the list of all samples in the directory
SAMPLE_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${sample_list})

echo -e '\n-------------------- WORKING WITH '${SAMPLE_ID}' SAMPLE --------------------\n'

mkdir -p ${sample_dir}/${SAMPLE_ID}/virome_discovery ${sample_dir}/${SAMPLE_ID}/virome_discovery/VIBRANT

echo -e '---- PREDICTING ORFs with parallel PRODIGAL-GV ----\n'

# Clean environment, load modules 
module purge; ml prodigal-gv/2.11.0-GCCcore-12.2.0 Python/3.11.3-GCCcore-12.3.0; module list

mkdir -p ${sample_dir}/${SAMPLE_ID}/virome_discovery/VIBRANT/Prodigal-gv
# https://raw.githubusercontent.com/apcamargo/prodigal-gv/master/parallel-prodigal-gv.py
python parallel-prodigal-gv.py \
	-t ${SLURM_CPUS_PER_TASK} \
	-q \
	-i ${sample_dir}/${SAMPLE_ID}_metaspades_scaffolds.fa \
	-a ${sample_dir}/${SAMPLE_ID}/virome_discovery/VIBRANT/Prodigal-gv/${SAMPLE_ID}.faa \
	-o ${sample_dir}/${SAMPLE_ID}/virome_discovery/VIBRANT/Prodigal-gv/${SAMPLE_ID}_prodigal.out

echo -e '---- RUNNING VIBRANT ----\n'

# Clean environment, and activate conda environment
module purge; ml Anaconda3; module list
source activate /scratch/hb-llnext/conda_envs/Vibrant_env; conda list

# Run VIBRANT
/scratch/hb-llnext/conda_envs/Vibrant_env/bin/VIBRANT_run.py \
	-i ${sample_dir}/${SAMPLE_ID}/virome_discovery/VIBRANT/Prodigal-gv/${SAMPLE_ID}.faa \
	-folder ${sample_dir}/${SAMPLE_ID}/virome_discovery/VIBRANT/ \
	-f prot \
	-t ${SLURM_CPUS_PER_TASK} \
	-l 10000 \
	-no_plot
	
# Select only viral contigs with >= 10,000 bp (-l frag does not seem to work)
awk -F"_| " '$4 >= 10000' ${sample_dir}/${SAMPLE_ID}/virome_discovery/VIBRANT/VIBRANT_${SAMPLE_ID}/VIBRANT_phages_${SAMPLE_ID}/${SAMPLE_ID}.phages_combined.txt > \
${sample_dir}/${SAMPLE_ID}/virome_discovery/VIBRANT/VIBRANT_${SAMPLE_ID}/VIBRANT_phages_${SAMPLE_ID}/${SAMPLE_ID}.phages_combined_10000.txt

# Remove intermediate files
rm -r ${sample_dir}/${SAMPLE_ID}/virome_discovery/VIBRANT/VIBRANT_${SAMPLE_ID}/VIBRANT_HMM_tables_unformatted_${SAMPLE_ID}

conda deactivate
