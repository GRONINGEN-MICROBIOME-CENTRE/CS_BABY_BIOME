#!/bin/bash
#SBATCH --job-name=Viral_detection_CT3
#SBATCH --output=Viral_detection_CT3_%A_%a.out
#SBATCH --mem=50gb
#SBATCH --time=00:59:00
#SBATCH --cpus-per-task=24
#SBATCH --export=NONE
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

sample_dir=$1 #directory with the metaSPAdes contigs/scaffolds
sample_list=${sample_dir}/$2 #file with the list of all samples in the directory
SAMPLE_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${sample_list})

echo -e '\n-------------------- WORKING WITH '${SAMPLE_ID}' SAMPLE --------------------\n'

mkdir -p ${sample_dir}/${SAMPLE_ID}/virome_discovery ${sample_dir}/${SAMPLE_ID}/virome_discovery/CT3

# --- SWITCH TO TMP (large number of intermediate files) --- 
mkdir -p ${TMPDIR}/${SAMPLE_ID}; cd ${TMPDIR}/${SAMPLE_ID}

echo -e '---- RUNNING CENOTE-TAKER3 ----\n'

# Load modules  and activate environment
# pip install pyrodigal and pyhmmer (if not working)
module purge; ml Anaconda3; module list
source activate /scratch/hb-llnext/conda_envs/CT3_env; conda list

# Run Cenote-Taker3
cenotetaker3 \
	-c ${sample_dir}/${SAMPLE_ID}_metaspades_scaffolds.fa \
	-r CenoteTaker3 \
	-p F \
	-db virion \
	--minimum_length_circular 10000 \
	--minimum_length_linear 10000 \
	--caller prodigal-gv \
	--cenote-dbs /scratch/hb-llnext/databases/CenoteTaker3_DB \
	-t ${SLURM_CPUS_PER_TASK}

echo -e '---- MOVING RESULTS TO /SCRATCH -----\n'

rsync -av ./CenoteTaker3/CenoteTaker3_cenotetaker.log ${sample_dir}/${SAMPLE_ID}/virome_discovery/CT3
rsync -av ./CenoteTaker3/CenoteTaker3_virus_summary.tsv ${sample_dir}/${SAMPLE_ID}/virome_discovery/CT3
rsync -av ./CenoteTaker3/final_genes_to_contigs_annotation_summary.tsv ${sample_dir}/${SAMPLE_ID}/virome_discovery/CT3

# Remove intermeidate files
rm -r ./CenoteTaker3/ct_processing # intermediate
rm -r ./CenoteTaker3/sequin_and_genome_maps # we will re-create it afterwards for selected contigs

conda deactivate
