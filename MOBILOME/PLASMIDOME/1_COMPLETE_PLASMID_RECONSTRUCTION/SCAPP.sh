#!/bin/bash
#SBATCH --job-name=SCAPP
#SBATCH --output=SCAPP_%A_%a.out
#SBATCH --mem=80gb
#SBATCH --time=2-0
#SBATCH --cpus-per-task=24
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --partition=himem

sample_dir=$1 #directory with the FASTQ files with clean reads
SPAdes_graph_dir=$2 #directory with the assembly graphs 
sample_list=${sample_dir}/$3 #file with the list of all samples in the directory
SAMPLE_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${sample_list})

echo -e '\n-------------------- WORKING WITH '${SAMPLE_ID}' SAMPLE --------------------'

echo -e '\n---- RUNNING SCAPP ON '${SAMPLE_ID}' SAMPLE ----'

mkdir -p SCAPP_RESULTS/
mkdir -p SCAPP/CONTIGS  SCAPP/LOG_files SCAPP/OUTPUT_files

# Clean environment, load conda environment 
module purge; ml Anaconda3; module list
source activate /scratch/hb-llnext/conda_envs/SCAPP_env; conda list

scapp \
    -g ${SPAdes_graph_dir}/${SAMPLE_ID}_graph.fastg \
    -o SCAPP_RESULTS/${SAMPLE_ID} \
    -r1 ${sample_dir}/${SAMPLE_ID}_kneaddata_paired_1.fastq.gz \
    -r2 ${sample_dir}/${SAMPLE_ID}_kneaddata_paired_2.fastq.gz \
    -p ${SLURM_CPUS_PER_TASK}

echo -e '\n---- Generating folders with CONTIGS and LOG FILES----'

# Add sample name to output files
mv SCAPP_RESULTS/${SAMPLE_ID}/*confident_cycs.fasta SCAPP_RESULTS/${SAMPLE_ID}/${SAMPLE_ID}_SCAPP_contigs.fa 
mv SCAPP_RESULTS/${SAMPLE_ID}/logs/scapp.log SCAPP_RESULTS/${SAMPLE_ID}/${SAMPLE_ID}_scapp.log

# Transfer them to new folder
rsync -av $(find SCAPP_RESULTS/${SAMPLE_ID} -name "${SAMPLE_ID}_SCAPP_contigs.fa" -type f) SCAPP/CONTIGS
rsync -av $(find SCAPP_RESULTS/${SAMPLE_ID} -name "${SAMPLE_ID}_scapp.log" -type f) SCAPP/LOG_files

echo -e '\n---- Checking errors ----'

mkdir -p SCAPP/SUMMARY_RESULTS

if ! grep -Eq 'Error|FAILED' SCAPP/LOG_files/${SAMPLE_ID}_scapp.log; then
    echo -e "SCAPP for sample ${SAMPLE_ID} completed successfully.\n" > SCAPP/SUMMARY_RESULTS/${SAMPLE_ID}_summary.txt
fi
echo -e "\nErrors checked for sample ${SAMPLE_ID}."

echo -e '\n---- Removing intermediate results ----'
rm -r SCAPP_RESULTS/${SAMPLE_ID}

echo -e "\nIntermediate results for sample ${SAMPLE_ID} have been removed."

echo -e '\n---- SCAPP step DONE ----'
