#!/bin/bash
#SBATCH --job-name=metaplasmidSPAdes
#SBATCH --output=metaplasmidSPAdes_%A_%a.out
#SBATCH --mem=60gb
#SBATCH --time=06:59:00
#SBATCH --cpus-per-task=24
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --partition=regular

sample_dir=$1 #directory with the FASTQ files with clean reads
unmatched_sample_dir=$2 #directory with the FASTQ files of unmatched reads
sample_list=${sample_dir}/$3 #file with the list of all samples in the directory
SAMPLE_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${sample_list})

echo -e '\n-------------------- WORKING WITH '${SAMPLE_ID}' SAMPLE --------------------'

echo -e '\n---- RUNNING metaplasmidSPAdes ON '${SAMPLE_ID}' SAMPLE ----'

mkdir -p MSP_RESULTS/
mkdir -p metaPlasmidSPAdes/CONTIGS metaPlasmidSPAdes/LOG_files metaPlasmidSPAdes/OUTPUT_files

# Clean environment, load modules 
module purge; module load SPAdes; module list

metaplasmidspades.py \
    -1 ${sample_dir}/${SAMPLE_ID}_kneaddata_paired_1.fastq.gz \
    -2 ${sample_dir}/${SAMPLE_ID}_kneaddata_paired_2.fastq.gz \
    -s ${unmatched_sample_dir}/${SAMPLE_ID}_kneaddata_unmatched.fastq.gz \
    -o MSP_RESULTS/${SAMPLE_ID} \
    --threads ${SLURM_CPUS_PER_TASK}

echo -e '\n---- Generating folders with CONTIGS, PATHS and LOG FILES----'

# Add sample name to output files
mv MSP_RESULTS/${SAMPLE_ID}/contigs.fasta MSP_RESULTS/${SAMPLE_ID}/${SAMPLE_ID}_metaplasmidspades_contigs.fa 
mv MSP_RESULTS/${SAMPLE_ID}/spades.log MSP_RESULTS/${SAMPLE_ID}/${SAMPLE_ID}_spades.log

# Transfer them to new folder
rsync -av $(find MSP_RESULTS/${SAMPLE_ID} -name "${SAMPLE_ID}_metaplasmidspades_contigs.fa" -type f) metaPlasmidSPAdes/CONTIGS
rsync -av $(find MSP_RESULTS/${SAMPLE_ID} -name "${SAMPLE_ID}_spades.log" -type f) metaPlasmidSPAdes/LOG_files

echo -e '\n---- Checking errors ----'

mkdir -p metaPlasmidSPAdes/SUMMARY_RESULTS

if ! grep -Eq 'Error|FAILED' metaPlasmidSPAdes/LOG_files/${SAMPLE_ID}_spades.log; then
    echo -e "metaPlasmidSPAdes assembly for sample ${SAMPLE_ID} completed successfully.\n" > metaPlasmidSPAdes/SUMMARY_RESULTS/${SAMPLE_ID}_summary.txt
fi
echo -e "\nErrors checked for sample ${SAMPLE_ID}."

echo -e '\n---- Removing intermediate results ----'
rm -r MSP_RESULTS/${SAMPLE_ID}

echo -e "\nIntermediate results for sample ${SAMPLE_ID} have been removed."

echo -e '\n---- metaplasmidSPAdes step DONE ----'
