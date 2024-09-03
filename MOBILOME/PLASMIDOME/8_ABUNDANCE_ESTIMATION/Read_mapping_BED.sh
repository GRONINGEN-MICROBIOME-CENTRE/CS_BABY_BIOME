#!/bin/bash
#SBATCH --job-name=Read_mapping_BED
#SBATCH --output=Read_mapping_BED_%A_%a.out
#SBATCH --mem=30gb
#SBATCH --time=00:59:00
#SBATCH --cpus-per-task=16
#SBATCH --export=NONE
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

BATCH=$1 #path to the batch folder with the quality trimmed reads
sample_list=${BATCH}/$2 #file with the list of samples
SAMPLE_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${sample_list})
Bowtie_DB=$3 # Path to Bowtie index generated from the plasmid sequences
plasmids=$4 # FASTA file with all the plasmid sequences

# Create directories if they don't exist
if [[ ! -d ${BATCH}/Alignment_results ]]; then
    mkdir ${BATCH}/Alignment_results
fi

if [[ ! -d ${BATCH}/Breadth_coverage_results ]]; then
    mkdir ${BATCH}/Breadth_coverage_results
fi

echo '-------------------- WORKING WITH '${SAMPLE_ID}' SAMPLE --------------------'

# Clean environment, load modules
module purge; ml Python/3.10.8-GCCcore-12.2.0 Bowtie2 SAMtools BEDTools; module list

# Map the reads
bowtie2 \
	--very-sensitive \
	-x $Bowtie_DB/plasmids_DB \
	-1 ${BATCH}/${SAMPLE_ID}_kneaddata_paired_1.fastq.gz \
	-2 ${BATCH}/${SAMPLE_ID}_kneaddata_paired_2.fastq.gz \
  -U ${BATCH}/${SAMPLE_ID}_kneaddata_unmatched.fastq.gz \
	--no-unal \
	--threads ${SLURM_CPUS_PER_TASK} \
	-S ${BATCH}/${SAMPLE_ID}_all_plasmid_alignments.sam

samtools view \
	-S ${BATCH}/${SAMPLE_ID}_all_plasmid_alignments.sam \
	-b \
	-o ${BATCH}/${SAMPLE_ID}_all_plasmid_alignments.bam

samtools sort \
	${BATCH}/${SAMPLE_ID}_all_plasmid_alignments.bam \
	-o ${BATCH}/Alignment_results/${SAMPLE_ID}.sorted.bam \
	-@ $((${SLURM_CPUS_PER_TASK}-1))

samtools index \
	-@ $((${SLURM_CPUS_PER_TASK}-1)) \
	${BATCH}/Alignment_results/${SAMPLE_ID}.sorted.bam

# Get coverage final tables 
bedtools coverage \
	-a $Bowtie_DB/plasmids_DB.bed \
	-b ${BATCH}/Alignment_results/${SAMPLE_ID}.sorted.bam \
	> ${BATCH}/Breadth_coverage_results/${SAMPLE_ID}.coverage.txt

# Remove intermediate files
rm ${BATCH}/${SAMPLE_ID}_all_plasmid_alignments.sam ${BATCH}/${SAMPLE_ID}_all_plasmid_alignments.bam
