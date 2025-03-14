#!/bin/bash
#SBATCH --job-name=summary_Read_mapping
#SBATCH --output=summary_Read_mapping.out
#SBATCH --mem=8gb
#SBATCH --time=00:29:00
#SBATCH --cpus-per-task=16
#SBATCH --export=NONE
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

plasmid_genomes=$1 # file with the names of plasmid genomes
samples=$2 # file with sample names
bed_coverage_output=$3 # path to BED coverage output files *.coverage.txt
kneadData_log=$4 # directory with KneadData LOG files

# Generate a file with the read alignment rate
for file in OUTPUT_files/*.out; do
    sample_name=$(grep -oP 'WORKING WITH \K[^ ]+' "$file")
    alignment_rate=$(grep -oP '\d+\.\d+% overall alignment rate' "$file" | grep -oP '\d+\.\d+')
    echo "$sample_name $alignment_rate" >> Bowtie_alignment_rate.txt
done

# Generate a table with sample names and number of clean reads
echo -e "sample\tclean_reads" > CS_Baby_Biome_nreads_195.txt
# Loop through all stats files 
for file in ${kneadData_log}/*.log; do
  # Extract Sample name and number of clean reads 
  sample_name=$(echo "$file" | cut -f8 -d "/" | cut -f1 -d "_")
  clean_reads_1=$(cat $file | grep "final pair1" | cut -f7 -d ":" | sed 's/^ *//' | cut -f1 -d ".")
  clean_reads_2=$(cat $file | grep "final pair2" | cut -f7 -d ":" | sed 's/^ *//' | cut -f1 -d ".")
  unmatched_reads_1=$(cat $file | grep "final orphan1" | cut -f7 -d ":" | sed 's/^ *//' | cut -f1 -d ".")
  unmatched_reads_2=$(cat $file | grep "final orphan2" | cut -f7 -d ":" | sed 's/^ *//' | cut -f1 -d ".")
  final_reads=$(( $clean_reads_1 + $clean_reads_2 + $unmatched_reads_1 + $unmatched_reads_2 ))
  # Append the data to the result file
  echo -e "$sample_name\t$final_reads" >> CS_Baby_Biome_nreads_195.txt
  echo -e "$sample_name\t$clean_reads_2" >> CS_Baby_Biome_nreads_195_FQ2.txt
done

# Clean environment, load modules
module purge; module load R; module list

# Execute the R script
Rscript summary_read_mapping.R $plasmid_genomes $samples $bed_coverage_output CS_Baby_Biome_nreads_195.txt
