#!/bin/bash
#SBATCH --job-name=check_Read_mapping
#SBATCH --output=check_Read_mapping.out
#SBATCH --mem=4gb
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=4
#SBATCH --export=NONE
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

directory=$1 #directory with output files of read mapping step

# Search for output files with "error" word and extract sample name
grep -r -E 'error' "$directory" |
while IFS= read -r file; do
	file_name=$(echo "$file" | sed 's/:.*//')
  sample_name=$(cat $file_name | grep -a '^-------------------- WORKING WITH' | sed -E 's/^\-+ WORKING WITH //' | sed -E 's/ SAMPLE \-+$//')
  echo "$sample_name contains 'error' word in the output file"
  echo "$sample_name" >> list_files_with_error.txt
done

# Search for output files without "COMPLETED" word and extract sample name
grep -rL -E 'COMPLETED' "$directory" |
while IFS= read -r file; do
	sample_name=$(cat $file | grep -a '^-------------------- WORKING WITH' | sed -E 's/^\-+ WORKING WITH //' | sed -E 's/ SAMPLE \-+$//')
  echo "$sample_name does not contain 'COMPLETED' word in the output file"
  echo "$sample_name" >> list_files_rerun.txt
done

