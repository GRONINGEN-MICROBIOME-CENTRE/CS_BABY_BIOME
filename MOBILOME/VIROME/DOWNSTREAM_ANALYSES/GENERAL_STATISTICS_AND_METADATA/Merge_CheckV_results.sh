#!/bin/bash
#SBATCH --job-name=merge_CheckV_results
#SBATCH --output=merge_CheckV_results.out
#SBATCH --mem=4gb
#SBATCH --time=00:09:59
#SBATCH --cpus-per-task=4
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --partition=regular

CheckV_output=$1 # CheckV output folder
mkdir -p CheckV_merged_results # directory to store results

echo -e '\n-------------------- MERGING CHECKV OUTPUT --------------------'

# Concatenate quality_summary.tsv, contamination.tsv, completeness.tsv and complete_genomes.tsv
# Important to use 'NR>1' to avoid concatenating headers
find $CheckV_output -name "quality_summary.tsv" -type f -exec head -1 {} \; -quit > CheckV_merged_results/quality_summary.tsv
find $CheckV_output -type f -name "quality_summary.tsv" -exec sed -e 1d '{}' \; >> CheckV_merged_results/quality_summary.tsv
find $CheckV_output -name "contamination.tsv" -type f -exec head -1 {} \; -quit > CheckV_merged_results/contamination_final.tsv
find $CheckV_output -type f -name "contamination.tsv" -exec sed -e 1d '{}' \; >> CheckV_merged_results/contamination_final.tsv
find $CheckV_output -name "completeness.tsv" -type f -exec head -1 {} \; -quit > CheckV_merged_results/completeness.tsv
find $CheckV_output -type f -name "completeness.tsv" -exec sed -e 1d '{}' \; >> CheckV_merged_results/completeness.tsv
find $CheckV_output -name "complete_genomes.tsv" -type f -exec head -1 {} \; -quit > CheckV_merged_results/complete_genomes.tsv
find $CheckV_output -type f -name "complete_genomes.tsv" -exec sed -e 1d '{}' \; >> CheckV_merged_results/complete_genomes.tsv
