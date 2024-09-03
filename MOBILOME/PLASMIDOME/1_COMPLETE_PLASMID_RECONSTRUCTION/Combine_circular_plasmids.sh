#!/bin/bash
#SBATCH --job-name=Combine_circular_plasmids
#SBATCH --output=Combine_circular_plasmids.out
#SBATCH --mem=2gb
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=24
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --partition=regular

metaPlasmidSPAdes_dir=$1 #directory with the metaPlasmidSPAdes contigs
SCAPP_dir=$2 #directory with the SCAPP contigs

ls $metaPlasmidSPAdes_dir | cut -f1 -d "_" | sort | uniq > list_samples.txt

echo -e '\n---- MERGING metaPlasmidSPAdes and SCAPP plasmids ----'

mkdir -p COMBINED_RESULTS/
mkdir -p COMBINED_RESULTS/CONTIGS COMBINED_RESULTS/RESULT_files COMBINED_RESULTS/OUTPUT_files

# For each sample, select metaPlasmidSPAdes and SCAPP contigs and create a combined file
while read i
do
	sample=$i
	grep ">" $metaPlasmidSPAdes_dir/${sample}_metaplasmidspades_contigs.fa | cut -f2 -d ">"  > COMBINED_RESULTS/RESULT_files/${sample}_MSP_circular_contigs.txt
	grep ">" $SCAPP_dir/${sample}_SCAPP_contigs.fa | cut -f2 -d ">"  > COMBINED_RESULTS/RESULT_files/${sample}_SCAPP_circular_contigs.txt
	cat COMBINED_RESULTS/RESULT_files/${sample}_MSP_circular_contigs.txt COMBINED_RESULTS/RESULT_files/${sample}_SCAPP_circular_contigs.txt | sort | uniq > COMBINED_RESULTS/RESULT_files/${sample}_combined_circular_contigs.txt
  
# Generate the FASTA file with the combined contigs. Add sample name to contig_IDs.
  cat $metaPlasmidSPAdes_dir/${sample}_metaplasmidspades_contigs.fa $SCAPP_dir/${sample}_SCAPP_contigs.fa > COMBINED_RESULTS/CONTIGS/${sample}_combined_circular_contigs.fa

done	< list_samples.txt

# Generate the TXT files and FASTA with the combined contigs. Add sample name to contig_IDs.
sed "s/^>/>${sample}_/" COMBINED_RESULTS/CONTIGS/*_combined_circular_contigs.fa > COMBINED_RESULTS/all_combined_circular_contigs.fa
sed "s/^/${sample}_/" COMBINED_RESULTS/RESULT_files/*_MSP_circular_contigs.txt > COMBINED_RESULTS/all_MSP_circular_contigs.txt
sed "s/^/${sample}_/" COMBINED_RESULTS/RESULT_files/*_SCAPP_circular_contigs.txt > COMBINED_RESULTS/all_SCAPP_circular_contigs.txt
sed "s/^/${sample}_/" COMBINED_RESULTS/RESULT_files/*_combined_circular_contigs.txt > COMBINED_RESULTS/all_combined_circular_contigs.txt

rm list_samples.txt

echo -e '\n---- MERGING step DONE ----'
