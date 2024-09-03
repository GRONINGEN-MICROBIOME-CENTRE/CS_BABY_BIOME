#!/bin/bash
#SBATCH --job-name=geNomad_plasmids_complete
#SBATCH --output=geNomad_plasmids_complete.out
#SBATCH --mem=20gb
#SBATCH --time=00:29:00
#SBATCH --cpus-per-task=24
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

contig_file=$1 #path to FASTA file with contigs
contig_file_name="$(basename "${contig_file}")" #extract filename
contig_file_name="${contig_file_name%.*}" #extract filename without the extension
output_dir=$2 #path to output directory

# Clean environment, load modules and activate conda environment
module purge; ml Anaconda3; module list
source activate /home2/p304845/Conda_envs/geNomad_conda/; conda list 

# Run geNomad
genomad end-to-end \
        --enable-score-calibration \
	--disable-find-proviruses \
        $contig_file \
        $output_dir \
        --cleanup \
	--threads ${SLURM_CPUS_PER_TASK} \
        /scratch/hb-llnext/databases/geNomad_db/

conda deactivate

mv ${output_dir}/${contig_file_name}_summary/${contig_file_name}_plasmid_summary.tsv geNomad_plasmid_summary.tsv
rsync -av ${output_dir}/${contig_file_name}_summary/${contig_file_name}_plasmid.fna .
mv ${contig_file_name}_plasmid.fna geNomad_plasmid.fna

# Clean environment, load modules 
module purge; ml SeqKit; module list

# Process geNomad output files to select:

# A) Plasmids with >= 1000bp and FDR < 0.05 (plasmid origin)
awk -F'\t' 'NR == 1 || ($2 >= 1000 && $7 < 0.05) {print}' geNomad_plasmid_summary.tsv > geNomad_selected_plasmids_summary.tsv
awk 'NR > 1 {print $1}' geNomad_selected_plasmids_summary.tsv | cut -f1 > geNomad_plasmid_list.txt
seqkit grep -f geNomad_plasmid_list.txt geNomad_plasmid.fna > geNomad_plasmid_sequences.fna

# Remove intermediate files
rm -r ${output_dir}/${contig_file_name}_*

exec > geNomad_summary_info.txt
echo -e "\n" 

echo -e "################################### SUMMARY STATS - geNomad ###################################\n"
n_plasmids_geNomad=$(awk 'NR > 1 {print $1}' geNomad_plasmid_summary.tsv | wc -l)
n_plasmids_FDR_1kb=$(cat geNomad_plasmid_list.txt | wc -l)
n_plasmids_FDR_1kb_with_DTR=$(cat geNomad_plasmid_summary.tsv | grep DTR | wc -l)
echo "The number of contigs classified as plasmids by geNomad is: $n_plasmids_geNomad"
echo "The number of contigs classified as plasmids by geNomad with FDR < 0.05 and length > 1kb is: $n_plasmids_FDR_1kb"
echo -e "The number of contigs classified as plasmids by geNomad with FDR < 0.05 and length > 1kb that have DTRs is: $n_plasmids_FDR_1kb_with_DTR\n"

echo -e "The summary stats of the plasmids identified by geNomad (FDR < 0.05 and length > 1kb):\n"
seqkit stats -a geNomad_plasmid_sequences.fna
echo -e "\n"

echo -e "The file geNomad_plasmid_sequences.fna contains the selected plasmids predicted by geNomad.\n"
echo -e "###################################################### END ######################################################\n"
