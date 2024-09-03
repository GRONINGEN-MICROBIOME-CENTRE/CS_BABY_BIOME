#!/bin/bash
#SBATCH --job-name=geNomad_plasmids
#SBATCH --output=geNomad_plasmids.out
#SBATCH --mem=40gb
#SBATCH --time=00:59:00
#SBATCH --cpus-per-task=24
#SBATCH --open-mode=truncate
#SBATCH --partition=himem

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
        --relaxed \
        --cleanup \
	--threads ${SLURM_CPUS_PER_TASK} \
        /scratch/hb-llnext/databases/geNomad_db/

conda deactivate

mv ${output_dir}/${contig_file_name}_summary/${contig_file_name}_plasmid_summary.tsv geNomad_plasmid_summary.tsv

