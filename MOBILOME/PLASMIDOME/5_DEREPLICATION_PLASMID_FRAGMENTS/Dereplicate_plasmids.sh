#!/bin/bash
#SBATCH --job-name=Dereplicate_plasmids
#SBATCH --output=Dereplicate_plasmids.out
#SBATCH --mem=20gb
#SBATCH --time=04:59:00
#SBATCH --cpus-per-task=24
#SBATCH --export=NONE
#SBATCH --open-mode=truncate
#SBATCH --partition=himem

plasmids=$1 #path to FASTA file with all the plasmids
plasmids_filename="$(basename "${plasmids}")" #extract filename
plasmids_filename="${plasmids_filename%.*}" #extract filename without the extension
blast_db=$2 #path to blast+ database generated from the sequences
min_ani=$3 #Minimum average nucleotide identity
min_tcov=$4 #Minimum alignment coverage of shorter sequence
min_qcov=$5 #Minimum alignment coverage of longer sequence

# Clean environment, load modules
module purge; ml Python/3.10.8-GCCcore-12.2.0 CheckV; module list

# Run dereplication for contigs with CheckV scripts

#Using the blast+ database as input, calculate pairwise ANI by combining local alignments between sequence pairs:
python /scratch/p304845/CS_Baby_Biome_prueba/PLASMIDS/4_DEREPLICATION_CIRCULAR_PLASMIDS/anicalc.py \
    -i $blast_db  \
    -o ${plasmids_filename}_ani.tsv

#Finally, perform UCLUST-like clustering:
python /scratch/p304845/CS_Baby_Biome_prueba/PLASMIDS/4_DEREPLICATION_CIRCULAR_PLASMIDS/aniclust.py \
    --fna $plasmids \
    --ani ${plasmids_filename}_ani.tsv \
    --out ${plasmids_filename}_clusters.tsv \
    --min_ani $min_ani \
    --min_tcov $min_tcov \
    --min_qcov $min_qcov
