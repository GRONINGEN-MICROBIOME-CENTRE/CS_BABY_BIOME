#!/bin/bash
#SBATCH --job-name=CheckV_derep_CS_viruses
#SBATCH --output=CheckV_derep_CS_viruses.out
#SBATCH --mem=20gb
#SBATCH --time=06:59:00
#SBATCH --cpus-per-task=24
#SBATCH --export=NONE
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

contig_file=$1 #file with all the renamed viral genomes used as input for the dereplication
blast_db=$2 #path to blast+ database generated from all viral sequences

mkdir -p CS_DEREPLICATION

# Load SeqKit
module purge; ml SeqKit Python/3.10.8-GCCcore-12.2.0 CheckV seqtk; ml list

# Extract only CS genomes from the FASTA file with the renamed viral sequences
seqkit grep -r -p "CS_" $contig_file > CS_viral_genomes.fa

# Extract only comparisons between CS genomes from the blast results
cat $blast_db | grep -Ev "MGV|GPD|Gulyaeva|Benler|Shah|IMG_VR|ELGV|RefSeq|Guerin|Yutin" > BLASTN_RESULTS/CS_viruses_blast.tsv 

echo -e '\n-------------------- DEDUPLICATING CS VIRAL GENOMES --------------------'

# Run deduplication for contigs with CheckV scripts

# Using the blast+ database as input, calculate pairwise ANI by combining local alignments between sequence pairs:
python /scratch/p304845/CS_Baby_Biome/6_DEREPLICATION_vOTUs/anicalc.py -i BLASTN_RESULTS/CS_viruses_blast.tsv  -o CS_DEREPLICATION/CS_viruses_ani.tsv

# Perform UCLUST-like clustering (99% ANI + 95% AF):
python /scratch/p304845/CS_Baby_Biome/6_DEREPLICATION_vOTUs/aniclust.py --fna CS_viral_genomes.fa --ani CS_DEREPLICATION/CS_viruses_ani.tsv --out CS_DEREPLICATION/CS_viral_genome_clusters.tsv --min_ani 99 --min_tcov 95 --min_qcov 0

# Extract the representatives viral genomes (deduplicated genomes)
cat CS_DEREPLICATION/CS_viral_genome_clusters.tsv  | cut -f1 > CS_dedup_viral_genomes.txt
seqtk subseq CS_viral_genomes.fa CS_dedup_viral_genomes.txt > CS_dedup_viral_genomes.fa

echo -e '\n-------------------- DEREPLICATING CS VIRAL GENOMES INTO vOTUs --------------------'

# Perform UCLUST-like clustering using the MIUVIG recommended-parameters (95% ANI + 85% AF):
# Use deduplicated genomes as input
python /scratch/p304845/CS_Baby_Biome/6_DEREPLICATION_vOTUs/aniclust.py --fna CS_dedup_viral_genomes.fa --ani CS_DEREPLICATION/CS_viruses_ani.tsv --out CS_DEREPLICATION/CS_vOTU_clusters.tsv --min_ani 95 --min_tcov 85 --min_qcov 0

# Generate a TXT and a FASTA file with the representative sequences of each CS vOTU (only CS)
cat CS_DEREPLICATION/CS_vOTU_clusters.tsv | cut -f1 > rep_seqs_CS_vOTUs.txt
seqtk subseq CS_viral_genomes.fa rep_seqs_CS_vOTUs.txt > rep_seqs_CS_vOTUs.fa

# Set permissions 
chmod 440 rep_seqs_CS_vOTUs.txt rep_seqs_CS_vOTUs.fa
