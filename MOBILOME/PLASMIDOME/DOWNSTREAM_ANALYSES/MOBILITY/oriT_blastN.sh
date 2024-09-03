#!/bin/bash
#SBATCH --job-name=oriT_blastN
#SBATCH --output=oriT_blastN.out
#SBATCH --mem=20gb
#SBATCH --time=00:09:00
#SBATCH --cpus-per-task=32
#SBATCH --export=NONE
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

contig_file=$1 #path to FASTA file with predicted plasmid contigs
oriT_sequences=$2 #path to FASTA with oriT sequences

mkdir -p BLASTN_RESULTS

# Load BLAST
module purge; ml BLAST+/2.13.0-gompi-2022a; ml list

#First, create a blast+ database:
makeblastdb \
	-in $contig_file \
	-dbtype nucl \
	-out BLASTN_RESULTS/oriTs

#Next, use blastn to perform the comparison of plasmid DB vs oriT:
blastn \
	-query $oriT_sequences \
	-db BLASTN_RESULTS/oriTs \
	-task blastn-short \
	-outfmt '6 std qlen slen qseq sseq' \
	-out BLASTN_RESULTS/oriT_blast.tsv \
	-dust no \
	-num_threads ${SLURM_CPUS_PER_TASK}


# Compute the query (oriT) coverage
awk 'NR>1 {print $0"\t"($4/$13)*100}' BLASTN_RESULTS/oriT_blast.tsv > BLASTN_RESULTS/oriT_blast_with_coverage.tsv

# Add header	
header="qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen\tqseq\tsseq\tqcov"
sed -i "1i$header" BLASTN_RESULTS/oriT_blast_with_coverage.tsv 


