#!/bin/bash
#SBATCH --job-name=Plasmid_mobility_CONJScan
#SBATCH --output=Plasmid_mobility_CONJScan_%A_%a.out
#SBATCH --mem=20gb
#SBATCH --time=00:29:00
#SBATCH --cpus-per-task=16
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

plasmid_dir=$1 #directory with the individual FASTA files for each plasmid
plasmid_list=${plasmid_dir}/$2 #file with the list of samples
PLASMID_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${plasmid_list})


echo -e '-------------------- WORKING WITH '${PLASMID_ID}' PLASMID --------------------\n'

echo -e '-------------------- PREDICTING PLASMIDS PROTEINS WITH PRODIGAL  --------------------\n'

mkdir -p ORF_PREDICTION

# Clean environment, load modules 
module purge; ml prodigal HMMER; module list

# Run prodigal to get the protein predictions
prodigal \
    -f gff \
    -o ORF_PREDICTION/${PLASMID_ID}_proteins.gff \
    -a ORF_PREDICTION/${PLASMID_ID}_proteins.faa \
    -i ${plasmid_dir}/${PLASMID_ID}.fa  \
    -p meta

echo -e '\n-------------------- RUNNING CONJScan  --------------------'

mkdir -p CONJSCAN_RESULTS OUTPUT_files

echo -e '\n---- Copying PRODIGAL results to tmpdir ----'

mkdir -p ${TMPDIR}/CONJSCAN_RESULTS ${TMPDIR}/CONJSCAN_RESULTS/${PLASMID_ID}
rsync -av ORF_PREDICTION/${PLASMID_ID}_proteins.faa ${TMPDIR}/ORF_PREDICTION/

# Run CONJScan
macsyfinder \
    --db-type ordered_replicon \
    -o ${TMPDIR}/CONJSCAN_RESULTS/${PLASMID_ID} \
	  --sequence-db ORF_PREDICTION/${PLASMID_ID}_proteins.faa \
	  --models CONJScan/Plasmids all 	
	 
echo -e '\n---- Copying CONJScan results from tmpdir ----'

rm -r ${TMPDIR}/CONJSCAN_RESULTS/${PLASMID_ID}/hmmer_results
rm ${TMPDIR}/CONJSCAN_RESULTS/${PLASMID_ID}/best_solution_loners.tsv ${TMPDIR}/CONJSCAN_RESULTS/${PLASMID_ID}/best_solution_multisystems.tsv
rm ${TMPDIR}/CONJSCAN_RESULTS/${PLASMID_ID}/*conf ${TMPDIR}/CONJSCAN_RESULTS/${PLASMID_ID}/*log
rm ${TMPDIR}/CONJSCAN_RESULTS/${PLASMID_ID}/*rejected*
rsync -av ${TMPDIR}/CONJSCAN_RESULTS/${PLASMID_ID} CONJSCAN_RESULTS

