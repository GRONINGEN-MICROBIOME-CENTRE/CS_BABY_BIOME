#!/bin/bash
#SBATCH --job-name=CONJscan_results_processing
#SBATCH --output=CONJscan_results_processing.out
#SBATCH --mem=8gb
#SBATCH --time=00:29:00
#SBATCH --cpus-per-task=24
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

CONJScan_results=$1 #path to CONJScan results folder 

# Create an empty TSV to store results
echo -e "Plasmid_ID\tCONJ_MOB_elements" > Plasmid_CONJ_MOB_elements.tsv

# Loop through each plasmid in the CONJScan results
for plasmid in $(ls $CONJScan_results); do
  # Extract the plasmid ID and the conjugative systems from the all_systems.txt
  plasmid_id=$(echo $plasmid)
  conj_systems=$(grep "model =" ${CONJScan_results}/${plasmid}/all_systems.txt | cut -f3 -d "/")
  # Convert spaces to commas and print plasmid ID and conjugative systems (tab-separated)
  conj_systems=$(echo $conj_systems | tr " " ",")
  # Append the plasmid ID and conjugative systems to the TSV file
  echo -e "$plasmid_id\t$conj_systems" >> Plasmid_CONJ_MOB_elements.tsv
done

