#!/bin/bash
#SBATCH --job-name=inStrain_compare
#SBATCH --output=inStrain_compare.out
#SBATCH --mem=200gb
#SBATCH --time=06:59:00
#SBATCH --cpus-per-task=24
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --partition=regular

IS_profiles=$1 #path to directory with the inStrain profile results

echo '-------------------- RUNNING inStrain COMPARE  --------------------'

# Activate conda environment
source /clusterfs/jgi/scratch/science/metagen/afernandezpato/Tools/miniconda3/bin/activate \
    /clusterfs/jgi/scratch/science/metagen/afernandezpato/Tools/inStrain; conda list

mkdir -p OUTPUT_files

#Run inStrain compare (all but breadth 0.75 are default parameters)
inStrain compare \
	  -i $IS_profiles/* \
	  -o 3_INSTRAIN_COMPARISONS \
	  --breadth 0.75 \
	  --min_cov 5 \
	  --min_freq 0.05 \
	  --fdr 1e-06 \
	  -p ${SLURM_CPUS_PER_TASK} 

