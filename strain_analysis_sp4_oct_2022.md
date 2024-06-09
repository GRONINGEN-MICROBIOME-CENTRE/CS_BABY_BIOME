# Strainphlan 4.0 analysis 

Adapted from Biobakery (StrainPhlAn 4.0). 
https://github.com/biobakery/MetaPhlAn/wiki/StrainPhlAn-4

Authors: Trishla Sinha.
Description: The script shows how strain profiling was performed using Strainphlan 4 
Languages: Bash and R.   

## Step 1: Reconstruct all species strains

https://github.com/GRONINGEN-MICROBIOME-CENTRE/gmc-mgs-pipeline/blob/main/GMH_pipe.py 

## Step 2: Profile the clades present in the samples (profileClades.sh)

```
#!/bin/bash

#SBATCH --mem=14gb
#SBATCH --time=6-23:30:00
#SBATCH --cpus-per-task=4
#SBATCH --open-mode=truncate
#SBATCH --job-name=SP4pr
#SBATCH --error=__SP4_profile.err
#SBATCH --output=__SP4_profile.out

# NOTES:
# script profiles all clades in the dataset in given folder ($1)
# puts results in the current folder!
# Adding --mutation_rates will give a mutation rates table for each of the alignes markers and a summary table for the concatenated MSA
# Removing the --print_clades only will actually run it 

# PARAMS
N=10 # --marker_in_n_samples
S=10 # --sample_with_n_markers
DB=/scratch/hb-tifn/condas/conda_biobakery4/lib/python3.9/site-packages/metaphlan/metaphlan_databases/mpa_vOct22_CHOCOPhlAnSGB_202212.pkl

# purge modules
module purge
# load conda
ml Anaconda3/2022.05
# load conda env
source activate /scratch/hb-tifn/condas/conda_biobakery4/
# run clade profiling
strainphlan -s *.pkl --database /scratch/hb-tifn/condas/conda_biobakery4/lib/python3.9/site-packages/metaphlan/metaphlan_databases/mpa_vOct22_CHOCOPhlAnSGB_202212.pkl --marker_in_n_samples ${N} --sample_with_n_markers ${S} --print_clades_only --phylophlan_mode accurate --output_dir . > strainphlan4_clades_${N}.txt

```
### Execution 

```
sbatch ./profileClades.sh 

```
Thu Jan 18 22:10:03 2024: Start StrainPhlAn 4.0.6 execution
Thu Jan 18 22:10:03 2024: Loading MetaPhlAn mpa_vOct22_CHOCOPhlAnSGB_202212 database...
Thu Jan 18 22:10:24 2024: Done.
Thu Jan 18 22:10:27 2024: Detecting clades...
Thu Jan 18 23:42:37 2024: Done.
Thu Jan 18 23:42:38 2024: Detected clades: 
Thu Jan 18 23:42:38 2024:       t__SGB10068: in 168 samples.
Thu Jan 18 23:42:38 2024:       t__SGB17248: in 163 samples.
Thu Jan 18 23:42:38 2024:       t__SGB6939: in 112 samples.
Thu Jan 18 23:42:38 2024:       t__SGB17247: in 101 samples.
Thu Jan 18 23:42:38 2024:       t__SGB6952: in 95 samples.
Thu Jan 18 23:42:38 2024:       t__SGB6936: in 94 samples.
Thu Jan 18 23:42:38 2024:       t__SGB17256: in 87 samples.
Thu Jan 18 23:42:38 2024:       t__SGB7962: in 75 samples.
Thu Jan 18 23:42:38 2024:       t__SGB1814: in 74 samples.
Thu Jan 18 23:42:38 2024:       t__SGB17244_group: in 72 samples.
Thu Jan 18 23:42:38 2024:       t__SGB1836: in 71 samples.
Thu Jan 18 23:42:38 2024:       t__SGB1934: in 69 samples.
Thu Jan 18 23:42:38 2024:       t__SGB17234: in 69 samples.
Thu Jan 18 23:42:38 2024:       t__SGB14535_group: in 67 samples.
Thu Jan 18 23:42:38 2024:       t__SGB8007_group: in 65 samples.
Thu Jan 18 23:42:38 2024:       t__SGB10119: in 65 samples.
Thu Jan 18 23:42:38 2024:       t__SGB9663_group: in 61 samples.
Thu Jan 18 23:42:38 2024:       t__SGB4933: in 57 samples.



We next process this output file to select only the clade names

```
cat  strainphlan4_clades_10.txt | grep t__ | cut -f 2 | cut -f 1 -d ':' > CS_Baby_Biome_clade_names_Oct_2022_db.txt

```

This will give us the names of each species found: 
t__SGB10068
t__SGB17248
t__SGB6939
t__SGB17247
t__SGB6952
t__SGB6936
t__SGB17256
t__SGB7962
t__SGB1814
t__SGB17244_group
t__SGB1836
t__SGB1934
t__SGB17234
t__SGB14535_group
t__SGB8007_group
t__SGB10119
t__SGB9663_group
t__SGB4933
t__SGB10120
t__SGB1877
t__SGB15316
t__SGB15132
t__SGB4837
t__SGB2318
t__SGB2303
t__SGB15254
t__SGB4874



## Step 3: Build the multiple sequence alignment (sp4_runMarkerComparison.sh)


```
#!/bin/bash
#SBATCH --mem=24gb
#SBATCH --time=0-09:59:00
#SBATCH --cpus-per-task=8
#SBATCH --open-mode=truncate

# NOTES:
# > $1 is input folder
# > $2 is clade name

echo "Invoking runMarkerComparison.sh"
echo "CLs: ${1} ${2}"

# HELP / USE
if [[ $# -ne 2 ]]; 
    then 
    echo "ERROR: script requires two command line parameters:"
    echo " <input folder with pkl files> <clade name>"
    exit 2
fi

# PARAMS
# ===================
N=4 # --marker_in_n_samples
S=10 # --sample_with_n_markers 
MODE=accurate # {accurate,fast}
CONDA=/scratch/hb-tifn/condas/conda_biobakery4/
#CM= # clade markers

# purge modules
module purge

# load conda
ml Anaconda3/2022.05   
# load conda env
source deactivate
source activate ${CONDA}

# prep results folder (where clade result goes)
mkdir ${2}
# run strainphlan for that clade
echo "strainphlan -s ${1}/*.pkl --output_dir ./${2} --clade ${2} --marker_in_n_samples ${N} --sample_with_n_markers ${S} --nprocs 8 --phylophlan_mode ${MODE}"
strainphlan -s ${1}/*.pkl --database /scratch/hb-tifn/condas/conda_biobakery4/lib/python3.9/site-packages/metaphlan/metaphlan_databases/mpa_vOct22_CHOCOPhlAnSGB_202212.pkl --output_dir ./${2} --clade ${2} --marker_in_n_samples ${N} --sample_with_n_markers ${S} --$


```

### Execution

```
for i in $(cat CS_Baby_Biome_clade_names_Oct_2022_db.txt); do sbatch sp4_runMarkerComparison.sh  /scratch/p280306/CS_BABY_BIOME_AMS_NEXT/biobakery_oct_2022_results_CS_Baby_Biome_2024/strainphlan4 $i; done
```
Where CS_BABY_BIOME_clades_names.txt contains a list of names of all (sub) species identified in the previous step 
This will perform MSA and create .tre files and .aln files for each of the (sub)species you feed it in 
