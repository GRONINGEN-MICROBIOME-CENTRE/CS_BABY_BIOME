################################################################################
##### CS Baby Biome: vOTU AMR analysis - associations
### Author(s): Asier Fernández-Pato
### Last updated: 28th June, 2024
################################################################################

# Set working directory 
setwd("~/Desktop/PhD/Projects/CS Baby Biome/Final_analysis/VIRUSES/")

# Read abundance table and metadata tables
Virus_abundance_infants <- read.delim("Abundance_table/CS_Abundance_Table_INFANTS_17052024.txt")
Sample_metadata_infants <- read.delim("Metadata_CS/CS_Baby_Biome_Sample_Metadata_Infants_17052024.txt")
Virus_metadata <- read.delim("Metadata_CS/CS_Baby_Biome_Viral_Metadata_17052024.txt")


##************************************************************************
# 0. Phenotype selection and AMR data processing: Infants
##************************************************************************

# Read non-transformed AMR abundance tables 
AM_abundance <- read.delim("Abundance_table/CS_Baby_Biome_AM_Abundance_Table.txt", check.names = F)
AMR_Gene_abundance <- read.delim("Abundance_table/CS_Baby_Biome_AMR_Gene_Abundance_Table.txt", check.names = F)
AMR_Class_abundance <- read.delim("Abundance_table/CS_Baby_Biome_AMR_Class_Abundance_Table.txt", check.names = F)
AMR_Mechanism_abundance <- read.delim("Abundance_table/CS_Baby_Biome_AMR_Mechanism_Abundance_Table.txt", check.names = F)

# Read AMR count data
AM_presence <- read.delim("Abundance_table/CS_Baby_Biome_AM_Count_Table.txt", check.names = F)
AMR_Gene_presence <- read.delim("Abundance_table/CS_Baby_Biome_AMR_Gene_Count_Table.txt", check.names = F)
AMR_Class_presence <- read.delim("Abundance_table/CS_Baby_Biome_AMR_Class_Count_Table.txt", check.names = F)
AMR_Mechanism_presence <- read.delim("Abundance_table/CS_Baby_Biome_AMR_Mechanism_Count_Table.txt", check.names = F)

# Subset only samples in metadata
AM_abundance <- AM_abundance[rownames(AM_abundance) %in% Sample_metadata_infants$bioSampleId,]
AMR_Gene_abundance <- AMR_Gene_abundance[rownames(AMR_Gene_abundance) %in% Sample_metadata_infants$bioSampleId,]
AMR_Class_abundance <- AMR_Class_abundance[rownames(AMR_Class_abundance) %in% Sample_metadata_infants$bioSampleId,]
AMR_Mechanism_abundance <- AMR_Mechanism_abundance[rownames(AMR_Mechanism_abundance) %in% Sample_metadata_infants$bioSampleId,]

AM_presence <- AM_presence[rownames(AM_presence) %in% Sample_metadata_infants$bioSampleId,]
AMR_Gene_presence <- AMR_Gene_presence[rownames(AMR_Gene_presence) %in% Sample_metadata_infants$bioSampleId,]
AMR_Class_presence <- AMR_Class_presence[rownames(AMR_Class_presence) %in% Sample_metadata_infants$bioSampleId,]
AMR_Mechanism_presence <- AMR_Mechanism_presence[rownames(AMR_Mechanism_presence) %in% Sample_metadata_infants$bioSampleId,]


# Select only phenotypes of interest from metadata
Sample_metadata_infants <- Sample_metadata_infants[,c("DNA_concentration_ng_ul","read_depth", "bioSampleId",
                                                      "CS_BABY_BIOME_ID", "Timepoint_categorical", "Timepoint_numeric", 
                                                      "preg_gest_age","pre_preg_bmi_mother", "infant_birthweight","infant_sex",
                                                      "feeding_mode_pragmatic", "living_situation", "cats_dogs", "rand_AB","richness",
                                                      "shannon", "Virome.composition.infant", "Viral_relab_temperates", "AMR_gene_counts"), ]
# Convert character variables to factors
Sample_metadata_infants[sapply(Sample_metadata_infants, is.character)] <- lapply(Sample_metadata_infants[sapply(Sample_metadata_infants, is.character)],  #convert character columns to factors
                                                                                 as.factor)

##################################
# Process AMR tables: Log-transform abundances and filter by prevalence
##################################

# Generate log transformed abundance tables
## CLR transformation not suitable here: data hardly compositional after read mapping
## Log transformation can help with data normality necessary for linear models
pseudocount_AM_abundance <- min(AM_abundance[AM_abundance > 0]) /2
pseudocount_AMR_Gene_abundance <- min(AMR_Gene_abundance[AMR_Gene_abundance > 0]) /2
pseudocount_AMR_Class_abundance <- min(AMR_Class_abundance[AMR_Class_abundance > 0]) /2
pseudocount_AMR_Mechanism_abundance <- min(AMR_Mechanism_abundance[AMR_Mechanism_abundance > 0]) /2

AM_abundance_log <- AM_abundance
AM_abundance_log[AM_abundance_log  == 0] <- pseudocount_AM_abundance
AMR_Gene_abundance_log <- AMR_Gene_abundance
AMR_Gene_abundance_log[AMR_Gene_abundance_log  == 0] <- pseudocount_AMR_Gene_abundance
AMR_Class_abundance_log <- AMR_Class_abundance
AMR_Class_abundance_log[AMR_Class_abundance_log  == 0] <- pseudocount_AMR_Class_abundance
AMR_Mechanism_abundance_log <- AMR_Mechanism_abundance
AMR_Mechanism_abundance_log[AMR_Mechanism_abundance_log  == 0] <- pseudocount_AMR_Mechanism_abundance

AM_abundance_log <- log(AM_abundance_log)
AMR_Gene_abundance_log <- log(AMR_Gene_abundance_log)
AMR_Class_abundance_log <- log(AMR_Class_abundance_log)
AMR_Mechanism_abundance_log <- log(AMR_Mechanism_abundance_log)

# Select only features with a prevalence > 10%
# We observe no features with more than 10% prevalence across samples
# We conclude that AMRs in phages from early life are rare and we do not have the power to explore associations
AM_abundance_prev_list <- colnames(AM_abundance)[colSums(AM_abundance != 0) > 159*0.1] 
AMR_Gene_abundance_prev_list <- colnames(AMR_Gene_abundance)[colSums(AMR_Gene_abundance != 0) > 159*0.1] 
AMR_Gene_presence_prev_list <- colnames(AMR_Gene_presence)[colSums(AMR_Gene_presence != 0) > 159*0.1] # to include "AMR_Gene_counts" 
AMR_Class_abundance_prev_list <- colnames(AMR_Class_abundance)[colSums(AMR_Class_abundance != 0) > 159*0.1] 
AMR_Mechanism_abundance_prev_list <- colnames(AMR_Mechanism_abundance)[colSums(AMR_Mechanism_abundance != 0) > 159*0.1] 


##************************************************************************
# 1. General statistics
##************************************************************************

# Estimate number of vOTUs with AMRs in infants (W1-W6)
Virus_metadata_infants <- Virus_metadata[Virus_metadata$Virus_ID %in% rownames(Virus_abundance_infants),]
length(which(Virus_metadata_infants$AMR_genes>0))

# Estimate total number of AMRs in infant vOTUs (W1-W6)
sum(Virus_metadata_infants$AMR_genes)

# Estimate number of different AMR genes found in infant vOTUs
length(colSums(AMR_Gene_abundance)[colSums(AMR_Gene_abundance) != 0])

# Estimate AMR gene content (proportion) infant gut phages
sum(Virus_metadata_infants$AMR_genes) / sum(Virus_metadata_infants$n_genes)

# Compare with vOTU and PTU plasmid gene content in infants
# Phages: 49278 total genes (11 AMR, 49267 NO AMR)
# Plasmids: 41051 total genes (132 AMR, 40919 NO AMR)
Count_AMR_table <- matrix(c(132, 40919,
                       11, 49267),
                     nrow = 2, byrow = TRUE,
                     dimnames = list(Genome_Type = c("Plasmids", "Phages"),
                                     ARGs = c("Yes", "No")))

chisq_result <- chisq.test(Count_AMR_table)
print(chisq_result)

#****************
# Save output
#****************
# Save Virus metadata for infants
write.table(Virus_metadata_infants,"Metadata_CS/CS_Baby_Biome_Viral_Metadata_Infants_17052024.txt", sep = "\t", 
            row.names = T, col.names = T, quote = FALSE)



