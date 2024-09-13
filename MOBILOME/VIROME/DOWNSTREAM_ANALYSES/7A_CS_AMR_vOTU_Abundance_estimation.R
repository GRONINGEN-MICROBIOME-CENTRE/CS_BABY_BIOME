################################################################################
##### CS Baby Biome: vOTU AMR abundance estimation
### Author(s): Asier Fernández-Pato
### Last updated: 28th June, 2024
################################################################################

# Set working directory 
setwd("~/Desktop/PhD/Projects/CS Baby Biome/Final_analysis/VIRUSES/")

#****************
# Load libraries
#****************
library(dplyr)
library(tidyverse)

#****************
# Define functions
#****************

# Function to estimate the abundance of AMR features (genes, classes, mechanisms). Input:
# vOTU_abundance: vOTU abundance table (vOTUs in rows and samples in columns)
# AM_count: AM feature count table (vOTUs in rows and AM feature in columns)
calculate_AMR_abundance <- function(vOTU_abundance, AM_count) {
  # Initialize an empty data frame for the gene abundance table
  gene_abundance_table <- data.frame(matrix(ncol = ncol(vOTU_abundance), nrow = ncol(AM_count)))
  colnames(gene_abundance_table) <- colnames(vOTU_abundance)
  rownames(gene_abundance_table) <- colnames(AM_count)
  
  # Generate a table with merged vOTU abundances and AMR counts
  merged_data <- merge(vOTU_abundance, AM_count, by = 0, all.x = TRUE)
  merged_data <- replace(merged_data, is.na(merged_data), 0)
  rownames(merged_data) <- merged_data[, 1]
  merged_data <- merged_data[, -1]
  
  # Loop through samples and AMRs to estimate the abundances
  for (sample in colnames(merged_data)[1:(ncol(vOTU_abundance))]) {
    vOTUs_present <- rownames(merged_data)[merged_data[, sample] != 0] # vOTUs present in the sample
    abundance_vOTUs <- merged_data[merged_data[, sample] != 0, sample] # the abundance of vOTUs present in the sample
    
    for (AMR in colnames(merged_data)[(ncol(vOTU_abundance) + 1):ncol(merged_data)]) {
      vOTUs_AMR <- rownames(merged_data)[merged_data[, AMR] != 0] # vOTUs that have the AMR
      counts_AMR_vOTUs <-  merged_data[merged_data[, AMR] != 0, AMR] # counts of AMRs the vOTUs
      
      # Select vOTUs present in the sample that have the AMR gene / class
      vOTUs <- intersect(vOTUs_AMR, vOTUs_present) # vOTUs with the AMR present in the sample
      counts <- as.numeric(counts_AMR_vOTUs[vOTUs_AMR %in% vOTUs_present]) # the AMR counts of the vOTUs present in the sample
      abundance <- as.numeric(abundance_vOTUs[match(vOTUs, vOTUs_present)]) # the abundance of the vOTUs with the AMR present in the sample
      gene_abundance <- ifelse(length(counts) > 0 && length(abundance) > 0, sum(counts * abundance), 0)
      
      # Add AMR abundance value to the corresponding sample column
      gene_abundance_table[AMR, sample] <- gene_abundance
    }
  }
  return(gene_abundance_table)
}

calculate_AMR_counts <- function(vOTU_abundance, AM_count) {
  # Initialize an empty data frame for the count table
  count_table <- data.frame(matrix(ncol = ncol(vOTU_abundance), nrow = ncol(AM_count)))
  colnames(count_table) <- colnames(vOTU_abundance)
  rownames(count_table) <- colnames(AM_count)
  
  # Generate a table with merged vOTU abundances and AMR counts
  merged_data <- merge(vOTU_abundance, AM_count, by = 0, all.x = TRUE)
  merged_data <- replace(merged_data, is.na(merged_data), 0)
  rownames(merged_data) <- merged_data[, 1]
  merged_data <- merged_data[, -1]
  
  # Loop through samples and AMRs to estimate the abundances
  for (sample in colnames(merged_data)[1:(ncol(vOTU_abundance))]) {
    vOTUs_present <- rownames(merged_data)[merged_data[, sample] != 0] # vOTUs present in the sample
    abundance_vOTUs <- merged_data[merged_data[, sample] != 0, sample] # the abundance of vOTUs present in the sample
    
    for (AMR in colnames(merged_data)[(ncol(vOTU_abundance) + 1):ncol(merged_data)]) {
      vOTUs_AMR <- rownames(merged_data)[merged_data[, AMR] != 0] # vOTUs that have the AMR
      counts_AMR_vOTUs <-  merged_data[merged_data[, AMR] != 0, AMR] # counts of AMRs the vOTUs
      
      # Select vOTUs present in the sample that have the AMR gene / class
      vOTUs <- intersect(vOTUs_AMR, vOTUs_present) # vOTUs with the AMR present in the sample
      counts <- as.numeric(counts_AMR_vOTUs[vOTUs_AMR %in% vOTUs_present]) # the AMR counts of the vOTUs present in the sample
      abundance <- as.numeric(abundance_vOTUs[match(vOTUs, vOTUs_present)]) # the abundance of the vOTUs with the AMR present in the sample
      feature_count <- ifelse(length(counts) > 0 && length(abundance) > 0, sum(counts), 0)
      
      # Add AMR feature count value to the corresponding sample column
      count_table[AMR, sample] <- feature_count
    }
  }
  return(count_table)
}

Sample_metadata <- read.delim("Metadata_CS/CS_Baby_Biome_Sample_Metadata_17052024.txt")
Sample_metadata_infants <- read.delim("Metadata_CS/CS_Baby_Biome_Sample_Metadata_Infants_17052024.txt")
Virus_metadata <- read.delim("Metadata_CS//CS_Baby_Biome_Viral_Metadata_17052024.txt")


#*******************************
# Processing of CARD annotation with RGI
#*******************************

# Read CARD AMR annotation table
CARD_annotation <- read.delim("7_AMR_ANALYSIS/RGI_CARD.txt", sep = "\t", header = T, check.names = F)
CARD_annotation$Virus_ID <- sub(" #.*", "", CARD_annotation$ORF_ID)
CARD_annotation$Virus_ID <- sub("_[^_]*$", "", CARD_annotation$Virus_ID)
CARD_annotation$Contig <-NULL

# Generate tables with the counts of (per Virus):
# A) AMR Gene
# B) AMR Classes
# C) AMR Mechanisms
# D) Antibiotics/Antimicrobials (group those closely related: e.g. Colistin A and Colistin B)
AMR_Gene_count <- CARD_annotation %>%
  count(Virus_ID, Best_Hit_ARO) %>%
  pivot_wider(names_from = Best_Hit_ARO, values_from = n, values_fill = 0) %>%
  column_to_rownames(var = "Virus_ID")

AMR_Class_count <- CARD_annotation %>%
  filter(`Drug Class` != "") %>%
  separate_rows(`Drug Class`, sep = "; ") %>%
  group_by(`Drug Class`, Virus_ID) %>%
  summarise(n = sum(!is.na(`Drug Class`))) %>%
  pivot_wider(names_from = `Drug Class`, values_from = n, values_fill = 0) %>%
  column_to_rownames(var = "Virus_ID")

AMR_Class_count$`streptogramin A antibiotic` <- NULL # a more general class already present
AMR_Class_count$`streptogramin B antibiotic` <- NULL

AMR_Mechanism_count <- CARD_annotation %>%
  filter(`Resistance Mechanism` != "") %>%
  separate_rows(`Resistance Mechanism`, sep = "; ") %>%
  group_by(`Resistance Mechanism`, Virus_ID) %>%
  summarise(n = sum(!is.na(`Resistance Mechanism`))) %>%
  pivot_wider(names_from = `Resistance Mechanism`, values_from = n, values_fill = 0) %>%
  column_to_rownames(var = "Virus_ID")

AM_count <- CARD_annotation %>%
  filter(Antibiotic != "") %>%
  separate_rows(Antibiotic, sep = "; ") %>%
  group_by(Antibiotic, Virus_ID) %>%
  summarise(n = sum(!is.na(Antibiotic))) %>%
  pivot_wider(names_from = Antibiotic, values_from = n, values_fill = 0) %>%
  column_to_rownames(var = "Virus_ID")

#********************************************
# Estimate AM feature abundances and counts
#********************************************

# Load Virus contigs abundance table
Virus_abundance <- read.delim("Abundance_table/CS_Abundance_Table_17052024.txt", sep = "\t", header = T, check.names = F)
Virus_abundance_present <- Virus_abundance[rowSums(Virus_abundance) != 0, , drop = FALSE]

# Estimate AMR abundances
AM_abundance <- data.frame(t(calculate_AMR_abundance(Virus_abundance, AM_count)),check.names = F)
AMR_Gene_abundance <- data.frame(t(calculate_AMR_abundance(Virus_abundance, AMR_Gene_count)),check.names = F)
AMR_Class_abundance <- data.frame(t(calculate_AMR_abundance(Virus_abundance, AMR_Class_count)),check.names = F)
AMR_Mechanism_abundance <- data.frame(t(calculate_AMR_abundance(Virus_abundance, AMR_Mechanism_count)),check.names = F)

# Estimate the AMR gene counts per sample
AM_presence <- data.frame(t(calculate_AMR_counts(Virus_abundance, AM_count)),check.names = F) 
AMR_Gene_presence <- data.frame(t(calculate_AMR_counts(Virus_abundance, AMR_Gene_count)),check.names = F) 
AMR_Class_presence <- data.frame(t(calculate_AMR_counts(Virus_abundance, AMR_Class_count)),check.names = F) 
AMR_Mechanism_presence <- data.frame(t(calculate_AMR_counts(Virus_abundance, AMR_Mechanism_count)),check.names = F) 

#********************************************
# Summarizing results (adding to metadata tables)
#********************************************

# Order results according to Sample_metadata (and filter infant samples from M6 and M12)
AM_abundance <- AM_abundance[rownames(AM_abundance) %in% Sample_metadata$bioSampleId,]
AMR_Gene_abundance <- AMR_Gene_abundance[rownames(AMR_Gene_abundance) %in% Sample_metadata$bioSampleId,]
AMR_Class_abundance <- AMR_Class_abundance[rownames(AMR_Class_abundance) %in% Sample_metadata$bioSampleId,]
AMR_Mechanism_abundance <- AMR_Mechanism_abundance[rownames(AMR_Mechanism_abundance) %in% Sample_metadata$bioSampleId,]

AM_presence <- AM_presence[rownames(AM_presence) %in% Sample_metadata$bioSampleId,]
AMR_Gene_presence <- AMR_Gene_presence[rownames(AMR_Gene_presence) %in% Sample_metadata$bioSampleId,]
AMR_Class_presence <- AMR_Class_presence[rownames(AMR_Class_presence) %in% Sample_metadata$bioSampleId,]
AMR_Mechanism_presence <- AMR_Mechanism_presence[rownames(AMR_Mechanism_presence) %in% Sample_metadata$bioSampleId,]

# Add total AMR gene abundance and total AMR gene count variables to Sample_metadata
AMR_Gene_abundance$AMR_gene_abundance <- rowSums(AMR_Gene_abundance)
AMR_Gene_presence$AMR_gene_counts <- rowSums(AMR_Gene_presence)
Sample_metadata$AMR_gene_abundance <- AMR_Gene_abundance$AMR_gene_abundance
Sample_metadata$AMR_gene_counts <- AMR_Gene_presence$AMR_gene_counts

# Add variables with total AMR gene abundance and presence also to Sample_metadata_infants
AMR_Gene_abundance_infants <- AMR_Gene_abundance[rownames(AMR_Gene_abundance) %in% Sample_metadata_infants$bioSampleId,]
AMR_Gene_presence_infants <- AMR_Gene_presence[rownames(AMR_Gene_presence) %in% Sample_metadata_infants$bioSampleId,]
Sample_metadata_infants$AMR_gene_abundance <- AMR_Gene_abundance_infants$AMR_gene_abundance
Sample_metadata_infants$AMR_gene_counts <- AMR_Gene_presence_infants$AMR_gene_counts

# Add AMR counts to Viral metadata
AMR_genes <- data.frame(cbind("Virus_ID"= rownames(AMR_Gene_count), "AMR_genes"= rowSums(AMR_Gene_count)))
Virus_metadata <- left_join(Virus_metadata, AMR_genes, by="Virus_ID")
Virus_metadata$AMR_genes[is.na(Virus_metadata$AMR_genes)] <- 0


#****************
# Save output
#****************
# Save Sample metadata
write.table(Sample_metadata,"Metadata_CS/CS_Baby_Biome_Sample_Metadata_17052024.txt", sep = "\t", 
            row.names = T, col.names = T, quote = FALSE)
write.table(Sample_metadata_infants,"Metadata_CS/CS_Baby_Biome_Sample_Metadata_Infants_17052024.txt", sep = "\t", 
            row.names = T, col.names = T, quote = FALSE)
write.table(Virus_metadata,"Metadata_CS/CS_Baby_Biome_Viral_Metadata_17052024.txt", sep = "\t", 
            row.names = T, col.names = T, quote = FALSE)

# Save AMR abundance and count tables
write.table(AM_abundance,"Abundance_table/CS_Baby_Biome_AM_Abundance_Table.txt", sep = "\t", 
            row.names = T, col.names = T, quote = FALSE)
write.table(AMR_Gene_abundance,"Abundance_table/CS_Baby_Biome_AMR_Gene_Abundance_Table.txt", sep = "\t", 
            row.names = T, col.names = T, quote = FALSE)
write.table(AMR_Class_abundance,"Abundance_table/CS_Baby_Biome_AMR_Class_Abundance_Table.txt", sep = "\t", 
            row.names = T, col.names = T, quote = FALSE)
write.table(AMR_Mechanism_abundance,"Abundance_table/CS_Baby_Biome_AMR_Mechanism_Abundance_Table.txt", sep = "\t", 
            row.names = T, col.names = T, quote = FALSE)


write.table(AM_presence,"Abundance_table/CS_Baby_Biome_AM_Count_Table.txt", sep = "\t", 
            row.names = T, col.names = T, quote = FALSE)
write.table(AMR_Gene_presence,"Abundance_table/CS_Baby_Biome_AMR_Gene_Count_Table.txt", sep = "\t", 
            row.names = T, col.names = T, quote = FALSE)
write.table(AMR_Class_presence,"Abundance_table/CS_Baby_Biome_AMR_Class_Count_Table.txt", sep = "\t", 
            row.names = T, col.names = T, quote = FALSE)
write.table(AMR_Mechanism_presence,"Abundance_table/CS_Baby_Biome_AMR_Mechanism_Count_Table.txt", sep = "\t", 
            row.names = T, col.names = T, quote = FALSE)
