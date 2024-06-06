################################################################################
##### CS Baby Biome: AMR gene abundance analysis in MAGs RGI-CARD2023
### Author(s):Asier Fernández-Pato
### Last updated: 1st February, 2024
################################################################################

#****************
# Load modules
#****************
library(tidyverse)
library(dplyr)
library(ggplot2)


# Set working directory 
setwd("~/Desktop/PhD/Projects/CS Baby Biome/New_analysis/")

#****************
# Define functions
#****************

# Function to estimate the abundance of AMR features (genes, classes, mechanisms). Input:
# Abundance table: MAG abundance table (MAGs in rows and samples in columns)
# AM_count: AM feature count table (MAGs in rows and AM feature in columns)
calculate_AMR_abundance <- function(Bin_abundance, AM_count) {
  # Initialize an empty data frame for the gene abundance table
  abundance_table <- data.frame(matrix(ncol = ncol(Bin_abundance), nrow = ncol(AM_count)))
  colnames(abundance_table) <- colnames(Bin_abundance)
  rownames(abundance_table) <- colnames(AM_count)
  
  # Generate a table with merged Bin abundances and AMR counts
  merged_data <- merge(Bin_abundance, AM_count, by = 0, all.x = TRUE)
  merged_data <- replace(merged_data, is.na(merged_data), 0)
  rownames(merged_data) <- merged_data[, 1]
  merged_data <- merged_data[, -1]
  
  # Loop through samples and AMRs to estimate the abundances
  for (sample in colnames(merged_data)[1:(ncol(Bin_abundance))]) {
    Bins_present <- rownames(merged_data)[merged_data[, sample] != 0] # Bins present in the sample
    abundance_Bins <- merged_data[merged_data[, sample] != 0, sample] # the abundance of Bins present in the sample
    
    for (AMR in colnames(merged_data)[(ncol(Bin_abundance) + 1):ncol(merged_data)]) {
      Bins_AMR <- rownames(merged_data)[merged_data[, AMR] != 0] # Bins that have the AMR
      counts_AMR_Bins <-  merged_data[merged_data[, AMR] != 0, AMR] # counts of AMRs the Bins
      
      # Select Bins present in the sample that have the AMR gene / class
      Bins <- intersect(Bins_AMR, Bins_present) # Bins with the AMR present in the sample
      counts <- as.numeric(counts_AMR_Bins[Bins_AMR %in% Bins_present]) # the AMR counts of the Bins present in the sample
      abundance <- as.numeric(abundance_Bins[match(Bins, Bins_present)]) # the abundance of the Bins with the AMR present in the sample
      gene_abundance <- ifelse(length(counts) > 0 && length(abundance) > 0, sum(counts * abundance), 0)
      
      # Add AMR abundance value to the corresponding sample column
      abundance_table[AMR, sample] <- gene_abundance
    }
  }
  return(abundance_table)
}

calculate_AMR_counts <- function(Bin_abundance, AM_count) {
  # Initialize an empty data frame for the count table
  count_table <- data.frame(matrix(ncol = ncol(Bin_abundance), nrow = ncol(AM_count)))
  colnames(count_table) <- colnames(Bin_abundance)
  rownames(count_table) <- colnames(AM_count)
  
  # Generate a table with merged Bin abundances and AMR counts
  merged_data <- merge(Bin_abundance, AM_count, by = 0, all.x = TRUE)
  merged_data <- replace(merged_data, is.na(merged_data), 0)
  rownames(merged_data) <- merged_data[, 1]
  merged_data <- merged_data[, -1]
  
  # Loop through samples and AMRs to estimate the abundances
  for (sample in colnames(merged_data)[1:(ncol(Bin_abundance))]) {
    Bins_present <- rownames(merged_data)[merged_data[, sample] != 0] # Bins present in the sample
    abundance_Bins <- merged_data[merged_data[, sample] != 0, sample] # the abundance of Bins present in the sample
    
    for (AMR in colnames(merged_data)[(ncol(Bin_abundance) + 1):ncol(merged_data)]) {
      Bins_AMR <- rownames(merged_data)[merged_data[, AMR] != 0] # Bins that have the AMR
      counts_AMR_Bins <-  merged_data[merged_data[, AMR] != 0, AMR] # counts of AMRs the Bins
      
      # Select Bins present in the sample that have the AMR gene / class
      Bins <- intersect(Bins_AMR, Bins_present) # Bins with the AMR present in the sample
      counts <- as.numeric(counts_AMR_Bins[Bins_AMR %in% Bins_present]) # the AMR counts of the Bins present in the sample
      abundance <- as.numeric(abundance_Bins[match(Bins, Bins_present)]) # the abundance of the Bins with the AMR present in the sample
      feature_count <- ifelse(length(counts) > 0 && length(abundance) > 0, sum(counts), 0)
      
      # Add AMR feature count value to the corresponding sample column
      count_table[AMR, sample] <- feature_count
    }
  }
  return(count_table)
}

#****************
# Processing of CARD annotation with RGI
#****************

# Read CARD AMR annotation table
CARD_annotation <- read.delim("AMR_BACTERIA/B_APPROACH/Merged_AMR_annotation_strict.txt", sep = "\t", header = T, check.names = F)
CARD_annotation$Contig <- result <- sub("_[^_]+$", "", CARD_annotation$Contig)

# Generate tables with the counts of (per MAG):
# A) AMR Gene
# B) AMR Classes
# C) AMR Mechanisms
# D) Antibiotics/Antimicrobials (group those closely related: e.g. Colistin A and Colistin B)
AMR_Gene_count <- CARD_annotation %>%
  count(Contig, Best_Hit_ARO) %>%
  pivot_wider(names_from = Best_Hit_ARO, values_from = n, values_fill = 0) %>%
  column_to_rownames(var = "Contig")

AMR_Class_count <- CARD_annotation %>%
  filter(`Drug Class` != "") %>%
  separate_rows(`Drug Class`, sep = "; ") %>%
  group_by(`Drug Class`, Contig) %>%
  summarise(n = sum(!is.na(`Drug Class`))) %>%
  pivot_wider(names_from = `Drug Class`, values_from = n, values_fill = 0) %>%
  column_to_rownames(var = "Contig")

AMR_Class_count$`streptogramin A antibiotic` <- NULL # a more general class already present
AMR_Class_count$`streptogramin B antibiotic` <- NULL

AMR_Mechanism_count <- CARD_annotation %>%
  filter(`Resistance Mechanism` != "") %>%
  separate_rows(`Resistance Mechanism`, sep = "; ") %>%
  group_by(`Resistance Mechanism`, Contig) %>%
  summarise(n = sum(!is.na(`Resistance Mechanism`))) %>%
  pivot_wider(names_from = `Resistance Mechanism`, values_from = n, values_fill = 0) %>%
  column_to_rownames(var = "Contig")

AM_count <- CARD_annotation %>%
  filter(Antibiotic != "") %>%
  separate_rows(Antibiotic, sep = "; ") %>%
  group_by(Antibiotic, Contig) %>%
  summarise(n = sum(!is.na(Antibiotic))) %>%
  pivot_wider(names_from = Antibiotic, values_from = n, values_fill = 0) %>%
  column_to_rownames(var = "Contig")

# Remove columns with a more general AB already present (bacitracin, gentamicin, polymyxin)
AM_count[c("bacitracin A", "bacitracin B", "bacitracin F", 
           'gentamicin B', 'gentamicin C', 'polymyxin B1', 
           'polymyxin B2', 'polymyxin B3', 'polymyxin B4')] <- NULL
  
#**************************************
# Estimate AM feature abundances and count tables
#**************************************

# Load bin abundance table and process it
Bin_abundance <- read.delim("AMR_BACTERIA/B_APPROACH/Merged_bin_abundances.txt", sep = "\t", header = T, check.names = F)
rownames(Bin_abundance) <- Bin_abundance$Genome
Bin_abundance$Genome <- NULL
colnames(Bin_abundance) <- sub("^(.*?)_.*", "\\1", colnames(Bin_abundance))
Bin_abundance_present <- Bin_abundance[rowSums(Bin_abundance) != 0, , drop = FALSE] # All 516 are present

# Estimate AMR abundances
AM_abundance <- calculate_AMR_abundance(Bin_abundance, AM_count) #117
AMR_Gene_abundance <- calculate_AMR_abundance(Bin_abundance, AMR_Gene_count) #188
AMR_Class_abundance <- calculate_AMR_abundance(Bin_abundance, AMR_Class_count) #30
AMR_Mechanism_abundance <- calculate_AMR_abundance(Bin_abundance, AMR_Mechanism_count) #6

# Estimate the AMR gene counts per sample
AM_presence <- calculate_AMR_counts(Bin_abundance, AM_count) #117
AMR_Gene_presence <- calculate_AMR_counts(Bin_abundance, AMR_Gene_count) #188
AMR_Class_presence <- calculate_AMR_counts(Bin_abundance, AMR_Class_count) #30
AMR_Mechanism_presence <- calculate_AMR_counts(Bin_abundance, AMR_Mechanism_count) #6

# Select only features with a prevalence > 5%
AM_abundance_prev_list <- rownames(AM_abundance)[rowSums(AM_abundance != 0) > 195*0.05] #111
AMR_Gene_abundance_prev_list <- rownames(AMR_Gene_abundance)[rowSums(AMR_Gene_abundance != 0) > 195*0.05] #163
AMR_Class_abundance_prev_list <- rownames(AMR_Class_abundance)[rowSums(AMR_Class_abundance != 0) > 195*0.05] #30
AMR_Mechanism_abundance_prev_list <- rownames(AMR_Mechanism_abundance)[rowSums(AMR_Mechanism_abundance != 0) > 195*0.05] #6

AM_abundance_prev <- AM_abundance[AM_abundance_prev_list,]
AMR_Gene_abundance_prev <- AMR_Gene_abundance[AMR_Gene_abundance_prev_list,]
AMR_Class_abundance_prev  <- AMR_Class_abundance[AMR_Class_abundance_prev_list,]
AMR_Mechanism_abundance_prev <- AMR_Mechanism_abundance[AMR_Mechanism_abundance_prev_list,]

AM_presence_prev <- AM_presence[AM_abundance_prev_list,]
AMR_Gene_presence_prev <- AMR_Gene_presence[AMR_Gene_abundance_prev_list,]
AMR_Class_presence_prev  <- AMR_Class_presence[AMR_Class_abundance_prev_list,]
AMR_Mechanism_presence_prev <- AMR_Mechanism_presence[AMR_Mechanism_abundance_prev_list,]


#****************
# Save output
#****************

# Save AMR abundance and count tables
write.table(AM_abundance,"AMR_BACTERIA/B_APPROACH/ABUNDANCE TABLES/CS_Baby_Biome_AM_Abundance_Table.txt", sep = "\t", 
            row.names = T, col.names = T, quote = FALSE)
write.table(AMR_Gene_abundance,"AMR_BACTERIA/B_APPROACH/ABUNDANCE TABLES/CS_Baby_Biome_AMR_Gene_Abundance_Table.txt", sep = "\t", 
            row.names = T, col.names = T, quote = FALSE)
write.table(AMR_Class_abundance,"AMR_BACTERIA/B_APPROACH/ABUNDANCE TABLES/CS_Baby_Biome_AMR_Class_Abundance_Table.txt", sep = "\t", 
            row.names = T, col.names = T, quote = FALSE)
write.table(AMR_Mechanism_abundance,"AMR_BACTERIA/B_APPROACH/ABUNDANCE TABLES/CS_Baby_Biome_AMR_Mechanism_Abundance_Table.txt", sep = "\t", 
            row.names = T, col.names = T, quote = FALSE)

write.table(AM_presence,"AMR_BACTERIA/B_APPROACH/COUNT TABLES/CS_Baby_Biome_AM_Count_Table.txt", sep = "\t", 
            row.names = T, col.names = T, quote = FALSE)
write.table(AMR_Gene_presence,"AMR_BACTERIA/B_APPROACH/COUNT TABLES/CS_Baby_Biome_AMR_Gene_Count_Table.txt", sep = "\t", 
            row.names = T, col.names = T, quote = FALSE)
write.table(AMR_Class_presence,"AMR_BACTERIA/B_APPROACH/COUNT TABLES/CS_Baby_Biome_AMR_Class_Count_Table.txt", sep = "\t", 
            row.names = T, col.names = T, quote = FALSE)
write.table(AMR_Mechanism_presence,"AMR_BACTERIA/B_APPROACH/COUNT TABLES/CS_Baby_Biome_AMR_Mechanism_Count_Table.txt", sep = "\t", 
            row.names = T, col.names = T, quote = FALSE)

