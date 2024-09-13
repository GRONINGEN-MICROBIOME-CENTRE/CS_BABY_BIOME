################################################################################
##### CS Baby Biome: AMR abundance estimation in plasmids
### Author(s): Asier Fernández-Pato
### Last updated: 22nd June, 2024
################################################################################

# Set working directory 
setwd("~/Desktop/PhD/Projects/CS Baby Biome/Final_analysis/PLASMIDS/")

#****************
# Load libraries
#****************
library(dplyr)
library(tidyverse)
library(lme4)
library(lmerTest)
library(ggplot2)

#****************
# Define functions
#****************

# Function to estimate the abundance of AMR features (genes, classes, mechanisms). Input:
# Plasmid_abundance: Plasmid bundance table (plasmids in rows and samples in columns)
# AM_count: AM feature count table (plasmids in rows and AM feature in columns)
calculate_AMR_abundance <- function(Plasmid_abundance, AM_count) {
  # Initialize an empty data frame for the gene abundance table
  gene_abundance_table <- data.frame(matrix(ncol = ncol(Plasmid_abundance), nrow = ncol(AM_count)))
  colnames(gene_abundance_table) <- colnames(Plasmid_abundance)
  rownames(gene_abundance_table) <- colnames(AM_count)
  
  # Generate a table with merged plasmid abundances and AMR counts
  merged_data <- merge(Plasmid_abundance, AM_count, by = 0, all.x = TRUE)
  merged_data <- replace(merged_data, is.na(merged_data), 0)
  rownames(merged_data) <- merged_data[, 1]
  merged_data <- merged_data[, -1]
  
  # Loop through samples and AMRs to estimate the abundances
  for (sample in colnames(merged_data)[1:(ncol(Plasmid_abundance))]) {
    plasmids_present <- rownames(merged_data)[merged_data[, sample] != 0] # plasmids present in the sample
    abundance_plasmids <- merged_data[merged_data[, sample] != 0, sample] # the abundance of plasmids present in the sample
    
    for (AMR in colnames(merged_data)[(ncol(Plasmid_abundance) + 1):ncol(merged_data)]) {
      plasmids_AMR <- rownames(merged_data)[merged_data[, AMR] != 0] # plasmids that have the AMR
      counts_AMR_plasmids <-  merged_data[merged_data[, AMR] != 0, AMR] # counts of AMRs the plasmids
      
      # Select plasmids present in the sample that have the AMR gene / class
      plasmids <- intersect(plasmids_AMR, plasmids_present) # plasmids with the AMR present in the sample
      counts <- as.numeric(counts_AMR_plasmids[plasmids_AMR %in% plasmids_present]) # the AMR counts of the plasmids present in the sample
      abundance <- as.numeric(abundance_plasmids[match(plasmids, plasmids_present)]) # the abundance of the plasmids with the AMR present in the sample
      gene_abundance <- ifelse(length(counts) > 0 && length(abundance) > 0, sum(counts * abundance), 0)
      
      # Add AMR abundance value to the corresponding sample column
      gene_abundance_table[AMR, sample] <- gene_abundance
    }
  }
  return(gene_abundance_table)
}

calculate_AMR_counts <- function(Plasmid_abundance, AM_count) {
  # Initialize an empty data frame for the count table
  count_table <- data.frame(matrix(ncol = ncol(Plasmid_abundance), nrow = ncol(AM_count)))
  colnames(count_table) <- colnames(Plasmid_abundance)
  rownames(count_table) <- colnames(AM_count)
  
  # Generate a table with merged plasmid abundances and AMR counts
  merged_data <- merge(Plasmid_abundance, AM_count, by = 0, all.x = TRUE)
  merged_data <- replace(merged_data, is.na(merged_data), 0)
  rownames(merged_data) <- merged_data[, 1]
  merged_data <- merged_data[, -1]
  
  # Loop through samples and AMRs to estimate the abundances
  for (sample in colnames(merged_data)[1:(ncol(Plasmid_abundance))]) {
    plasmids_present <- rownames(merged_data)[merged_data[, sample] != 0] # plasmids present in the sample
    abundance_plasmids <- merged_data[merged_data[, sample] != 0, sample] # the abundance of plasmids present in the sample
    
    for (AMR in colnames(merged_data)[(ncol(Plasmid_abundance) + 1):ncol(merged_data)]) {
      plasmids_AMR <- rownames(merged_data)[merged_data[, AMR] != 0] # plasmids that have the AMR
      counts_AMR_plasmids <-  merged_data[merged_data[, AMR] != 0, AMR] # counts of AMRs the plasmids
      
      # Select plasmids present in the sample that have the AMR gene / class
      plasmids <- intersect(plasmids_AMR, plasmids_present) # plasmids with the AMR present in the sample
      counts <- as.numeric(counts_AMR_plasmids[plasmids_AMR %in% plasmids_present]) # the AMR counts of the plasmids present in the sample
      abundance <- as.numeric(abundance_plasmids[match(plasmids, plasmids_present)]) # the abundance of the plasmids with the AMR present in the sample
      feature_count <- ifelse(length(counts) > 0 && length(abundance) > 0, sum(counts), 0)
      
      # Add AMR feature count value to the corresponding sample column
      count_table[AMR, sample] <- feature_count
    }
  }
  return(count_table)
}

Sample_metadata <- read.delim("METADATA_TABLES/CS_Baby_Biome_Sample_Metadata_25062024.txt")
Sample_metadata_infants <- read.delim("METADATA_TABLES/CS_Baby_Biome_Sample_Metadata_Infants_25062024.txt")
Plasmid_metadata <- read.delim("METADATA_TABLES/CS_Baby_Biome_Plasmid_Metadata_25062024.txt")


#*******************************
# Processing of CARD annotation with RGI
#*******************************

# Read CARD AMR annotation table
CARD_annotation <- read.delim("5_AMR_ANALYSIS/RGI_CARD_annotation.txt", sep = "\t", header = T, check.names = F)
CARD_annotation$Plasmid_ID <- sub(" #.*", "", CARD_annotation$ORF_ID)
CARD_annotation$Plasmid_ID <- sub("_[^_]*$", "", CARD_annotation$Plasmid_ID)
CARD_annotation$Contig <-NULL

# Generate tables with the counts of (per plasmid):
# A) AMR Gene
# B) AMR Classes
# C) AMR Mechanisms
# D) Antibiotics/Antimicrobials (group those closely related: e.g. Colistin A and Colistin B)
AMR_Gene_count <- CARD_annotation %>%
  count(Plasmid_ID, Best_Hit_ARO) %>%
  pivot_wider(names_from = Best_Hit_ARO, values_from = n, values_fill = 0) %>%
  column_to_rownames(var = "Plasmid_ID")

AMR_Class_count <- CARD_annotation %>%
  filter(`Drug Class` != "") %>%
  separate_rows(`Drug Class`, sep = "; ") %>%
  group_by(`Drug Class`, Plasmid_ID) %>%
  summarise(n = sum(!is.na(`Drug Class`))) %>%
  pivot_wider(names_from = `Drug Class`, values_from = n, values_fill = 0) %>%
  column_to_rownames(var = "Plasmid_ID")

AMR_Class_count$`streptogramin A antibiotic` <- NULL # a more general class already present
AMR_Class_count$`streptogramin B antibiotic` <- NULL

AMR_Mechanism_count <- CARD_annotation %>%
  filter(`Resistance Mechanism` != "") %>%
  separate_rows(`Resistance Mechanism`, sep = "; ") %>%
  group_by(`Resistance Mechanism`, Plasmid_ID) %>%
  summarise(n = sum(!is.na(`Resistance Mechanism`))) %>%
  pivot_wider(names_from = `Resistance Mechanism`, values_from = n, values_fill = 0) %>%
  column_to_rownames(var = "Plasmid_ID")

AM_count <- CARD_annotation %>%
  filter(Antibiotic != "") %>%
  separate_rows(Antibiotic, sep = "; ") %>%
  group_by(Antibiotic, Plasmid_ID) %>%
  summarise(n = sum(!is.na(Antibiotic))) %>%
  pivot_wider(names_from = Antibiotic, values_from = n, values_fill = 0) %>%
  column_to_rownames(var = "Plasmid_ID")

# Remove columns with a more general AB already present (gentamicin and polymyxin)
AM_count[c('gentamicin A', 'gentamicin B', 'gentamicin C', 'polymyxin B1', 
           'polymyxin B2', 'polymyxin B3', 'polymyxin B4')] <- NULL


#********************************************
# Estimate AM feature abundances and counts
#********************************************

# Load plasmid contigs abundance table
Plasmid_abundance <- read.delim("ABUNDANCE_TABLES/CS_Baby_Biome_Plasmid_Abundance_Table.txt", sep = "\t", header = T, check.names = F)
Plasmid_abundance_present <- Plasmid_abundance[rowSums(Plasmid_abundance) != 0, , drop = FALSE]

# Estimate AMR abundances
AM_abundance <- data.frame(t(calculate_AMR_abundance(Plasmid_abundance, AM_count)),check.names = F)
AMR_Gene_abundance <- data.frame(t(calculate_AMR_abundance(Plasmid_abundance, AMR_Gene_count)),check.names = F)
AMR_Class_abundance <- data.frame(t(calculate_AMR_abundance(Plasmid_abundance, AMR_Class_count)),check.names = F)
AMR_Mechanism_abundance <- data.frame(t(calculate_AMR_abundance(Plasmid_abundance, AMR_Mechanism_count)),check.names = F)

# Estimate the AMR gene counts per sample
AM_presence <- data.frame(t(calculate_AMR_counts(Plasmid_abundance, AM_count)),check.names = F) 
AMR_Gene_presence <- data.frame(t(calculate_AMR_counts(Plasmid_abundance, AMR_Gene_count)),check.names = F) 
AMR_Class_presence <- data.frame(t(calculate_AMR_counts(Plasmid_abundance, AMR_Class_count)),check.names = F) 
AMR_Mechanism_presence <- data.frame(t(calculate_AMR_counts(Plasmid_abundance, AMR_Mechanism_count)),check.names = F) 

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

# AMR counts already added to Plasmid metadata

#*******************************
# Plotting maternal results
#*******************************

# As association analysis will be focused only on infants, here we:
#-- Test if ARG counts are different between mothers and infants
#-- Generate a plot for to compare ARG counts between mothers and infants

MM_type_ARG_counts <- lmer(AMR_gene_counts ~ Type + DNA_concentration_ng_ul + read_depth +  (1|CS_BABY_BIOME_ID),
                           REML = F, data = Sample_metadata)

summary(MM_type_ARG_counts)


pdf('5_AMR_ANALYSIS/PLOTS/Mother_total_AMR_count.pdf', width=3.9, height=3.5)
ggplot(Sample_metadata, 
       aes(x=Type, y=AMR_gene_counts, fill = Type)) +
  geom_boxplot(aes(fill=Type), alpha=0.7, width=0.8) +
  geom_point(aes(color=Type), size=1.2, position=position_jitterdodge(), alpha=0.7) +
  labs(x = "Sample Type", y = "AMR Gene Count", fill = "Type") +
  scale_fill_manual(values=c("#008080", "#4B0082")) + 
  scale_color_manual(values=c("#008080", "#4B0082")) +
  theme_classic()  +
  theme(plot.title = element_text(size=18, hjust = 0.5, face="bold"), 
        axis.text = element_text(size=14),
        axis.title = element_text(size=16), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.position = "bottom",
        legend.text = element_text(size=12),
        legend.title = element_text(size=13)) +
  guides(color = FALSE)
dev.off()


#****************
# Save output
#****************
# Save Sample metadata
write.table(Sample_metadata,"METADATA_TABLES/CS_Baby_Biome_Sample_Metadata_25062024.txt", sep = "\t", 
            row.names = T, col.names = T, quote = FALSE)
write.table(Sample_metadata_infants,"METADATA_TABLES/CS_Baby_Biome_Sample_Metadata_Infants_25062024.txt", sep = "\t", 
            row.names = T, col.names = T, quote = FALSE)

# Save AMR abundance and count tables
write.table(AM_abundance,"ABUNDANCE_TABLES/CS_Baby_Biome_AM_Abundance_Table.txt", sep = "\t", 
            row.names = T, col.names = T, quote = FALSE)
write.table(AMR_Gene_abundance,"ABUNDANCE_TABLES/CS_Baby_Biome_AMR_Gene_Abundance_Table.txt", sep = "\t", 
            row.names = T, col.names = T, quote = FALSE)
write.table(AMR_Class_abundance,"ABUNDANCE_TABLES/CS_Baby_Biome_AMR_Class_Abundance_Table.txt", sep = "\t", 
            row.names = T, col.names = T, quote = FALSE)
write.table(AMR_Mechanism_abundance,"ABUNDANCE_TABLES/CS_Baby_Biome_AMR_Mechanism_Abundance_Table.txt", sep = "\t", 
            row.names = T, col.names = T, quote = FALSE)


write.table(AM_presence,"ABUNDANCE_TABLES/CS_Baby_Biome_AM_Count_Table.txt", sep = "\t", 
            row.names = T, col.names = T, quote = FALSE)
write.table(AMR_Gene_presence,"ABUNDANCE_TABLES/CS_Baby_Biome_AMR_Gene_Count_Table.txt", sep = "\t", 
            row.names = T, col.names = T, quote = FALSE)
write.table(AMR_Class_presence,"ABUNDANCE_TABLES/CS_Baby_Biome_AMR_Class_Count_Table.txt", sep = "\t", 
            row.names = T, col.names = T, quote = FALSE)
write.table(AMR_Mechanism_presence,"ABUNDANCE_TABLES/CS_Baby_Biome_AMR_Mechanism_Count_Table.txt", sep = "\t", 
            row.names = T, col.names = T, quote = FALSE)
