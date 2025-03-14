################################################################################
##### CS Baby Biome: Plasmid AMR analysis - associations
### Author(s): Asier Fernández-Pato
### Last updated: 23th February, 2025
################################################################################

#****************
# Load libraries
#****************
library(tidyverse)
library(dplyr)

#****************
# Define functions
#****************
# Function to check ARG presence in a cluster
check_ARG_in_cluster <- function(cluster) {
  cluster_members <- na.omit(cluster[-1])  # Remove representative col and NAs
  total_plasmids <- length(cluster_members)
  
  has_ARG <- cluster_members %in% names(ARG_presence)  # Check for ARG presence
  ARG_count <- sum(has_ARG)
  ARG_proportion <- ifelse(total_plasmids > 0, ARG_count / total_plasmids, NA)
  
  rep_has_ARG <- cluster[1] %in% names(ARG_presence)  # Check if representative has ARG
  
  return(data.frame(
    Representative = cluster[1],
    Total_Plasmids = total_plasmids,
    ARG_Count = ARG_count,
    ARG_Proportion = ARG_proportion,
    Rep_HAS_ARG = rep_has_ARG
  ))
}

# Set working directory 
setwd("~/Desktop/PhD/Projects/CS Baby Biome/Final_analysis/PLASMIDS/")

##************************************************************************
# 1. Circular plasmids: Read Data
##************************************************************************

# Read information about circular plasmids and their clusters
Circular_plasmids <- read.delim("5_AMR_ANALYSIS/CHECK_CIRCULAR_PLASMID_GENOMES/plasmid_sequence_ids.txt",
                                header = F, col.names = "Plasmid_ID")
Circular_plasmids_clusters <- read.delim("5_AMR_ANALYSIS/CHECK_CIRCULAR_PLASMID_GENOMES/all_plasmid_complete_sequences_clusters.tsv",
                                header = F)
Circular_plasmids_clusters_processed <-  data.frame(Circular_plasmids_clusters[,-2], 
                                               tstrsplit(as.character(Circular_plasmids_clusters[,2]), ",", fixed = TRUE))
colnames(vOTU_clusters_present_processed) [1] <- c("rep_seq")

# Read ARG annotations
ARG_annotation <- read.delim("5_AMR_ANALYSIS/CHECK_CIRCULAR_PLASMID_GENOMES/RGI_CARD_annotation.txt",
                                header = T)
ARG_annotation$Plasmid_ID <- sub("_[^_]+$", "", ARG_annotation$Contig)

##************************************************************************
# 1B. Check percentage of AMR detection concordance per cluster
##************************************************************************

# Create a lookup table for ARG presence
ARG_presence <- setNames(rep(TRUE, nrow(ARG_annotation)), ARG_annotation$Plasmid_ID)

# Apply function to all rows in Circular_plasmids_clusters_processed
ARG_summary <- do.call(rbind, apply(Circular_plasmids_clusters_processed, 1, check_ARG_in_cluster))

##************************************************************************
# 1C. Summary statistics
##************************************************************************

# Number of clusters with ARGs
length(which(ARG_summary$ARG_Count>0)) #113

# How many of those have a representative genome without ARG
ARG_summary_present <- ARG_summary[ARG_summary$Rep_HAS_ARG == TRUE, ]
ARG_summary_absent<- ARG_summary[ARG_summary$Rep_HAS_ARG == FALSE, ]
table(ARG_summary_present$Rep_HAS_ARG) #1

# Check median and 1st and 3rd quartile for concondance of ARG presence across clusters with ARG presence
quantile(ARG_summary_present$ARG_Proportion, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)

# Only 7/112 clusters with ARGs in the representative did not show perfectly consistent ARG presence across all members in the cluster
# So median, 1st quartile and 3rd quartile =1 in the consistency
# Only 1/4222 clusters with no ARG in their representative had ARGs detected in any genome within the cluster

# Check how many of the circular PTUs present in early infant samples showed consistent results

# Read abundance table and metadata tables
Plasmid_abundance_infants <- read.delim("ABUNDANCE_TABLES/CS_Abundance_Table_INFANTS_25062024.txt")
Sample_metadata_infants <- read.delim("METADATA_TABLES/CS_Baby_Biome_Sample_Metadata_Infants_25062024.txt")
Renamed_plasmids <- read.delim("METADATA_TABLES/Renamed_plasmids.txt", header = F)
Plasmid_metadata <- read.delim("METADATA_TABLES/CS_Baby_Biome_Plasmid_Metadata_Infants_25062024.txt")

#Select only samples from W01-W06
Plasmid_abundance_infants <-Plasmid_abundance_infants[, colnames(Plasmid_abundance_infants) 
                                                      %in% Sample_metadata_infants$NG_ID]

# Get only circular plasmids present in infants (W01-W06) (n=848)
rename_map <- setNames(Renamed_plasmids$V3, Renamed_plasmids$V1)
ARG_summary$Representative <- rename_map[ARG_summary$Representative]
ARG_summary_W01_W06 <- ARG_summary[ARG_summary$Representative %in% rownames(Plasmid_abundance_infants),]


# Get summary stats for infant circular plasmids
ARG_summary_W01_W06_present <- ARG_summary_W01_W06[ARG_summary_W01_W06$Rep_HAS_ARG == TRUE, ]
ARG_summary_W01_W06_absent <- ARG_summary_W01_W06[ARG_summary_W01_W06$Rep_HAS_ARG == FALSE, ]

# Only 3/27 infant clusters with ARGs in the representative did not show perfectly consistent ARG presence across all members in the cluster
# None (0/821) infant clusters with no ARG in their representative had ARGs detected in any genome within the cluster

##************************************************************************
# 2. Fragmented plasmids: Read Data
##************************************************************************
# Read information about fragmented plasmids and their clusters
Fragmented_plasmids <- read.delim("5_AMR_ANALYSIS/CHECK_FRAGMENTED_PLASMID_GENOMES/sequence_ids_fragments.txt",
                                header = F, col.names = "Plasmid_ID")
Fragmented_plasmids_clusters <- read.delim("5_AMR_ANALYSIS/CHECK_FRAGMENTED_PLASMID_GENOMES/geNomad_plasmid_fragment_sequences_clusters.tsv",
                                         header = F)
Fragmented_plasmids_clusters_processed <-  data.frame(Fragmented_plasmids_clusters[,-2], 
                                                    tstrsplit(as.character(Fragmented_plasmids_clusters[,2]), ",", fixed = TRUE))
colnames(vOTU_clusters_present_processed) [1] <- c("rep_seq")

# Read ARG annotations
ARG_annotation_fragments <- read.delim("5_AMR_ANALYSIS/CHECK_FRAGMENTED_PLASMID_GENOMES/RGI_CARD_annotation_fragments.txt",
                             header = T)
ARG_annotation_fragments$Plasmid_ID <- sub("_[^_]+$", "", ARG_annotation_fragments$Contig)

##************************************************************************
# 2B. Check percentage of AMR detection concordance per cluster
##************************************************************************

# Create a lookup table for ARG presence
ARG_presence <- setNames(rep(TRUE, nrow(ARG_annotation_fragments)), ARG_annotation_fragments$Plasmid_ID)

# Apply function to all rows in Circular_plasmids_clusters_processed
ARG_summary_fragments <- do.call(rbind, apply(Fragmented_plasmids_clusters_processed, 1, check_ARG_in_cluster))

##************************************************************************
# 2C. Summary statistics
##************************************************************************
# Number of clusters with ARGs
length(which(ARG_summary_fragments$ARG_Count>0)) #75

# How many of those have a representative genome without ARG
ARG_summary_fragments_present <- ARG_summary_fragments[ARG_summary_fragments$Rep_HAS_ARG == TRUE, ]
ARG_summary_fragments_absent<- ARG_summary_fragments[ARG_summary_fragments$Rep_HAS_ARG == FALSE, ]
table(ARG_summary_fragments_present$Rep_HAS_ARG) #1

# Only 1/2452 clusters with no ARG in their representative had ARGs detected in any genome within the cluster

##************************************************************************
# Save  output files
##************************************************************************
# Add the renamed representative column
ARG_summary_fragments_absent$Renamed_representative <- rename_map[ARG_summary_fragments_absent$Representative]
ARG_summary_absent$Renamed_representative <- rename_map[ARG_summary_absent$Representative]

write.table(ARG_summary_fragments_absent,"5_AMR_ANALYSIS/ARG_summary_per_fragmented_PTU_no_ARG_representative.txt", sep = "\t", row.names = T, quote = FALSE)
write.table(ARG_summary_absent,"5_AMR_ANALYSIS/ARG_summary_per_circular_PTU_no_ARG_representative.txt", sep = "\t", row.names = T, quote = FALSE)
