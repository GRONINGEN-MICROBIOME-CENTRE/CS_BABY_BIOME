################################################################################
##### CS Baby Biome: Host Assignment Analysis - Validation
### Author(s):Asier Fernández-Pato
### Last updated: 18th July, 2024
################################################################################


#****************
# Define functions
#****************
# Function to compute Jaccard index
jaccard_index <- function(x, y) {
  intersection <- sum(x & y)
  union <- sum(x | y)
  return(intersection / union)
}

# Set working directory 
setwd("~/Desktop/PhD/Projects/CS Baby Biome/Final_analysis/VIRUSES/")

# Read abundance table and metadata tables
Abundance_table_infants <- read.delim("Abundance_table/CS_Abundance_Table_INFANTS_17052024.txt")
Sample_metadata_infants <- read.delim("Metadata_CS/CS_Baby_Biome_Sample_Metadata_Infants_17052024.txt")
Virus_metadata <- read.delim("Metadata_CS/CS_Baby_Biome_Viral_Metadata_17052024.txt")

# Load bacterial abundances
bacterial_abundances <- read.delim("4_BACTERIAL_ANALYSIS/Metaphlan_4_genus_21062024.txt")
bacterial_abundances_infants <- bacterial_abundances[rownames(bacterial_abundances) %in% Sample_metadata_infants$bioSampleId,]

##************************************************************************
# 1. Preprocessing
##************************************************************************

##################################
# Preprocess vOTU abundance table
##################################

Abundance_table_infants <- t(Abundance_table_infants)

# Calculate the prevalence of each vOTU
vOTU_nonZero  <- colSums(Abundance_table_infants>0)/nrow(Abundance_table_infants) 

# Get the list of vOTUs present in more than 5% samples
vOTU_keep <- colnames(Abundance_table_infants)[vOTU_nonZero > 0.05]

# Remove the vOTUs that are not prevalent
Abundance_table_infants_filtered <- as.data.frame(Abundance_table_infants[,vOTU_keep])

# Order abundance table according to Sample metadata
Abundance_table_infants_filtered <- Abundance_table_infants_filtered[match(Sample_metadata_infants$bioSampleId,
                                                                           rownames(Abundance_table_infants_filtered)), ]

##################################
# Preprocess host abundance table
##################################

# Select bacterial genera that have been assigned to the prevalent vOTUs
genus_keep <- unique(Virus_metadata[Virus_metadata$Virus_ID %in% vOTU_keep, "Bacterial_genus_host"])
genus_keep_present <- genus_keep[genus_keep %in% colnames(bacterial_abundances_infants)]
bacterial_abundances_infants <- bacterial_abundances_infants[,genus_keep_present ]


##************************************************************************
# 2. Estimating Pearson correlations between vOTUs and their hosts
##************************************************************************

Correlation_results <- data.frame(vOTU = character(),
                      host = character(),
                      rho_Coef = numeric(),
                      p_value = numeric(),
                      p_value_adjusted = numeric(),
                      stringsAsFactors = FALSE)

for (i in 1:ncol(Abundance_table_infants_filtered)) {
  # Get the vOTU and host names
  vOTU <- colnames(Abundance_table_infants_filtered)[i]
  host <- Virus_metadata[Virus_metadata$Virus_ID %in% vOTU, "Bacterial_genus_host"]
  
  # If no host is found, we skip the iteration
  if (!host  %in% colnames(bacterial_abundances_infants)) next
  
  # Get the abundance data for the vOTU and its host
  vOTU_abundance <- Abundance_table_infants_filtered[, i]
  host_abundance <- bacterial_abundances_infants[, host]
  
  # Perform Spearman correlation test
  spearman_test <- cor.test(vOTU_abundance, host_abundance, method = "spearman", use = "complete.obs")
  spearman_corr <- spearman_test$estimate
  p_value <- spearman_test$p.value
  
  # Perform FDR adjustment
  p_value_adjusted <- p.adjust(p_value, method = "BH")
  
  # Get the results
  Correlation_results <- rbind(Correlation_results, data.frame(vOTU = vOTU,
                                       host = host,
                                       rho_Coef = spearman_corr,
                                       p_value = p_value,
                                       p_value_adjusted = p_value_adjusted,
                                       stringsAsFactors = FALSE))
}

Correlation_results <- Correlation_results[order(Correlation_results$p_value_adjusted), ]


##************************************************************************
# 3. Estimating Jaccard Index representing the fraction of metagenomes that contain both the vOTU and the predicted host
##************************************************************************
# *Only Pearson correlation was considered for the manuscript

Jaccard_results <- data.frame(vOTU = character(),
                              host = character(),
                              jaccard_index = numeric(),
                              stringsAsFactors = FALSE)

for (i in 1:ncol(Abundance_table_infants_filtered)) {
  # Get the vOTU and host names
  vOTU <- colnames(Abundance_table_infants_filtered)[i]
  host <- Virus_metadata[Virus_metadata$Virus_ID %in% vOTU, "Bacterial_genus_host"]
  
  # If no host is found, skip to the next iteration
  if (!host %in% colnames(bacterial_abundances_infants)) next
  
  # Get the presence/absence information for the vOTU and its host
  vOTU_presence <- Abundance_table_infants_filtered[, i] > 0
  host_presence <- bacterial_abundances_infants[, host] > 0
  
  # Compute Jaccard index
  jaccard_value <- jaccard_index(vOTU_presence, host_presence)
  
  # Get the results
  Jaccard_results <- rbind(Jaccard_results, data.frame(vOTU = vOTU,
                                                       host = host,
                                                       jaccard_index = jaccard_value,
                                                       stringsAsFactors = FALSE))
}


##*************
# Save output
#**************
write.table(Correlation_results, "5_HOST_ASSIGNMENT/CS_vOTU_Host_Pearson_scorrelation.txt", sep = "\t", row.names = F, quote = FALSE)
write.table(Jaccard_results, "5_HOST_ASSIGNMENT/CS_vOTU_Host_Jaccard_Index.txt", sep = "\t", row.names = F, quote = FALSE)

