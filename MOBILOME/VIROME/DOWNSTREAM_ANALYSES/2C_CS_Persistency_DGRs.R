################################################################################
##### CS Baby Biome: Viral persistence - DGR analysis and co-ocurrence with hosts
### Author(s):Asier Fernández-Pato
### Last updated: 18th February, 2025
################################################################################

#****************
# Load modules
#****************
library(ggplot2)

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
Virus_metadata_infants <- Virus_metadata[Virus_metadata$Virus_ID %in% rownames(Abundance_table_infants),]

# Generate the metadata and abundance tables for mothers 
Sample_metadata_mothers <- Sample_metadata[Sample_metadata$Type == "Mother",]
Abundance_table_mothers <- Abundance_table[,Sample_metadata$Type == "Mother"]
Abundance_table_mothers <- Abundance_table_mothers[rowSums(Abundance_table_mothers)>0,] 

# Load bacterial abundances
bacterial_abundances <- read.delim("4_BACTERIAL_ANALYSIS/Metaphlan_4_genus_21062024.txt")
bacterial_abundances_infants <- bacterial_abundances[rownames(bacterial_abundances) %in% Sample_metadata_infants$bioSampleId,]

##***************************************
# 1. DGR analysis in mothers and infants
##***************************************

# First, subselect only phages present in infants (W01-W06) and mothers
Virus_metadata_infants <- Virus_metadata[Virus_metadata$Virus_ID %in% rownames(Abundance_table_infants),]
Virus_metadata_mothers <- Virus_metadata[Virus_metadata$Virus_ID %in% rownames(Abundance_table_mothers),]
Virus_metadata$Origin <- ifelse(Virus_metadata$Virus_ID %in% Virus_metadata_infants$Virus_ID & 
                                  Virus_metadata$Virus_ID %in%  Virus_metadata_mothers$Virus_ID, "Both",
                                ifelse(Virus_metadata$Virus_ID %in% Virus_metadata_infants$Virus_ID, "Infant", "Maternal"))

# Compare proportion of infant (W1-W6) and maternal vOTUs with DGRs
contingency_table <- table(Virus_metadata$DGR, Virus_metadata$Origin)
proportions <- prop.table(contingency_table, margin = 2)

# Perform chi-squared test (we now have a 2x3 table)
chi_squared_test_result <- chisq.test(contingency_table)


# Plot of DGR frequency in mothers and infants
contingency_table_plot <- contingency_table[, -1] # Remove the "Both" column
contingency_table_plot[, "Infant"] <- contingency_table_plot[, "Infant"] + contingency_table[, "Both"]
contingency_table_plot[, "Maternal"] <- contingency_table_plot[, "Maternal"] + contingency_table[, "Both"]
prop_table_plot <- prop.table(contingency_table_plot, margin = 1)  

# Convert to data frame for ggplot
DGR_origin_ggplot <- as.data.frame.matrix(prop_table_plot)
DGR_origin_ggplot$Type <- colnames(DGR_origin_ggplot)
DGR_origin_ggplot$DGR <- rownames(DGR_origin_ggplot)
DGR_origin_ggplot<- DGR_origin_ggplot %>%
  pivot_longer(cols = c("Infant", "Maternal"), 
               names_to = "Group", 
               values_to = "Proportion")
DGR_origin_ggplot$Group <- factor(DGR_origin_ggplot$Group)

pdf('2_DIVERSITY/Plots/DGRs_mother_infant.pdf', width=4.2, height=3.2)
ggplot(DGR_origin_ggplot, aes(x = DGR, y = Proportion, fill = Group)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "DGR Presence", y = "Proportion", fill = "Group") +
  scale_fill_manual(name = "Origin", values = c(alpha("#008080", 0.7),alpha("#4B0082", 0.8))) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 14),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10))
  )
dev.off()


##***************************************
# 2. DGRs and viral persistency
##***************************************
# Then, check if persistent or highly-persistent vOTUs are more likely to have DGRs

# Generate a contingency table 
contingency_table_persist <- table(Virus_metadata_infants$DGR, Virus_metadata_infants$Persistent)
contingency_table_high_persist <- table(Virus_metadata_infants$DGR, Virus_metadata_infants$Highly_persistent)

# Perform Fisher's exact test
fisher_test_result_persist <- fisher.test(contingency_table_persist)
fisher_test_result_high_persist <- fisher.test(contingency_table_high_persist)

##***************************************
# 3. Co-ocurrence of persistent vOTUs with viral hosts
##***************************************

#Select persistent vOTUs
Persistent_vOTUs <- Virus_metadata_infants$Virus_ID[Virus_metadata_infants$Persistent =="Yes"]

# Select bacterial genera that have been assigned to the persistent vOTUs
genus_keep <- unique(Virus_metadata[Virus_metadata$Virus_ID %in% Persistent_vOTUs, "Bacterial_genus_host"])
genus_keep_present <- genus_keep[genus_keep %in% colnames(bacterial_abundances_infants)]
bacterial_abundances_infants <- bacterial_abundances_infants[,genus_keep_present ]

# Process vOTU abundance table to remove vOTus that are not persistent
Abundance_table_infants <- data.frame(t(Abundance_table_infants))
Abundance_table_infants_filtered <- as.data.frame(Abundance_table_infants[,Persistent_vOTUs])

# Estimate the Jaccard Index representing the fraction of metagenomes that contain both the vOTU and the predicted host
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
  
  # Subset the host_presence to only include samples where vOTU is present
  vOTU_present_indices <- which(vOTU_presence)  
  host_presence_filtered <- host_presence[vOTU_present_indices] 
  vOTU_presence_filtered <- vOTU_presence[vOTU_present_indices]  
  
  # Compute Jaccard index only for samples where vOTU is present
  jaccard_value <- jaccard_index(vOTU_presence_filtered, host_presence_filtered)
  
  # Get the results
  Jaccard_results <- rbind(Jaccard_results, data.frame(vOTU = vOTU,
                                                       host = host,
                                                       jaccard_index = jaccard_value,
                                                       stringsAsFactors = FALSE))
}



##*************
# Save output
#**************
write.table(Jaccard_results, "CS_persistent_vOTU_Host_Jaccard_Index.txt", sep = "\t", row.names = F, quote = FALSE)

