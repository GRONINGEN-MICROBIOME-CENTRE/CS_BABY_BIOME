################################################################################
##### CS Baby Biome: Plasmid summary stats
### Author(s): Asier Fernández-Pato
### Last updated: 25th June, 2024
################################################################################

# Set working directory 
setwd("~/Desktop/PhD/Projects/CS Baby Biome/Final_analysis/PLASMIDS/")

#****************
# Load libraries
#****************
library(dplyr)
library(data.table)
library(ggplot2)
library(waffle)
library(dunn.test)

#****************
# Define functions
#****************

# Function to estimate summary stats of numeric variables
summary_stats <- function(vector) {
  summary <- c(
    mean = mean(vector),
    median = median(vector),
    q1 = quantile(vector, 0.25),
    q3 = quantile(vector, 0.75),
    min = min(vector),
    max = max(vector)
  )
  return(summary)
}

# Function to perform Kruskal-Wallis test and Dunn test if results is significant
Kruskal_wallis_dunn_test <- function(response, group) {
  kruskal_test <- kruskal.test(response ~ group)
  print(kruskal_test)
  
  if (kruskal_test$p.value < 0.05) {
    dunn_test_result <- dunn.test(response, group, method = "bh")
    print(dunn_test_result)
  } else {
    cat("Kruskal-Wallis test is not significant; no post-hoc Dunn test performed.\n")
  }
}

##************************************************************************
# 1. Load metadata table for the 195 CS Baby Biome samples 
#*************************************************************************
# Load the metadata table 
Sample_metadata <- read.delim("../Metadata_CS/Metadata_EGA_CS_BABY_BIOME.txt") 

# Update clean read numbers
Sample_metadata$NG_ID <- Sample_metadata$bioSampleId
Clean_reads <- read.delim("../Metadata_CS/CS_Baby_Biome_nreads_195.txt", header = T) 
Sample_metadata$read_depth <- Clean_reads$clean_reads

# Add variable to distinguish between maternal and infant samples
Sample_metadata$Type <- ifelse(Sample_metadata$Timepoint_categorical == "MOM", "Mother", "Infant")

# Exclude infants of M6 and M12 from the metadata
#Sample_metadata <- Sample_metadata[!Sample_metadata$Timepoint_categorical %in% c("M06", "M12"),]
#Abundance_table <- Abundance_table[, colnames(Abundance_table) %in% Sample_metadata$NG_ID]

# Create only infant metadata
Sample_metadata_infants <- Sample_metadata[Sample_metadata$Type == "Infant", ]

##************************************************************************
# 2. Process final abundance table  
#*************************************************************************
# Read abundance table
Abundance_table <- read.delim("ABUNDANCE_TABLES/CS_Baby_Biome_Plasmid_Abundance_Table.txt") 

# Select only the samples from the metadata
Abundance_table <- Abundance_table[, colnames(Abundance_table) %in% Sample_metadata$NG_ID]

# Remove from the table those plasmids not present in any sample
Absent_plasmids <- rownames(Abundance_table[rowSums(Abundance_table)==0, ]) #3000
Abundance_table <- Abundance_table[!rownames(Abundance_table) %in% Absent_plasmids,] #3,477
Present_plasmids <- rownames(Abundance_table) #3,477

# Reorder abundance table to match Final metadata
Abundance_table <- Abundance_table[, Sample_metadata$NG_ID]

# Repeat processing for infant samples
Abundance_table_infants <- Abundance_table[, colnames(Abundance_table) %in% Sample_metadata_infants$NG_ID]
Absent_plasmids_infants <- rownames(Abundance_table_infants[rowSums(Abundance_table_infants)==0, ]) #1,478
Abundance_table_infants <- Abundance_table_infants[!rownames(Abundance_table_infants) %in% Absent_plasmids_infants,] 
Present_plasmids_infants <- rownames(Abundance_table_infants) #1,999
Abundance_table_infants <- Abundance_table_infants[, Sample_metadata_infants$NG_ID] #1,999

#****************
# Write results
#****************
# Abundance table and list of present and excluded plasmids (with no presence in CS samples)
write.table(Abundance_table,"ABUNDANCE_TABLES/CS_Abundance_Table_25062024.txt", sep = "\t", row.names = T, quote = FALSE)
write.table(Absent_plasmids,"ABUNDANCE_TABLES/plasmids_with_no_presence.txt", sep = "\n", row.names = F, col.names = F, quote = FALSE)
write.table(Present_plasmids,"ABUNDANCE_TABLES/Present_plasmids.txt", sep = "\n", row.names = F, col.names = F, quote = FALSE)
write.table(Abundance_table_infants,"ABUNDANCE_TABLES/CS_Abundance_Table_INFANTS_25062024.txt", sep = "\t", row.names = T, quote = FALSE)
write.table(Absent_plasmids_infants,"ABUNDANCE_TABLES/plasmids_with_no_presence_INFANTS.txt", sep = "\n", row.names = F, col.names = F, quote = FALSE)
write.table(Present_plasmids_infants,"ABUNDANCE_TABLES/Present_plasmids_INFANTS.txt", sep = "\n", row.names = F, col.names = F, quote = FALSE)

##************************************************************************
# 3. Generation of metadata table for each plasmid
#*************************************************************************

Plasmid_metadata <- data.frame(matrix(nrow=nrow(Abundance_table),ncol=0))
Plasmid_metadata$Plasmid_ID <- rownames(Abundance_table)

# Add information from geNomad output
geNomad_table <- read.delim("1_GENERAL_STATS/geNomad_plasmid_summary.tsv") 
geNomad_table <- geNomad_table[,c("seq_name", "length", "n_genes","topology", "n_hallmarks")]
colnames(geNomad_table)[1] <- "Plasmid_ID"
Plasmid_metadata <- left_join(Plasmid_metadata, geNomad_table, by = "Plasmid_ID")
Plasmid_metadata[Plasmid_metadata$Plasmid_ID=="IMG_circular_plasmid_795",2] <- 19770 # add length for missing plasmid in geNomad output

# Add information from Mobility analysis
mobility_result <- read.delim("3_MOBILITY//CS_Plasmid_Mobility.txt")
colnames(mobility_result)[c(3,4)] <- c("oriTs", "Mobility")
Plasmid_metadata <- left_join(Plasmid_metadata, mobility_result, by = "Plasmid_ID")

# Add information regarding the plasmid topology (circular, linear)
# Circular: plasmids generated with metaPlasmidSPAdes, SCAPP, complete IMG plasmids, and geNomad plasmids (>1kb, FDR < 0.05, DTR)
# Linear: plasmids identified by geNomad (>1kb, FDR < 0.05, no DTR, >=1 plasmid hallmark)
Plasmid_metadata$topology_simple <- NA
Plasmid_metadata <- Plasmid_metadata %>%
  mutate(topology_simple = case_when(
    grepl("circular", Plasmid_ID, ignore.case = TRUE) ~ "Circular",
    grepl("fragment", Plasmid_ID, ignore.case = TRUE) & topology == "DTR" ~ "Circular",
    is.na(topology_simple) ~ "Linear",
    TRUE ~ topology_simple
  ))

# Add CARD annotation results
# Generate a AMR count variable and a binary variable (Yes/No)
CARD_annotation <- read.delim("5_AMR_ANALYSIS/RGI_CARD_annotation.txt", sep = "\t", header = T, check.names = F)
CARD_annotation$Plasmid_ID <- sub(" #.*", "", CARD_annotation$ORF_ID)
CARD_annotation$Plasmid_ID <- sub("_[^_]*$", "", CARD_annotation$Plasmid_ID)
CARD_annotation$Contig <-NULL

AMR_gene_counts <- CARD_annotation %>%
  group_by(Plasmid_ID) %>%
  summarise(AMR_genes = n_distinct(Best_Hit_ARO))  

Plasmid_metadata <- Plasmid_metadata %>%
  left_join(AMR_gene_counts, by = "Plasmid_ID") %>%
  mutate(AMR_genes = ifelse(is.na(AMR_genes), 0, AMR_genes))  # Replace NA with 0 if there are no AMR genes

Plasmid_metadata <- Plasmid_metadata %>%
  mutate(AMR_genes_binary = ifelse(AMR_genes > 0, "Yes", "No"))

# Estimate the proportion of genes annotated as AMR
Plasmid_metadata$AMR_gene_proportion <- 100*(Plasmid_metadata$AMR_genes / Plasmid_metadata$n_genes)

# Generate a metadata table only for circular plasmids
Plasmid_metadata_circular <- Plasmid_metadata[Plasmid_metadata$topology_simple == "Circular",]
Plasmid_metadata_linear <- Plasmid_metadata[Plasmid_metadata$topology_simple == "Linear",]


##************************************************************************
# 4. Statistical comparisons
#*************************************************************************
Plasmid_metadata$topology_simple <- factor(Plasmid_metadata$topology_simple, levels = c("Circular", "Linear"))
Plasmid_metadata$Mobility <- factor(Plasmid_metadata$Mobility, levels = c("Conjugative", "Mobilizable", "Non-mobilizable"))

# Estimate mean abundance and prevalence and add to Plasmid_metadata
Plasmid_metadata$Mean_abundance <- rowMeans(Abundance_table)
Plasmid_metadata$Prevalence <- 100*(rowSums(Abundance_table >0) / ncol(Abundance_table))

# Generate a metadata table only for plasmids present in infants
Plasmid_metadata_infants <- Plasmid_metadata[Plasmid_metadata$Plasmid_ID %in% Present_plasmids_infants,]

#A) Compare the length between circular and fragment PTUs
wilcox_test_length <- wilcox.test(length ~ topology_simple, 
                                     data = Plasmid_metadata, 
                                     alternative = "greater")


#B) Check if Topology/Mobility are associated with prevalence / mean abundance
wilcox_test_topology_abundance <- wilcox.test(Mean_abundance ~ topology_simple, 
                                              data = Plasmid_metadata, 
                                              alternative = "greater")

wilcox_test_topology_prevalence <- wilcox.test(Prevalence ~ topology_simple, 
                                              data = Plasmid_metadata, 
                                              alternative = "greater")

test_mobility_abundance <- Kruskal_wallis_dunn_test(Plasmid_metadata$Mean_abundance, Plasmid_metadata$Mobility)
test_mobility_prevalence <- Kruskal_wallis_dunn_test(Plasmid_metadata$Prevalence, Plasmid_metadata$Mobility)


##************************************************************************
# 5. Calculation of summary statistics of plasmids
#*************************************************************************

#Estimate summary stats for all plasmids
table(Plasmid_metadata$topology_simple)  #Topology (circular or linear)
table(Plasmid_metadata$Mobility)  #Mobility
summary_stats(Plasmid_metadata$length) #Summary stats of the length of the genomes
summary_stats(Plasmid_metadata$n_genes[!is.na(Plasmid_metadata$n_genes)]) # Summary stats of number of predicted genes 
summary_stats(Plasmid_metadata$n_hallmarks[!is.na(Plasmid_metadata$n_hallmarks)]) #Summary stats of the number of plasmid hallmark genes

#Estimate summary stats for plasmids present in infants
table(Plasmid_metadata_infants$topology_simple)  #Topology (circular or linear)
table(Plasmid_metadata_infants$Mobility)  #Mobility
summary_stats(Plasmid_metadata_infants$length) #Summary stats of the length of the genomes
summary_stats(Plasmid_metadata_infants$n_hallmarks) #Summary stats of the number of plasmid hallmark genes
length(which(Plasmid_metadata_infants$AMR_genes > 0)) # Count number of plasmids with AMR genes
sum(Plasmid_metadata_infants$AMR_genes) # Count number of AMR genes in infant plasmids
100*(sum(Plasmid_metadata_infants$AMR_genes, na.rm = T)/sum(Plasmid_metadata_infants$n_genes, na.rm = T)) # % of genes annotated as AMR

#Estimate summary stats for circular plasmids 
table(Plasmid_metadata_circular$topology_simple)  #Topology 
table(Plasmid_metadata_circular$Mobility)  #Mobility
summary_stats(Plasmid_metadata_circular$length) #Summary stats of the length of the genomes
length(which(Plasmid_metadata_circular$AMR_genes > 0)) # Count number of plasmids with AMR genes
sum(Plasmid_metadata_circular$AMR_genes) # Count number of AMR genes in infant circular plasmids

#Estimate summary stats for linear plasmids 
table(Plasmid_metadata_linear$topology_simple)  #Topology
table(Plasmid_metadata_linear$Mobility)  #Mobility
summary_stats(Plasmid_metadata_linear$length) #Summary stats of the length of the genomes
sum(Plasmid_metadata_linear$AMR_genes) # Count number of AMR genes in linear infant plasmids

##************************************************************************
# 6. Generation of plots
#*************************************************************************

#############################
# Circularity of plasmids and mobility
#############################

Circularity_prop <- c(Circular = 37,  Linear = 63)
Mobility_prop<- c(Non_mobilizable = 53,  Mobilizable = 40, Conjugative = 7)

Circularity_waffle <- waffle(Circularity_prop,legend_pos = "bottom", size = 1, )
Mobility_waffle <- waffle(Mobility_prop,legend_pos = "bottom", size = 1, colors = c("#FFEDA0", "#FEB24C", "#FC4E2A"))

pdf('1_GENERAL_STATS/PLOTS/Stats_plasmids.pdf', width=3, height=5)
iron (Circularity_waffle,Mobility_waffle)
dev.off()

#############################
# Association of topology with length and mean abundance
#############################

pdf('1_GENERAL_STATS/PLOTS/Topology_length.pdf', width=2.6, height=3.2)
ggplot(Plasmid_metadata[!is.na(Plasmid_metadata$topology_simple), ], aes(x = topology_simple, y = log(length))) +
  geom_violin(aes(fill=topology_simple), alpha=0.7, width=0.8, trim = F) +
  geom_jitter(aes(color=topology_simple), alpha=0.1) + 
  geom_boxplot(width=0.2, fill="white", color="black", outliers = FALSE) + 
  labs(x = 'Plasmid circularity', y = "Genome Length (log10(bp))") + 
  scale_fill_manual(values = c("#6BBFA4", "#5FA4DA")) + 
  scale_color_manual(values = c("#6BBFA4", "#5FA4DA")) + 
  theme_classic() +
  theme(plot.title = element_text(size = 18, hjust = 0.5, face = "bold"), 
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 15), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.position = "none")
dev.off()

pdf('1_GENERAL_STATS/PLOTS/Topology_abundance.pdf', width=2.6, height=3.2)
ggplot(Plasmid_metadata[!is.na(Plasmid_metadata$topology_simple), ], aes(x = topology_simple, y = log(Mean_abundance))) +
  geom_violin(aes(fill=topology_simple), alpha=0.7, width=0.8, trim = F) +
  geom_jitter(aes(color=topology_simple), alpha=0.1) + 
  geom_boxplot(width=0.2, fill="white", color="black", outliers = FALSE) + 
  labs(x = 'Plasmid circularity', y = "log10(Mean Abundance)") + 
  scale_fill_manual(values = c("#6BBFA4", "#5FA4DA")) + 
  scale_color_manual(values = c("#6BBFA4", "#5FA4DA")) + 
  theme_classic() +
  theme(plot.title = element_text(size = 18, hjust = 0.5, face = "bold"), 
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 15), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.position = "none")
dev.off()

#############################
# Association of mobility with prevalence and mean abundance
#############################

pdf('1_GENERAL_STATS/PLOTS/Mobility_abundance.pdf', width=3.2, height=3.2)
ggplot(Plasmid_metadata[!is.na(Plasmid_metadata$Mobility), ], aes(x = Mobility, y = log(Mean_abundance))) +
  geom_violin(aes(fill=Mobility), alpha=0.7, width=0.8, trim = F) +
  geom_jitter(aes(color=Mobility), alpha=0.1) + 
  geom_boxplot(width=0.2, fill="white", color="black", outliers = FALSE) + 
  labs(x = 'Plasmid mobility', y = "log10(Mean Abundance)") + 
  scale_fill_manual(values = c("#EA5133", "#F9B14F", "#FFEDA0")) + 
  scale_color_manual(values = c("#EA5133", "#F9B14F", "#FFEDA0")) + 
  theme_classic() +
  theme(plot.title = element_text(size = 18, hjust = 0.5, face = "bold"), 
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 15), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.position = "none")
dev.off()

pdf('1_GENERAL_STATS/PLOTS/Mobility_prevalence.pdf', width=3.2, height=3.2)
ggplot(Plasmid_metadata[!is.na(Plasmid_metadata$Mobility), ], aes(x = Mobility, y = Prevalence)) +
  geom_violin(aes(fill=Mobility), alpha=0.7, width=0.8, trim = F) +
  geom_jitter(aes(color=Mobility), alpha=0.1) + 
  geom_boxplot(width=0.2, fill="white", color="black", outliers = FALSE) + 
  labs(x = 'Plasmid mobility', y = "Prevalence") + 
  scale_fill_manual(values = c("#EA5133", "#F9B14F", "#FFEDA0")) + 
  scale_color_manual(values = c("#EA5133", "#F9B14F", "#FFEDA0")) + 
  theme_classic() +
  theme(plot.title = element_text(size = 18, hjust = 0.5, face = "bold"), 
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 15), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.position = "none")
dev.off()

#****************
# Write results
#****************
write.table(Plasmid_metadata,"METADATA_TABLES/CS_Baby_Biome_Plasmid_Metadata_25062024.txt", sep = "\t", row.names = F, quote = FALSE) # Plasmid metadata
write.table(Plasmid_metadata_infants,"METADATA_TABLES/CS_Baby_Biome_Plasmid_Metadata_Infants_25062024.txt", sep = "\t", row.names = F, quote = FALSE) # Plasmid metadata
write.table(Sample_metadata, "METADATA_TABLES/CS_Baby_Biome_Sample_Metadata_25062024.txt", sep = "\t", row.names = F, quote = FALSE)
write.table(Sample_metadata_infants, "METADATA_TABLES/CS_Baby_Biome_Sample_Metadata_Infants_25062024.txt", sep = "\t", row.names = F, quote = FALSE)
