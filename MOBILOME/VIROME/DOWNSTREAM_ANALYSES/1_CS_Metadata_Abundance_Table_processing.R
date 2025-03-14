################################################################################
##### CS Baby Biome: Metadata and abundance table processing
### Author(s): Asier Fernández-Pato
### Last updated: 19th July, 2024
################################################################################

# Set working directory 
setwd("~/Desktop/PhD/Projects/CS Baby Biome/Final_analysis/VIRUSES/")

#****************
# Load libraries
#****************
library(dplyr)
library(data.table)
library(Biostrings)
library(UpSetR)
library(ggplot2)

#****************
# Define functions
#****************
# Function to count the number of viruses with each completeness group in a row
count_viruses <- function(row) {
  viruses <- row[!is.na(row)]
  complete_quality_count <- sum(viruses %in% CheckV_quality_all$virus_ID[CheckV_quality_all$quality == "Complete"])
  high_quality_count <- sum(viruses %in% CheckV_quality_all$virus_ID[CheckV_quality_all$quality == "High-quality"])
  medium_quality_count <- sum(viruses %in% CheckV_quality_all$virus_ID[CheckV_quality_all$quality == "Medium-quality"])
  return(data.frame(complete_quality_count, high_quality_count, medium_quality_count))
}

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

# Function to estimate the origin of genomes from each vOTU at a given taxonomic level (genus or family)
estimate_origin <- function(tax_level, clusters, databases) {
  results <- data.frame(matrix(ncol=1, nrow=(nrow(clusters))))
  colnames(results) <- paste0("vOTU_", toupper(tax_level), "_composition_extended")
  
  for (database in databases) {
    match_rows <- apply(clusters, 1, function(row) any(grepl(database, row, ignore.case = TRUE)))
    results[, paste0("vOTU_", toupper(tax_level), "_composition_extended")][match_rows] <-
      paste0(results[, paste0("vOTU_", toupper(tax_level), "_composition_extended")][match_rows], "&", database)
  }
  results[, paste0("vOTU_", toupper(tax_level), "_composition_extended")] <-
    substr(results[, paste0("vOTU_", toupper(tax_level), "_composition_extended")], 4, nchar(
      results[, paste0("vOTU_", toupper(tax_level), "_composition_extended")]
    ))
  # Add vOTU representative
  results$Virus_ID <- clusters$V1
  
  return(results)
}

# Function to find the cluster representative for a given Virus_ID
find_cluster <- function(virus_id, cluster_matrix) {
  row_index <- which(cluster_matrix == virus_id, arr.ind = TRUE)
  if (nrow(row_index) > 0) {
    return(cluster_matrix[row_index[1, 1], 1])  
  } else {
    return(NA)  
  }
}

# Function to generate combination of DBs for Upset plot
generate_DB_combinations <- function(databases, prefix = "", index = 1) {
  combinations <- c()
  if (index == length(databases)) {
    return(paste0(prefix, databases[index]))
  }
  combinations <- c(combinations, paste0(prefix, databases[index]))
  for (i in (index + 1):length(databases)) {
    sub_combinations <- generate_combinations(databases, paste0(prefix, databases[index], "&"), i)
    combinations <- c(combinations, sub_combinations)
  }
  return(combinations)
}


##************************************************************************
# 1. Load metadata table for the 195 CS Baby Biome samples 
#*************************************************************************
# Load the metadata table 
Sample_metadata <- read.delim("Metadata_CS/Metadata_EGA_CS_BABY_BIOME.txt") 

# Update clean read numbers
Sample_metadata$NG_ID <- Sample_metadata$bioSampleId
Clean_reads <- read.delim("Metadata_CS/CS_Baby_Biome_nreads_195.txt", header = T) 
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
Abundance_table <- read.delim("Abundance_table/CS_Baby_Biome_Abundance_Table_RPKM.txt") 

# Select only the samples from the metadata
#Abundance_table <- Abundance_table[, colnames(Abundance_table) %in% Sample_metadata$NG_ID]

# Remove from the table those viruses not present in any sample
Absent_viruses <- rownames(Abundance_table[rowSums(Abundance_table)==0, ]) #100,250
Abundance_table <- Abundance_table[!rownames(Abundance_table) %in% Absent_viruses,] #2,263
Present_viruses <- rownames(Abundance_table) #2,263

# Reorder abundance table to match Final metadata
Abundance_table <- Abundance_table[, Sample_metadata$NG_ID]

# Repeat processing for infant samples
Abundance_table_infants <- Abundance_table[, colnames(Abundance_table) %in% Sample_metadata_infants$NG_ID]
Absent_viruses_infants <- rownames(Abundance_table_infants[rowSums(Abundance_table_infants)==0, ]) #1,107
Abundance_table_infants <- Abundance_table_infants[!rownames(Abundance_table_infants) %in% Absent_viruses_infants,] #1,156
Present_viruses_infants <- rownames(Abundance_table_infants) #1,156
Abundance_table_infants <- Abundance_table_infants[, Sample_metadata_infants$NG_ID]

#****************
# Write results
#****************
# Abundance table and list of present and excluded viruses (with no presence in CS samples)
write.table(Abundance_table,"Abundance_table/CS_Abundance_Table_17052024.txt", sep = "\t", row.names = T, quote = FALSE)
write.table(Absent_viruses,"Abundance_table/Viruses_with_no_presence.txt", sep = "\n", row.names = F, col.names = F, quote = FALSE)
write.table(Present_viruses,"Abundance_table/Present_viruses.txt", sep = "\n", row.names = F, col.names = F, quote = FALSE)
write.table(Abundance_table_infants,"Abundance_table/CS_Abundance_Table_INFANTS_17052024.txt", sep = "\t", row.names = T, quote = FALSE)
write.table(Absent_viruses_infants,"Abundance_table/Viruses_with_no_presence_INFANTS.txt", sep = "\n", row.names = F, col.names = F, quote = FALSE)
write.table(Present_viruses_infants,"Abundance_table/Present_viruses_INFANTS.txt", sep = "\n", row.names = F, col.names = F, quote = FALSE)

##************************************************************************
# 3. Generation of metadata table for each vOTU representative
#*************************************************************************
# From here on, decide whether to use all vOTUs or only infant vOTUs.

#**********************************
#A. Add DB of origin and virus IDs
#**********************************
Virus_metadata <- read.delim("Metadata_CS/Dereplication_input_sequences_nodup_DB_origin.txt", header=F, col.names = c("Virus_ID_original", "DB", "Virus_ID"))
Virus_metadata <- Virus_metadata %>%
  filter(Virus_ID %in% Present_viruses)
rownames(Virus_metadata) <- Virus_metadata$Virus_ID

#**********************************
#B. Add sequence length and GC content
#**********************************
# Read file and merge it with viral metadata
Length_GC <- read.delim("1_GENERAL_STATISTICS_AND_METADATA/All viruses/Length_GC_rep_seqs.txt", header=F, col.names = c("Virus_ID", "Length", "GC")) 
Virus_metadata <- left_join(Virus_metadata, Length_GC, by = "Virus_ID")

#**********************************
#C. Add CheckV quality
#**********************************
# Estimate the CheckV quality for all viral genomes used as input for dereplication. Reason for this:
#- Quality/Completeness is missing for viruses from external DBs
# Load file and add quality information to non-CS viruses
CheckV_quality_all <- read.delim("1_GENERAL_STATISTICS_AND_METADATA/All viruses/quality_summary.tsv", header=T)
CheckV_quality_all <- CheckV_quality_all[,c("contig_id", "viral_genes", "checkv_quality", "completeness_method")]
colnames(CheckV_quality_all) <- c("Virus_ID", "Viral_genes", "Quality", "Completeness_method")

Virus_metadata <- left_join(Virus_metadata, CheckV_quality_all, by = "Virus_ID")

#**********************************
#D. Add number of genes predicted by prodigal-gv
#**********************************\

vOTU_sequences <- readAAStringSet("1_GENERAL_STATISTICS_AND_METADATA/All viruses/viruses_present_proteins.faa")
Virus_ids <- character()

for (header in names(vOTU_sequences)) {
  Virus_id <- gsub(" #.*", "", header) 
  Virus_id <- gsub("_[^_]*$", "", Virus_id, perl = TRUE)  
  Virus_ids <- c(Virus_ids, Virus_id)
}

Virus_proteins_counts <- table(Virus_ids)
Virus_proteins <- data.frame("Virus_ID"=names(Virus_proteins_counts), "n_genes"= as.numeric(Virus_proteins_counts))
Virus_metadata <- left_join(Virus_metadata, Virus_proteins, by = "Virus_ID")

#**********************************
#E. Add the number of viral genomes in the vOTU
#**********************************
# Read file with vOTU clusters
vOTU_clusters <- read.delim("1_GENERAL_STATISTICS_AND_METADATA/All viruses/viral_clusters.tsv", header = F)
vOTU_clusters_present <- vOTU_clusters[vOTU_clusters$V1 %in% Present_viruses,]
vOTU_clusters_present_processed <-  data.frame(vOTU_clusters_present[,-2], 
                                               tstrsplit(as.character(vOTU_clusters_present[,2]), ",", fixed = TRUE))
colnames(vOTU_clusters_present_processed) [1] <- c("rep_seq")

# Get number of viral genomes per vOTU
Number_viruses_vOTUs <- data.frame("Virus_ID"= vOTU_clusters_present_processed$rep_seq,
                                   "Number_genomes_vOTU" = rowSums(!is.na(vOTU_clusters_present_processed)) -1)
  
# Merge it with viral metadata
Virus_metadata <- left_join(Virus_metadata, Number_viruses_vOTUs, by = "Virus_ID")

# Add also a factor variable 
Virus_metadata <- Virus_metadata %>%
  mutate(Number_genomes_vOTU_factor = cut(Number_genomes_vOTU,
                                          breaks = c(0, 2, 10, 100, Inf),
                                          labels = c("Singleton", "< 10 genomes", "< 100 genomes", ">= 100 genomes"),
                                          include.lowest = TRUE,
                                          right = FALSE))
Virus_metadata$Number_genomes_vOTU_factor <- factor(Virus_metadata$Number_genomes_vOTU_factor)

#**********************************
#E. Add the DB composition of the vOTUs
#**********************************

#############################
#F.1. DB composition of vOTUs 
#############################
databases <- c("CS","GPD", "MGV", "IMG_VR", "ELGV","Shah","Benler","RefSeq", "Gulyaeva", "Guerin","Yutin")

# Estimate the origin of genomes from each vOTU
results <- data.frame(matrix(ncol=1, nrow=(nrow(vOTU_clusters_present_processed))))
colnames(results) <- "vOTU_DB_composition_extended"

for (database in databases) {
  match_rows <- apply(vOTU_clusters_present_processed, 1, function(row) any(grepl(database, row, ignore.case = TRUE)))
  results$vOTU_DB_composition_extended[match_rows] <- paste0(results$vOTU_DB_composition_extended[match_rows], "&", database)
}
results$vOTU_DB_composition_extended <- substr(results$vOTU_DB_composition_extended, 4, nchar(results$vOTU_DB_composition_extended))
results$Virus_ID <- vOTU_clusters_present_processed$rep_seq

# Merge it with viral metadata
Virus_metadata <- left_join(Virus_metadata, results, by = "Virus_ID")


#############################
#F.2. vOTU exclusively composed of CS Baby Biome sequences
#############################
# Estimate which vOTUs have viral genomes that do not match "CS" pattern (from other DBs)
vOTUs_with_external_genomes <- apply(vOTU_clusters_present_processed, 1, function(row) {
  any(!is.na(row) & !grepl("CS", row, ignore.case = TRUE))
})

# Summarize results
vOTU_composition <- data.frame("Virus_ID" = vOTU_clusters_present_processed$rep_seq,
                                  "vOTU_DB_composition" = ifelse(vOTUs_with_external_genomes, "OTHER_DBs", "CS"))

# Merge it with viral metadata
Virus_metadata <- left_join(Virus_metadata, vOTU_composition, by = "Virus_ID")

#**********************************
#G. Taxonomic assignment
#**********************************
###################
# geNomad taxonomy
###################
# Load geNomad taxonomic assignment
taxonomy_viruses <- read.delim("4_TAXONOMY/viruses_present_taxonomy.tsv", 
                               sep = "\t", header = T, col.names = c("Virus_ID", "Number_genes_with_taxonomy", 
                                                                     "Agreement", "Taxid", "Lineage"))                           
taxonomy_viruses <- taxonomy_viruses[,c("Virus_ID", "Lineage")]

# Split viral taxonomy
tax_levels <- strsplit(as.character(taxonomy_viruses$Lineage), ";")

# Determine the maximum number of levels in the taxonomy column and fill empty values with NA
max_levels <- max(sapply(tax_levels, length))
tax_levels <- lapply(tax_levels, function(x) c(x, rep(NA, max_levels - length(x))))

# Bind the taxonomy vectors to the original data frame
taxonomy_viruses  <- cbind(taxonomy_viruses, do.call("rbind", tax_levels))
taxonomy_viruses$Lineage <- NULL

# Rename the new columns
colnames(taxonomy_viruses)[2:8] <- c("Superkingdom", "Other", "Kingdom", "Phylum", "Class", "Order", "Family")

# As there are some viruses without order assigned, some family-level taxonomy is assigned as class
# Move this taxonomic level to the correct column
assign_to_family <- grep("viridae", taxonomy_viruses$Order, value=T)
taxonomy_viruses$Family[grep("viridae", taxonomy_viruses$Order)] <- assign_to_family
taxonomy_viruses$Order[grep("viridae", taxonomy_viruses$Order)] <- NA

# Merge it with viral metadata
Virus_metadata <- left_join(Virus_metadata, taxonomy_viruses, by = "Virus_ID")

#**********************************
#H. Lifestyle assignment
#**********************************
# Read BACPHLIP lifestyle prediction output
lifestyle_viruses <- read.delim("3_LIFESTYLE/viruses_present.fa.bacphlip", 
                                col.names = c("Virus_ID", "Virulent", "Temperate"))

# Remove predictions for eukaryotic viruses: Bamfordvirae and Shotokuvirae kingdoms
eukaryotic_viruses <- Virus_metadata[Virus_metadata$Kingdom %in%
                                                c("Bamfordvirae", "Shotokuvirae"),
                                     "Virus_ID"]
lifestyle_viruses <- lifestyle_viruses[!lifestyle_viruses$Virus_ID %in% eukaryotic_viruses,]

# Generate final lifestyle prediction:
## A) Viruses with Temperate > 0.5 --> Temperate
## B) Viruses with Virulent > 0.5:
## B.1) If Completeness > 90% (CheckV) --> Virulent
## B.2) If Completeness < 90% (CheckV) --> NA
quality <- Virus_metadata[,c("Virus_ID", "Quality")]
quality_lifestyle <- inner_join(lifestyle_viruses, quality, by = "Virus_ID")

quality_lifestyle <- quality_lifestyle %>%
  mutate(Lifestyle = case_when(
    Temperate > 0.5 ~ "Temperate",
    Virulent > 0.5 & (Quality == "High-quality" | Quality == "Complete") ~ "Virulent",
  )) 
quality_lifestyle <- quality_lifestyle[, c("Virus_ID", "Lifestyle")]

# Merge it with viral metadata
Virus_metadata <- left_join(Virus_metadata, quality_lifestyle, by = "Virus_ID")

#**********************************
#I. DGR identification
#**********************************
# Read DGR identification results
DGR_results <- read.delim("DGR_ANALYSIS/All_RT_1k_DGR_detection_filtered_075A.tsv")

Virus_metadata$DGR <- ifelse(Virus_metadata$Virus_ID %in% DGR_results$Contig, "Yes", "No")

#**********************************
#J. Genus and family-level clustering
#**********************************

# Load the genus_clusters.txt and family_clusters.txt files
genus_clusters <- read.delim("1_GENERAL_STATISTICS_AND_METADATA/All viruses/genus_clusters.txt", header = F)
family_clusters <- read.delim("1_GENERAL_STATISTICS_AND_METADATA/All viruses/family_clusters.txt", header = F)

# A) Add genus and family-level cluster information to metadata
family_clusters_matrix <- as.matrix(family_clusters)
genus_clusters_matrix <- as.matrix(genus_clusters)

Virus_metadata$Family_cluster <- sapply(Virus_metadata$Virus_ID, function(id) find_cluster(id, family_clusters_matrix))
Virus_metadata$Genus_cluster <- sapply(Virus_metadata$Virus_ID, function(id) find_cluster(id, genus_clusters_matrix))

# B) Estimate the origin of genomes from each vOTU at genus and family-level (vOTU_DB_composition extended)
results_genus_level <- estimate_origin("genus", genus_clusters, databases)
results_family_level <- estimate_origin("family", family_clusters, databases)

# Add a variable vOTU_DB_composition indicating vOTUs specific to CS
vOTUs_genus_with_external_genomes <- apply(genus_clusters, 1, function(row) {
  any(nchar(row) > 0 & !grepl("CS", row, ignore.case = TRUE))
})
vOTUs_family_with_external_genomes <- apply(family_clusters, 1, function(row) {
  any(nchar(row) > 0 & !grepl("CS", row, ignore.case = TRUE))
})

results_genus_level$vOTU_composition_genus <- ifelse(vOTUs_genus_with_external_genomes, "OTHER_DBs", "CS")
results_family_level$vOTU_composition_family <- ifelse(vOTUs_family_with_external_genomes, "OTHER_DBs", "CS")

# Add the number of viral genomes per vOTU at each taxonomic level
results_genus_level$Number_genomes_vOTU <- apply(genus_clusters, 1, function(row) {
  sum(nchar(row) > 0)
})
results_family_level$Number_genomes_vOTU <- apply(family_clusters, 1, function(row) {
  sum(nchar(row) > 0)
})


##************************************************************************
# 4. Calculation of summary statistics of vOTUs
#*************************************************************************

#Estimate summary stats for all vOTUs
table(Virus_metadata$DB) #DB distribution of vOTU representatives
table(Virus_metadata$Quality)  #Quality distribution of vOTU representatives
summary_stats(Virus_metadata$Length) #Summary stats of the length of all vOTU representatives
summary_stats(Virus_metadata$GC) #Summary stats of the %GC of all vOTU representatives
summary_stats(Virus_metadata$Viral_genes) # Summary stats of number of predicted viral genes 
summary_stats(Virus_metadata$Number_genomes_vOTU) #Summary stats of the number of genomes per vOTU

# Estimate the summary stats for vOTUs with a CS genome as representative 
table(Virus_metadata$Quality[Virus_metadata$DB=="CS"])
summary_stats(Virus_metadata$Length[Virus_metadata$DB=="CS"])
summary_stats(Virus_metadata$GC[Virus_metadata$DB=="CS"])
summary_stats(Virus_metadata$Number_genomes_vOTU[Virus_metadata$DB=="CS"])

# Estimate the summary stats for vOTUs with only CS genomes 
table(Virus_metadata$Quality[Virus_metadata$vOTU_DB_composition=="CS"])
summary_stats(Virus_metadata$Length[Virus_metadata$vOTU_DB_composition=="CS"])
summary_stats(Virus_metadata$GC[Virus_metadata$vOTU_DB_composition=="CS"])
summary_stats(Virus_metadata$Number_genomes_vOTU[Virus_metadata$vOTU_DB_composition=="CS"])

# Estimate the number of genus and family-level vOTUs with only CS genomes
# (vOTUs with only CS representatives at genus or family level and only CS genomes at species level)
length(which(results_genus_level$Virus_ID[results_genus_level$vOTU_composition_genus == "CS"] %in%
        Virus_metadata$Virus_ID[Virus_metadata$vOTU_DB_composition == "CS"]))
length(which(results_family_level$Virus_ID[results_family_level$vOTU_composition_family == "CS"] %in%
               Virus_metadata$Virus_ID[Virus_metadata$vOTU_DB_composition == "CS"]))

##************************************************************************
# 5. Generation of summary plots
#*************************************************************************

#############################
# Origin of vOTU representatives
#############################

# Generate UpSet plot at species, genus and family levels
DB_distribution_species <- table(Virus_metadata$vOTU_DB_composition_extended)
DB_distribution_genus <- table(results_genus_level$vOTU_GENUS_composition_extended)
DB_distribution_family <- table(results_family_level$vOTU_FAMILY_composition_extended)

pdf('1_GENERAL_STATISTICS_AND_METADATA/Plots/Upset_DBs_species_vOTUs.pdf', width=23, height=11)
upset(fromExpression(DB_distribution_species), 
      nintersects = 30, 
      nsets = 9, 
      order.by = "freq", 
      decreasing = T, 
      mb.ratio = c(0.45, 0.55),
      text.scale = 3.9, 
      point.size = 5, 
      line.size = 1, 
      sets.bar.color = "darkblue",
      shade.color = "grey",
      shade.alpha = 0.15, show.numbers = F,
      sets.x.label = "Number vOTUs",
)
dev.off()

# Generate a barplot with DB distribution of representatives
DB_distribution <- table(Virus_metadata$DB)
DB_distribution <- data.frame(
  DB = names(DB_distribution),
  Count = as.numeric(DB_distribution) 
)
DB_distribution$Description <- rep ("vOTU database distribution", length(DB_distribution$DB))
DB_distribution <- DB_distribution %>%
  mutate(Proportion = 100*(Count / sum(Count)))
DB_distribution$DB <- c("Benler et al", "CS", "ELGV","GPD", "Gulyaeva et al",
                        "IMG_VR","MGV", "RefSeq", "Shah et al")

colors <- c("#8DD3C7", "#FFDAB9", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", 
            "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

pdf('1_GENERAL_STATISTICS_AND_METADATA/Plots/DB_distribution.pdf', width=6.3, height=2)
ggplot(DB_distribution, aes(x = Description, y = Proportion, fill = DB)) +
  geom_col() +  
  coord_flip() +
  scale_fill_manual(values = colors) + 
  theme_classic() +
  theme(axis.text=element_text(size=14),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 13),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "top")
dev.off()

#############################
# Length/GC%/Viral genes distribution
#############################
# Generate histogram for length distribution of vOTU representatives
pdf('1_GENERAL_STATISTICS_AND_METADATA/Plots/Length_distribution_log.pdf', width=6, height=1.5)
ggplot(Virus_metadata, aes(x = Length)) +
  geom_histogram(bins = 40, fill = alpha("steelblue", 0.7)) +
  labs(x = "Length (kbp)", y = "Frequency") +  # Modify the axis labels as needed
  theme_classic() +
  scale_x_log10(labels = scales::number_format(scale = 0.001)) +  
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 15),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 13),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    legend.position = "top"
  )
dev.off()

# Generate histogram for %GC distribution of vOTU representatives
pdf('1_GENERAL_STATISTICS_AND_METADATA/Plots/GC_distribution.pdf', width=6, height=1.5)
ggplot(Virus_metadata, aes(x = GC)) +
  geom_histogram(bins = 40, fill = alpha("coral", 0.7)) +
  labs(x = "GC content (%)", y = "Frequency") +  # Modify the axis labels as needed
  theme_classic() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 15),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 13),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    legend.position = "top"
  )
dev.off()

# Generate histogram for viral genes distribution of vOTU representatives
pdf('1_GENERAL_STATISTICS_AND_METADATA/Plots/Viral genes_distribution.pdf', width=6, height=1.5)
ggplot(Virus_metadata, aes(x = (Viral_genes))) +
  geom_histogram(bins = 40, fill = alpha("plum", 0.7)) +
  labs(x = "Number of predicted viral genes", y = "Frequency") +  # Modify the axis labels as needed
  theme_classic() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 15),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 13),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    legend.position = "top"
  )
dev.off()

#############################
# Distribution of viral taxonomy and lifestyle
#############################
# Generate a barplot with the viral taxonomy at class level

# Change all values other than "Caudoviricetes" to "Other"
Viral_class <- Virus_metadata$Class[!is.na(Virus_metadata$Class)]
Viral_class[Viral_class != "Caudoviricetes"] <- "Other"
Viral_class_distribution <- table(Viral_class)

# Generate a dataframe with the counts and proportions
Viral_class_distribution <- data.frame(
  Class = names(Viral_class_distribution),
  Count = as.numeric(Viral_class_distribution) 
)
Viral_class_distribution$Description <- rep ("Viral class distribution", length(Viral_class_distribution$Class))
Viral_class_distribution <- Viral_class_distribution %>%
  mutate(Proportion = 100*(Count / sum(Count)))


pdf('1_GENERAL_STATISTICS_AND_METADATA/Plots/All/Viral_classes_barplot.pdf', width=4, height=1.7)
ggplot(Viral_class_distribution, aes(x = Description, y = Proportion, fill = Class)) +
  geom_col() +  
  coord_flip() +
  scale_fill_manual(values = c("#ADD8E6", "#E6E6FA")) + 
  theme_classic() +
  theme(axis.text=element_text(size=13),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 13.5),
        legend.text = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "top")
dev.off()

# Generate a barplot with the predicted viral lyfestyle
Viral_lifestyle_distribution <- table(Virus_metadata$Lifestyle, useNA = "ifany")

# Generate a dataframe with the counts and proportions
Viral_lifestyle_distribution <- data.frame(
  Lifestyle = names(Viral_lifestyle_distribution),
  Count = as.numeric(Viral_lifestyle_distribution) 
)
Viral_lifestyle_distribution$Description <- rep ("Viral lifestyle distribution", length(Viral_lifestyle_distribution$Lifestyle))
Viral_lifestyle_distribution <- Viral_lifestyle_distribution %>%
  mutate(Proportion = 100*(Count / sum(Count)))


pdf('1_GENERAL_STATISTICS_AND_METADATA/Plots/All/Viral_lifestyle_barplot.pdf', width=4, height=1.7)
ggplot(Viral_lifestyle_distribution, aes(x = Description, y = Proportion, fill = Lifestyle)) +
  geom_col() +  
  coord_flip() +
  scale_fill_manual(values = c("#E0C9A8", "#C8D9BF")) + 
  theme_classic() +
  theme(axis.text=element_text(size=13),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 13.5),
        legend.text = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "top")
dev.off()


#############################
# Distribution of number of viral genomes
#############################
# Generate pie chart for the distrubution of genomes in vOTUs
palette <- c(alpha("steelblue", 0.7), alpha("palegreen3", 0.7), 
             alpha("coral", 0.7), alpha("orchid", 0.7))

pdf('1_GENERAL_STATISTICS_AND_METADATA/Plots/Number_genomes_piechart.pdf', width=5, height=2.5)
ggplot(Virus_metadata, aes(x = "", fill = Number_genomes_vOTU_factor)) +
  geom_bar(size=0.6) +
  coord_polar(theta = "y") +
  labs(x = "Number of genomes vOTU", y = NULL, fill = "vOTU genome distribution") +
  theme_void() +
  scale_fill_manual(values = palette) +
  theme(legend.title = element_text(size = 18),
    legend.text = element_text(size = 15),
    legend.position = "right",
    legend.key.height = unit(1.4, "lines")  
  )
dev.off()


  
#############################
# Completeness distribution
#############################

Virus_metadata$Quality <- factor(Virus_metadata$Quality, levels = c("Complete", "High-quality",
                                                                                "Medium-quality","Low-quality",
                                                                                "Not-determined"))
# Generate a pie chart
palette <- c("#D8BFD8", "#C0C0C0", "#B0E0E6", "#E6E6FA", "#F0E68C")

pdf('1_GENERAL_STATISTICS_AND_METADATA/Plots/Completeness_piechart.pdf', width=5, height=2.5)
ggplot(Virus_metadata, aes(x = "", fill = Quality)) +
  geom_bar(size=0.6) +
  coord_polar(theta = "y") +
  labs(x = "Number of genomes vOTU", y = NULL, fill = "Completeness distribution") +
  theme_void() +
  scale_fill_manual(values = palette, labels = c("Complete", "High quality",
                                                 "Medium quality","Low quality",
                                                 "Not determined")) +
  theme(legend.title = element_text(size = 18),
        legend.text = element_text(size = 15),
        legend.position = "right",
        legend.key.height = unit(1.4, "lines")  
  )
dev.off()

#****************
# Write results
#****************
# Sample metadata
write.table(Sample_metadata,"Metadata_CS/CS_Baby_Biome_Sample_Metadata_17052024.txt", sep = "\t", row.names = F, quote = FALSE) 
# Virus metadata
write.table(Virus_metadata,"Metadata_CS/CS_Baby_Biome_Viral_Metadata_17052024.txt", sep = "\t", row.names = F, quote = FALSE) 
