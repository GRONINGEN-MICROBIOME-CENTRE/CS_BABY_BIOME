################################################################################
##### CS Baby Biome: Host Assignment Analysis
### Author(s):Asier Fernández-Pato
### Last updated: 24th July, 2024
################################################################################

#****************
# Load modules
#****************
library(dplyr)
library(tidyverse)
library(grDevices)
library(phyloseq)
library(microbiome)
library(ggplot2)

# Set working directory 
setwd("~/Desktop/PhD/Projects/CS Baby Biome/Final_analysis/VIRUSES/")

# Read abundance table and metadata tables
Abundance_table <- read.delim("Abundance_table/CS_Abundance_Table_17052024.txt")
Sample_metadata <- read.delim("Metadata_CS/CS_Baby_Biome_Sample_Metadata_17052024.txt")
Virus_metadata <- read.delim("Metadata_CS/CS_Baby_Biome_Viral_Metadata_17052024.txt")

##************************************************************************
# 1. Host Assignment Analysis: iPHOP results preprocessing
##************************************************************************

#A) Read iPHOP prediction output at genus level with default DB
viral_hosts_df <- read.delim("6_HOST_ASSIGNMENT/Host_prediction_to_genus_m90_default_db.csv", 
                                sep = ",", header = T)

# Filter iPHOP predictions to get only the top hit per virus
viral_hosts_df <- viral_hosts_df %>%
  group_by(Virus) %>%
  dplyr::slice(which.max(Confidence.score))

viral_hosts_df <- viral_hosts_df[,c(1,3,4)]
colnames(viral_hosts_df) <- c("Virus_ID", "Bacterial_host", "Score") #predictions for 2092/2263 viruses (92.4%)

#B) Read iPHOP prediction output at genus level with enriched DB
viral_hosts_enr <- read.delim("6_HOST_ASSIGNMENT/Host_prediction_to_genus_m90_enriched_db.csv", 
                          sep = ",", header = T)

# Filter iPHOP predictions to get only the top hit per virus
viral_hosts_enr <- viral_hosts_enr %>%
  group_by(Virus) %>%
  dplyr::slice(which.max(Confidence.score))

viral_hosts_enr <- viral_hosts_enr[,c(1,3,4)]
colnames(viral_hosts_enr) <- c("Virus_ID", "Bacterial_host", "Score") #predictions for 2076/2263 viruses (91.7%)


# Select prediction with highest score for the previous 2 prediction approaches
# In total, we get confident predictions for 2132/2263 viruses (94.2%) (1952/2078 after removing M06-M12)
viral_hosts_combined <- rbind(viral_hosts_df, viral_hosts_enr)

viral_hosts_combined <- viral_hosts_combined %>%
  group_by(Virus_ID) %>%
  filter(Score == max(Score)) %>%
  dplyr::slice(1) %>%  # Keep the first occurrence
  ungroup()
  
viral_hosts_combined$Score <- NULL

# Merge with Virus metadata
Virus_metadata <- left_join(Virus_metadata, viral_hosts_combined, by = "Virus_ID")

# Extract family and genus-level hosts
Virus_metadata$Bacterial_family_host <- str_extract(Virus_metadata$Bacterial_host, "(?<=f__)[^;]+")
Virus_metadata$Bacterial_genus_host <- str_extract(Virus_metadata$Bacterial_host, "(?<=g__)[^;]+")
Virus_metadata$Bacterial_host <- NULL

# Set host columns as a factor (ordering them by occurence)
Virus_metadata$Bacterial_family_host <- factor(Virus_metadata$Bacterial_family_host)
Virus_metadata$Bacterial_genus_host <- factor(Virus_metadata$Bacterial_genus_host)


##################################
# Generation of Phyloseq object
##################################

# Split host taxonomy
tax_levels <- strsplit(as.character(viral_hosts_combined$Bacterial_host), ";")
tax_levels <- lapply(tax_levels, function(x) substring(x, 4))

# Determine the maximum number of levels in the taxonomy column and fill empty values with NA
max_levels <- max(sapply(tax_levels, length))
tax_levels <- lapply(tax_levels, function(x) c(x, rep(NA, max_levels - length(x))))

# Bind the taxonomy vectors to the original data frame
taxonomy_hosts  <- cbind(viral_hosts_combined, do.call("rbind", tax_levels))

# Rename the new columns
colnames(taxonomy_hosts)[3:8] <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
taxonomy_hosts$Bacterial_host <- NULL

# As some viruses do not have taxonomy, add them manually to the taxonomy table
# Create a vector of Virus_IDs with no match in taxonomy_viruses
taxonomy_hosts <- taxonomy_hosts[taxonomy_hosts$Virus_ID %in% rownames(Abundance_table), ] 
missing_IDs <- setdiff(rownames(Abundance_table), taxonomy_hosts$Virus_ID)

# Create a new dataframe with NA values for missing Virus_IDs
missing_hosts_viruses <- data.frame(Virus_ID = missing_IDs)
missing_hosts_viruses[, names(taxonomy_hosts)[2:7]] <- NA

# Combine taxonomy_viruses and missing_viruses
taxonomy_hosts <- rbind(taxonomy_hosts, missing_hosts_viruses)

# Reorder the rows based on Virus_ID matching row_names
taxonomy_hosts <- data.frame(taxonomy_hosts[match(rownames(Abundance_table), taxonomy_hosts$Virus_ID), ])

# Set virus IDs as rownames 
rownames(taxonomy_hosts) <- taxonomy_hosts$Virus_ID
taxonomy_hosts$Virus_ID <- NULL

#**Due to differences in taxonomy assignment of MAGs in enriched DB some viruses are assigned = "Genus" but != "Phylum" (changed in GTDB tax)
# I manually adapt taxonomy_hosts table it to make it uniform (and avoid problems for visualization with phyloseq)
taxonomy_hosts <- read.delim("6_HOST_ASSIGNMENT/taxonomy_hosts.tsv", row.names = 1)
colnames(taxonomy_hosts) <- c("Domain","Phylum","Class","Order","Family","Genus")
taxonomy_hosts <- taxonomy_hosts[match(rownames(Abundance_table), rownames(taxonomy_hosts)), ]

# Transform dataframes into tables (for Phyloseq)
Abundance_table_matrix <- as.matrix(Abundance_table)
taxonomy_hosts_matrix <- as.matrix(taxonomy_hosts)

# Generate phyloseq object
vOTU <- otu_table(Abundance_table_matrix, taxa_are_rows = TRUE)
TAX <- tax_table(taxonomy_hosts_matrix)
samples <- sample_data(Sample_metadata)
sample_names(samples) <- Sample_metadata$NG_ID

Phyloseq_hosts <- phyloseq(vOTU, TAX, samples)


##************************************************************************
# 2. Host Taxonomic Analysis in Infants
##************************************************************************

# Generate Virus metadata for those present in infant samples
Abundance_table_infants <- Abundance_table[,Sample_metadata$Type == "Infant"]
viruses_present_infants <- rownames(Abundance_table_infants)[rowSums(Abundance_table_infants)>0]
Virus_metadata_infants <- Virus_metadata[Virus_metadata$Virus_ID %in% viruses_present_infants, ]

# Generate a barplot of the host assignment of vOTUs at family level
Viral_family_host_counts_infants <- Virus_metadata_infants %>%
  count(Bacterial_family_host) %>%
  na.omit()

Viral_family_host_counts_infants <- head(arrange(Viral_family_host_counts_infants, desc(n)), 8) # select top 8

family_colors <-  c("#E6AB02", "#1B9E77", "#D95F02", "#8DD3C7", "#A6CEE3", "#F768A1", "#FDB462", "#B3DE69")
family_colors <- adjustcolor(family_colors, alpha.f = 0.6)
family_names <- unique(Viral_family_host_counts_infants$Bacterial_family_host)
color_mapping <- setNames(family_colors, family_names)

pdf('6_HOST_ASSIGNMENT/Plots/Family_host_barplot_infants.pdf', width=4.5, height=3.1)
ggplot(Viral_family_host_counts_infants, aes(x = reorder(Bacterial_family_host, -n), y = n)) +
  geom_bar(stat = "identity", color = NA, aes(fill = Bacterial_family_host)) +
  coord_flip() +
  ylab("Number of vOTUs") +
  scale_fill_manual(values = color_mapping) +
  theme_classic() +
  theme(axis.text = element_text(size = 14, colour = "black"),
        axis.title = element_text(size = 15),
        axis.text.y = element_text(face = "italic"),
        axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "none")
dev.off()

# Generate a barplot of the host assignment of vOTUs at genus level
Viral_genus_host_counts_infants <- Virus_metadata_infants %>%
  count(Bacterial_genus_host) %>%
  na.omit()

Viral_genus_host_counts_infants <- head(arrange(Viral_genus_host_counts_infants, desc(n)), 8) # select top 8

genus_colors <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5")
genus_colors <- adjustcolor(genus_colors, alpha.f = 0.6)
genus_names <- unique(Viral_genus_host_counts_infants$Bacterial_genus_host)
color_mapping <- setNames(genus_colors, genus_names)


pdf('6_HOST_ASSIGNMENT/Plots/Genus_host_barplot_infants.pdf', width=4.5, height=3.1)
ggplot(Viral_genus_host_counts_infants, aes(x = reorder(Bacterial_genus_host, -n), y = n)) +
  geom_bar(stat = "identity", color = NA, aes(fill = Bacterial_genus_host)) +
  coord_flip() +
  ylab("Number of vOTUs") +
  scale_fill_manual(values = color_mapping) +
  theme_classic() +
  theme(axis.text = element_text(size = 14, colour = "black"),
        axis.title = element_text(size = 15),
        axis.text.y = element_text(face = "italic"),
        axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "none") 
dev.off()

##################################
# Host abundance analysis
##################################

# Select only infant samples
Phyloseq_hosts_infants <- subset_samples(Phyloseq_hosts, Type =="Infant")

#********
# Family-level
#********
# Generate family-level abundances
# if unknown category wanted keep: subset_taxa(!is.na(Family))
# if other category wanted: aggregate_rare(level="Family", detection = 0.05, prevalence = 0.05) 
# Select top 10: 
#topf <- names(sort(taxa_sums(Phyloseq_hosts_infants_family), TRUE)[1:10])
#Phyloseq_hosts_infants_family <- prune_taxa(topf, Phyloseq_hosts_infants_family)
# Make abundances compositional: microbiome::transform(Phyloseq_hosts_infants_family, transform = "compositional")
Phyloseq_hosts_infants_family <- Phyloseq_hosts_infants %>% 
  subset_taxa(!is.na(Family)) %>% 
  tax_glom(taxrank = "Family")  

top_families <- names(sort(taxa_sums(Phyloseq_hosts_infants_family), TRUE)[1:12])
Phyloseq_hosts_infants_family_filt <- prune_taxa(top_families, Phyloseq_hosts_infants_family)
Phyloseq_hosts_infants_family_filt <-  microbiome::transform(Phyloseq_hosts_infants_family_filt, transform = "compositional")

top_family_names <- tax_table(Phyloseq_hosts_infants_family_filt)@.Data[,5]
taxa_names(Phyloseq_hosts_infants_family_filt) <- top_family_names

# Generate abundance plot at family and genus level
pdf('6_HOST_ASSIGNMENT/Plots/Infant_Family_Abundance_barplot.pdf', width = 6, height = 3.1)
Phyloseq_hosts_infants_family_filt  %>%
  plot_composition(average_by = "Timepoint_categorical") +
  #scale_y_continuous(labels = function(x) paste0(x * 100), limits = c(0, 1)) +
  scale_fill_manual(values = c("#8DD3C7", "#AEC7E8", "#BEBADA", "#FB8072", "#FFFFB3",
                               "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9",
                               "#BC80BD", "#CCEBC5", "#FFED6F", "#E31A1C", "#FD8D3C")) +
  theme_classic() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 13, face = "italic"),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        panel.background = element_blank(),
        plot.background = element_blank()) +
  ylab("Relative Abundance (%)") +
  xlab("Timepoint") +
  guides(fill = guide_legend(title = "Family"))
dev.off()

#********
# Genus-level
#********

# Generate vOTU abundances grouped by bacterial host at genus level and save abundance table
Phyloseq_hosts_infants_genus <- Phyloseq_hosts_infants %>% 
  subset_taxa(!is.na(Genus)) %>%
  tax_glom(taxrank = "Genus") 

genus_names <- tax_table(Phyloseq_hosts_infants_genus)@.Data[,6]
taxa_names(Phyloseq_hosts_infants_genus) <- genus_names

# Save viral abundance table grouped by bacterial host at genus level
Bacterial_genus_host_vOTU_abundance_infants <- data.frame(t(otu_table(Phyloseq_hosts_infants_genus)))
write.table(Bacterial_genus_host_vOTU_abundance_infants, "6_HOST_ASSIGNMENT/Bacterial_genus_host_vOTU_abundance_infants_17052024.txt",
            sep = "\t", row.names = T, quote = FALSE)

top_genera <- names(sort(taxa_sums(Phyloseq_hosts_infants_genus), TRUE)[1:12])
Phyloseq_hosts_infants_genus_filt <- prune_taxa(top_genera, Phyloseq_hosts_infants_genus)
Phyloseq_hosts_infants_genus_filt <-  microbiome::transform(Phyloseq_hosts_infants_genus_filt, transform = "compositional")

pdf('6_HOST_ASSIGNMENT/Plots/Infant_Genus_Abundance_barplot.pdf', width = 6, height = 3.1)
Phyloseq_hosts_infants_genus_filt  %>%
  plot_composition(average_by = "Timepoint_categorical") +
  scale_y_continuous(labels = function(x) paste0(x * 100), limits = c(0, 1)) +
  scale_fill_manual(values = c("#8DD3C7", "#AEC7E8", "#BEBADA", "#FB8072", "#FFFFB3",
                               "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9",
                               "#BC80BD", "#CCEBC5", "#FFED6F", "#E31A1C", "#FD8D3C")) +
  theme_classic() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 13,face = "italic"),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        panel.background = element_blank(),
        plot.background = element_blank()) +
  ylab("Relative Abundance (%)") +
  xlab("Timepoint") +
  guides(fill = guide_legend(title = "Genus"))
dev.off()


##*************
# Save output
#**************
# Save virus metadata
write.table(Virus_metadata, "Metadata_CS/CS_Baby_Biome_Viral_Metadata_17052024.txt", sep = "\t", row.names = F, quote = FALSE)

