################################################################################
##### CS Baby Biome: Bacterial abundance analysis
### Author(s):Asier Fernández-Pato
### Last updated: 20th June, 2024
################################################################################

#****************
# Define functions
#****************

# Function to process Metaphlan4 abundance table
Clean_MP4 = function(Microbial){
  unknown_rel = Microbial$UNCLASSIFIED
  Microbial %>% select(-UNCLASSIFIED) -> M2
  lapply(colnames(M2)[2:dim(M2)[2]], FUN= function(x){ str_split(x,"\\.")[[1]] -> y ;  return(y[length(y)]) }) %>% as_vector() -> New_names
  colnames(M2) = c("ID", New_names)
  #to go back to 100% without unknown we need to apply Back_to_100 in each taxonomic level 
  lapply(str_split(colnames(select(M2, -ID)), "__"), FUN=function(x){ x[1] }) %>% as_vector() ->Taxonomies ; Taxonomies %>% unique() -> Taxonomy_Levels
  Bugs = colnames(select(M2, -ID)) ; tibble(Bug = Bugs, Level = Taxonomies) -> Bug_array
  M3 = tibble(ID = M2$ID)
  for (T_L in Taxonomy_Levels){
    filter(Bug_array, Level == T_L) -> Bugs
    Back_to_100(select(M2, Bugs$Bug)) -> Normalized
    M3 = left_join(M3, mutate(Normalized, ID = M2$ID ),  by= "ID" )
  }
  return(list(unknown_rel, M2, M3, Bug_array))
}

# Function to recalculate the relative abundance after removing the unclassified
Back_to_100 = function(Microbial){
  #We can go back to 100%, but this function should be applied per taxonomic level
  apply(Microbial,1, FUN= function(x){ sum(x) -> y ;  return(y) }) %>% as_vector() -> Total_per_sample
  Microbial/Total_per_sample -> Normalized
  return(as_tibble(Normalized))
}

#****************
# Load modules
#****************
library(dplyr)
library(purrr)
library(stringr)
library(phyloseq)
library(microbiome)
library(ggplot2)

#*# Set working directory 
setwd("~/Desktop/PhD/Projects/CS Baby Biome/Final_analysis/VIRUSES/")

# Read abundance table and metadata tables
Abundance_table <- read.delim("Abundance_table/CS_Abundance_Table_17052024.txt")
Sample_metadata <- read.delim("Metadata_CS/CS_Baby_Biome_Sample_Metadata_17052024.txt")
Virus_metadata <- read.delim("Metadata_CS/CS_Baby_Biome_Viral_Metadata_17052024.txt")

##************************************************************************
# 1. Bacterial Abundance Analysis: Metaphlan results preprocessing
##************************************************************************
# Read Metaphlan4 output table
bacterial_abundances <- read.delim("4_BACTERIAL_ANALYSIS/taxa_clean_mp4_oct_2022_01_04_2024.txt", 
                          sep = "\t", header = T)

# Select only those samples present in Sample metadata
bacterial_abundances <- bacterial_abundances[rownames(bacterial_abundances) %in% Sample_metadata$NG_ID,]

# Process the table to use it as input of Clean_MP4 function
tax <- bacterial_abundances
tax$ID <- rownames(bacterial_abundances)
tax$ID=make.unique(tax$ID, sep = '_')
tax <- tax[,c(ncol(tax), 1:(ncol(tax) - 1))] # bringing ID to the 1st column
tax$NG_ID <- NULL

# Process Metaphlan4 table
List_microbial = Clean_MP4(tax) # This removes the unclassified column and shortens the taxa name. 
Unknown_fraction = List_microbial[[1]] # unknown fractions
Microbial_abundance = List_microbial[[2]] # relative abundances without re-assessment
Microbial_abundance_reassigned = List_microbial[[3]] # relative abundances after re-assessment
Taxonomy_level = List_microbial[[4]] # Microbes and their specific taxonomic levels

# Get species level abundance table
microbes <-as.data.frame(Microbial_abundance)
row.names(microbes)<-microbes$ID
microbes$ID=NULL
species_abundance <- microbes[,grep("s__",colnames(microbes))]
species_abundance <- species_abundance[,-which(colSums(species_abundance)==0)] #remove those species not present
colnames(species_abundance) <- sub("^s__", "", colnames(species_abundance))  

# Get genus level abundance table
genus_abundance <- microbes[,grep("g__",colnames(microbes))]
colnames(genus_abundance) <- str_extract(colnames(genus_abundance), "(?<=__).*")
genus_abundance <- genus_abundance[,-which(colSums(genus_abundance)==0)] #remove those genera not present

# Create the taxonomy table and format it as desired
Tax_table <- as.data.frame(t(species_abundance))
clade <- grep("s__", colnames(tax), value = T)
clade <- clade[grep("t__", clade, invert = TRUE)]
names_clade <- sub(".*\\.", "", clade)
names_clade <- sub("^s__", "", names_clade)
names_clade_exclude <- which(!names_clade %in% colnames(species_abundance))
clade <- clade [-names_clade_exclude]

Tax_table$Clade <- clade
Tax_table <- separate(data = Tax_table, col = Clade, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = "\\.")
Tax_table <- Tax_table[,(ncol(Tax_table)-6):ncol(Tax_table)]
Tax_table <- data.frame(lapply(Tax_table, function(x) str_extract(x, "(?<=__).*")))
rownames(Tax_table) <- Tax_table$Species

##################################
# Generation of Phyloseq object
##################################

# Reorder samples according to Sample metadata
species_abundance <- data.frame(species_abundance[match(Sample_metadata$NG_ID,rownames(species_abundance)),])

# Transform dataframes into tables (for Phyloseq)
Abundance_table_matrix <- as.matrix(t(species_abundance))
taxonomy_hosts_matrix <- as.matrix(Tax_table)

# Generate phyloseq object
vOTU <- otu_table(Abundance_table_matrix, taxa_are_rows = TRUE)
TAX <- tax_table(taxonomy_hosts_matrix)
samples <- sample_data(Sample_metadata)
sample_names(samples) <- Sample_metadata$bioSampleId

Phyloseq_bacteria <- phyloseq(vOTU, TAX, samples)

##************************************************************************
# 2. Bacterial Taxonomic Analysis in Infants
##************************************************************************

##################################
# Host abundance analysis
##################################

# Select only infant samples
Phyloseq_bacteria_infants <- subset_samples(Phyloseq_bacteria, Type =="Infant")

Phyloseq_bacteria_infants@sam_data$Timepoint_categorical <- factor(Phyloseq_bacteria_infants@sam_data$Timepoint_categorical,
                                                                    levels = c("W01","W02", "W03","W04","W05","W06"))

# Generate genus-level abundances
Phyloseq_bacteria_infants_genus <- Phyloseq_bacteria_infants %>% 
  subset_taxa(!is.na(Genus)) %>%
  aggregate_taxa(level = "Genus") 

top_genera <- names(sort(taxa_sums(Phyloseq_bacteria_infants_genus), TRUE)[1:10])
Phyloseq_bacteria_infants_genus <- prune_taxa(top_genera, Phyloseq_bacteria_infants_genus)
Phyloseq_bacteria_infants_genus <-  microbiome::transform(Phyloseq_bacteria_infants_genus, transform = "compositional")

# Generate abundance plot at genus level
pdf('5_BACTERIAL_ANALYSIS/Plots/Infant_Genus_Abundance_barplot.pdf', width = 6, height = 3.1)
Phyloseq_bacteria_infants_genus  %>%
  plot_composition(average_by = "Timepoint_categorical") +
  scale_y_continuous(labels = function(x) paste0(x * 100)) +
  scale_fill_manual(values = c("#8DD3C7", "#AEC7E8", "#BEBADA", "#FB8072", "#FFFFB3",
                               "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9",
                               "#CCEBC5", "#BC80BD", "#FFED6F", "#E31A1C", "#FD8D3C")) +
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
write.table(species_abundance, "4_BACTERIAL_ANALYSIS/Metaphlan_4_species_21062024.txt", sep = "\t", row.names = F, quote = FALSE)
write.table(genus_abundance, "4_BACTERIAL_ANALYSIS/Metaphlan_4_genus_21062024.txt", sep = "\t", row.names = T, quote = FALSE)

