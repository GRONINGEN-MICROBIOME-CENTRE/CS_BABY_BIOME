################################################################################
##### CS Baby Biome: Plasmid host prediction analysis: CRISPR spacers
### Author(s): Asier Fernández-Pato
### Last updated: 20th February, 2025
################################################################################

# Set working directory 
setwd("~/Desktop/PhD/Projects/CS Baby Biome/Final_analysis/PLASMIDS/")

#****************
# Load libraries
#****************
library(dplyr)
library(tidyr)

#****************
# Define functions
#****************
# Function to assign majority taxonomy for plasmids with multiple matches
# It assigns the more specific taxon that represents >70% of the CRISPR spacer hits 
get_majority_taxonomy <- function(df) {
  tax_levels <- c("Species", "Genus", "Family", "Order", "Class", "Phylum", "Domain")
  
  for (level in tax_levels) {
    # Exclude NAs, then count occurrences of each taxon at the given level
    tax_count <- df %>%
      filter(!is.na(!!sym(level))) %>%
      count(!!sym(level)) %>%
      mutate(percentage = n / sum(n) * 100) %>%
      filter(percentage >= 70) %>%
      arrange(desc(n)) %>%
      slice_head(n = 1)
    
    if (nrow(tax_count) > 0) {
      return(data.frame(Plasmid = df$Plasmid[1], Level_Assigned = level, Taxon = tax_count[[level]]))
    }
  }
  
  return(data.frame(Plasmid = df$Plasmid[1], Level_Assigned = "None", Taxon = NA))
}


##************************************************************************
# 1. Load results from Blastn matches to CRISPR spacers
#*************************************************************************
# Load blastn results
col_names <- c("Plasmid", "Spacer_ID", "Identity", "Alignment_Length", "Mismatches",
               "Plasmid_Start", "Plasmid_End", "Spacer_Start", "Spacer_End")

blast_results <- read.delim("7_HOST_PREDICTION/PTU_crispr_results.tsv", 
                            header = FALSE, col.names = col_names)

# Load length of CRISPR spacers
CRISPR_spacer_length <- read.delim("7_HOST_PREDICTION/CRISPR_length.tsv", 
                                   header = FALSE, col.names = c("Spacer_ID", "Spacer_Length"))

# Merge the blast results with the spacer length
blast_results <- merge(blast_results, CRISPR_spacer_length, by = "Spacer_ID", all.x = TRUE)

# Filter and add Spacer_Coverage column
filtered_blast_results <- blast_results %>%
  mutate(Spacer_Coverage = Alignment_Length / Spacer_Length) %>%
  filter(Alignment_Length >= 25,             
         Spacer_Coverage >= 0.95,  
         Mismatches <= 1)

##************************************************************************
# 2. Associate bacterial hosts to each PTU based on matches
#*************************************************************************

# Read metadata for CRISPR spacers
CRISPR_spacer_metadata <- read.delim("7_HOST_PREDICTION/crispr_spacers_filtered_clustered.tsv")
CRISPR_spacer_metadata <- CRISPR_spacer_metadata[CRISPR_spacer_metadata$spacer_id %in% 
                                                   filtered_blast_results$Spacer_ID,]
# Add taxonomy to blast results
CRISPR_spacer_metadata <- CRISPR_spacer_metadata %>%
  separate(gtdb_lineage, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", remove = FALSE) %>%
  mutate(across(c(Domain, Phylum, Class, Order, Family, Genus, Species), ~ gsub("^[a-zA-Z]+__", "", .)))

filtered_blast_results <- merge(filtered_blast_results, 
                                CRISPR_spacer_metadata[, c("spacer_id", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")], 
                                by.x = "Spacer_ID", by.y = "spacer_id", all.x = TRUE)


################################################################################
# Select as host the most specific taxon that represents >70% of the CRISPR spacer hits
################################################################################

# Count how many CRISPR spacers matched each plasmid
plasmid_counts <- filtered_blast_results %>%
  group_by(Plasmid) %>%
  summarise(num_spacers = n(), .groups = 'drop')

# Separate plasmids with one vs. multiple spacer matches
single_match_plasmids <- filtered_blast_results %>%
  filter(Plasmid %in% plasmid_counts$Plasmid[plasmid_counts$num_spacers == 1])

multi_match_plasmids <- filtered_blast_results %>%
  filter(Plasmid %in% plasmid_counts$Plasmid[plasmid_counts$num_spacers > 1])

# Convert blank values to NA
single_match_plasmids[single_match_plasmids == ""] <- NA
multi_match_plasmids[multi_match_plasmids == ""] <- NA

# Directly assign taxonomy for single-matched plasmids
single_match_taxonomy <- single_match_plasmids %>%
  dplyr::mutate(
    Level_Assigned = case_when(
      !is.na(Species) ~ "Species", !is.na(Genus) ~ "Genus", !is.na(Family) ~ "Family",
      !is.na(Order) ~ "Order", !is.na(Class) ~ "Class", !is.na(Phylum) ~ "Phylum",
      !is.na(Domain) ~ "Domain",TRUE ~ NA_character_
    ),
    Taxon = coalesce(Species, Genus, Family, Order, Class, Phylum, Domain)
  ) %>%
  dplyr::select(Plasmid, Level_Assigned, Taxon)

# Apply the majority taxonomy function to multi-matched plasmids
multi_match_taxonomy <- multi_match_plasmids %>%
  group_split(Plasmid) %>%
  lapply(get_majority_taxonomy) %>%
  bind_rows()

# Combine single and multi-matched taxonomy results, removing NAs
final_plasmid_taxonomy <- bind_rows(single_match_taxonomy, multi_match_taxonomy) %>%
  filter(!is.na(Taxon))

##************************************************************************
# 3. Summary statistics
#*************************************************************************

# Number of plasmids with CRISPR spacer matches (min 1) meeting the criteria
length(table(filtered_blast_results$Plasmid)) #2,079

# Number of plasmids with assigned taxonomy
length(table(final_plasmid_taxonomy$Plasmid)) #2,074 (2,072 if excluding 2 assigned at Domain level to unclassified bacteria)
length(which(final_plasmid_taxonomy$Level_Assigned == "Species")) #847
length(which(final_plasmid_taxonomy$Level_Assigned == "Genus")) #339
length(which(final_plasmid_taxonomy$Level_Assigned == "Family")) #313

# Count the number of plasmids per bacterial species
Bacterial_host_counts_plasmids <- final_plasmid_taxonomy %>%
  filter(Level_Assigned == "Species") %>%  
  group_by(Taxon) %>%
  summarise(n = n_distinct(Plasmid), .groups = 'drop')

# Select the top 10 species with the most plasmids
Bacterial_host_counts_plasmids_top10 <- Bacterial_host_counts_plasmids %>%
  top_n(10, n) %>%
  mutate(Taxon = paste0(substr(Taxon, 1, 1), ".", gsub("^[^ ]+ ", "", Taxon)))

# Generate a barplot with the counts of plasmids per bacterial species (top 10)
pdf('7_HOST_PREDICTION/Plots/Bacterial_host_counts_plasmids_INFANTS.pdf', width = 3, height = 3.5)
ggplot(Bacterial_host_counts_plasmids_top10, aes(x = reorder(Taxon, -n), y = n)) +
  geom_bar(stat = "identity", width = 0.7, fill = "#D3D3D3", color = "#707070") +
  coord_flip() +
  ylab("Number of PTUs per host") + 
  theme_classic() +
  theme(
    axis.text = element_text(size = 16),
    axis.text.y =  element_text(size = 14, face = "italic"),
    axis.title = element_text(size = 17),
    axis.title.y = element_blank(),
    legend.position = "none"
  )
dev.off()

# Write results
#****************
write.table(final_plasmid_taxonomy,"7_HOST_PREDICTION//CS_Plasmid_Host_prediction.txt", sep = "\t", row.names = F, quote = FALSE) 


