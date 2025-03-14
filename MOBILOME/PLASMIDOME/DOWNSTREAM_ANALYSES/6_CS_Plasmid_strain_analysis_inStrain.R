################################################################################
##### CS Baby Biome: Plasmid strain transmission 
### Author(s): Asier Fernández-Pato
### Last updated: 2nd July, 2024
################################################################################


#****************
# Load libraries
#****************
library(dplyr)
library(data.table)
library(ggplot2)
library(ggmosaic)

#****************
# Define functions
#****************
pairwise_chisq <- function(cont_table, levels) {
  results <- matrix(NA, nrow = length(levels), ncol = length(levels))
  # Get column names 
  col_names <- colnames(cont_table)
  for (i in 1:(length(levels) - 1)) {
    for (j in (i + 1):length(levels)) {
      # Subset rows and columns 
      sub_table <- cont_table[c(levels[i], levels[j]), ]
      # Perform chi-square test
      result <- chisq.test(sub_table)
      results[i, j] <- results[j, i] <- result$p.value
    }
  }
  rownames(results) <- colnames(results) <- levels
  return(results)
}

# Function to extract COG terms for BAKTA annotation
extract_COG <- function(x) {
  terms <- trimws(strsplit(x, ",")[[1]]) # Split the string by commas and trim whitespace
  cog_terms <- grep("COG:", terms, value = TRUE) # Find terms that contain "COG"
  cog_counts <- sapply(cog_terms, function(term) sum(gregexpr("COG", term)[[1]] > 0)) # Check how many times "COG" appears in each term
  single_cog_term <- cog_terms[cog_counts == 1] # Find the term with "COG" exactly once
  if (length(single_cog_term) == 1) { # If there is only one term, extract the part after the semicolon
    result <- strsplit(single_cog_term, ":")[[1]][2]
  } else {
    result <- NA
  }
  return(result)
}

# Set working directory 
setwd("~/Desktop/PhD/Projects/CS Baby Biome/Final_analysis/")


##************************************************************************
# 1. Load metadata and abundance table for the 195 CS Baby Biome samples 
#*************************************************************************
# Read metadata tables (with all infants including M06 and M12)
Sample_metadata <- read.delim("VIRUSES/Metadata_CS/Metadata_EGA_CS_BABY_BIOME.txt")

# Update clean read numbers
Sample_metadata$NG_ID <- Sample_metadata$bioSampleId
Clean_reads <- read.delim("VIRUSES/Metadata_CS/CS_Baby_Biome_nreads_195.txt", header = T) 
Sample_metadata$read_depth <- Clean_reads$clean_reads

# Add variable to distinguish between maternal and infant samples
Sample_metadata$Type <- ifelse(Sample_metadata$Timepoint_categorical == "MOM", "Mother", "Infant")
Sample_metadata_infants <- Sample_metadata[Sample_metadata$Type == "Infant", ]
Sample_metadata_mothers <- Sample_metadata[Sample_metadata$Timepoint_categorical=="MOM",]

# Read abundance table
Abundance_table <- read.delim("PLASMIDS/ABUNDANCE_TABLES/CS_Abundance_Table_25062024.txt")
Abundance_table_infants <- Abundance_table[, colnames(Abundance_table) %in% Sample_metadata_infants$NG_ID]
Abundance_table_infants <- Abundance_table_infants[rowSums(Abundance_table_infants) != 0, ]

Plasmid_metadata <- read.delim("PLASMIDS/METADATA_TABLES/CS_Baby_Biome_Plasmid_Metadata_All_25062024.txt")
Plasmid_metadata_infants <- Plasmid_metadata[Plasmid_metadata$Plasmid_ID %in% rownames(Abundance_table_infants),]


##************************************************************************
# 2. Select plasmids present in at least 2 samples
#*************************************************************************

# Select only those plasmids present in at least 2 samples
PTUs_inStrain <- rownames(Abundance_table)[rowSums(Abundance_table != 0) >= 2]

# Save list of vOTUs to be included for strain analysis
write.table(PTUs_inStrain,"PLASMIDS/STRAIN_TRANSMISION_inSTRAIN/list_PTUs_2samples.txt", sep = "\t", 
            row.names = F, col.names = F, quote = FALSE)

##************************************************************************
# 3. Metadata and inStrain results formatting/prcessing
#*************************************************************************
# After mapping the vOTUs_inStrain genomes to the 195 CS samples and running inStrain, we load the results
inStrain_results <- fread("PLASMIDS/6_STRAIN_TRANSMISION_inStrain/3_INSTRAIN_COMPARISONS_comparisonsTable.tsv")

#Extract comparisons where the genomes were present in both samples with >75% of the bases covered in both
inStrain_results_filtered <- inStrain_results[inStrain_results$percent_genome_compared >= 0.75 ,]
inStrain_results_filtered$name1 <- gsub("\\.sorted\\.bam$", "", inStrain_results_filtered$name1)
inStrain_results_filtered$name2 <- gsub("\\.sorted\\.bam$", "", inStrain_results_filtered$name2)
inStrain_results_filtered <- inStrain_results_filtered[, c("scaffold","name1","name2","popANI")]
colnames(inStrain_results_filtered) <- c("Plasmid_ID", "Sample1", "Sample2", "popANI")
inStrain_results <- inStrain_results_filtered

# Add metadata to the table:
## Relatedness: mother-infant pair vs non-pair (infant vs infant should be NA) 
## Within vs between individual comparison
## MOM-Week1 transmission vs NO (to explore vertical transmission in the first week of life)
inStrain_results$Mobility <- Plasmid_metadata$Mobility[match(inStrain_results$Plasmid_ID, Plasmid_metadata$Plasmid_ID)]
inStrain_results$Topology <- Plasmid_metadata$topology_simple[match(inStrain_results$Plasmid_ID, Plasmid_metadata$Plasmid_ID)]
inStrain_results$FAM_ID1 <- Sample_metadata$CS_BABY_BIOME_ID[match(inStrain_results$Sample1, Sample_metadata$bioSampleId)]
inStrain_results$FAM_ID2 <- Sample_metadata$CS_BABY_BIOME_ID[match(inStrain_results$Sample2, Sample_metadata$bioSampleId)]
inStrain_results$Timepoint_1 <- Sample_metadata$Timepoint_categorical[match(inStrain_results$Sample1, Sample_metadata$bioSampleId)]
inStrain_results$Timepoint_2 <- Sample_metadata$Timepoint_categorical[match(inStrain_results$Sample2, Sample_metadata$bioSampleId)]
inStrain_results$Type1 <- Sample_metadata$Type[match(inStrain_results$Sample1, Sample_metadata$bioSampleId)]
inStrain_results$Type2 <- Sample_metadata$Type[match(inStrain_results$Sample2, Sample_metadata$bioSampleId)]
inStrain_results$Feeding_mode1 <- Sample_metadata$feeding_mode_pragmatic[match(inStrain_results$Sample1, Sample_metadata$bioSampleId)]
inStrain_results$Feeding_mode2 <- Sample_metadata$feeding_mode_pragmatic[match(inStrain_results$Sample2, Sample_metadata$bioSampleId)]
inStrain_results$Living_env1 <- Sample_metadata$living_situation[match(inStrain_results$Sample1, Sample_metadata$bioSampleId)]
inStrain_results$Living_env2 <- Sample_metadata$living_situation[match(inStrain_results$Sample2, Sample_metadata$bioSampleId)]
inStrain_results$Comparison <- ifelse(
  inStrain_results$Type1 == "Mother" & inStrain_results$Type2 == "Mother", "Mother-Mother",
  ifelse(
    inStrain_results$Type1 == "Infant" & inStrain_results$Type2 == "Infant","Infant-Infant",
         "Mother-Infant")
)

inStrain_results$Individual_ID1 <- ifelse(
  inStrain_results$Timepoint_1 == "MOM",
  paste(inStrain_results$FAM_ID1, inStrain_results$Timepoint_1, sep = "_"),
  paste(inStrain_results$FAM_ID1, "INFANT", sep = "_")
)
inStrain_results$Individual_ID2 <- ifelse(
  inStrain_results$Timepoint_2 == "MOM",
  paste(inStrain_results$FAM_ID2, inStrain_results$Timepoint_2, sep = "_"),
  paste(inStrain_results$FAM_ID2, "INFANT", sep = "_")
)

inStrain_results$Within_Between <- ifelse(inStrain_results$Individual_ID1 == inStrain_results$Individual_ID2,
                                          "Within", "Between")

inStrain_results$Mom_Infant_pair <- ifelse(
  inStrain_results$Comparison == "Mother-Infant" & inStrain_results$FAM_ID1 == inStrain_results$FAM_ID2, "Pair",
  ifelse(
    inStrain_results$Comparison == "Mother-Infant" & inStrain_results$FAM_ID1 != inStrain_results$FAM_ID2,
    "Not pair",
    NA
  )
)

# Add ID to identify individual mother-infant events
inStrain_results$Mother_infant_comparison_ID <- ifelse(
  inStrain_results$Comparison == "Mother-Infant",
  ifelse(grepl("MOM", inStrain_results$Individual_ID1),
         paste(inStrain_results$Plasmid_ID, inStrain_results$Individual_ID1, inStrain_results$Individual_ID2, sep = "_"),
         paste(inStrain_results$Plasmid_ID, inStrain_results$Individual_ID2, inStrain_results$Individual_ID1, sep = "_")),
  NA
)

inStrain_results$Mother_infant_transmission <- ifelse(
  inStrain_results$Comparison == "Mother-Infant" & inStrain_results$popANI > 0.99999 , "Yes", "No")


# Add living environment, pets and feeding information to the comparison (only infants)
# Note that comparisons with M06, M12 and MOM samples are NA for these phenotypes (also within W01-W06 with different feeding/living values)
# We exclude Pets for now as some W1-W6 samples have also NAs for it
inStrain_results$Feeding <- ifelse(
  inStrain_results$Feeding_mode1 == "formula_feeding" & inStrain_results$Feeding_mode2 == "formula_feeding", "Formula-fed",
  ifelse(
    inStrain_results$Feeding_mode1 == "breast_feeding" & inStrain_results$Feeding_mode2 == "breast_feeding", "Breastfed",
    ifelse(
      inStrain_results$Feeding_mode1 == "mixed_feeding" & inStrain_results$Feeding_mode2 == "mixed_feeding", "Mix-fed",NA
      )
    )
  )

inStrain_results$Living_env <- ifelse(
  inStrain_results$Living_env1 == "city" & inStrain_results$Living_env2 == "city", "City",
  ifelse(
    inStrain_results$Living_env1 == "farm" & inStrain_results$Living_env2 == "farm", "Farm",
    ifelse(
      inStrain_results$Living_env1 == "village" & inStrain_results$Living_env2 == "village", "Village",NA
    )
  )
)


##************************************************************************
# 4. Exploring plasmid strain distances
#*************************************************************************

# Explore distances in related vs unrelated mother-infant pairs
inStrain_results_mother_infant <- inStrain_results[inStrain_results$Comparison == "Mother-Infant",]

# Perform statistical test
inStrain_results_mother_infant$Mom_Infant_pair <- factor(inStrain_results_mother_infant$Mom_Infant_pair,
                                                         levels = c("Pair", "Not pair"))

# Perform statistical test
# Due to lack of independence (repeated measures within infants), we used a permutation-based test
# Extract the observed test statistic: difference in medians between pairs and non-pairs
observed_stat <- median(inStrain_results_mother_infant$popANI[inStrain_results_mother_infant$Mom_Infant_pair == "Pair"]) - 
  median(inStrain_results_mother_infant$popANI[inStrain_results_mother_infant$Mom_Infant_pair == "Not pair"])

# Number of permutations
num_permutations <- 1000
permuted_stats <- numeric(num_permutations)

# Perform the permutation test
set.seed(123)
for (i in 1:num_permutations) {
  permuted_labels <- sample(inStrain_results_mother_infant$Mom_Infant_pair)
  permuted_stat <- median(inStrain_results_mother_infant$popANI[permuted_labels == "Pair"]) - 
    median(inStrain_results_mother_infant$popANI[permuted_labels == "Not pair"])
  permuted_stats[i] <- permuted_stat
}

# Calculate the p-value: proportion of permuted stats greater than or equal to observed (one-sided hypothesis)
p_value <- mean(permuted_stats >= observed_stat) 

# Generate plot
pdf('PLASMIDS/6_STRAIN_TRANSMISION_inStrain/Mother_infant_distances_relatedness.pdf', width=3.2, height=3.2)
ggplot(inStrain_results_mother_infant, 
       aes(x = Mom_Infant_pair, y = popANI, fill = Mom_Infant_pair)) +
  geom_jitter(aes(color=Mom_Infant_pair), alpha=0.7,) + 
  geom_boxplot(width=0.5, fill="white", color="black") + 
  scale_fill_manual(values=c("#800020", "#D3D3D3")) + 
  scale_color_manual(values=c("#800020", "#D3D3D3")) +
  labs(x = "Relatedness", y = "popANI", fill = "Relatedness") +
  theme_classic() +
  theme(plot.title = element_text(size=18, hjust = 0.5, face="bold"), 
        axis.text = element_text(size=14),
        axis.title = element_text(size=15), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.position = "none") 
dev.off()

# Explore distances according to feeding mode in related mother-infant pairs
inStrain_results_related_mother_infant <- inStrain_results[inStrain_results$Comparison == "Mother-Infant" &
                                                          inStrain_results$Mom_Infant_pair == "Pair",]

# Add infant feeding mode information (for MOM-W01 to MOM-W06 comparisons)
inStrain_results_related_mother_infant$Feeding <- NULL
inStrain_results_related_mother_infant$Feeding <- ifelse(!is.na(inStrain_results_related_mother_infant$Feeding_mode1),
                                                         inStrain_results_related_mother_infant$Feeding_mode1,
                                                         ifelse(!is.na(inStrain_results_related_mother_infant$Feeding_mode2),
                                                                inStrain_results_related_mother_infant$Feeding_mode2,
                                                                NA))


# Perform statistical test (permutation-based test)
inStrain_results_related_mother_infant_noNA <- inStrain_results_related_mother_infant[!is.na(inStrain_results_related_mother_infant$Feeding), ]

# Calculate the observed statistic: difference in medians between feeding groups
observed_stat <- median(inStrain_results_related_mother_infant_noNA$popANI[inStrain_results_related_mother_infant_noNA$Feeding == "formula_feeding"]) - 
  median(inStrain_results_related_mother_infant_noNA$popANI[inStrain_results_related_mother_infant_noNA$Feeding == "mixed_feeding"])

# Set the number of permutations
num_permutations <- 1000
permuted_stats <- numeric(num_permutations)

# Perform the permutation test
set.seed(123)
for (i in 1:num_permutations) {
  permuted_labels <- sample(inStrain_results_related_mother_infant_noNA$Feeding)  
  permuted_stat <- median(inStrain_results_related_mother_infant_noNA$popANI[permuted_labels == "formula_feeding"]) - 
    median(inStrain_results_related_mother_infant_noNA$popANI[permuted_labels == "mixed_feeding"])
  permuted_stats[i] <- permuted_stat
}

# Calculate the p-value: proportion of permuted stats greater than or equal to observed (one-sided hypothesis)
p_value <- mean(permuted_stats >= observed_stat)

# Generate plot
pdf('PLASMIDS/6_STRAIN_TRANSMISION_inStrain/Mother_infant_distances_related_feeding.pdf', width=3.2, height=3.2)
ggplot(inStrain_results_related_mother_infant[!is.na(inStrain_results_related_mother_infant$Feeding),], 
       aes(x = Feeding, y = popANI, fill = Feeding)) +
  geom_violin(aes(fill=Feeding), alpha=0.7, width=0.8, trim = F) +
  geom_jitter(aes(color=Feeding), alpha=0.7,) + 
  geom_boxplot(width=0.2, fill="white", color="black") + 
  scale_fill_manual(values=c("#C7C7FF", "#DFF0D8")) + 
  scale_color_manual(values=c("#C7C7FF", "#DFF0D8")) +
  labs(x = "Relatedness", y = "popANI", fill = "Relatedness") +
  theme_classic() +
  theme(plot.title = element_text(size=18, hjust = 0.5, face="bold"), 
        axis.text = element_text(size=14),
        axis.title = element_text(size=15), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.position = "none") 
dev.off()


##************************************************************************
# 5. Exploring plasmid transmission
#*************************************************************************

# For A and B, we retain only 1 value for each Mother_infant_comparison_ID (each mother-infant-PTU comparison)
# This will allow to count only once each transmission event (and not multiple times if it happens in  multiple timepoints)
# The transmission value will be "Yes" if we detect potential transmission at any timepoint

# A) Including mother-infant comparisons at all timepoints
#-----------------------------------------------------------

# Count potential transmission events in related/unrelated pairs
# Build contingency table and test significance
mother_infant_transmission_data <- na.omit(inStrain_results %>%
  group_by(Mother_infant_comparison_ID) %>%
  summarise(Mother_infant_transmission = ifelse(any(Mother_infant_transmission == "Yes"), "Yes", "No")))

#3,358 mother-infant comparisons (including comparisons of the same plasmid in the same pair at different timepoint)
#2,034 mother-infant comparisons (excluding comparisons of the same plasmid in the same pair at different timepoint)
Pairedness <- na.omit(inStrain_results[,c("Mother_infant_comparison_ID", "Mom_Infant_pair")])
Pairedness <- unique(Pairedness) # Unique comparisons
mother_infant_transmission_data <- left_join(mother_infant_transmission_data,
                                             Pairedness, by = "Mother_infant_comparison_ID")
contingency_table <- table(mother_infant_transmission_data$Mom_Infant_pair,
                           mother_infant_transmission_data$Mother_infant_transmission)

# Perform chi-square test
chi2_test <- chisq.test(contingency_table)
chi2_test


# B) Including mother-infant comparisons at W01-W06 (exclude M06/M12)
#-----------------------------------------------------------

# Filter inStrain results to exclude M06/M12 comparisons
inStrain_results_W01_W06 <- inStrain_results %>%
  filter(!grepl("M06|M12", Timepoint_1) & !grepl("M06|M12", Timepoint_2)) 

# Count potential transmission events in related/unrelated pairs
mother_infant_transmission_data_W01_W06 <- inStrain_results_W01_W06 %>%
  group_by(Mother_infant_comparison_ID) %>%
  summarise(Mother_infant_transmission = ifelse(any(Mother_infant_transmission == "Yes"), "Yes", "No"))

#585 non-M06/M12 mother-infant comparisons (excluding comparisons of the same plasmid in the same pair at different timepoint)
Pairedness_W01_W06 <- na.omit(inStrain_results_W01_W06[,c("Mother_infant_comparison_ID", "Mom_Infant_pair")])
Pairedness_W01_W06 <- unique(Pairedness_W01_W06) # Unique comparisons
mother_infant_transmission_data_W01_W06 <- left_join(mother_infant_transmission_data_W01_W06,
                                                     Pairedness_W01_W06, by = "Mother_infant_comparison_ID")
contingency_table_W01_W06 <- table(mother_infant_transmission_data_W01_W06$Mom_Infant_pair,
                                   mother_infant_transmission_data_W01_W06$Mother_infant_transmission)

# Perform Fisher test (lower sample-size)
fisher_test_W01_W06 <- fisher.test(contingency_table_W01_W06)
fisher_test_W01_W06

#################################
# Barplots transmission by relatedness
#################################
prop_table <- prop.table(contingency_table, margin = 1)  
prop_table_W01_W06 <- prop.table(contingency_table_W01_W06, margin = 1)  

# Convert to data frame for ggplot
transmission_relatedness_ggplot <- as.data.frame.matrix(prop_table)
transmission_relatedness_ggplot$Mother_infant_pair <- rownames(transmission_relatedness_ggplot)
transmission_relatedness_ggplot <- tidyr::gather(transmission_relatedness_ggplot, key = "Mother_infant_transmission",
                                                 value = "Proportion", -Mother_infant_pair)
transmission_relatedness_ggplot$Mother_infant_pair <- factor(transmission_relatedness_ggplot$Mother_infant_pair,
                                                             levels = c("Pair", "Not pair"))

transmission_relatedness_ggplot_W01_W06 <- as.data.frame.matrix(prop_table_W01_W06)
transmission_relatedness_ggplot_W01_W06$Mother_infant_pair <- rownames(transmission_relatedness_ggplot_W01_W06)
transmission_relatedness_ggplot_W01_W06 <- tidyr::gather(transmission_relatedness_ggplot_W01_W06,
                                                         key = "Mother_infant_transmission", value = "Proportion", -Mother_infant_pair)
transmission_relatedness_ggplot_W01_W06$Mother_infant_pair <- factor(transmission_relatedness_ggplot_W01_W06$Mother_infant_pair,
                                                                     levels = c("Pair", "Not pair"))

# Generate plots
pdf('PLASMIDS/6_STRAIN_TRANSMISION_inStrain/Transmission_events_relatedness.pdf', width=4.2, height=3.2)
ggplot(transmission_relatedness_ggplot, aes(x = product(Mother_infant_pair), fill = Mother_infant_transmission, weight = Proportion)) +
  geom_mosaic(aes(x = product(Mother_infant_pair), fill = Mother_infant_transmission, weight = Proportion)) +
  labs(x = "Mother-infant Pair", y = "Proportion", fill = "Transmission") +
  scale_fill_manual(name = "Transmission", values = c("No" = "#D3D3D3", "Yes" = "#1F77B4")) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 14),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
dev.off()

pdf('PLASMIDS/6_STRAIN_TRANSMISION_inStrain/Transmission_events_relatedness_W01_W06.pdf', width=4.2, height=3.2)
ggplot(transmission_relatedness_ggplot_W01_W06, aes(x = product(Mother_infant_pair), fill = Mother_infant_transmission, weight = Proportion)) +
  geom_mosaic(aes(x = product(Mother_infant_pair), fill = Mother_infant_transmission, weight = Proportion)) +
  labs(x = "Mother-infant Pair", y = "Proportion", fill = "Transmission") +
  scale_fill_manual(name = "Transmission", values = c("No" = "#D3D3D3", "Yes" = "#1F77B4")) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 14),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
dev.off()


# C) Enrichment of PTU potential transmission at early or late timepoints 
#-----------------------------------------------------------

# Considering the results, we test the enrichment of comparisons above the popANI threshold in later timepoints (M06 and M12).
# Here, unlike for mother_infant_transmission_data, we keep all the comparisons instead of one value per mother-infant-plasmid (so we can check the timepoint)
enrichment_transmission_data_timepoint <- inStrain_results[inStrain_results$Comparison== "Mother-Infant" & 
                                                             inStrain_results$Mom_Infant_pair == "Pair",]

early_timepoints <- c("W01", "W02", "W03", "W04", "W05", "W06")

enrichment_transmission_data_timepoint$Early_late_timepoint <- ifelse(
  enrichment_transmission_data_timepoint$Timepoint_1 %in% early_timepoints | 
    enrichment_transmission_data_timepoint$Timepoint_2 %in% early_timepoints,
  "Early","Late")

contingency_table_enrichment_timepoint <- table(enrichment_transmission_data_timepoint$Early_late_timepoint, 
                           enrichment_transmission_data_timepoint$Mother_infant_transmission)

# Perform chi-square test
chi2_test <- chisq.test(contingency_table_enrichment_timepoint)
chi2_test

# Generate the plot
prop_table_enrichment_timepoint <- prop.table(contingency_table_enrichment_timepoint, margin = 1)  

# Convert to data frame for ggplot
transmission_enrichment_timepoint <- as.data.frame.matrix(prop_table_enrichment_timepoint)
transmission_enrichment_timepoint$Timepoint <- rownames(transmission_enrichment_timepoint)
transmission_enrichment_timepoint <- tidyr::gather(transmission_enrichment_timepoint, 
                                                 key = "Mother_infant_transmission", value = "Proportion", -Timepoint)

pdf('PLASMIDS/6_STRAIN_TRANSMISION_inStrain/Transmission_early_late_timepoints.pdf', width=4.2, height=3.2)
ggplot(transmission_enrichment_timepoint, aes(x = product(Timepoint), fill = Mother_infant_transmission, weight = Proportion)) +
  geom_mosaic(aes(x = product(Timepoint), fill = Mother_infant_transmission, weight = Proportion)) +
  labs(x = "Time point", y = "Proportion", fill = "Transmission") +
  scale_fill_manual(name = "Transmission", values = c("No" = "#D3D3D3", "Yes" = "#1F77B4")) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 14),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
dev.off()


# In addition, check plasmid mobility in transmitted/shared plasmids
# Generate a stacked barplot showing mobility of tranmsitted palsmdis per timepoint 
enrichment_transmission_data_timepoint_transmission <- enrichment_transmission_data_timepoint[enrichment_transmission_data_timepoint$Mother_infant_transmission == "Yes",]
mobility_transmitted_plasmids <- data.frame(cbind(enrichment_transmission_data_timepoint_transmission[,c("Plasmid_ID", "Mobility",
                                                                                                         "Timepoint_1", "Timepoint_2" )]))
mobility_transmitted_plasmids <- mobility_transmitted_plasmids %>%
  mutate(Timepoint = case_when(
    Timepoint_1 != "MOM" ~ Timepoint_1,
    Timepoint_2 != "MOM" ~ Timepoint_2,
    TRUE ~ NA_character_
  ))

mobility_transmitted_plasmids$Mobility <- factor(
  mobility_transmitted_plasmids$Mobility,
  levels = c("Non-mobilizable", "Mobilizable", "Conjugative")
)

mobility_transmitted_plasmids$Timepoint <- factor(
  mobility_transmitted_plasmids$Timepoint,
  levels = c("W01", "W02", "W03", "W04", "W05", "W06", "M06", "M12")
)

# Create the stacked bar plot
pdf('PLASMIDS/6_STRAIN_TRANSMISION_inStrain/Transmission_mobily_timepoint.pdf', width=6.7, height=3.2)
ggplot(mobility_transmitted_plasmids, aes(x = Timepoint, fill = Mobility)) +
  geom_bar(position = "fill") +
  labs(x = "Time point",
       y = "Proportion (%)") +
  scale_fill_manual(
    name = "Mobility",
    values = c("Non-mobilizable" = "#FFEDA0CC", 
               "Mobilizable" = "#F9B14FCC", 
               "Conjugative" = "#EA5133CC")) +
  scale_y_continuous(labels = scales::label_percent()) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 16),
    strip.background = element_rect(fill = "white")
  )
dev.off()


#################################
# Analysis transmission by feeding mode
#################################

# For this, consider only transmissions between related pairs
# Also, for feeding mode, 2 analysis will be performed for:
# A) All timepoints
# B) W01_W06 (excluding comparisons with M06 and M12)
 
# A) Related mother-infant comparisons at any timepoint
#-------------------------------------------------------
# First, add feeding information to the transmission table
# Extract the infant ID from comparison ID and check if that infant ID is associated with a feeding value in Sample_metadata
Feeding_info <- inStrain_results[,c("Mother_infant_comparison_ID", "Feeding", "Individual_ID1"             
                                  ,"Individual_ID2")]

Feeding_info <- Feeding_info %>%
  filter(!is.na(Mother_infant_comparison_ID)) %>%
  distinct(Mother_infant_comparison_ID, .keep_all = TRUE)  %>%
  mutate(
    INFANT_ID = ifelse(grepl("_INFANT$", Individual_ID1), sub("_INFANT$", "", Individual_ID1), sub("_INFANT$", "", Individual_ID2))
  )

Sample_metadata_feeding <-Sample_metadata[,c("CS_BABY_BIOME_ID", "feeding_mode_pragmatic"),]

Feeding_info_extended <- Feeding_info %>%
  left_join(Sample_metadata_feeding, by = c("INFANT_ID" = "CS_BABY_BIOME_ID")) %>%
  group_by(INFANT_ID) 

for (Comparison_ID in unique(Feeding_info$Mother_infant_comparison_ID)) {
  Feeding_mode <- unique(na.omit(Feeding_info_extended[Feeding_info_extended$Mother_infant_comparison_ID == Comparison_ID,
                                              "feeding_mode_pragmatic"]))
  Feeding_info[Feeding_info$Mother_infant_comparison_ID == Comparison_ID, "Feeding"] <- Feeding_mode
}

Feeding_info <- Feeding_info[,c("Mother_infant_comparison_ID", "Feeding")]

mother_infant_transmission_data <- left_join(mother_infant_transmission_data, Feeding_info,
                                             by = "Mother_infant_comparison_ID")

#Extract only comparisons between related pairs
mother_infant_transmission_data_related <- mother_infant_transmission_data[mother_infant_transmission_data$Mom_Infant_pair == "Pair",]

#-------------------------
# Perform chi-square tests and plots
#-------------------------

contingency_table_feeding_all <- table(mother_infant_transmission_data_related$Feeding,
                                    mother_infant_transmission_data_related$Mother_infant_transmission)

chi2_test_feeding_all <- chisq.test(contingency_table_feeding_all)
chi2_test_feeding_all 

# Generate plot 
prop_table_feeding_all <- prop.table(contingency_table_feeding_all, margin = 1)  

# Convert to data frame for ggplot
transmission_feeding_ggplot_all <- as.data.frame.matrix(prop_table_feeding_all)
transmission_feeding_ggplot_all$feeding <- rownames(transmission_feeding_ggplot_all)
transmission_feeding_ggplot_all <- tidyr::gather(transmission_feeding_ggplot_all, 
                                              key = "Mother_infant_transmission", value = "Proportion", -feeding)

transmission_feeding_ggplot_all$feeding <- factor(transmission_feeding_ggplot_all$feeding,
                                                  levels = c("breast_feeding", "mixed_feeding", "formula_feeding"),
                                                  labels = c("Breastfed", "Mix-fed", "Formula-fed"))

pdf('PLASMIDS/6_STRAIN_TRANSMISION_inStrain/Transmission_events_feeding.pdf', width=4.2, height=3.2)
ggplot(transmission_feeding_ggplot_all, aes(x = product(feeding), fill = Mother_infant_transmission, weight = Proportion)) +
  geom_mosaic(aes(x = product(feeding), fill = Mother_infant_transmission, weight = Proportion)) +
  labs(x = "Feeding mode", y = "Proportion", fill = "Transmission") +
  scale_fill_manual(name = "Transmission", values = c("No" = "#D3D3D3", "Yes" = "#1F77B4")) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 14),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
dev.off()

#B) Related mother-infant comparisons at W01-W06 
#-------------------------------------------------------

#Extract only comparisons between related pairs
mother_infant_transmission_data_W01_W06_related <- na.omit(mother_infant_transmission_data_W01_W06[mother_infant_transmission_data_W01_W06$Mom_Infant_pair == "Pair",])

# Add feeding information
mother_infant_transmission_data_W01_W06_related <- left_join(mother_infant_transmission_data_W01_W06_related, Feeding_info,
                                             by = "Mother_infant_comparison_ID")

#-------------------------
# Perform chi-square tests and plots
#-------------------------

contingency_table_feeding_W01_W06 <- table(mother_infant_transmission_data_W01_W06_related$Feeding,
                                       mother_infant_transmission_data_W01_W06_related$Mother_infant_transmission)

fisher_test_feeding_W01_W06 <- fisher.test(contingency_table_feeding_W01_W06)
fisher_test_feeding_W01_W06 #1 (low number of comparisons)

# Generate plot 
prop_table_feeding_W01_W06 <- prop.table(contingency_table_feeding_W01_W06, margin = 1)  

# Convert to data frame for ggplot
transmission_feeding_ggplot_W01_W06 <- as.data.frame.matrix(prop_table_feeding_W01_W06)
transmission_feeding_ggplot_W01_W06$feeding <- rownames(transmission_feeding_ggplot_W01_W06)
transmission_feeding_ggplot_W01_W06 <- tidyr::gather(transmission_feeding_ggplot_W01_W06, 
                                                 key = "Mother_infant_transmission", value = "Proportion", -feeding)

transmission_feeding_ggplot_W01_W06$feeding <- factor(transmission_feeding_ggplot_W01_W06$feeding,
                                                  levels = c("breast_feeding", "mixed_feeding", "formula_feeding"),
                                                  labels = c("Breastfed", "Mix-fed", "Formula-fed"))


pdf('PLASMIDS/6_STRAIN_TRANSMISION_inStrain/Transmission_events_feeding_W01W06.pdf', width=4.2, height=3.2)
ggplot(transmission_feeding_ggplot_W01_W06, aes(x = product(feeding), fill = Mother_infant_transmission, weight = Proportion)) +
  geom_mosaic(aes(x = product(feeding), fill = Mother_infant_transmission, weight = Proportion)) +
  labs(x = "Feeding mode", y = "Proportion", fill = "Transmission") +
  scale_fill_manual(name = "Transmission", values = c("No" = "#D3D3D3", "Yes" = "#1F77B4")) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 14),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
dev.off()


#################################
# Analysis transmission by plasmid topology and mobility (related and non-related transmission at all timepoints)
#################################
Topology <- na.omit(inStrain_results[,c("Mother_infant_comparison_ID", "Topology")])
Topology <- unique(Topology) # Unique comparisons
Mobility <- na.omit(inStrain_results[,c("Mother_infant_comparison_ID", "Mobility")])
Mobility <- unique(Mobility) # Unique comparisons

mother_infant_transmission_data_related <- left_join(mother_infant_transmission_data_related, Topology,
                                             by = "Mother_infant_comparison_ID")
mother_infant_transmission_data_related  <- left_join(mother_infant_transmission_data_related, Mobility,
                                             by = "Mother_infant_comparison_ID")
#-------------------------
# Perform chi-square tests and plots
#-------------------------

# A) Topology
#-------------
contingency_table_topology <- table(mother_infant_transmission_data_related$Topology,
                                    mother_infant_transmission_data_related$Mother_infant_transmission)

chi2_test_topology <- chisq.test(contingency_table_topology)
chi2_test_topology

# Generate plot 
prop_table_topology <- prop.table(contingency_table_topology, margin = 1)  

# Convert to data frame for ggplot
transmission_topology_ggplot <- as.data.frame.matrix(prop_table_topology)
transmission_topology_ggplot$Topology <- rownames(transmission_topology_ggplot)
transmission_topology_ggplot <- tidyr::gather(transmission_topology_ggplot, 
                                              key = "Mother_infant_transmission", value = "Proportion", -Topology)

pdf('PLASMIDS/6_STRAIN_TRANSMISION_inStrain/Transmission_events_topology.pdf', width=4.2, height=3.2)
ggplot(transmission_topology_ggplot, aes(x = product(Topology), fill = Mother_infant_transmission, weight = Proportion)) +
  geom_mosaic(aes(x = product(Topology), fill = Mother_infant_transmission, weight = Proportion)) +
  labs(x = "Topology", y = "Proportion", fill = "Transmission") +
  scale_fill_manual(name = "Transmission", values = c("No" = "#D3D3D3", "Yes" = "#1F77B4")) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 14),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.text.y = element_blank(),   
    axis.ticks.y = element_blank() 
  )
dev.off()

#B) Mobility
#-------------
contingency_table_mobility <- table(mother_infant_transmission_data_related$Mobility,
                                    mother_infant_transmission_data_related$Mother_infant_transmission)

chi2_test_mobility <- chisq.test(contingency_table_mobility)
chi2_test_mobility

# As we have 3 levels of mobility, check which levels show significant results
levels <- rownames(contingency_table_mobility)
pairwise_results <- pairwise_chisq(contingency_table_mobility, levels)
pairwise_results

# Generate plot 
prop_table_mobility <- prop.table(contingency_table_mobility, margin = 1)  

# Convert to data frame for ggplot
transmission_mobility_ggplot <- as.data.frame.matrix(prop_table_mobility)
transmission_mobility_ggplot$mobility <- rownames(transmission_mobility_ggplot)
transmission_mobility_ggplot <- tidyr::gather(transmission_mobility_ggplot, 
                                              key = "Mother_infant_transmission", value = "Proportion", -mobility)

# Generate plots
pdf('PLASMIDS/6_STRAIN_TRANSMISION_inStrain/Transmission_events_mobility.pdf', width=4.2, height=3.2)
ggplot(transmission_mobility_ggplot, aes(x = product(Mobility), fill = Mother_infant_transmission, weight = Proportion)) +
  geom_mosaic(aes(x = product(mobility), fill = Mother_infant_transmission, weight = Proportion)) +
  labs(x = "Mobility", y = "Proportion", fill = "Transmission") +
  scale_fill_manual(name = "Transmission", values = c("No" = "#D3D3D3", "Yes" = "#1F77B4")) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 14),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.text.y = element_blank(),   
    axis.ticks.y = element_blank() 
  )
dev.off()


##************************************************************************
# 5. Characterization and functional enrichment of transmitted pTUs
#*************************************************************************
transmitted_plasmids <- mother_infant_transmission_data_related %>%
  filter(Mom_Infant_pair == "Pair" & Mother_infant_transmission == "Yes") %>%
  mutate(Plasmid_ID = sub("_D98.*", "", Mother_infant_comparison_ID))

transmitted_plasmids$Feeding <- NULL

# Add AMR presence
AMR_data <- Plasmid_metadata[,c("Plasmid_ID", "AMR_genes_binary")]
transmitted_plasmids <- left_join(transmitted_plasmids, AMR_data, by = "Plasmid_ID")

# Functional enrichment analysis
# Load functional analysis results (BAKTA)
Functional_annotation_PTUs <- read.delim("PLASMIDS/6B_FUNCTIONAL_ANNOTATION/bakta_annotations_all_plasmids_merged.tsv")
colnames(Functional_annotation_PTUs)[1] <- "Plasmid_ID"
Functional_annotation_transmitted_PTUs <- Functional_annotation_PTUs[Functional_annotation_PTUs$Plasmid_ID %in%
                                                                       transmitted_plasmids$Plasmid_ID,] 

all_PTUs_COG <- sapply(Functional_annotation_PTUs$DbXrefs, extract_COG)
all_PTUs_COG <- as.character(all_PTUs_COG[!is.na(all_PTUs_COG)])
transmitted_PTUs_COG <- sapply(Functional_annotation_transmitted_PTUs$DbXrefs, extract_COG)
transmitted_PTUs_COG <- transmitted_PTUs_COG[!is.na(transmitted_PTUs_COG)]

# Generate a count table with the number of proteins belonging to each functional category
split_categories_all <- unlist(strsplit(all_PTUs_COG, ""))
all_PTUs_COG_counts <- table(split_categories_all)
split_categories_transmitted <- unlist(strsplit(transmitted_PTUs_COG, ""))
transmitted_PTUs_COG_counts <- table(split_categories_transmitted)

# Select only COG categories present in transmitted plasmids to test their enrichment
transmitted_categories <- names(transmitted_PTUs_COG_counts)

# Perform Fisher's Exact Test for each category in transmitted PTUs

# Create a data frame to store results
results_enrichment <- data.frame(
  Category = character(),
  PValue = numeric(),
  OddsRatio = numeric(),
  stringsAsFactors = FALSE
)

for (category in transmitted_categories) {
  # Create contingency table
  contingency_table <- matrix(c(
    transmitted_PTUs_COG_counts[category], sum(transmitted_PTUs_COG_counts) - transmitted_PTUs_COG_counts[category],
    all_PTUs_COG_counts[category], sum(all_PTUs_COG_counts) - all_PTUs_COG_counts[category]
  ), nrow = 2, byrow = TRUE)
  
  # Perform Fisher's Exact Test
  fisher_test_result <- fisher.test(contingency_table)
  
  # Store result
  results_enrichment <- rbind(results_enrichment, data.frame(
    Category = category,
    PValue = fisher_test_result$p.value,
    OddsRatio = fisher_test_result$estimate
  ))
}

# Add FDR adjustment to the results data frame
results_enrichment$FDR <- p.adjust(results_enrichment$PValue, method = "BH")


# Generate plots for significant categories
contingency_table_L <- matrix(c(transmitted_PTUs_COG_counts["L"], sum(transmitted_PTUs_COG_counts) - transmitted_PTUs_COG_counts["L"],
                              all_PTUs_COG_counts["L"], sum(all_PTUs_COG_counts) - all_PTUs_COG_counts["L"]), nrow = 2, byrow = TRUE)
contingency_table_N <- matrix(c(transmitted_PTUs_COG_counts["N"], sum(transmitted_PTUs_COG_counts) - transmitted_PTUs_COG_counts["N"],
                                all_PTUs_COG_counts["N"], sum(all_PTUs_COG_counts) - all_PTUs_COG_counts["N"]), nrow = 2, byrow = TRUE)
contingency_table_M <- matrix(c(transmitted_PTUs_COG_counts["M"], sum(transmitted_PTUs_COG_counts) - transmitted_PTUs_COG_counts["M"],
                                all_PTUs_COG_counts["M"], sum(all_PTUs_COG_counts) - all_PTUs_COG_counts["M"]), nrow = 2, byrow = TRUE)

prop_table_L <- prop.table(contingency_table_L, margin = 1)
category_L_enrichment <- data.frame(prop_table_L)
rownames(category_L_enrichment) <- c("Transmitted", "All")
colnames(category_L_enrichment) <- c("L", "Other category")
category_L_enrichment$Transmission <- rownames(category_L_enrichment)

prop_table_N <- prop.table(contingency_table_N, margin = 1)
category_N_enrichment <- data.frame(prop_table_N)
rownames(category_N_enrichment) <- c("Transmitted", "All")
colnames(category_N_enrichment) <- c("N", "Other category")
category_N_enrichment$Transmission <- rownames(category_N_enrichment)

prop_table_M <- prop.table(contingency_table_M, margin = 1)
category_M_enrichment <- data.frame(prop_table_M)
rownames(category_M_enrichment) <- c("Transmitted", "All")
colnames(category_M_enrichment) <- c("M", "Other category")
category_M_enrichment$Transmission <- rownames(category_M_enrichment)

combined_enrichment <- bind_rows(
  category_L_enrichment %>% mutate(Category = "L"),
  category_N_enrichment %>% mutate(Category = "N"),
  category_M_enrichment %>% mutate(Category = "M")
)

combined_enrichment_long <- combined_enrichment %>%
  pivot_longer(cols = c("L", "N", "M", "Other category"), names_to = "CategoryType", values_to = "Proportion") %>%
  mutate(Category = factor(Category, levels = c("L", "N", "M")),
         CategoryType = factor(CategoryType, levels = c("L", "N", "M", "Other category")),
         Transmission = factor(Transmission, levels = c("Transmitted", "All")))

# Create the stacked bar plot
pdf('PLASMIDS/6_STRAIN_TRANSMISION_inStrain/Functional_categories_proportions_stacked.pdf', width=6.8, height=3.6)
ggplot(combined_enrichment_long, aes(x = Transmission, y = Proportion, fill = CategoryType)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~Category) +
  labs(x = "Transmission", y = "Proportion", fill = "Category Type") +
  scale_fill_manual(name = "Category Type", values = c("L" = "#0055AA", "N" = "#00BFC4", "M" = "#F8766D", "Other category" = "#D3D3D3")) +
  scale_y_continuous(labels = scales::label_percent()) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 16),
    strip.background = element_rect(fill = "white"),
  )
dev.off()

#****************
# Write results
#****************
write.table(results_enrichment,"PLASMIDS/6_FUNCTIONAL_ANNOTATION/Functional_enrichment_transmitted_PTUs.txt", sep = "\t", 
            row.names = T, col.names = T, quote = FALSE)
