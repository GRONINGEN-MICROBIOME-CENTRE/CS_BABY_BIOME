################################################################################
##### CS Baby Biome: Viral strain transmission
### Author(s): Asier Fernández-Pato
### Last updated: 30th June, 2024
################################################################################


#****************
# Load libraries
#****************
library(dplyr)
library(data.table)
library(ggplot2)


# Set working directory 
setwd("~/Desktop/PhD/Projects/CS Baby Biome/Final_analysis/VIRUSES/")


##************************************************************************
# 1. Load metadata and abundance table for the 195 CS Baby Biome samples 
#*************************************************************************
# Read metadata tables (with all infants including M06 and M12)
Sample_metadata <- read.delim("Metadata_CS/Metadata_EGA_CS_BABY_BIOME.txt")

# Update clean read numbers
Sample_metadata$NG_ID <- Sample_metadata$bioSampleId
Clean_reads <- read.delim("Metadata_CS/CS_Baby_Biome_nreads_195.txt", header = T) 
Sample_metadata$read_depth <- Clean_reads$clean_reads

# Add variable to distinguish between maternal and infant samples
Sample_metadata$Type <- ifelse(Sample_metadata$Timepoint_categorical == "MOM", "Mother", "Infant")
Sample_metadata_infants <- Sample_metadata[Sample_metadata$Type == "Infant", ]
Sample_metadata_mothers <- Sample_metadata[Sample_metadata$Timepoint_categorical=="MOM",]

# Read abundance table
Abundance_table <- read.delim("Abundance_table/CS_Baby_Biome_Abundance_Table_RPKM.txt")
Abundance_table_infants <- Abundance_table[, colnames(Abundance_table) %in% Sample_metadata_infants$NG_ID]

Virus_metadata <- read.delim("Metadata_CS/CS_Baby_Biome_Viral_Metadata_17052024.txt")
Virus_metadata_infants <- Virus_metadata[Virus_metadata$Virus_ID %in% rownames(Abundance_table_infants),]


##************************************************************************
# 2. Select viruses present in at least 2 samples
#*************************************************************************

# Remove from the table those viruses not present in any sample
Absent_viruses <- rownames(Abundance_table[rowSums(Abundance_table)==0, ]) 
Abundance_table <- Abundance_table[!rownames(Abundance_table) %in% Absent_viruses,] 
Present_viruses <- rownames(Abundance_table) #2,263

# Select only those viruses present in at least 2 samples
vOTUs_inStrain <- rownames(Abundance_table)[rowSums(Abundance_table != 0) >= 2]

# Save list of vOTUs to be included for strain analysis
write.table(vOTUs_inStrain,"8_STRAIN_TRANSMISION_inSTRAIN/list_vOTU_2samples.txt", sep = "\t", 
            row.names = F, col.names = F, quote = FALSE)

##************************************************************************
# 3. Process in inStrain results
#*************************************************************************
# After mapping the vOTUs_inStrain genomes to the 195 CS samples and running inStrain, we load the results
inStrain_results <- fread("8_STRAIN_TRANSMISION_inStrain/3_INSTRAIN_COMPARISONS_comparisonsTable.tsv")

#Extract comparisons where the genomes were present in both samples with >75% of the bases covered in both
inStrain_results_filtered <- inStrain_results[inStrain_results$percent_genome_compared >= 0.75 ,]
inStrain_results_filtered$name1 <- gsub("\\.sorted\\.bam$", "", inStrain_results_filtered$name1)
inStrain_results_filtered$name2 <- gsub("\\.sorted\\.bam$", "", inStrain_results_filtered$name2)
inStrain_results_filtered <- inStrain_results_filtered[, c("scaffold","name1","name2","popANI")]
colnames(inStrain_results_filtered) <- c("Virus_ID", "Sample1", "Sample2", "popANI")
inStrain_results <- inStrain_results_filtered

# Add metadata to the table:
## Relatedness: mother-infant pair vs non-pair (infant vs infant should be NA) 
## Within vs between individual comparison
## MOM-Week1 transmission vs NO (to explore vertical transmission in the first week of life)
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
inStrain_results$Host <- Virus_metadata$Bacterial_genus_host[match(inStrain_results$Virus_ID, Virus_metadata$Virus_ID)]
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
         paste(inStrain_results$Virus_ID, inStrain_results$Individual_ID1, inStrain_results$Individual_ID2, sep = "_"),
         paste(inStrain_results$Virus_ID, inStrain_results$Individual_ID2, inStrain_results$Individual_ID1, sep = "_")),
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
# 4. Exploring viral strain distances
#*************************************************************************

# Explore distances in related vs unrelated mother-infant pairs
inStrain_results_mother_infant <- inStrain_results[inStrain_results$Comparison == "Mother-Infant",]

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

# Set seed for reproducibility
set.seed(123)

# Perform the permutation test
for (i in 1:num_permutations) {
  permuted_labels <- sample(inStrain_results_mother_infant$Mom_Infant_pair)
  permuted_stat <- median(inStrain_results_mother_infant$popANI[permuted_labels == "Pair"]) - 
    median(inStrain_results_mother_infant$popANI[permuted_labels == "Not pair"])
  permuted_stats[i] <- permuted_stat
}

# Calculate the p-value: proportion of permuted stats greater than or equal to observed (one-sided hypothesis)
p_value <- mean(permuted_stats >= observed_stat) 


# Generate plot
pdf('STRAIN_TRANSMISION_inSTRAIN/PLOTS/Mother_infant_distances_relatedness.pdf', width=3.2, height=3.2)
ggplot(inStrain_results_mother_infant, 
       aes(x = Mom_Infant_pair, y = popANI, fill = Mom_Infant_pair)) +
  geom_jitter(aes(color = Mom_Infant_pair), alpha = 0.7) + 
  geom_boxplot(width = 0.5, fill = "white", color = "black") + 
  scale_fill_manual(values = c("#800020", "#D3D3D3")) + 
  scale_color_manual(values = c("#800020", "#D3D3D3")) +
  labs(x = "Relatedness", y = "popANI", fill = "Relatedness") +
  scale_y_continuous(breaks = seq(0.96, 1.00, by = 0.01)) +  # Set specific y-axis breaks
  theme_classic() +
  theme(
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"), 
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 15), 
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.position = "none"
  ) 
dev.off()


##************************************************************************
# 5. Exploring viral transmission
#*************************************************************************

# A) Including mother-infant comparisons at all timepoints
#-----------------------------------------------------------

# Count potential transmission events in related/unrelated pairs
# Build contingency table and test significance
mother_infant_transmission_data <- na.omit(inStrain_results %>%
                                             group_by(Mother_infant_comparison_ID) %>%
                                             summarise(Mother_infant_transmission = ifelse(any(Mother_infant_transmission == "Yes"), "Yes", "No")))

#165 mother-infant comparisons (including comparisons of the same virus in the same pair at different timepoint)
#98 mother-infant comparisons (excluding comparisons of the same virus in the same pair at different timepoint)
Pairedness <- na.omit(inStrain_results[,c("Mother_infant_comparison_ID", "Mom_Infant_pair")])
Pairedness <- unique(Pairedness) # Unique comparisons
mother_infant_transmission_data <- left_join(mother_infant_transmission_data, Pairedness, by = "Mother_infant_comparison_ID")
contingency_table <- table(mother_infant_transmission_data$Mom_Infant_pair, mother_infant_transmission_data$Mother_infant_transmission)

# Perform chi-square test
chi2_test <- chisq.test(contingency_table)
chi2_test


# B) Including mother-infant comparisons at W01-W06 (exclude M06/M12)
#-----------------------------------------------------------

# Filter inStrain results to exclude M06/M12 comparisons
inStrain_results_W01_W06 <- inStrain_results %>%
  filter(!grepl("M06|M12", Timepoint_1) & !grepl("M06|M12", Timepoint_2)) 

# Count potential transmission events in related/unrelated pairs
mother_infant_transmission_data_W01_W06 <- na.omit(inStrain_results_W01_W06 %>%
  group_by(Mother_infant_comparison_ID) %>%
  summarise(Mother_infant_transmission = ifelse(any(Mother_infant_transmission == "Yes"), "Yes", "No")))

#30 non-M06/M12 mother-infant comparisons (excluding comparisons of the same plasmid in the same pair at different timepoint)
Pairedness_W01_W06 <- na.omit(inStrain_results_W01_W06[,c("Mother_infant_comparison_ID", "Mom_Infant_pair")])
Pairedness_W01_W06 <- unique(Pairedness_W01_W06) # Unique comparisons
mother_infant_transmission_data_W01_W06 <- left_join(mother_infant_transmission_data_W01_W06, Pairedness_W01_W06, by = "Mother_infant_comparison_ID")
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
transmission_relatedness_ggplot$Mother_infant_pair <- factor(transmission_relatedness_ggplot$Mother_infant_pair, levels = c("Pair", "Not pair"))

transmission_relatedness_ggplot_W01_W06 <- as.data.frame.matrix(prop_table_W01_W06)
transmission_relatedness_ggplot_W01_W06$Mother_infant_pair <- rownames(transmission_relatedness_ggplot_W01_W06)
transmission_relatedness_ggplot_W01_W06 <- tidyr::gather(transmission_relatedness_ggplot_W01_W06, key = "Mother_infant_transmission",
                                                         value = "Proportion", -Mother_infant_pair)
transmission_relatedness_ggplot_W01_W06$Mother_infant_pair <- factor(transmission_relatedness_ggplot_W01_W06$Mother_infant_pair,
                                                                     levels = c("Pair", "Not pair"))


# Generate plots
pdf('STRAIN_TRANSMISION_inStrain/PLOTS/Transmission_events_relatedness.pdf', width=4.2, height=3.2)
ggplot(transmission_relatedness_ggplot, aes(x = Mother_infant_pair, y = Proportion, fill = Mother_infant_transmission)) +
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

pdf('STRAIN_TRANSMISION_inStrain/PLOTS/Transmission_events_relatedness_W01_W06.pdf', width=4.2, height=3.2)
ggplot(transmission_relatedness_ggplot_W01_W06, aes(x = Mother_infant_pair, y = Proportion, fill = Mother_infant_transmission)) +
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

# No analysis of transmission by feeding mode will be made due to low numbers (and only happening in M12)


##************************************************************************
# 5. Characterization of transmitted vOTUs (related pairs)
#*************************************************************************
# Select transmitted viruses
transmitted_viruses <- mother_infant_transmission_data %>%
  filter(Mom_Infant_pair == "Pair" & Mother_infant_transmission == "Yes") %>%
  mutate(Virus_ID = sub("_D98.*", "", Mother_infant_comparison_ID))


# Add host information
mother_infant_transmission_data_host <- inStrain_results[,c("Virus_ID", "Host")]
transmitted_viruses  <- left_join(transmitted_viruses,
                                          mother_infant_transmission_data_host,
                                          by= "Virus_ID")

transmitted_viruses  <- unique(transmitted_viruses)


##************************************************************************
# 6. Other plots
#*************************************************************************
pdf('STRAIN_TRANSMISION_inSTRAIN/PLOTS/Within_between_distances.pdf', width=3.2, height=3.2)
ggplot(inStrain_results[!is.na(inStrain_results$Within_Between),], 
       aes(x = Within_Between, y = popANI, fill = Within_Between)) +
  geom_jitter(aes(color=Within_Between), alpha=0.7,) + 
  geom_boxplot(width=0.5, fill="white", color="black") + 
  scale_fill_manual(values=c("#0055AA", "#C40003")) + 
  scale_color_manual(values=c("#0055AA", "#C40003")) +
  labs(x = "Relatedness", y = "popANI", fill = "Relatedness") +
  theme_classic() +
  theme(plot.title = element_text(size=18, hjust = 0.5, face="bold"), 
        axis.text = element_text(size=14),
        axis.title = element_text(size=15), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.position = "none") 
dev.off()

# Select only unrelated distances infant-infant and mother-mother (strains are more similar between infants or mothers (unrelated))
pdf('STRAIN_TRANSMISION_inSTRAIN/PLOTS/Unrelated_distances.pdf', width=3.2, height=3.2)
ggplot(inStrain_results[!is.na(inStrain_results$Comparison) & inStrain_results$Within_Between == "Between",], 
       aes(x = Comparison, y = popANI, fill = Comparison)) +
  geom_jitter(aes(color=Comparison), alpha=0.7,) + 
  geom_boxplot(width=0.5, fill="white", color="black") + 
  #scale_fill_manual(values=c("#0055AA", "#C40003")) + 
  #scale_color_manual(values=c("#0055AA", "#C40003")) +
  labs(x = "Relatedness", y = "popANI", fill = "Comparison") +
  theme_classic() +
  theme(plot.title = element_text(size=18, hjust = 0.5, face="bold"), 
        axis.text = element_text(size=14),
        axis.title = element_text(size=15), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.position = "none") 
dev.off()

#########
# Feeding mode
#########
# Check if infants with specific feeding mode have more similar strains
pdf('STRAIN_TRANSMISION_inSTRAIN/PLOTS/Feeding_mode_distances.pdf', width=3.2, height=3.2)
ggplot(inStrain_results[!is.na(inStrain_results$Feeding) & inStrain_results$Within_Between == "Between",], 
       aes(x = Feeding, y = popANI, fill = Feeding)) +
  geom_jitter(aes(color=Feeding), alpha=0.7,) + 
  geom_violin() + 
  geom_boxplot(width=0.5, fill="white", color="black") + 
  #scale_fill_manual(values=c("#0055AA", "#C40003")) + 
  #scale_color_manual(values=c("#0055AA", "#C40003")) +
  labs(x = "Feeding mode", y = "popANI", fill = "Feeding") +
  theme_classic() +
  theme(plot.title = element_text(size=18, hjust = 0.5, face="bold"), 
        axis.text = element_text(size=14),
        axis.title = element_text(size=15), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.position = "none") 
dev.off()


# Check strain similarity by bacterial host
pdf('STRAIN_TRANSMISION_inSTRAIN/PLOTS/Feeding_mode_distances_by_host.pdf', width=3.2, height=3.2)
ggplot(inStrain_results_related_mother_infant[!is.na(inStrain_results_related_mother_infant$Feeding) &
                          !is.na(inStrain_results_related_mother_infant$Host) &
                            inStrain_results_related_mother_infant$Within_Between == "Between",],
       aes(x = Feeding, y = popANI, fill = Feeding)) +
  geom_jitter(aes(color = Feeding), alpha = 0.7, position = position_jitterdodge(jitter.width = 0.2)) + 
  geom_violin(alpha = 0.5, position = position_dodge(width = 0.75)) + 
  geom_boxplot(width = 0.3, fill = "white", color = "black", position = position_dodge(width = 0.75)) + 
  facet_wrap(~ Host, scales = "free_x") + 
  labs(x = "Feeding Level", y = "popANI", fill = "Feeding") +
  theme_classic() +
  theme(plot.title = element_text(size = 18, hjust = 0.5, face = "bold"), 
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 15), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.position = "none")
dev.off()


#########
# Living environment and pets
#########
# Check if infants with specific living environment have more related strains
pdf('STRAIN_TRANSMISION_inSTRAIN/PLOTS/Living_environment_distances_infants.pdf', width=3.2, height=3.2)
ggplot(inStrain_results[!is.na(inStrain_results$Living_env) 
                        & inStrain_results$Comparison == "Infant-Infant"
                        & inStrain_results$Within_Between == "Between",], 
       aes(x = Living_env, y = popANI, fill = Living_env)) +
  geom_jitter(aes(color=Living_env), alpha=0.7,) + 
  geom_violin() + 
  geom_boxplot(width=0.5, fill="white", color="black") + 
  #scale_fill_manual(values=c("#0055AA", "#C40003")) + 
  #scale_color_manual(values=c("#0055AA", "#C40003")) +
  labs(x = "Living environment", y = "popANI", fill = "Living_env") +
  theme_classic() +
  theme(plot.title = element_text(size=18, hjust = 0.5, face="bold"), 
        axis.text = element_text(size=14),
        axis.title = element_text(size=15), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.position = "none") 
dev.off()

