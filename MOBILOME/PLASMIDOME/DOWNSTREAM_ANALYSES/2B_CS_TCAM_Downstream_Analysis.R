################################################################################
##### CS Baby Biome: Plasmid Diversity analysis: TCAM downstream analysis
### Author(s): Asier Fernández-Pato, Trishla Sinha
### Last updated: 26th June, 2024
################################################################################

#****************
# Load libraries
#****************
library(vegan)
library(dplyr)
library(ggplot2)


#****************
# Define functions
#****************

permanova_analysis <- function(Distance, selection) {
  # Initialize results data frame
  n_phenotypes <- ncol(selection) - 3  # Exclude the first 3 columns (covariates)
  adonis_results_raw <- data.frame(Df = numeric(n_phenotypes),
                                   F = numeric(n_phenotypes),
                                   R2 = numeric(n_phenotypes),
                                   p_value = numeric(n_phenotypes))
  rownames(adonis_results_raw) <- colnames(selection)[4:ncol(selection)]
  
  # Loop through phenotypes (excluding the first 3)
  for (i in 1:n_phenotypes) {
    # Subset data with complete observations
    na_rows <- which(is.na(selection[, i + 3]))  # Adjust index for excluded columns
    distmat_cleaned <- Distance[!(rownames(Distance) %in% row.names(selection)[na_rows]),
                                !(colnames(Distance) %in% row.names(selection)[na_rows])]
    phenos2 <- selection[row.names(selection) %in% row.names(distmat_cleaned), ]
    
    # Check if any data remains after subsetting
    if (nrow(distmat_cleaned) == 0) {
      warning(paste("No data remaining for phenotype", colnames(selection)[i + 3], sep = " "))
      # Handle empty data case (e.g., set results to NA or skip iteration)
      next  # Skip to the next iteration
    }
    
    # Run Adonis test
    ad1 <- adonis2(distmat_cleaned ~ phenos2[[i + 3]],  # Adjust index for excluded columns
                   permutations = 10000,
                   parallel = 8,
                   na.action = na.fail,
                   by = "margin")
    
    # Save test results
    adonis_results_raw[i, ] <- c(ad1$Df[1], ad1$F[1], ad1$R2[1], ad1$"Pr(>F)"[1])
  }
  
  # Calculate FDR-adjusted p-values
  adonis_results_raw$FDR <- p.adjust(adonis_results_raw$p_value, method = "BH")
  
  # Return results data frame
  return(adonis_results_raw)
}


permanova_analysis_birthweight <- function(Distance, selection) {
  # Identify birthweight column index
  birthweight_col <- which(colnames(selection) == "infant_birthweight")

  # Initialize results data frame
  n_phenotypes <- ncol(selection) - 3  # Exclude the first 3 columns (covariates) and birthweight
  adonis_results <- data.frame(Df = numeric(n_phenotypes),
                                   F = numeric(n_phenotypes),
                                   R2 = numeric(n_phenotypes),
                                   p_value = numeric(n_phenotypes))
  rownames(adonis_results) <- colnames(selection)[4:ncol(selection)]
  
  # Loop through phenotypes (excluding the first 3 and birthweight)
  for (i in 1:n_phenotypes) {
    # Subset data with complete observations
    na_rows <- which(is.na(selection[, c(i + 3, birthweight_col)]))  # Adjust index for excluded columns
    distmat_cleaned <- Distance[!(rownames(Distance) %in% row.names(selection)[na_rows]),
                                !(colnames(Distance) %in% row.names(selection)[na_rows])]
    phenos2 <- selection[row.names(selection) %in% row.names(distmat_cleaned), ]
    birthweight_cleaned <- selection[row.names(selection) %in% row.names(distmat_cleaned), birthweight_col]
    
    # Check if any data remains after subsetting
    if (nrow(distmat_cleaned) == 0) {
      warning(paste("No data remaining for phenotype", colnames(selection)[i + 3], sep = " "))
      next  # Skip to the next iteration
    }
    
    # Run Adonis test with birthweight as covariate (using formula)
    ad1 <- adonis2(distmat_cleaned ~ birthweight_cleaned + phenos2[[i + 3]],
                   permutations = 10000,
                   parallel = 8,
                   na.action = na.fail,
                   by = "margin")
    
    # Save test results
    adonis_results[i, ] <- c(ad1$Df[2], ad1$F[2], ad1$R2[2], ad1$"Pr(>F)"[2])
  }
  
  # Remove rows corresponding to covariates from final table 
  adonis_results <- adonis_results[!grepl("infant_birthweight", rownames(adonis_results)), ]
  
  # Calculate FDR-adjusted p-values
  adonis_results$FDR <- p.adjust(adonis_results$p_value, method = "BH")
  
  # Return results data frame
  return(adonis_results)
}


# Set working directory 
setwd("~/Desktop/PhD/Projects/CS Baby Biome/Final_analysis/PLASMIDS/")

# Read abundance table and metadata tables
Sample_metadata <- read.delim("METADATA_TABLES/CS_Baby_Biome_Sample_Metadata_25062024.txt") 
Sample_metadata_infants <- read.delim("METADATA_TABLES/CS_Baby_Biome_Sample_Metadata_Infants_25062024.txt")
Plasmid_metadata <- read.delim("METADATA_TABLES/CS_Baby_Biome_Plasmid_Metadata_25062024.txt")
Plasmid_metadata_infants <- read.delim("METADATA_TABLES/CS_Baby_Biome_Plasmid_Metadata_Infants_25062024.txt")
Plasmid_abundance <- read.delim("ABUNDANCE_TABLES/CS_Baby_Biome_Plasmid_Abundance_Table.txt")
Plasmid_abundance_infants <- read.delim("ABUNDANCE_TABLES/CS_Baby_Biome_Plasmid_Infants_Abundance_Table.txt")

#************************************************
# 1. Preprocessing and selection of phenotypes
#************************************************

# After running TCAM (see TCAM.ipynb), read the outout tables
TCAM_all_plasmids <- read.delim("2_DIVERSITY/TCAM_All_plasmids_CS_BABY_BIOME.txt", sep = ",")
rownames(TCAM_all_plasmids) <- TCAM_all_plasmids [,1]
TCAM_all_plasmids  <- subset(TCAM_all_plasmids , select = -c(X, CS_BABY_BIOME_ID))

TCAM_circular_plasmids <- read.delim("2_DIVERSITY/TCAM_Circular_plasmids_CS_BABY_BIOME.txt", sep = ",")
rownames(TCAM_circular_plasmids) <- TCAM_circular_plasmids[,1]
TCAM_circular_plasmids <- subset(TCAM_circular_plasmids, select = -c(X))

TCAM_conjugative_plasmids <- read.delim("2_DIVERSITY/TCAM_Conjugative_plasmids_CS_BABY_BIOME.txt", sep = ",")
rownames(TCAM_conjugative_plasmids) <- TCAM_conjugative_plasmids[,1]
TCAM_conjugative_plasmids <- subset(TCAM_conjugative_plasmids, select = -c(X))


# Preparing phenotypes for associations
selection <- Sample_metadata_infants[, c("CS_BABY_BIOME_ID", "read_depth", "DNA_concentration_ng_ul",
                                         "rand_AB", "preg_gest_age","pre_preg_bmi_mother", "Timepoint_numeric",
                                         "infant_birthweight","infant_sex", "feeding_mode_pragmatic",
                                         "living_situation", "cats_dogs")]
selection <- selection[!duplicated(selection$CS_BABY_BIOME_ID), ]
rownames(selection)<-selection$CS_BABY_BIOME_ID


#************************************************
# 2. PERMANOVA analysis
#************************************************

# Generate the Aitchison distance matrices  
Distance_all <- vegdist(TCAM_all_plasmids, method = "euclidean" ) 
Distance_all <- data.frame(as.matrix((Distance_all)))
Distance_circular <- vegdist(TCAM_circular_plasmids, method = "euclidean" ) 
Distance_circular <- data.frame(as.matrix((Distance_circular)))
Distance_conjugative <- vegdist(TCAM_conjugative_plasmids, method = "euclidean" ) 
Distance_conjugative <- data.frame(as.matrix((Distance_conjugative)))

# Run PERMANOVA to estimate the effect of phenotypes on overall plasmid composition
# A) All plasmids
adonis_results_raw <- permanova_analysis(Distance_all, selection) #Birthweight nominally significant
adonis_results_corrected <- permanova_analysis_birthweight(Distance_all, selection) #AB nominally significant
# B) Circular plasmids
adonis_results_circular_raw <- permanova_analysis(Distance_circular, selection) #Birthweight nominally significant
adonis_results_circular_corrected <- permanova_analysis_birthweight(Distance_circular, selection) 
# B) Conjugative plasmids
adonis_results_conjugative_raw <- permanova_analysis(Distance_conjugative, selection)


#************************************************
# 3. PLOTTING RESULTS
#************************************************

# Generate plots of the 1st and 2nd factors of the TCAM analysis grouping samples by AB group, feeding mode and birthweight
all<-merge(selection, TCAM_all_plasmids, by="row.names")
all$feeding_mode_pragmatic<-factor(all$feeding_mode, levels = c("breast_feeding", "mixed_feeding", "formula_feeding"))
all$feeding_mode_pragmatic <- factor(all$feeding_mode_pragmatic, 
                                     levels = c("breast_feeding", "mixed_feeding", "formula_feeding"), 
                                     labels = c("Breastmilk", "Mixed", "Formula"))

# Add infant birthweight as a factor variable
median_birthweight <- median(all$infant_birthweight)
birthweight_category <- ifelse(all$infant_birthweight <= median_birthweight, "lower", "higher")
all$infant_birthweight_factor <- birthweight_category

pdf('2_DIVERSITY/PLOTS/TCAM_Plasmid_INFANTS_AB.pdf', width=3.6, height=3.7)
TCAM_plot<-ggplot(all ,aes(X0,X1, color = rand_AB))+
  geom_point(size = 1.5, alpha=0.7) +
  stat_ellipse(geom = "polygon", alpha = 0.0, aes(group = rand_AB, color=rand_AB), linetype = 2, linewidth = 0.8)+
  xlab("Factor 1 (8.91%)")+
  ylab("Factor 2 (5.43%)")+
  scale_color_manual(values=c("#0055AA","#C40003")) +
  scale_fill_manual(values=c("#0055AA","#C40003")) +
  theme_bw()+
  theme(axis.text=element_text(size=14,  family = "Helvetica"),
        axis.title = element_text(size = 15,  family = "Helvetica"), 
        panel.background = element_blank(), 
        panel.grid = element_blank(),  
        panel.border = element_rect(fill = NA, linewidth =1.2, colour = "grey30"), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 7)),
        legend.title = element_blank(), 
        legend.text = element_text(size = 15, colour = "grey30"),
        legend.position = "bottom") +
  guides(fill = "none")  
ggMarginal(TCAM_plot, type= "boxplot", groupFill=T)
dev.off()

pdf('2_DIVERSITY/PLOTS/TCAM_Plasmid_INFANTS_Feeding.pdf', width=3.6, height=3.7)
TCAM_plot<-ggplot(all ,aes(X0,X1, color = feeding_mode_pragmatic))+
  geom_point(size = 1.5, alpha=0.7) +
  stat_ellipse(geom = "polygon", alpha = 0.0, aes(group = feeding_mode_pragmatic, color=feeding_mode_pragmatic), linetype = 2, linewidth = 0.8)+
  xlab("Factor 1 (8.91%)")+
  ylab("Factor 2 (5.43%)")+
  scale_color_manual(values=c("#FFA255", "#BFE2C5", "#9A9AFF")) +
  scale_fill_manual(values=c("#FFA255", "#BFE2C5", "#9A9AFF")) +
  theme_bw()+
  theme(axis.text=element_text(size=13,  family = "Helvetica"),
        axis.title = element_text(size = 14,  family = "Helvetica"), 
        panel.background = element_blank(), 
        panel.grid = element_blank(),  
        panel.border = element_rect(fill = NA, linewidth =1.2, colour = "grey30"), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.title = element_blank(), 
        legend.text = element_text(size = 14, colour = "grey30"),
        legend.position = "bottom") +
  guides(fill = "none")  
ggMarginal(TCAM_plot, type= "densigram", groupFill=T)
dev.off()

pdf('2_DIVERSITY/PLOTS/TCAM_Plasmid_INFANTS_Birthweight.pdf', width=3.6, height=3.7)
TCAM_plot<-ggplot(all ,aes(X0,X1, color = infant_birthweight_factor))+
  geom_point(size = 1.5, alpha=0.7) +
  stat_ellipse(geom = "polygon", alpha = 0.0, aes(group = infant_birthweight_factor, color=infant_birthweight_factor), linetype = 2, linewidth = 0.8)+
  xlab("Factor 1 (8.91%)")+
  ylab("Factor 2 (5.43%)")+
  scale_color_manual(values = c("#008080", "#FF7F50")) + 
  scale_fill_manual(values = c("#008080", "#FF7F50")) + 
  theme_bw()+
  theme(axis.text=element_text(size=13,  family = "Helvetica"),
        axis.title = element_text(size = 14,  family = "Helvetica"), 
        panel.background = element_blank(), 
        panel.grid = element_blank(),  
        panel.border = element_rect(fill = NA, linewidth =1.2, colour = "grey30"), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.title = element_blank(), 
        legend.text = element_text(size = 14, colour = "grey30"),
        legend.position = "bottom") +
  guides(fill = "none")  
ggMarginal(TCAM_plot, type= "densigram", groupFill=T)
dev.off()

#****************
# Save output
#****************
write.table(adonis_results_raw,"2_DIVERSITY/PERMANOVA_Plasmid_TCAM_Infants_RAW.txt", sep = "\t", 
            row.names = T, col.names = T, quote = FALSE)
write.table(adonis_results_corrected,"2_DIVERSITY/PERMANOVA_Plasmid_TCAM_Infants_CORRECT_BIRTHWEIGHT.txt", sep = "\t", 
            row.names = T, col.names = T, quote = FALSE)




