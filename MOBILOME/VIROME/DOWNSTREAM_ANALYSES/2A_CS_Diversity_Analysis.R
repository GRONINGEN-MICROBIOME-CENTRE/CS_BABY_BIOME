################################################################################
##### CS Baby Biome: vOTU diversity, specificity and persistence analysis
### Author(s):Asier Fernández-Pato
### Last updated: 6th June, 2024
################################################################################

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

# Function to estimate post-hoc pairwise significance for timepoints 
Timepoint_significance_analysis <- function(data, response_variable, time_variable, group_variable) { 
  # Generate a linear mixed model
  model <- lmer(paste(response_variable, "~", time_variable, "+ DNA_concentration_ng_ul + read_depth + (1|", group_variable, ")"),
                data = data, REML = FALSE)
  # Perform post-hoc pairwise comparisons with Tukey test (multicomp package)
  pairwise_results <- glht(model, linfct = mcp(Timepoint_categorical = "Tukey"))
  pairwise_results <- summary(pairwise_results)
  # Return the overall results and pairwise comparison results
  result <- list(pairwise_results = pairwise_results)
  return(result)
}

#Function to combine timepoints in order for distances comparison 
combine_timepoints <- function(x, y) { # Define function to create the combined timepoint
  paste0(min(x, y), "-", max(x, y))
}

#****************
# Load modules
#****************
library(vegan)
library(lme4)
library(lmerTest)
library(multcomp)
library(dplyr)
library(tidyr)
library(RLRsim)
library(ggplot2)
library(ggrepel)
library(ggExtra)
library(ggpubr)
library(reshape2)


# Set working directory 
setwd("~/Desktop/PhD/Projects/CS Baby Biome/Final_analysis/VIRUSES/")

# Read abundance table and metadata tables
Abundance_table <- read.delim("Abundance_table/CS_Abundance_Table_17052024.txt")
Sample_metadata <- read.delim("Metadata_CS/CS_Baby_Biome_Sample_Metadata_17052024.txt")
Virus_metadata <- read.delim("Metadata_CS/CS_Baby_Biome_Viral_Metadata_17052024.txt")

##************************************************************************
# 0. Preprocessing
#*************************************************************************

# Exclude infants of M6 and M12 from the metadata
Sample_metadata <- Sample_metadata[!Sample_metadata$Timepoint_categorical %in% c("M06", "M12"),]
Abundance_table <- Abundance_table[, colnames(Abundance_table) %in% Sample_metadata$NG_ID]

# Remove viruses that are not present after sample removal from Abundance table and virus metadata
Abundance_table <- Abundance_table[rowSums(Abundance_table != 0) > 0, ] 
rownames(Virus_metadata) <- Virus_metadata$Virus_ID
Virus_metadata <- Virus_metadata[rownames(Abundance_table), ]
  
# Set Timepoint and Type variable from metadata as factors
Sample_metadata$Timepoint_categorical <- factor(Sample_metadata$Timepoint_categorical, 
                                                levels=c("MOM", "W01", "W02", "W03","W04",
                                                         "W05", "W06"))
Sample_metadata$Type <- factor(Sample_metadata$Type)


##************************************************************************
# 1. Exploratory analyses
#*************************************************************************

# Generate the metadata and abundance tables for mothers and infants
# Estimate the number of vOTU present in maternal and infant samples
Sample_metadata_infants <- Sample_metadata[Sample_metadata$Type == "Infant",]
Abundance_table_infants <- Abundance_table[,Sample_metadata$Type == "Infant"]
Abundance_table_infants <- Abundance_table_infants[rowSums(Abundance_table_infants)>0,] #933 viruses

Sample_metadata_mothers <- Sample_metadata[Sample_metadata$Type == "Mother",]
Abundance_table_mothers <- Abundance_table[,Sample_metadata$Type == "Mother"]
Abundance_table_mothers <- Abundance_table_mothers[rowSums(Abundance_table_mothers)>0,] #1,174 viruses 

# Generate a variable with the total RPKM counts of the viral fraction per sample
# Estimate the summary statistics for mothers and infants
Sample_metadata$Viral_abundance <- colSums(Abundance_table)
summary_stats(Sample_metadata$Viral_abundance [Sample_metadata$Timepoint_categorical =="MOM"])
summary_stats(Sample_metadata$Viral_abundance [Sample_metadata$Timepoint_categorical !="MOM"])

# Add viral prevalence to metadata
Virus_metadata$Prevalence <- rowSums(Abundance_table != 0)


##************************************************************************
# 2. Overall diversity analyses: Richness, alpha and beta diversity
#*************************************************************************

##################################
# 2.1. Estimation of alpha diversity and richness 
##################################

# Estimate shannon diversity and richness
Sample_metadata$richness  <-specnumber(t(Abundance_table))                 
Sample_metadata$shannon <- vegan::diversity(t(Abundance_table), index = "shannon") 

# Test significance: different in mothers vs babies
# A) Richness
MM_type_richness <- lmer(richness ~ Type + DNA_concentration_ng_ul + read_depth + (1|CS_BABY_BIOME_ID), REML = F, data = Sample_metadata)
summary(MM_type_richness)

# B) Shannon Index
MM_type_shannon <- lmer(shannon ~ Type + DNA_concentration_ng_ul + read_depth + (1|CS_BABY_BIOME_ID), REML = F, data = Sample_metadata)
summary(MM_type_shannon) 

# Plot the richness and diversity in mothers and babies
pdf('2_DIVERSITY/Plots/Shannon_Index_comp.pdf', width=2.5, height=3.2)
ggplot(Sample_metadata[!is.na(Sample_metadata$Timepoint_categorical),], aes(x=Type, y=shannon)) +
  geom_violin(aes(fill=Type), alpha=0.5, width=0.8, trim = F) +
  geom_jitter(aes(color=Type), alpha=0.7) +  
  geom_boxplot(width=0.2, fill="white", color="black") + 
  labs(x = 'Sample type', y = 'Shannon Index') + 
  scale_fill_manual(values=c("#008080", "#4B0082")) +  
  scale_color_manual(values=c("#008080", "#4B0082")) +  
  theme_classic() +
  theme(
    plot.title = element_text(size=18, hjust = 0.5, face="bold"), 
    axis.text = element_text(size=14),
    axis.title = element_text(size=15), 
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.position = "none" 
  )
dev.off()

pdf('2_DIVERSITY/Plots/Richness_comp.pdf', width=2.9, height=3.2)
ggplot(Sample_metadata[!is.na(Sample_metadata$Timepoint_categorical),], aes(x=Type, y=richness)) +
  geom_violin(aes(fill=Type), alpha=0.5, width=0.8, trim = F) +
  geom_jitter(aes(color=Type), alpha=0.7) +  # Added color aesthetic
  geom_boxplot(width=0.2, fill="white", color="black") + 
  labs(x = 'Sample type', y = 'Viral Richness') + 
  scale_fill_manual(values=c("#008080", "#4B0082")) + 
  scale_color_manual(values=c("#008080", "#4B0082")) + 
  theme_classic() + 
  theme(plot.title = element_text(size=18, hjust = 0.5, face="bold"), 
        axis.text = element_text(size=14),
        axis.title = element_text(size=15), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.position = "none") 
dev.off()

# Plot the correlations with read depth and DNA concentration (change variables to be plotted)
pdf('2_DIVERSITY/Plots/Spearman_cor_Shannon_Read_depth_splitted.pdf', width=2.5, height=2.9)
ggscatter(Sample_metadata, x = "shannon", y = "read_depth",
          add = "reg.line", color = "Type" , conf.int = FALSE, shape = 19, alpha = 0.4) + # "#FFB900" "#0072B2"
  stat_cor(aes(color = Type), label.y.npc = 0.2, label.x.npc = 0.02, cor.coef.name = "rho", method = "spearman", size = 4.2) +
  labs(x = "Shannon Index", y = "Read depth (million)") +
  scale_color_manual(values=c("#008080", "#4B0082")) +
  scale_y_continuous(labels = scales::label_number(scale = 1e-6)) +
  facet_wrap(~.,) +
  #xlim(c(0,6.3)) +
  #ylim(0,60) +
  theme_classic() + 
  border(size = ) +
  theme(axis.text=element_text(size=11.5),
        axis.title = element_text(size = 12.5), 
        strip.background = element_rect(fill = "grey90", linewidth = 1),
        strip.text = element_text(size = 13),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.position = "none") 
dev.off()


##################################
# 2.2. Overall Beta diversity analysis (Mothers and Infants)
##################################

# Calculate dissimilarities in microbiome composition between samples
dist_vOTUs <- vegdist(t(Abundance_table), index = "bray")

# NMDS analysis
# Subset metadata to select phenotypes
NMDS_phenos <- Sample_metadata[Sample_metadata$NG_ID %in% colnames(Abundance_table),]
row.names(NMDS_phenos) <- NMDS_phenos$NG_ID
NMDS_phenos <- NMDS_phenos[,c("Type", "read_depth", "DNA_concentration_ng_ul", 
                              "shannon", "Timepoint_categorical")]

# Run NMDS in vegan (metaMDS)
NMS <- metaMDS(t(Abundance_table), distance = "bray", k=2)
en_all<- envfit(NMS, NMDS_phenos, permutations = 999, na.rm = TRUE)
en_all$factors 

centroids_all <- as.data.frame(scores(en_all, "factors"))
centroids_all$Type <- c(gsub('Type', '', row.names(centroids_all)))
centroids_all$Type[-c(1,2)] <- NA
centroids_all$Timepoint <- c(gsub('Timepoint_categorical', '', row.names(centroids_all)))
centroids_all$Timepoint[c(1,2)] <- NA

data.scores_all = as.data.frame(scores(NMS, "sites"))
data.scores_all <- merge(data.scores_all, NMDS_phenos, by='row.names')
row.names(data.scores_all) <- data.scores_all$Row.names
data.scores_all$Row.names <- NULL

# Calculate centroid and distances for each timepoint
data.scores_all <- data.scores_all %>%
  group_by(Timepoint_categorical) %>%
  mutate(centroid_NMDS1 = mean(NMDS1),
         centroid_NMDS2 = mean(NMDS2),
         distance_from_centroid = sqrt((NMDS1 - centroid_NMDS1)^2 + (NMDS2 - centroid_NMDS2)^2))

# Use IQR method to identify and remove outliers within each timepoint
data.scores_all <- data.scores_all %>%
  group_by(Timepoint_categorical) %>%
  mutate(Q1 = quantile(distance_from_centroid, 0.25),
         Q3 = quantile(distance_from_centroid, 0.75),
         IQR = Q3 - Q1,
         lower_bound = Q1 - 1.5 * IQR,
         upper_bound = Q3 + 1.5 * IQR) %>%
  filter(distance_from_centroid >= lower_bound & distance_from_centroid <= upper_bound) %>%
  ungroup()

# Outliers identified based on beta diversity analysis(based on IQR method):
outliers <- Sample_metadata$bioSampleId[which(!Sample_metadata$read_depth %in% data.scores_all$read_depth)]

# Generate NMDS plot with maternal (blue) and infant (green) samples
pdf('2_DIVERSITY/Plots/Viral_vOTUs_ALL_Bray_NMDS_without_outliers.pdf', width=5.5, height=4)
NMDS_plot_all <- ggplot(data = data.scores_all, aes(x = NMDS1, y = NMDS2, color=Type)) + 
  geom_point(size = 1.5, alpha=0.7) + 
  geom_point(data=centroids_all, aes(fill=Type),shape=NA, size=4, color='black') + 
  stat_ellipse(geom = "polygon", alpha = 0.0, aes(group = Type, color=Type), linetype = 2, linewidth = 1)+
  theme_bw()+
  scale_color_manual(values=c("#008080", "#4B0082")) +
  theme(axis.text=element_text(size=16),
        axis.title = element_text(size = 17), 
        panel.background = element_blank(), 
        panel.grid = element_blank(), 
        panel.border = element_rect(fill = NA, linewidth =1.2, colour = "grey30"), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.title = element_blank(), 
        legend.text = element_text(size = 16, colour = "grey30"),
        legend.position = "bottom") +
  guides(fill = "none")  
ggMarginal(NMDS_plot_all, type="densigram", groupFill=T)
dev.off()

##################################
# 2.3. NMDS of 1 dimension
##################################
# PERMANOVA is not suitable for longitudinal data
# Generate a NMDS that represents the virome composition in one dimension
NMS1 <- metaMDS(t(Abundance_table), distance = "bray", k=1)
data.scores1 <- as.data.frame(scores(NMS1 , "sites"))
data.scores1$NG_ID <- rownames(data.scores1)
colnames(data.scores1) <- c("Virome composition", "NG_ID")
rownames(data.scores1) <- NULL

# Add to the Sample_metadata as a variable "Virome composition"
Sample_metadata <- left_join(Sample_metadata, data.scores1 , by = "NG_ID")

# Check differences in virome composition between mothers and infants
MM_type_composition <- lmer(`Virome composition` ~ Type + DNA_concentration_ng_ul + read_depth + (1|CS_BABY_BIOME_ID),
                            REML = F, data = Sample_metadata)
summary(MM_type_composition)


##************************************************************************
# 3. Diversity analyses in Infants: Richness, alpha and beta diversity
#*************************************************************************

# Generate again sample metadata for infants (after estimating virome composition)
Sample_metadata_infants <- Sample_metadata[Sample_metadata$Type == "Infant",]

##################################
# 3.1. Richness and alpha diversity 
##################################

# Check statistical significance over time (time as continuous variable)
# A) Richness
MM_time_richness_infant <- lmer(richness ~ Timepoint_numeric + DNA_concentration_ng_ul + read_depth + (1|CS_BABY_BIOME_ID),
                                REML = F, data = Sample_metadata_infants)
summary(MM_time_richness_infant) 

# B) Shannon Index
MM_time_shannon_infant <- lmer(shannon ~ Timepoint_numeric + DNA_concentration_ng_ul + read_depth + (1|CS_BABY_BIOME_ID),
                               REML = F, data = Sample_metadata_infants)
summary(MM_time_shannon_infant) 

# Check statistical significance over time (per timepoint - time as discrete variable)
result_shannon_infant <- Timepoint_significance_analysis(data = Sample_metadata_infants,
                                                         response_variable = "shannon", 
                                                         time_variable = "Timepoint_categorical", 
                                                         group_variable = "CS_BABY_BIOME_ID")
result_richness_infant <- Timepoint_significance_analysis(data =  Sample_metadata_infants,
                                                          response_variable = "richness", 
                                                          time_variable = "Timepoint_categorical", 
                                                          group_variable = "CS_BABY_BIOME_ID")

# Check correlation of richness and alpha diversity with sequencing depth and DNA concentration
spearman_cor_richness_CONC_infants <- cor.test(Sample_metadata$DNA_concentration_ng_ul[Sample_metadata$Type == "Infant"], 
                                               Sample_metadata$richness[Sample_metadata$Type == "Infant"],
                                               method = "spearman") 
spearman_cor_richness_Read_DEPTH_infants <- cor.test(Sample_metadata$read_depth[Sample_metadata$Type == "Infant"],
                                                     Sample_metadata$richness[Sample_metadata$Type == "Infant"],
                                                     method = "spearman") 
spearman_cor_shannon_CONC_infants <- cor.test(Sample_metadata$DNA_concentration_ng_ul[Sample_metadata$Type == "Infant"],
                                              Sample_metadata$shannon[Sample_metadata$Type == "Infant"],
                                              method = "spearman")
spearman_cor_shannon_Read_DEPTH_infants <- cor.test(Sample_metadata$read_depth[Sample_metadata$Type == "Infant"],
                                                    Sample_metadata$shannon[Sample_metadata$Type == "Infant"],
                                                    method = "spearman") 


# Generate plots for richness and shannon index 
pdf('2_DIVERSITY/Plots/Shannon_Index_Infants.pdf', width=3.6, height=3.1)
ggplot(Sample_metadata[Sample_metadata$Type=="Infant" & !is.na(Sample_metadata$Timepoint_categorical) &
                         Sample_metadata$Timepoint_categorical %in% c("W01", "W02", "W03", 
                                                                      "W04", "W05", "W06"),]
       , aes(x=Timepoint_categorical, y=shannon)) + 
  geom_violin(aes(fill=Type), alpha=0.5, width=0.8, trim = F) +
  geom_jitter(alpha = 0.7, color ="#008080") +
  geom_boxplot(width=0.2, fill="white", color="black") + 
  labs(x = 'Timepoint', y = 'Shannon Index') + 
  scale_fill_manual(values=c("#008080")) + #"#008080", "#4B0082"
  theme_classic() + 
  #ylim(0,6) +
  theme(plot.title = element_text(size=18, hjust = 0.5, face="bold"), 
        axis.text = element_text(size=14),
        axis.title = element_text(size=16), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.position = "none")
dev.off()

pdf('2_DIVERSITY/Plots/Richness_Infants.pdf', width=3.8, height=3.1)
ggplot(Sample_metadata[Sample_metadata$Type=="Infant" & !is.na(Sample_metadata$Timepoint_categorical) &
                         Sample_metadata$Timepoint_categorical %in% c("W01", "W02", "W03", 
                                                                      "W04", "W05", "W06"),]
       , aes(x=Timepoint_categorical, y=richness)) +
  geom_violin(aes(fill=Type), alpha=0.5, width=0.8, trim = F) +
  geom_jitter(alpha = 0.7, color ="#008080") +
  geom_boxplot(width=0.2, fill="white", color="black") + 
  labs(x = 'Timepoint', y = 'Viral richness') + 
  scale_fill_manual(values=c("#008080")) + 
  theme_classic() + 
  theme(plot.title = element_text(size=18, hjust = 0.5, face="bold"), 
        axis.text = element_text(size=14),
        axis.title = element_text(size=16), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.position = "none") 
dev.off()


##################################
# 3.2. Beta diversity 
##################################

# Calculate dissimilarities in microbiome composition between infant samples
dist_vOTUs_infants <- vegdist(t(Abundance_table_infants), index = "bray")

## NDMS analysis
# Preparing phenotypes
NMDS_phenos_infants <- Sample_metadata_infants [Sample_metadata_infants$NG_ID %in% colnames(Abundance_table_infants),]
NMDS_phenos_infants <- NMDS_phenos_infants[,c("NG_ID", "read_depth", "DNA_concentration_ng_ul", "Timepoint_categorical")]
row.names(NMDS_phenos_infants) <- NMDS_phenos_infants$NG_ID
NMDS_phenos_infants$NG_ID <- NULL

# Run NMDS in vegan (metaMDS)
NMS_infants <- metaMDS(t(Abundance_table_infants), distance = "bray", k=2)
en_infants = envfit(NMS_infants, NMDS_phenos_infants, permutations = 1000, na.rm = TRUE)
en_infants$factors 

centroids_infants <- as.data.frame(scores(en_infants, "factors"))
centroids_infants$Timepoint_categorical <- c(gsub('Timepoint_categorical', '', row.names(centroids_infants)))

data.scores_infants <- as.data.frame(scores(NMS_infants, "sites"))
data.scores_infants <- merge(data.scores_infants, NMDS_phenos_infants, by='row.names')
row.names(data.scores_infants) <- data.scores_infants$Row.names
data.scores_infants$Row.names <- NULL

# Calculate centroid and distances for each timepoint
data.scores_infants <- data.scores_infants %>%
  group_by(Timepoint_categorical) %>%
  mutate(centroid_NMDS1 = mean(NMDS1),
         centroid_NMDS2 = mean(NMDS2),
         distance_from_centroid = sqrt((NMDS1 - centroid_NMDS1)^2 + (NMDS2 - centroid_NMDS2)^2))

# Use IQR method to identify and remove outliers within each timepoint
data.scores_infants <- data.scores_infants %>%
  group_by(Timepoint_categorical) %>%
  mutate(Q1 = quantile(distance_from_centroid, 0.25),
         Q3 = quantile(distance_from_centroid, 0.75),
         IQR = Q3 - Q1,
         lower_bound = Q1 - 1.5 * IQR,
         upper_bound = Q3 + 1.5 * IQR) %>%
  filter(distance_from_centroid >= lower_bound & distance_from_centroid <= upper_bound) %>%
  ungroup()

# Generate NMDS plot
pdf('2_DIVERSITY/Plots/Viral_vOTU_INFANTS_Bray_NMDS.pdf', width=3.6, height=3.7)
NMDS_plot_infants <- ggplot(data = data.scores_infants, aes(x = NMDS1, y = NMDS2, color=Timepoint_categorical)) + 
  geom_point(size = 1.5, alpha=0.7) + 
  geom_point(data=centroids_infants, aes(fill=Timepoint_categorical), shape=NA, size=4, color='black') + 
  stat_ellipse(geom = "polygon", alpha = 0.0, aes(group = Timepoint_categorical, color=Timepoint_categorical), linetype = 2, linewidth = 0.8)+
  scale_color_manual(values=c("#0055AA","#007ED3","#7FD2FF","#EAC862", "#FF9D1E", "#C40003")) +
  theme_bw()+
  theme(axis.text=element_text(size=13),
        axis.title = element_text(size = 14), 
        panel.background = element_blank(), 
        panel.grid = element_blank(),  
        panel.border = element_rect(fill = NA, linewidth =1.2, colour = "grey30"), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.title = element_blank(), 
        legend.text = element_text(size = 14, colour = "grey30"),
        legend.position = "bottom") +
  guides(fill = "none")  
ggMarginal(NMDS_plot_infants, type= "boxplot", groupFill=T)
dev.off()

##################################
# 3.3. NMDS of 1 dimension
##################################

# Generate a NMDS that represents the virome composition in one dimension (specific for infants)
NMS1 <- metaMDS(t(Abundance_table_infants), distance = "bray", k=1)
data.scores1 <- as.data.frame(scores(NMS1 , "sites"))
data.scores1$NG_ID <- rownames(data.scores1)
colnames(data.scores1) <- c("Virome composition infant", "NG_ID")
rownames(data.scores1) <- NULL

# Add to the Sample_metadata as a variable "Virome composition"
Sample_metadata_infants <- left_join(Sample_metadata_infants, data.scores1 , by = "NG_ID")

# Check differences in virome composition infants between timepoints
MM_type_composition_infants <- lmer(`Virome composition` ~ Timepoint_numeric + DNA_concentration_ng_ul + read_depth + (1|CS_BABY_BIOME_ID),
                                    REML = F, data = Sample_metadata_infants)
summary(MM_type_composition_infants) 


##************************************************************************
# 4. Virome specificity / stability over time
##************************************************************************

##################################
# 4.1. Prevalence analysis
##################################

# Estimate the number of vOTUs present in mothers and infants at different prevalence thresholds
N_samples_mothers <- ncol(Abundance_table_mothers)
N_samples_infants <- ncol(Abundance_table_infants)
Viral_presence_mothers <- rowSums(Abundance_table_mothers > 0)
Viral_presence_infants <- rowSums(Abundance_table_infants > 0)
Prevalence_mothers <- 100*(Viral_presence_mothers / N_samples_mothers)
Prevalence_infants <- 100*(Viral_presence_infants / N_samples_infants)
prevalence_thresholds <- c(5, 10)

# Count the number of viruses at each prevalence threshold
for (threshold in prevalence_thresholds) {
  num_viruses_at_threshold_mothers <- sum(Prevalence_mothers >= threshold)
  num_viruses_at_threshold_infants <- sum(Prevalence_infants >= threshold)
  cat(paste("Number of viruses with prevalence ≥ ", threshold, "%: in mothers", num_viruses_at_threshold_mothers, "\n"))
  cat(paste("Number of viruses with prevalence ≥ ", threshold, "%: in infants", num_viruses_at_threshold_infants, "\n"))
}

# Estimate number of viruses present at least in 2 individuals
vOTUs_multiple_infants <- c()
vOTUs_single_infant <- c()

for (row in 1:nrow(Abundance_table_infants)) {
  vOTU_name <- rownames(Abundance_table_infants)[row]
  samples_present <- colnames(Abundance_table_infants)[which(Abundance_table_infants[row,] !=0)]
  CS_ID_present <- Sample_metadata$CS_BABY_BIOME_ID[Sample_metadata$NG_ID %in% samples_present]
  if (sum(table(CS_ID_present) >= 1) >= 2) {
    vOTUs_multiple_infants<- c(vOTUs_multiple_infants, vOTU_name)
  } 
  if (sum(table(CS_ID_present) >= 1) == 1) {
    vOTUs_single_infant <- c(vOTUs_single_infant, vOTU_name)
  } 
}

# Estimate number of persistent vOTUs
persistent_vOTUs_3_timepoints <- c()
persistent_vOTUs_all_timepoints <- c()

for (row in 1:nrow(Abundance_table_infants)) {
  vOTU_name <- rownames(Abundance_table_infants)[row]
  samples_present <- colnames(Abundance_table_infants)[which(Abundance_table_infants[row,] !=0)]
  CS_ID_present <- Sample_metadata$CS_BABY_BIOME_ID[Sample_metadata$NG_ID %in% samples_present]
  if (sum(table(CS_ID_present) >= 3) >= 1) {
    persistent_vOTUs_3_timepoints <- c(persistent_vOTUs_3_timepoints, vOTU_name)
  } 
  if (sum(table(CS_ID_present) == 6) >= 1) {
    persistent_vOTUs_all_timepoints <- c(persistent_vOTUs_all_timepoints, vOTU_name)
  } 
}

# Add persistency results to Viral_metadata
Virus_metadata$Persistent <- ifelse(Virus_metadata$Virus_ID %in% persistent_vOTUs_3_timepoints, "Yes", "No")
Virus_metadata$Highly_persistent <- ifelse(Virus_metadata$Virus_ID %in% persistent_vOTUs_all_timepoints, "Yes", "No")

# Generate barplots
number_infants <- data.frame(
  Category = c("Multiple infants", "Single infant"),
  Proportion = c(11, 89)  
)

persistency <- data.frame(
  Timepoint = c("≥3 timepoints", "All timepoints"),
  Proportion = c(52.4, 9.9)  
)

pdf('2_DIVERSITY/Plots/vOTU_presence_multiple_infants.pdf', width=2.2, height=1.5)
ggplot(number_infants, aes(x = "", y = Proportion, fill = Category)) +
  geom_col() +  
  coord_flip() +
  scale_fill_manual(values = c("#7B9F80", "#E5D1D0")) + 
  theme_classic() +
  theme(axis.text=element_text(size=11),
        axis.title = element_text(size = 12),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "top")
dev.off()

pdf('2_DIVERSITY/Plots/vOTU_persistency.pdf', width=2.6, height=2.2)
ggplot(persistency, aes(x = Timepoint, y = Proportion, fill = Timepoint)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  coord_flip() +
  scale_fill_manual(values = c("#E0C9A8", "#C8D9BF")) + 
  theme_classic() +
  theme(axis.text=element_text(size=13),
        axis.title = element_text(size = 14),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "top")
dev.off()

##################################
# 4.2 DGR analysis
##################################

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


pdf('2_DIVERSITY/Plots/DGRs_mother_infant.pdf', width=4.2, height=3.2)
ggplot(DGR_origin_ggplot, aes(x = DGR, y = Proportion, fill = Group)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "DGR Presence", y = "Proportion", fill = "Origin") +
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


# Then, check if persistent or highly-persistent vOTUs are more likely to have DGRs

# Generate a contingency table 
contingency_table_persist <- table(Virus_metadata_infants$DGR, Virus_metadata_infants$Persistent)
contingency_table_high_persist <- table(Virus_metadata_infants$DGR, Virus_metadata_infants$Highly_persistent)

# Perform Fisher's exact test
fisher_test_result_persist <- fisher.test(contingency_table_persist)
fisher_test_result_high_persist <- fisher.test(contingency_table_high_persist)
                                  
                                  
##################################
# 4.3 Intra vs Interindividual distance comparison
##################################

my_pseudocount_normal <- min(Abundance_table_infants[Abundance_table_infants!=0])/2
aitchison_distances_infants_raw <- vegdist(t(Abundance_table_infants), method = "aitchison", pseudocount = my_pseudocount_normal)

# Transform distance matrix into a long format
distances_infants <- reshape2::melt(as.matrix(aitchison_distances_infants_raw))
distances_infants <- distances_infants[distances_infants$Var1 != distances_infants$Var2, ] # remove distances between the same sample
colnames(distances_infants) <- c("Sample1", "Sample2", "Distance")

# Add metadata to the long-formatted table
distances_infants$Timepoint1 <- Sample_metadata_infants$Timepoint_categorical[match(distances_infants$Sample1, Sample_metadata_infants$NG_ID)]
distances_infants$Timepoint2 <-  Sample_metadata_infants$Timepoint_categorical[match(distances_infants$Sample2, Sample_metadata_infants$NG_ID)]
distances_infants$CS_ID1 <- Sample_metadata_infants$CS_BABY_BIOME_ID[match(distances_infants$Sample1, Sample_metadata_infants$NG_ID)]
distances_infants$CS_ID2 <- Sample_metadata_infants$CS_BABY_BIOME_ID[match(distances_infants$Sample2, Sample_metadata_infants$NG_ID)]
distances_infants$AB1 <- Sample_metadata_infants$rand_AB[match(distances_infants$Sample1, Sample_metadata_infants$NG_ID)]
distances_infants$AB2 <-  Sample_metadata_infants$rand_AB[match(distances_infants$Sample2, Sample_metadata_infants$NG_ID)]
distances_infants$Type <- rep("infant", nrow(distances_infants))
distances_infants$Within_between <- ifelse(distances_infants$CS_ID1 == distances_infants$CS_ID2, "Within", "Between")

# Add a combined_timepoint variable, avoiding different levels for the same timepoint comparison (e.g. P12-B and B-P12 always represented as 1-3)
mapping <- c(W01 = 1, W02 = 2, W03 = 3, W04 = 4, W05 = 5, W06 = 6) # Change timepoints to values from 1 to 7
distances_infants$Timepoint1 <- mapping[as.character(distances_infants$Timepoint1)]
distances_infants$Timepoint2 <- mapping[as.character(distances_infants$Timepoint2)]
distances_infants$Combined_timepoint <- mapply(combine_timepoints, distances_infants$Timepoint1, distances_infants$Timepoint2) # Apply the function

mapping_back <- c("1-1" = "W01-W01", "1-2" = "W01-W02", "1-3" = "W01-W03", 
                  "1-4" = "W01-W04", "1-5" = "W01-W05", "1-6" = "W01-W06",
                  "2-2" = "W02-W02", "2-3" = "W02-W03", "2-4" = "W02-W04",
                  "2-5" = "W02-W05", "2-6" = "W02-W06",
                  "3-3" = "W03-W03", "3-4" = "W03-W04", "3-5" = "W03-W05",
                  "3-6" = "W03-W06",
                  "4-4" = "W04-W04", "4-5" = "W04-W05", "4-6" = "W04-W06",
                  "5-5" = "W05-W05", "5-6" = "W05-W06",
                  "6-6" = "W06-W06")

# Map back the values of the combined timepoint
distances_infants$Combined_timepoint <- mapping_back[distances_infants$Combined_timepoint] 

# Add a combined AB treatment variable
distances_infants$Combined_AB <- factor(mapply(combine_timepoints, distances_infants$AB1, distances_infants$AB2)) # Apply the function


# Select: within- and between-infants from different timepoints (will be ploted together)
distances_infants_across_timepoints <- distances_infants[!distances_infants$Combined_timepoint %in% 
                                                           c("W01-W01","W02-W02", "W03-W03","W04-W04", 
                                                             "W05-W05", "W06-W06"),]

# Test if W01-W06 intra-individual distances are smaller than between-individual distances
# For this we use a permutation based test (10,000 permutations)
W01_W06_within_distances <- distances_infants_across_timepoints[
  distances_infants_across_timepoints$Within_between == "Within" & 
    distances_infants_across_timepoints$Combined_timepoint == "W01-W06", 
  "Distance"]

Between_distances <- distances_infants_across_timepoints[
  distances_infants_across_timepoints$Within_between == "Between", 
  "Distance"]

# Create a combined vector with all distances and a label vector (1 = Within, 0 = Between)
all_distances <- c(W01_W06_within_distances, Between_distances)
labels <- c(rep(1, length(W01_W06_within_distances)), rep(0, length(Between_distances)))
 
# Define the observed difference in medians
observed_stat <- median(W01_W06_within_distances) - median(Between_distances) 

# Set the number of permutationsand perform the test
num_permutations <- 1000
permuted_stats <- numeric(num_permutations)
set.seed(123) 
for (i in 1:num_permutations) {
  permuted_labels <- sample(labels) 
  within_permuted <- all_distances[permuted_labels == 1]
  between_permuted <- all_distances[permuted_labels == 0]
  permuted_stats[i] <- median(within_permuted) - median(between_permuted) 
}

# Calculate the p-value: proportion of permuted stats smaller than or equal to observed (one-sided hypothesis)
p_value <- mean(permuted_stats <= observed_stat) 

# Generate the summarized tables for the plots
distances_infants_across_timepoints_summ <- distances_infants_across_timepoints %>%
  drop_na(Combined_timepoint) %>%
  group_by(Combined_timepoint, Within_between) %>%
  summarize(
    median_Dissimilarity = median(Distance),
    q1 = quantile(Distance, 0.25),
    q3 = quantile(Distance, 0.75)
  )

# Generate lineplot of mean Aitchison distances of different individuals over time (per timepoint comparison)
distances_infants_across_timepoints_summ_plot <- distances_infants_across_timepoints_summ[grep("W01", distances_infants_across_timepoints_summ$Combined_timepoint),]
distances_infants_across_timepoints_summ_plot$Combined_timepoint <- substr(distances_infants_across_timepoints_summ_plot$Combined_timepoint, 5, 7)
distances_infants_across_timepoints_summ_plot$Combined_timepoint <- factor(distances_infants_across_timepoints_summ_plot$Combined_timepoint,
                                                                           levels = c("W01","W02", "W03","W04", 
                                                                                      "W05", "W06"))

pdf('2_DIVERSITY/Plots/W1_Aitchison_distances_across_timepoints_withinbetween_infants_W1-W6.pdf', width=4.2, height=3.6)
ggplot(distances_infants_across_timepoints_summ_plot,
       aes(x = Combined_timepoint, y = median_Dissimilarity, group = Within_between, color = Within_between)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = q1, ymax = q3), width = 0.2, linetype = 5) +
  geom_ribbon(aes(ymin = q1, ymax = q3, fill = Within_between), alpha = 0.2, linetype = 0) +
  labs(x = 'Timepoint comparison', y = 'Aitchison Distance') +
  scale_color_manual(values=c("#333333","#F0E442")) +
  scale_fill_manual(values=c("#333333","#F0E442")) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    axis.title.y = element_text(margin = margin(r = 15)),
    axis.title.x = element_text(margin = margin(t = 13)), 
    legend.text = element_text(size = 16, colour = "grey30"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top", 
    legend.title = element_blank(),
  )
dev.off()


 #****************
 # Write results
 #****************
# Sample metadata
write.table(Sample_metadata,"Metadata_CS/CS_Baby_Biome_Sample_Metadata_17052024.txt", sep = "\t", row.names = F, quote = FALSE) 
write.table(Sample_metadata_mothers,"Metadata_CS/CS_Baby_Biome_Sample_Metadata_Mothers_17052024.txt", sep = "\t", row.names = F, quote = FALSE) 
write.table(Sample_metadata_infants,"Metadata_CS/CS_Baby_Biome_Sample_Metadata_Infants_17052024.txt", sep = "\t", row.names = F, quote = FALSE) 
 
# Abundance table
write.table(Abundance_table,"Abundance_table/CS_Abundance_Table_17052024.txt", sep = "\t", row.names = T, quote = FALSE)
write.table(Abundance_table_infants,"Abundance_table/CS_Abundance_Table_INFANTS_17052024.txt", sep = "\t", row.names = T, quote = FALSE)

# Virus metadata
write.table(Virus_metadata,"Metadata_CS/CS_Baby_Biome_Viral_Metadata_17052024.txt", sep = "\t", row.names = F, quote = FALSE) 

