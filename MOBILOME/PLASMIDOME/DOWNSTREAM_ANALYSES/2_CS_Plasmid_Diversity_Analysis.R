################################################################################
##### CS Baby Biome: Plasmid Diversity analysis
### Author(s): Asier Fernández-Pato
### Last updated: 24th February, 2025
################################################################################

#****************
# Load libraries
#****************
library(vegan)
library(lme4)
library(lmerTest)
library(multcomp)
library(RLRsim)
library(ggplot2)
library(ggrepel)
library(ggExtra)
library(ggpubr)
library(reshape2)
library(tidyr)
library(nnet)


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

# Set working directory 
setwd("~/Desktop/PhD/Projects/CS Baby Biome/Final_analysis/PLASMIDS/")

# Read abundance table and metadata tables
Sample_metadata <- read.delim("METADATA_TABLES/CS_Baby_Biome_Sample_Metadata_25062024.txt") 
Sample_metadata_infants <- read.delim("METADATA_TABLES/CS_Baby_Biome_Sample_Metadata_Infants_25062024.txt")
Plasmid_metadata <- read.delim("METADATA_TABLES/CS_Baby_Biome_Plasmid_Metadata_25062024.txt")
Plasmid_metadata_infants <- read.delim("METADATA_TABLES/CS_Baby_Biome_Plasmid_Metadata_Infants_25062024.txt")
Plasmid_abundance <- read.delim("ABUNDANCE_TABLES/CS_Abundance_Table_25062024.txt")

##************************************************************************
# 0. Preprocessing
#*************************************************************************

# Exclude infants of M6 and M12 from the metadata
Sample_metadata <- Sample_metadata[!Sample_metadata$Timepoint_categorical %in% c("M06", "M12"),]
Plasmid_abundance <- Plasmid_abundance[, colnames(Plasmid_abundance) %in% Sample_metadata$NG_ID]

# Remove plasmids that are not present after sample removal from Abundance table and plasmid metadata
Plasmid_abundance <- Plasmid_abundance[rowSums(Plasmid_abundance != 0) > 0, ] #3,348
rownames(Plasmid_metadata) <- Plasmid_metadata$Plasmid_ID
Plasmid_metadata <- Plasmid_metadata[rownames(Plasmid_abundance), ]

# Set Timepoint and Type variable from metadata as factors
Sample_metadata$Timepoint_categorical <- factor(Sample_metadata$Timepoint_categorical, 
                                                levels=c("MOM", "W01", "W02", "W03","W04",
                                                         "W05", "W06"))
Sample_metadata$Type <- factor(Sample_metadata$Type)


# Add plasmid abundance as variable to the Sample metadata
Sample_metadata$Plasmid_ab <- colSums(Plasmid_abundance)

##************************************************************************
# 1. Exploratory analyses
#*************************************************************************

# Generate the metadata and abundance tables for mothers and infants
# Estimate the number of pTUs present in maternal and infant samples
Sample_metadata_infants <- Sample_metadata[Sample_metadata$Type == "Infant",]
Plasmid_abundance_infants <- Plasmid_abundance[,Sample_metadata$Type == "Infant"]
Plasmid_abundance_infants <- Plasmid_abundance_infants[rowSums(Plasmid_abundance_infants)>0,] #1,533 plasmids

Sample_metadata_mothers <- Sample_metadata[Sample_metadata$Type == "Mother",]
Plasmid_abundance_mothers <- Plasmid_abundance[,Sample_metadata$Type == "Mother"]
Plasmid_abundance_mothers <- Plasmid_abundance_mothers[rowSums(Plasmid_abundance_mothers)>0,] #2,135 plasmids

# Estimate the summary statistics for the plasmid abundance mothers and infants
summary_stats(Sample_metadata$Plasmid_ab [Sample_metadata$Timepoint_categorical =="MOM"])
summary_stats(Sample_metadata$Plasmid_ab [Sample_metadata$Timepoint_categorical !="MOM"])


##************************************************************************
# 2. Overall diversity analyses (mother vs infant)
#*************************************************************************

##################################
# 2.1. Estimation of alpha diversity and richness 
##################################

# Estimate shannon diversity and richness (also per mobility type and topology)
Sample_metadata$plasmid_richness  <-specnumber(t(Plasmid_abundance))
Sample_metadata$plasmid_shannon <- vegan::diversity(t(Plasmid_abundance), index = "shannon")

# Test significance: plasmid richness/diversity/abundance in mothers vs babies
# A) Richness
MM_type_richness <- lmer(plasmid_richness ~ Type + DNA_concentration_ng_ul + read_depth + (1|CS_BABY_BIOME_ID),
                         REML = F,data = Sample_metadata)
summary(MM_type_richness) 

# B) Shannon Index
MM_type_shannon <- lmer(plasmid_shannon ~ Type + DNA_concentration_ng_ul + read_depth + (1|CS_BABY_BIOME_ID),
                        REML = F, data = Sample_metadata)
summary(MM_type_shannon) 

# Plot the richness and diversity in mothers and babies
pdf('2_DIVERSITY/PLOTS/Plasmid_Shannon_Index_comp.pdf', width=2.5, height=3.2)
ggplot(Sample_metadata[!is.na(Sample_metadata$Timepoint_categorical),], aes(x=Type, y=plasmid_shannon)) +
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

pdf('2_DIVERSITY/PLOTS/Plasmid_Richness_comp.pdf', width=2.9, height=3.2)
ggplot(Sample_metadata[!is.na(Sample_metadata$Timepoint_categorical),], aes(x=Type, y=plasmid_richness)) +
  geom_violin(aes(fill=Type), alpha=0.5, width=0.8, trim = F) +
  geom_jitter(aes(color=Type), alpha=0.7) +  # Added color aesthetic
  geom_boxplot(width=0.2, fill="white", color="black") + 
  labs(x = 'Sample type', y = 'Plasmid Richness') + 
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

##************************************************************************
# 2. Diversity analyses in Infants: Richness, alpha and beta diversity
#*************************************************************************

# Generate again sample metadata for infants (after estimating diversity)
Sample_metadata_infants <- Sample_metadata[Sample_metadata$Type == "Infant",]

##################################
# A. Richness and alpha diversity 
##################################

# Check statistical significance over time (time as continuous variable)
# A) Richness
MM_time_richness_infant <- lmer(plasmid_richness ~ Timepoint_numeric + DNA_concentration_ng_ul + read_depth + (1|CS_BABY_BIOME_ID),
                                REML = F, data = Sample_metadata_infants)
summary(MM_time_richness_infant)  

# B) Shannon Index
MM_time_shannon_infant <- lmer(plasmid_shannon ~ Timepoint_numeric + DNA_concentration_ng_ul + read_depth + (1|CS_BABY_BIOME_ID),
                               REML = F, data = Sample_metadata_infants)
summary(MM_time_shannon_infant) 

# Generate plots for richness and shannon index 

# A) Per timepoint
pdf('2_DIVERSITY/PLOTS/Plasmid_Shannon_Infants.pdf', width=3.6, height=3.1)
ggplot(Sample_metadata_infants
       , aes(x=Timepoint_categorical, y=plasmid_shannon)) + 
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

pdf('2_DIVERSITY/PLOTS/Plasmid_Richness_Infants.pdf', width=3.6, height=3.1)
ggplot(Sample_metadata_infants
       , aes(x=Timepoint_categorical, y=plasmid_richness)) + 
  geom_violin(aes(fill=Type), alpha=0.5, width=0.8, trim = F) +
  geom_jitter(alpha = 0.7, color ="#008080") +
  geom_boxplot(width=0.2, fill="white", color="black") + 
  labs(x = 'Timepoint', y = 'Plasmid richness') + 
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


##################################
# B. Beta diversity 
##################################

# Calculate dissimilarities in plasmid composition between infant samples
dist_plasmids_infants <- vegdist(t(Plasmid_abundance_infants), index = "bray")

## NDMS analysis
# Preparing phenotypes
NMDS_phenos_infants <- Sample_metadata_infants [Sample_metadata_infants$NG_ID %in% colnames(Plasmid_abundance_infants),]
NMDS_phenos_infants <- NMDS_phenos_infants[,c("NG_ID", "read_depth", "DNA_concentration_ng_ul",
                                              "Timepoint_categorical", "rand_AB", "feeding_mode_pragmatic")]
row.names(NMDS_phenos_infants) <- NMDS_phenos_infants$NG_ID
NMDS_phenos_infants$NG_ID <- NULL

# Run NMDS in vegan (metaMDS)
NMS_infants <- metaMDS(t(Plasmid_abundance_infants), distance = "bray", k=2)
en_infants = envfit(NMS_infants, NMDS_phenos_infants, permutations = 1000, na.rm = TRUE)
en_infants$factors 

centroids_infants <- as.data.frame(scores(en_infants, "factors"))
centroids_infants$Timepoint_categorical <- c(gsub('Timepoint_categorical', '', row.names(centroids_infants)))
centroids_infants$Timepoint_categorical <- c(c("W01", "W02", "W03", "W04", "W05", "W06"), rep(NA,5))
centroids_infants$rand_AB <-  c(rep(NA,6), c("No", "Yes"), rep(NA,3))
centroids_infants$Feeding_mode <- c(rep(NA,8), c("Breastmilk", "Formula", "Mixed"))

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

# Generate NMDS plot: according to timepoint
pdf('2_DIVERSITY/PLOTS/Plasmid_INFANTS_Timepoint_Bray_NMDS.pdf', width=3.6, height=3.7)
NMDS_plot_infants <- ggplot(data = data.scores_infants, aes(x = NMDS1, y = NMDS2, color=Timepoint_categorical)) + 
  geom_point(size = 1.5, alpha=0.7) + 
  geom_point(data=centroids_infants, aes(fill=Timepoint_categorical), shape=NA, size=4, color='black') + 
  stat_ellipse(geom = "polygon", alpha = 0.0, aes(group = Timepoint_categorical, color=Timepoint_categorical), linetype = 2, linewidth = 0.8)+
  scale_color_manual(values=c("#0055AA","#007ED3","#7FD2FF","#EAC862", "#FF9D1E", "#C40003"))+
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

# Generate NMDS plot: according to AB group
pdf('2_DIVERSITY/PLOTS/Plasmid_INFANTS_AB_Bray_NMDS.pdf', width=3.6, height=3.7)
NMDS_plot_infants <- ggplot(data = data.scores_infants, aes(x = NMDS1, y = NMDS2, color=rand_AB)) + 
  geom_point(size = 1.5, alpha=0.7) + 
  geom_point(data=centroids_infants, aes(fill=rand_AB), shape=NA, size=4, color='black') + 
  stat_ellipse(geom = "polygon", alpha = 0.0, aes(group = rand_AB, color=rand_AB), linetype = 2, linewidth = 0.8)+
  scale_color_manual(values=c("#0055AA","#C40003","#EAC862", "#7FD2FF",
                              "#007ED3", "#FF9D1E", "#FFACAA", "#BADA55")) +
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

# Generate NMDS plot: according to feeding mode
pdf('2_DIVERSITY/PLOTS/Plasmid_INFANTS_Feeding_Bray_NMDS.pdf', width=3.6, height=3.7)
NMDS_plot_infants <- ggplot(data = data.scores_infants, aes(x = NMDS1, y = NMDS2, color=feeding_mode_pragmatic)) + 
  geom_point(size = 1.5, alpha=0.7) + 
  geom_point(data=centroids_infants, aes(fill=Feeding_mode), shape=NA, size=4, color='black') + 
  stat_ellipse(geom = "polygon", alpha = 0.0, aes(group = feeding_mode_pragmatic, color=feeding_mode_pragmatic),
               linetype = 2, linewidth = 0.8)+
  scale_color_manual(values=c("#0055AA","#C40003","#EAC862", "#7FD2FF",
                              "#007ED3", "#FF9D1E", "#FFACAA", "#BADA55")) +
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
# 1.3. NMDS of 1 dimension
##################################

# Generate a NMDS that represents the plasmid composition in one dimension (specific for infants)
NMS1 <- metaMDS(t(Plasmid_abundance_infants), distance = "bray", k=1)
data.scores1 <- as.data.frame(scores(NMS1 , "sites"))
data.scores1$NG_ID <- rownames(data.scores1)
colnames(data.scores1) <- c("Plasmid composition infant", "NG_ID")
rownames(data.scores1) <- NULL

# Add to the Sample_metadata as a variable "Plasmid composition"
Sample_metadata_infants <- left_join(Sample_metadata_infants, data.scores1 , by = "NG_ID")


##################################
# B.2. Beta diversity TCAM
##################################

# Filter the infant plasmid abundance table to retain those plasmids present in at least 2 infants
Plasmid_abundance_infants_TCAM <- Plasmid_abundance_infants[rowSums(Plasmid_abundance_infants != 0) >= 2, ] # From 1,533 to 1,287 plasmids


##************************************************************************
# 4. Plasmid specificity / stability over time
##************************************************************************

##################################
# 4.1. Prevalence analysis
##################################

# Estimate number of plasmids present at least in 2 individuals
PTUs_multiple_infants <- c()
PTUs_single_infant <- c()

for (row in 1:nrow(Plasmid_abundance_infants)) {
  PTU_name <- rownames(Plasmid_abundance_infants)[row]
  samples_present <- colnames(Plasmid_abundance_infants)[which(Plasmid_abundance_infants[row,] !=0)]
  CS_ID_present <- Sample_metadata$CS_BABY_BIOME_ID[Sample_metadata$NG_ID %in% samples_present]
  if (sum(table(CS_ID_present) >= 1) >= 2) {
    PTUs_multiple_infants<- c(PTUs_multiple_infants, PTU_name)
  } 
  if (sum(table(CS_ID_present) >= 1) == 1) {
    PTUs_single_infant <- c(PTUs_single_infant, PTU_name)
  } 
}

# Estimate number of persistent PTUs
persistent_PTUs_3_timepoints <- c()
persistent_PTUs_all_timepoints <- c()

for (row in 1:nrow(Plasmid_abundance_infants)) {
  PTU_name <- rownames(Plasmid_abundance_infants)[row]
  samples_present <- colnames(Plasmid_abundance_infants)[which(Plasmid_abundance_infants[row,] !=0)]
  CS_ID_present <- Sample_metadata$CS_BABY_BIOME_ID[Sample_metadata$NG_ID %in% samples_present]
  if (sum(table(CS_ID_present) >= 3) >= 1) {
    persistent_PTUs_3_timepoints <- c(persistent_PTUs_3_timepoints, PTU_name)
  } 
  if (sum(table(CS_ID_present) == 6) >= 1) {
    persistent_PTUs_all_timepoints <- c(persistent_PTUs_all_timepoints, PTU_name)
  } 
}

# Add persistency results to Viral_metadata
Plasmid_metadata$Persistent <- ifelse(Plasmid_metadata$Plasmid_ID %in% persistent_PTUs_3_timepoints, "Yes", "No")
Plasmid_metadata$Highly_persistent <- ifelse(Plasmid_metadata$Plasmid_ID %in% persistent_PTUs_all_timepoints, "Yes", "No")


#Test if PTU persistency is associated with mobility or circularity 
chi_square_topology_persistency <- chisq.test(table(Plasmid_metadata$Highly_persistent, Plasmid_metadata$topology_simple))
chi_square_mobility_persistency <- chisq.test(table(Plasmid_metadata$Highly_persistent, Plasmid_metadata$Mobility))

contingency_table_topology_persistency <- table(Plasmid_metadata$Highly_persistent, Plasmid_metadata$topology_simple)
prop_table_topology_persistency <- prop.table(contingency_table_topology_persistency, margin = 2)  

contingency_table_mobility_persistency <- table(Plasmid_metadata$Highly_persistent, Plasmid_metadata$Mobility)
prop_table_mobility_persistency <- prop.table(contingency_table_mobility_persistency, margin = 2)  


# Generate barplots
number_infants <- data.frame(
  Category = c("Multiple infants", "Single infant"),
  Proportion = c(36.5, 63.5)  
)

persistency <- data.frame(
  Timepoint = c("≥3 timepoints", "All timepoints"),
  Proportion = c(62.2, 15)  
)

pdf('2_DIVERSITY/Plots/PTUs_presence_multiple_infants.pdf', width=2.2, height=1.5)
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

pdf('2_DIVERSITY/Plots/PTU_persistency.pdf', width=2.6, height=2.2)
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
# 4.2. Association of persistance with plasmid traits
##################################

Plasmid_metadata_infants <- Plasmid_metadata_infants[,c("Plasmid_ID", "Persistent", "Mobility", "topology_simple")]

Plasmid_abundance_infants$mean_abundance <- rowMeans(Plasmid_abundance_infants, na.rm = TRUE)
Plasmid_abundance_infants$Plasmid_ID <- rownames(Plasmid_abundance_infants)
Plasmid_abundance_infants <- Plasmid_abundance_infants[,c("Plasmid_ID", "mean_abundance")]
data_merged <- merge(Plasmid_metadata_infants, Plasmid_abundance_infants, by = "Plasmid_ID")
data_merged$Persistent_binary <- ifelse(data_merged$Persistent == "Yes", 1, 0)

# Fit logistic regression models
model_mob <- multinom(Persistent ~ Mobility + mean_abundance, data = data_merged)
summary(model_mob)
z_values <- coef(model_mob) / summary(model_mob)$coefficients[, "Std. Error"]
p_values <- 2 * pnorm(-abs(z_values))
p_values

model_circ <- glm(Persistent_binary ~ topology_simple + mean_abundance, data = data_merged, family = binomial)
summary(model_circ)

#****************
# Save output
#****************

# Save the abundance tables for TCAM analysis
write.table(Plasmid_abundance_infants_TCAM,"2_DIVERSITY/Plasmid_abundance_min2infants.txt", sep = "\t", 
            row.names = T, col.names = T, quote = FALSE)

# Save metadata and abundance tables
write.table(Sample_metadata,"METADATA_TABLES/CS_Baby_Biome_Sample_Metadata_25062024.txt", sep = "\t", 
            row.names = T, col.names = T, quote = FALSE)
write.table(Sample_metadata_infants,"METADATA_TABLES/CS_Baby_Biome_Sample_Metadata_Infants_25062024.txt", sep = "\t", 
            row.names = T, col.names = T, quote = FALSE)
write.table(Plasmid_abundance_infants,"ABUNDANCE_TABLES/CS_Baby_Biome_Plasmid_Infants_Abundance_Table.txt", sep = "\t", 
            row.names = T, col.names = T, quote = FALSE)
write.table(Plasmid_metadata,"METADATA_TABLES/CS_Baby_Biome_Plasmid_Metadata_25062024.txt", sep = "\t", 
            row.names = T, col.names = T, quote = FALSE)
write.table(Plasmid_metadata_infants,"METADATA_TABLES/CS_Baby_Biome_Plasmid_Metadata_Infants_25062024.txt", sep = "\t", 
            row.names = T, col.names = T, quote = FALSE)




