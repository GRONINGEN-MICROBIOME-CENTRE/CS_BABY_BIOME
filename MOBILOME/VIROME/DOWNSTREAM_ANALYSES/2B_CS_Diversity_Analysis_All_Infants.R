################################################################################
##### CS Baby Biome: vOTU diversity Analysis (early vs late infancy)
### Author(s):Asier Fernández-Pato
### Last updated: 6th June, 2024
################################################################################


#****************
# Load modules
#****************
library(vegan)
library(ggplot2)
library(ggrepel)
library(ggExtra)

# Set working directory 
setwd("~/Desktop/PhD/Projects/CS Baby Biome/Final_analysis/VIRUSES/")

# Read abundance table and metadata tables
Abundance_table <- read.delim("Abundance_table/CS_Abundance_Table_17052024.txt")
Sample_metadata <- read.delim("Metadata_CS/CS_Baby_Biome_Sample_Metadata_17052024.txt")
Virus_metadata <- read.delim("Metadata_CS/CS_Baby_Biome_Viral_Metadata_17052024.txt")


##************************************************************************
# 0. Preprocessing
#*************************************************************************
# Set Timepoint and Type variable from metadata as factors
Sample_metadata$Timepoint_categorical <- factor(Sample_metadata$Timepoint_categorical, 
                                                levels=c("MOM", "W01", "W02", "W03","W04",
                                                         "W05", "W06", "M06", "M12"))
Sample_metadata$Type <- factor(Sample_metadata$Type)

# Add a variable grouping infant samples into early or late infancy
Sample_metadata$Early_late_infancy <- ifelse(
  Sample_metadata$Timepoint_categorical %in% 
    c("W01", "W02", "W03", "W04", "W05", "W06"),"Early",
  ifelse(Sample_metadata$Timepoint_categorical == "MOM", NA, "Late")
)

##************************************************************************
# 1. Beta diversity: PERMANOVA analysis
#*************************************************************************

###############
# Infancy
###############
# Generate the sample metadata  and abundance table for infants 
Sample_metadata_infants <- Sample_metadata[Sample_metadata$Type == "Infant",]
Abundance_table_infants <- Abundance_table[,Sample_metadata$Type == "Infant"]

# Calculate dissimilarities in microbiome composition between samples
dist_vOTUs_infants <- vegdist(t(Abundance_table_infants), index = "bray")

##************************************************************************
# 2. Beta diversity: NMDS analysis
#*************************************************************************

# Subset metadata to select phenotypes for NMDS
NMDS_phenos <- Sample_metadata[Sample_metadata$NG_ID %in% colnames(Abundance_table),]
row.names(NMDS_phenos) <- NMDS_phenos$NG_ID
NMDS_phenos <- NMDS_phenos[,c("read_depth", "DNA_concentration_ng_ul", "Type", "Timepoint_numeric",
                              "Timepoint_categorical", "Early_late_infancy")]

# Run NMDS (metaMDS)
NMS <- metaMDS(t(Abundance_table), distance = "bray", k=2)
en_all<- envfit(NMS, NMDS_phenos, permutations = 1000, na.rm = TRUE)

centroids_all <- as.data.frame(scores(en_all, "factors"))
centroids_all$Type <- c(gsub('Type', '', row.names(centroids_all)))
centroids_all$Type[-1] <- NA
centroids_all$Timepoint <- c(gsub('Timepoint_categorical', '', row.names(centroids_all)))
centroids_all$Timepoint[-c(2:9)] <- NA
centroids_all$Infancy <- c(gsub('Early_late_infancy', '', row.names(centroids_all)))
centroids_all$Infancy[-c(10:11)] <- NA

phenos.scrs_all <- as.data.frame(scores(en_all, display = "vectors"))
phenos.scrs_all <- cbind(phenos.scrs_all, phenos = rownames(phenos.scrs_all))
phenos.scrs_all$Phenos <- c("Number of reads", 'DNA concentration', 'Age (weeks) ')

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

# Generate NMDS plot with infant samples coloured per timepoint (maternal samples grouped in one color)
pdf('2_DIVERSITY/Plots/Viral_vOTUs_ALL_Bray_NMDS_without_outliers_timepoint_W1M12.pdf', width=4.5, height=4)
NMDS_plot_all_timepoint <- ggplot(data = data.scores_all, aes(x = NMDS1, y = NMDS2, color=Timepoint_categorical)) + 
  geom_point(size = 1.5, alpha=0.7) + 
  geom_point(data=centroids_all, aes(fill=Timepoint),shape=NA, size=4, color='black') + 
  stat_ellipse(geom = "polygon", alpha = 0.0, 
               aes(group = Timepoint_categorical, color=Timepoint_categorical), linetype = 2, linewidth = 1) +
  geom_segment(data = phenos.scrs_all["Timepoint_numeric",],
  aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
  arrow = arrow(length = unit(0.2, "cm")), colour = "grey30", size =0.7) +
  geom_label_repel(data = phenos.scrs_all["Timepoint_numeric",], 
                   aes(x = NMDS1, y = NMDS2, label = "Age (weeks)"), color='black', alpha = 0.8, size = 4.5) +
  theme_bw()+
  scale_color_manual(values=c("#4B0082", "#C40003", "#EAC862", "#7FD2FF", "#007ED3", "#FF9D1E", "#3366CC", "#FFACAA", "#BADA55")) +
  theme(axis.text=element_text(size=14.5),
        axis.title = element_text(size = 15.5), 
        panel.background = element_blank(), 
        panel.grid = element_blank(),  
        panel.border = element_rect(fill = NA, linewidth =1.2, colour = "grey30"), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.title = element_blank(), 
        legend.text = element_text(size = 15, colour = "grey30"),
        legend.position = "bottom") +
  guides(fill = "none")  
ggMarginal(NMDS_plot_all_timepoint, type="boxplot", groupFill=T)
dev.off()

# Generate NMDS plot with infant samples grouped by early or late infancy 
data.scores_all_filtered <- data.scores_all[!is.na(data.scores_all$Early_late_infancy),]

pdf('2_DIVERSITY/Plots/Viral_vOTUs_Infants_Bray_NMDS_without_outliers_Early_Late_Infancy.pdf', width=4.2, height=4)
NMDS_plot_infancy <- ggplot(data = data.scores_all_filtered, 
                            aes(x = NMDS1, y = NMDS2, color=Early_late_infancy)) + 
  geom_point(size = 1.5, alpha=0.7) + 
  geom_point(data=centroids_all, aes(fill=Infancy),shape=NA, size=4, color='black') + 
  stat_ellipse(geom = "polygon", alpha = 0.0, 
               aes(group = Early_late_infancy, color=Early_late_infancy), linetype = 2, linewidth = 1)+
  theme_bw()+
  scale_color_manual(values=c("#C40003", "#A3C55A")) +
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
ggMarginal(NMDS_plot_infancy, type="densigram", groupFill=T)
dev.off()

# Generate a NMDS that represents the virome composition in one dimension (specific for infants)
NMS1 <- metaMDS(t(Abundance_table_infants), distance = "bray", k=1)
data.scores1 <- as.data.frame(scores(NMS1 , "sites"))
data.scores1$NG_ID <- rownames(data.scores1)
colnames(data.scores1) <- c("Virome composition infant", "NG_ID")
rownames(data.scores1) <- NULL

# Add to the Sample_metadata as a variable "Virome composition"
Sample_metadata_infants <- left_join(Sample_metadata_infants, data.scores1 , by = "NG_ID")

# Check differences in virome composition infants between timepoints
MM_type_composition_infants_early_late <- lmer(`Virome composition infant` ~ Early_late_infancy + 
                                                 DNA_concentration_ng_ul + read_depth + (1|CS_BABY_BIOME_ID),
                                               REML = F, data = Sample_metadata_infants)
summary(MM_type_composition_infants_early_late) 

