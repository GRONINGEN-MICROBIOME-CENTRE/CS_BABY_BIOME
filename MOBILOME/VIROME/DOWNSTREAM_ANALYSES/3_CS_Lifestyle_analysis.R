################################################################################
##### CS Baby Biome: Viral Lifestyle Analysis
### Author(s):Asier Fernández-Pato
### Last updated: 7th June, 2024
################################################################################

#****************
# Load modules
#****************
library(dplyr)
library(vegan)
library(lme4)
library(lmerTest)
library(RLRsim)
library(ggplot2)
library(ggdist)

# Set working directory 
setwd("~/Desktop/PhD/Projects/CS Baby Biome/Final_analysis/VIRUSES/")

# Read abundance table and metadata tables
Abundance_table <- read.delim("Abundance_table/CS_Abundance_Table_17052024.txt")
Sample_metadata <- read.delim("Metadata_CS/CS_Baby_Biome_Sample_Metadata_17052024.txt")
Virus_metadata <- read.delim("Metadata_CS/CS_Baby_Biome_Viral_Metadata_17052024.txt")
Sample_metadata_infants <- read.delim("Metadata_CS/CS_Baby_Biome_Sample_Metadata_Infants_17052024.txt")

# Set timepoint variable from metadata as factor
Sample_metadata$Timepoint_categorical <- factor(Sample_metadata$Timepoint_categorical, 
                                                levels=c("MOM", "W01", "W02", "W03","W04",
                                                         "W05", "W06"))

##************************************************************************
# 1. Viral lifestyle analysis: richness, diversity and abundance of temperate phages
##************************************************************************

# Check alpha diversity and richness only of temperate phages
temperate_viruses <-Virus_metadata$Virus_ID[which(Virus_metadata$Lifestyle == "Temperate")]
Abundance_table_temperates <- Abundance_table[rownames(Abundance_table) %in% temperate_viruses,]

# Add to sample metadata (including proportion of temperate phages)
Sample_metadata$richness_temperates  <-specnumber(t(Abundance_table_temperates))                
Sample_metadata$shannon_temperates <- vegan::diversity(t(Abundance_table_temperates), index = "shannon") 
Sample_metadata$prop_temperates <- 100*(Sample_metadata$richness_temperates / Sample_metadata$richness)

# Estimate the total abundance of temperate phages per sample
Sample_metadata$Viral_abundance_temperates <- colSums(Abundance_table_temperates)

# Analyze the proportion of the viral abundance corresponding to temperate phages
Sample_metadata$Viral_relab_temperates <- 100*(colSums(Abundance_table_temperates) / Sample_metadata$Viral_abundance)

# Add the estimated variables to the infant metadata
lifestyle_variables <- Sample_metadata[,c("NG_ID", "richness_temperates", "shannon_temperates",
                                          "prop_temperates", "Viral_relab_temperates")]
Sample_metadata_infants <- left_join(Sample_metadata_infants, lifestyle_variables, by = "NG_ID")

##################################
# Check statistical significance of differences
##################################

# Check differences in total/relative abundance of temperates between mothers and babies
MM_type_temperate_relab <- lmer(Viral_relab_temperates ~ Type + DNA_concentration_ng_ul + read_depth + (1|CS_BABY_BIOME_ID),
                                REML = F, data = Sample_metadata)
summary(MM_type_temperate_relab) 

MM_type_temperate_ab <- lmer(Viral_abundance_temperates ~ Type + DNA_concentration_ng_ul + read_depth + (1|CS_BABY_BIOME_ID),
                             REML = F, data = Sample_metadata)
summary(MM_type_temperate_ab) 


# Check correlation of the proportion/abundance of temperate phages with sequencing depth and DNA concentration
spearman_cor_prop_temperate_CONC <- cor.test(Sample_metadata$DNA_concentration_ng_ul,
                                             Sample_metadata$prop_temperates, method = "spearman",) 
spearman_cor_prop_temperate_Read_DEPTH <- cor.test(Sample_metadata$read_depth,
                                                   Sample_metadata$prop_temperates, method = "spearman") 

spearman_cor_Viral_relab_temperates_CONC <- cor.test(Sample_metadata$DNA_concentration_ng_ul,
                                                     Sample_metadata$Viral_abundance_temperates, method = "spearman") 
spearman_cor_Viral_relab_temperates_Read_DEPTH <- cor.test(Sample_metadata$read_depth,
                                                           Sample_metadata$Viral_abundance_temperates, method = "spearman") 

# Plot the proportion and relative abundances in mothers and babies
pdf('3_LIFESTYLE/Plots/Prop_temperates_comp.pdf', width=2.8, height=3.2)
ggplot(Sample_metadata[!is.na(Sample_metadata$Timepoint_categorical),], aes(x=Type, y=prop_temperates)) +
  geom_violin(aes(fill=Type), alpha=0.7, width=0.8, trim = F) +
  geom_jitter(aes(color=Type), alpha=0.7) + 
  geom_boxplot(width=0.2, fill="white", color="black") + 
  labs(x = 'Sample type', y = 'Proportion (%)') + 
  scale_fill_manual(values=c("#008080", "#4B0082")) + 
  scale_color_manual(values=c("#008080", "#4B0082")) +
  ylim(0,120) +
  scale_y_continuous(breaks = c(25, 50, 75, 100)) +
  theme_classic() +
  theme(plot.title = element_text(size=18, hjust = 0.5, face="bold"), 
        axis.text = element_text(size=14),
        axis.title = element_text(size=15), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.position = "none")
dev.off()

pdf('3_LIFESTYLE/Plots/Relab_temperates_comp.pdf', width=2.8, height=3.2)
ggplot(Sample_metadata[!is.na(Sample_metadata$Timepoint_categorical),], aes(x=Type, y=Viral_relab_temperates)) +
  geom_violin(aes(fill=Type), alpha=0.7, width=0.8, trim = F) +
  geom_jitter(aes(color=Type), alpha=0.7) + 
  geom_boxplot(width=0.2, fill="white", color="black", outliers = FALSE) + 
  labs(x = 'Sample type', y = 'Relative Abundance (%)') + 
  scale_fill_manual(values=c("#008080", "#4B0082")) + 
  scale_color_manual(values=c("#008080", "#4B0082")) +
  ylim(0,120) +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100)) +
  theme_classic() +
  theme(plot.title = element_text(size=18, hjust = 0.5, face="bold"), 
        axis.text = element_text(size=14),
        axis.title = element_text(size=15), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.position = "none") 
dev.off()


# Plot the correlations (change variables to be plotted)
pdf('3_LIFESTYLE/Plots/Spearman_cor_viral_prop_read_depth_splitted.pdf', width=2.5, height=2.9)
ggscatter(Sample_metadata, x = "prop_temperates", y = "read_depth",
          add = "reg.line", color = "Type" , conf.int = FALSE, shape = 19, alpha = 0.4) + # "#FFB900" "#0072B2"
  stat_cor(aes(color = Type), label.y.npc = 0.7, label.x.npc = 0.03, cor.coef.name = "rho", method = "spearman", size = 4.2) +
  labs(x = "Proportion (%)", y = "Number of reads (million)") +
  scale_color_manual(values=c("#008080", "#4B0082")) +
  scale_y_continuous(labels = scales::label_number(scale = 1e-6)) +
  facet_wrap(~.,) +
  theme_classic() + 
  border(size = ) +
  theme(axis.text=element_text(size=12),
        axis.title = element_text(size = 13), 
        strip.background = element_rect(fill = "grey90", linewidth = 1),
        strip.text = element_text(size = 13),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.position = "none") 
dev.off()


##************************************************************************
# 2. Viral lifestyle analysis in Infants
#*************************************************************************

# Check statistical significance over time (time as continuous variable - age(months))

# A) Proportion of temperates
MM_time_prop_temperates_infant <- lmer(prop_temperates ~ Timepoint_numeric + DNA_concentration_ng_ul + read_depth + (1|CS_BABY_BIOME_ID),
                                       REML = F,  data = Sample_metadata[Sample_metadata$Type=="Infant" &
                                                                           Sample_metadata$Timepoint_categorical %in% c("W01", "W02", "W03", "W04", "W05", "W06"),])
summary(MM_time_prop_temperates_infant) 

# B) Relative abundance of temperates
MM_time_relab_temperates_infant <- lmer(Viral_relab_temperates ~ Timepoint_numeric + DNA_concentration_ng_ul + read_depth + (1|CS_BABY_BIOME_ID), REML = F,
                                        data = Sample_metadata[Sample_metadata$Type=="Infant" &
                                                                 Sample_metadata$Timepoint_categorical %in% c("W01", "W02", "W03", "W04", "W05", "W06"),])
summary(MM_time_relab_temperates_infant)   


# Generate plots for the proportion and cumulative relative abundance of temperate phages
pdf('3_LIFESTYLE/Plots/Prop_temperates_infants.pdf', width=3.8, height=3.1)
ggplot(Sample_metadata[Sample_metadata$Type=="Infant" & Sample_metadata$Timepoint_categorical %in% c("W01", "W02", "W03", "W04", "W05", "W06"),],
       aes(x=Timepoint_categorical, y=prop_temperates, fill = Type)) +
  geom_violin(aes(fill=Type), alpha=0.5, width=0.8, trim = F) +
  geom_jitter(alpha = 0.7, color ="#008080") +
  geom_boxplot(width=0.2, fill="white", color="black", outliers = F) +
  labs(y="Proportion (%)", x="Timepoint") +
  scale_fill_manual(values=c("#008080")) + 
  theme_classic()  +
  theme(plot.title = element_text(size=18, hjust = 0.5, face="bold"), 
        axis.text = element_text(size=14),
        axis.title = element_text(size=16), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.position = "none") 
dev.off()

pdf('3_LIFESTYLE/Plots/Relab_temperates_infants.pdf', width=3.8, height=3.1)
ggplot(Sample_metadata[Sample_metadata$Type=="Infant" & Sample_metadata$Timepoint_categorical %in% c("W01", "W02", "W03", "W04", "W05", "W06"),],
       aes(x=Timepoint_categorical, y=Viral_relab_temperates, fill = Type)) +
  geom_violin(aes(fill=Type), alpha=0.5, width=0.8, trim=F) +
  geom_jitter(alpha = 0.7, color ="#008080") +
  geom_boxplot(width=0.2, fill="white", color="black", outliers = F) +
  labs(y="Relative Abundance (%)", x="Timepoint") +
  scale_fill_manual(values=c("#008080")) +
  theme_classic()  +
  theme(plot.title = element_text(size=18, hjust = 0.5, face="bold"), 
        axis.text = element_text(size=14),
        axis.title = element_text(size=16), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.position = "none") 
dev.off()

pdf('3_LIFESTYLE/Plots/Prop_temperates_infants_alternative.pdf', width=3.9, height=3.1)
ggplot(Sample_metadata_infants, aes(x=Timepoint_categorical, y=prop_temperates, fill = Type)) +
  stat_halfeye(
    width = 0.8, position = "dodge", adjust = 0.3, alpha = 0.5,
    justification=-.2, .width = 0, point_colour = NA) + 
  geom_boxplot(width=.25, outlier.color = NA, alpha = 0.7) +
  stat_histinterval(slab_color = "gray45",slab_linewidth = 0.3, outline_bars = T,  width = 0.9,
                    position = "dodge",alpha = 0.5, fill = "white",
                    justification=-.2,  breaks = 60, .width = 0) +  
  labs(y="Proportion (%)", x="Timepoint") +
  scale_fill_manual(values=c("#008080")) +  # "#3366CC" mother , "#4B0082"
  theme_classic()  +
  theme(plot.title = element_text(size=18, hjust = 0.5, face="bold"), 
        axis.text = element_text(size=14),
        axis.title = element_text(size=16), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.position = "none") 
dev.off()

##*************
# Save output
#**************
write.table(Sample_metadata,"Metadata_CS/CS_Baby_Biome_Sample_Metadata_17052024.txt", sep = "\t", row.names = F, quote = FALSE) 
write.table(Sample_metadata_infants,"Metadata_CS/CS_Baby_Biome_Sample_Metadata_Infants_17052024.txt", sep = "\t", row.names = F, quote = FALSE) 

