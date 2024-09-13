################################################################################
##### CS Baby Biome: Pipeline intermediate results
### Author(s): Asier Fernández-Pato
### Last updated: 19th July, 2024
################################################################################

# Set working directory 
setwd("~/Desktop/PhD/Projects/CS Baby Biome/Final_analysis/")

#****************
# Load libraries
#****************
library(dplyr)
library(ggplot2)

#*******************
# STEP1: Viral identification
#*******************
# 15,165 viral contigs identified
# 8,998 contigs predicted by VS2  
# 3,119 contigs predicted by DVF
# 6,570 contigs predicted by geNomad
# 6,521 contigs predicted by VIBRANT
# 8,693 contigs predicted by CT3

#*******************
# STEP2: Re-assembly of viral contigs
#*******************
# 15,165 initial viral contigs
# 111 predicted viral contigs were extended into circular genomes 
# 4,641 viral contigs were partially extended. 
# 9,860 viral contigs were NOT extended in total
## 5,370 failed due to COBRA rules
## 4,206 failed due to orphan end
## 284 failed because they were already circular

#Create dataframe with COBRA results
df_COBRA <- data.frame(
  Genomes = c("Extended circular", "Self circular", "Extended partial", "Not extended"),
  value = c(111, 284, 4641, 9576)
)
df_COBRA$Genomes <- factor_var <- factor(df_COBRA$Genomes, 
                                         levels = c("Extended circular", 
                                                    "Self circular",
                                                    "Extended partial", 
                                                    "Not extended"))

# Create a stacked barplot 
pdf('0_PIPELINE_RESULTS/Plots/COBRA_results_barplot.pdf', width = 4, height = 2)
ggplot(df_COBRA, aes(x = "", y = value, fill = Genomes)) + 
  geom_col(aes(alpha = 0.7)) +  # Set the alpha to 0.6 for the bars
  coord_flip() +
  labs(y="Number of viral genomes") +
  scale_fill_manual(values = alpha(c("#B71C1C", "#EF5350", "#26A69A", "#B0BEC5"), alpha = 0.7)) +  
  theme_classic() + 
  theme(
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 13.5),
    legend.text = element_text(size = 13),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    legend.position = "top"
  )
dev.off()

#*******************
# STEP3: geNomad filtering
#*******************
# 14,612 sequences generated after COBRA as input
# 7,863 sequences identified as viral by geNomad

#Create dataframe with geNomad results
df_geNomad <- data.frame(
  Origin = c("Viral contigs", "Other"),
  value = c(7863, 6749)
)

df_geNomad$Origin <- factor(df_geNomad$Origin, 
                                         levels = c("Viral contigs",
                                                    "Other"))

# Estimate the proportions
df_geNomad <- df_geNomad %>% 
  arrange(desc(Origin)) %>%
  mutate(prop = value / sum(df_geNomad$value) *100) 

# Generate a piechart
colors <- c("#66C2A5", "#1F78B4")

pdf('0_PIPELINE_RESULTS/Plots/geNomad_filtering.pdf', width = 5, height = 2.5)
ggplot(df_geNomad, aes(x = "", y = prop, fill = Origin)) +
  geom_bar(stat = "identity", width = 1, color = NA) +
  scale_fill_manual(values = alpha(colors, alpha = 0.8)) +  # Add alpha to the colors
  coord_polar("y", start = 0) +
  theme_void() + 
  theme(legend.position = "right", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 14))
dev.off()

#*******************
# STEP4: Quality check
#*******************
# 7,863 viral contigs after geNomad filtering
# 3,614 viral sequences have > 50% completeness
# Complete genomes: 417 
# >90-99%: 1,611 
# >50-90%: 1,586 
# <50% is: 3,608 
# Undetermined:642 

#Create DF with CheckV results 
df_CheckV <- data.frame(
  Completeness = c("Undetermined","Low","Medium", "High", "Complete"),
  value = c(642, 3608, 1586, 1611, 417)
)
df_CheckV$Completeness <- factor(df_CheckV$Completeness, levels=c("Complete", "High", "Medium", "Low", "Undetermined"))

# Create a stacked barplot 
pdf('0_PIPELINE_RESULTS/Plots/CheckV_results_barplot.pdf', width = 4, height = 2)
ggplot(df_CheckV, aes(x = "", y = value, fill = Completeness)) + 
  geom_col(aes(alpha = 0.9)) +  
  coord_flip() +
  labs(y="Number of viral genomes") +
  scale_fill_manual(values = alpha(c("#283593", "#9FA8DA", "#BBDEFB","#FFCDD2","#B0BEC5"), alpha = 0.9)) +  
  theme_classic() + 
  theme(
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 13.5),
    legend.text = element_text(size = 13),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    legend.position = "top"
  )
dev.off()

#*******************
# STEP6: Dereplication
#*******************
# 3,446 viral contigs after CheckV step
# Public sequences: 189680 (MGV) + 84537 (GPD) + 637 (Gulyaeva) + 249 (Guerin) + 
# 648 (Yutin) + 5199 (Viral RefSeq) + Benler (1480) + Shah (4,627) + ELGV (21,295) + IMGVR (36,064) = 342,500
# Total= 345,946
# After dereplication: 100,250 vOTUs
# Representative sequences as follows:
#- 859 from CS Baby Biome
#- 31,905 from GPD DB
#- 35,470 from MGV DB
#- 3,398 from Viral Refseq
#- 1,572 from Shah et al
#- 496 from Benler et al
#- 11,380 from IMGVR et al
#- 15,034 from ELGV et al
#- 136 from crAss DBs (86,28,22)

#Create DF with input sequences for dereplication
df_STEP6_v0TUs <- data.frame(
  Database = c("CS","GPD","MGV", "RefSeq","Shah et al", "Benler et al",
               "CrAss DBs", "ELGV", "IMG_VR"),
  value = c(3446, 84537, 189680, 5199, 4627, 1480, 1534, 21295, 36064)
)

df_STEP6_v0TUs$Database <- factor(df_STEP6_v0TUs$Database, levels =c("Benler et al", "CS", "GPD",
                                                                     "CrAss DBs", "MGV", "RefSeq",
                                                                     "Shah et al", "ELGV", 
                                                                     "IMG_VR"))
# Estimate the proportions
df_STEP6_v0TUs <- df_STEP6_v0TUs %>% 
  arrange(desc(Database)) %>%
  mutate(prop = value / sum(df_STEP6_v0TUs$value) *100) 

# Generate a piechart
colors <- c("#8DD3C7", "#FFDAB9", "#BEBADA", "#FB8072", "#80B1D3",
            "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

pdf('0_PIPELINE_RESULTS/Plots/Dereplication_piechart.pdf', width = 5, height = 2.5)
ggplot(df_STEP6_v0TUs, aes(x = "", y = prop, fill = Database)) +
  geom_bar(stat = "identity", width = 1, color = NA) +
  scale_fill_manual(values = alpha(colors, alpha = 0.8)) +  # Add alpha to the colors
  coord_polar("y", start = 0) +
  theme_void() + 
  theme(legend.position = "right", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 14))
dev.off()

