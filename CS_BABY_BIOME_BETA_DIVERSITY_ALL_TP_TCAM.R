# Title: "CS_BABY_BIOME_BETA_DIVERSITY_ALL_TIMEPOINTS and TCAM on one "
# Author: "Trishla Sinha"
# Date: "7/03/2023"
# Last update: "01/04/2024"

#load packages 
library(tidyverse)
library(stringr)
library(vegan)
library(RColorBrewer)
library(wesanderson)
library(ggrepel)
library(mgcv) 
library(reshape2)
library(vtable)
library(dplyr)
library(ggpubr)
library(lmerTest)
library(reshape2)
library(plyr)
library(EnvStats)
library(foreach)
library(ggplot2)

# Loading data 
metadata<-read.delim("~/Desktop/CS_Baby_Biome/submission/2024_submission/analysis/Metadata_EGA_CS_BABY_BIOME.txt")
metadata[sapply(metadata, is.character)] <- lapply(metadata[sapply(metadata, is.character)],  #convert character columns to factors
                                                   as.factor)
metadata$Timepoint_categorical<-factor(metadata$Timepoint_categorical, levels = c("W01", "W02", "W03", "W04", "W05", "W06", "M06", "M12"))
metadata$growth_p_limited<-factor(metadata$growth_p_limited, levels = c("P<10", "P10-P50", "P51-90", ">P90"))
metadata$feeding_mode<-factor(metadata$feeding_mode, levels = c("breast_feeding", "mixed_feeding", "formula_feeding"))
row.names(metadata)<-metadata$bioSampleId

#FILTERING TAXA AT SUB SPECIES, SPECIES, GENUS LEVELS 
metadata_infants<-metadata %>% drop_na(Timepoint_numeric) # Only infant metadata
metadata_infants <- metadata_infants[!(metadata_infants$Timepoint_categorical %in% c("M06", "M12")), ]
taxa <-read.delim("~/Desktop/CS_Baby_Biome/submission/2024_submission/analysis/TAXA_MP4_2022/taxa_clean_mp4_oct_2022_01_04_2024.txt")
taxa=taxa[row.names(taxa)%in% rownames(metadata_infants),] # Only taxa present in infants 
setdiff (rownames(taxa), rownames(metadata_infants)) # all good 

taxa_sp=taxa
taxa_sp=taxa_sp[,grep("s__",colnames(taxa_sp))]
taxa_sp=taxa_sp[,-grep("t__",colnames(taxa_sp))]
taxa_sub_species<-taxa[,grep("t__",colnames(taxa))]

taxa_genus=taxa
taxa_genus=taxa_genus[,grep("g__",colnames(taxa_genus))]
taxa_genus=taxa_genus[,-grep("s__",colnames(taxa_genus))]

#Simplify names
colnames(taxa_sp)=gsub(".*s__","",colnames(taxa_sp))
colnames(taxa_sub_species)=gsub(".*t__","",colnames(taxa_sub_species))
colnames(taxa_genus)=gsub(".*g__","",colnames(taxa_genus))

# calculating alpha diversity per sample 
my_alpha_diversity=data.frame(Shannon=diversity(taxa_sp, index = "shannon"), Simpson=diversity(taxa_sp, index = "simpson"))
phenos_imputed_div=merge(my_alpha_diversity, metadata_infants, by="row.names")
row.names(phenos_imputed_div)=phenos_imputed_div$Row.names
phenos_imputed_div$Row.names=NULL
phenos_imputed_div$rand_AB=factor(phenos_imputed_div$rand_AB)
# calculating richness per sample 
richness <- ddply(taxa_sp,~rownames(taxa_sp),function(x) {
  data.frame(RICHNESS=sum(x[-1]>0, na.rm = T))
})
row.names(richness)<-richness$`rownames(taxa_sp)`
richness$`rownames(taxa_sp)`=NULL
phenos_imputed_div<-merge(richness, phenos_imputed_div, by="row.names")
row.names(phenos_imputed_div)<-phenos_imputed_div$Row.names

phenos_imputed_div= phenos_imputed_div %>%
  filter(!is.na(feeding_mode))
taxa_sp=taxa_sp[row.names(taxa_sp)%in% rownames(phenos_imputed_div),]

taxa_sp$CS_Baby_Biome_ID<-substr(row.names (taxa_sp), 1, 5)
unique_counts <- sapply(taxa_sp, function(x) length(unique(taxa_sp$CS_Baby_Biome_ID[x >0.05])))
taxa_sp_filt <- taxa_sp[, unique_counts >= 2]
taxa_sp_filt$CS_Baby_Biome_ID=NULL

#Transformations
my_pseudocount_normal=min(taxa_sp_filt[taxa_sp_filt!=0])/2# 

taxa_sp_filt<-taxa_sp_filt[match(row.names(taxa_sp_filt),row.names(phenos_imputed_div)),]
my_distance_CLR=vegdist(taxa_sp_filt, method = "aitchison", pseudocount=my_pseudocount_normal) 
taxa_sp_filt_CLR <- decostand(taxa_sp_filt, "clr", pseudocount=my_pseudocount_normal)
mypcoa_CLR=cmdscale(my_distance_CLR, k = 20, eig = T)
my_var_CLR=round(mypcoa_CLR$eig*100/sum(mypcoa_CLR$eig),2)[1:20]
barplot (my_var_CLR)
mypcoa_df_CLR=as.data.frame(mypcoa_CLR$points)
names(mypcoa_df_CLR) <- c('PC1','PC2','PC3','PC4','PC5', 'PC6','PC7','PC8','PC9','PC10',
                          'PC11','PC12','PC13','PC14','PC15', 'PC16','PC17','PC18','PC19','PC20')
phenos_imputed_div_pc_CLR=merge(mypcoa_df_CLR,phenos_imputed_div, by="row.names")
row.names (phenos_imputed_div_pc_CLR)<-phenos_imputed_div_pc_CLR$Row.names


names(phenos_imputed_div_pc_CLR)[names(phenos_imputed_div_pc_CLR) == "rand_AB"] <- "AB"

#setwd("/Users/trishlasinha/Desktop/CS_Baby_Biome/supplementary/figures/")
# Figure 1 D
#pdf('beta_diversity_AB_cs_baby_biome.pdf', width=8, height=5)

#ggplot(phenos_imputed_div_pc_CLR, aes(PC1, PC2, fill=AB))  +
  #scale_fill_manual(values = c("#54b7a6", "#fa2b13"))+
  #scale_color_manual(values = c("#54b7a6", "#fa2b13"))+
  #geom_point(size=3, pch=21) + theme_bw(base_size = 22) + xlab(paste("PC1",my_var_CLR[1], sep = " ")) + 
  # geom_text(data = phenos_imputed_div_pc, aes(label = row.names(phenos_imputed_div_pc)), 
  #position = position_dodge(width=0.9),  size=3)+
  #geom_line(aes(group = CS_BABY_BIOME_ID))+
  #ylab(paste("PC2",my_var_CLR[2], sep = " "))+ggtitle("Aitchison distance")#+ scale_fill_manual(values = c("firebrick", "lightblue"))

#dev.off()

setwd("/Users/trishlasinha/Desktop/CS_Baby_Biome/submission/2024_submission/figures/figures_individual/figure_1/")
pdf('beta_diversity_feeding_mode_all_tp_cs_baby_biome.pdf', width=4, height=5)

overall <- ggplot(phenos_imputed_div_pc_CLR , aes(PC1, PC2, color = feeding_mode)) +
  geom_point(size = 2, alpha = 0.8) +
  stat_ellipse(aes(group = feeding_mode, fill = feeding_mode, color = feeding_mode),
               type = "norm", linetype = 2, geom = "polygon", alpha = 0.05, show.legend = F) +
  xlab(paste("PCo1=", round(my_var_CLR[1], digits = 2), "%", sep = "")) +
  ylab(paste("PCo2=", round(my_var_CLR[2], digits = 2), "%", sep = "")) +
  scale_fill_manual(values = wes_palette("Moonrise3", n = 5)) +
  scale_color_manual(name=NULL,
                     breaks = c("breast_feeding", "mixed_feeding" ,"formula_feeding"),
                     labels = c("Breast", "Mixed", "Formula"),
                     values = wes_palette("Moonrise3", n = 5)) +
  theme(plot.subtitle = element_text(vjust = 1),
        plot.caption = element_text(vjust = 1),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(size = 18),
        legend.key = element_rect(fill = NA),
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        panel.background = element_rect(fill = NA))

overall <- ggExtra::ggMarginal(overall, type = "histogram", groupColour = F, groupFill = TRUE,
                               xparams = list(bins = 60, alpha = 0.8, position = 'identity', color = 'white'),
                               yparams = list(bins = 60, alpha = 0.8, position = 'identity', color = 'white'))

overall
dev.off()


######################## TENSOR FACTORIZATION ######################

# Condensing all the data points into one using TCAM ()
tcam<-read.delim("~/Desktop/CS_Baby_Biome/submission/2024_submission/analysis/TCAM_CS_BABY_BIOME_01_04_2024.txt", sep = ",")
row.names(tcam)<-tcam$X
tcam$X=NULL
tcam$CS_BABY_BIOME_ID=NULL
tcam$Cluster=NULL
str (tcam)
rownames(tcam) <- substr(rownames(tcam), nchar(rownames(tcam)) - 1, nchar(rownames(tcam)))

# Making distance matrix 
selection<-metadata_infants[,c("CS_BABY_BIOME_ID", "cefazoline_measurement_mg_L", "rand_AB", "pre_preg_bmi_mother",   "preg_gest_age", "preg_weight_gain", "gravida", "para", "infant_birthweight", "infant_sex", "mother_age_at_delivery","APGAR_1", "APGAR_5", "growth_p_limited", "feeding_mode_pragmatic", "living_situation", "cats_dogs"  )]
selection <- selection[!duplicated(selection$CS_BABY_BIOME_ID), ]
row.names(selection)<-selection$CS_BABY_BIOME_ID
summary (selection$feeding_mode_pragmatic)
rownames(selection) <- substr(rownames(selection), nchar(rownames(selection)) - 1, nchar(rownames(selection)))
selection<-selection[match(row.names(selection),row.names(tcam)),]

#breast_feeding formula_feeding   mixed_feeding 
#11              10               7 

# Making distance matrix 

Distance=vegdist(tcam,method = "euclidean" ) 
Distance=as.matrix(Distance)


phenotypes<-selection[,2:17]
# Running adonis for all other factors 
adonis_results_raw<- data.frame(matrix(ncol =4, nrow= ncol(phenotypes)))  
rownames(adonis_results_raw) <- colnames(phenotypes)
colnames(adonis_results_raw) <- c("Df", "F", "R2", "p-value")


adon<-foreach(i=1:ncol(phenotypes),.combine=rbind)%do%{
  #k<-which(!is.na(phenotypes[,i]))
  r<-which(is.na(phenotypes[,i]))
  r2 <- row.names(phenotypes)[r]
  distmat_cleaned <- Distance[!(rownames(Distance) %in% r2),
                              !(colnames(Distance)) %in% r2]  
  phenos2 <- phenotypes[rownames(phenotypes) %in% rownames(distmat_cleaned),]
  
  ad1<-adonis2(distmat_cleaned ~ phenos2[[i]],
               permutations=10000,parallel=8, na.action = na.fail,
               by="margin") 
  print(ad1)
  adonis_results_raw[i,] <- c(ad1$Df[1], ad1$F[1], ad1$R2[1],ad1$'Pr(>F)'[1])
}

setwd("/Users/trishlasinha/Desktop/CS_Baby_Biome/submission/2024_submission/analysis/RESULTS/")
write.table(adonis_results_raw, "adonis_results_tcam_raw_no_correction_01_04_2024.txt", sep = "\t", quote = F)

setwd("/Users/trishlasinha/Desktop/CS_Baby_Biome/supplementary/figures/")
pdf('beta_diversit_explained_variance_tcam_cs_baby_biome_raw_no_correction.pdf', width=8, height=5)

# Basic barplot
adonis_results_raw_2<-adonis_results_raw %>% drop_na()
p <- ggplot(data = adonis_results_raw_2, aes(x = reorder(row.names(adonis_results_raw_2), R2), y = R2)) +
  geom_bar(stat = "identity", aes(fill = ifelse(`p-value` < 0.05, "Significant", "Non-Significant"))) +
  coord_flip() +
  xlab("Phenotype") +
  theme_bw() +
  scale_fill_manual(values = c("Significant" = "blue", "Non-Significant" = "grey")) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.title = element_blank())


p

dev.off()


# Running adonis for all other factors correcting for feeding mode 
adonis_results<- data.frame(matrix(ncol =4, nrow= ncol(phenotypes)))  
rownames(adonis_results) <- colnames(phenotypes)
colnames(adonis_results) <- c("Df", "F", "R2", "p-value")


adon<-foreach(i=1:ncol(phenotypes),.combine=rbind)%do%{
  #k<-which(!is.na(phenotypes[,i]))
  r<-which(is.na(phenotypes[,i]))
  r2 <- row.names(phenotypes)[r]
  distmat_cleaned <- Distance[!(rownames(Distance) %in% r2),
                              !(colnames(Distance)) %in% r2]  
  phenos2 <- phenotypes[rownames(phenotypes) %in% rownames(distmat_cleaned),]
  
  ad1<-adonis2(distmat_cleaned ~ phenos2[["feeding_mode_pragmatic"]]+phenos2[[i]],
               permutations=10000,parallel=8, na.action = na.fail,
               by="margin") 
  print(ad1)
  adonis_results[i,] <- c(ad1$Df[2], ad1$F[2], ad1$R2[2],ad1$'Pr(>F)'[2])
}

adonis_results<-adonis_results %>% drop_na()
setwd("/Users/trishlasinha/Desktop/CS_Baby_Biome/submission/2024_submission/analysis/RESULTS/")
write.table (adonis_results, "adonis_results_tcam_correction_feeding_mode_01_04_2024.txt", row.names = T, quote = F, sep="\t")

setwd("/Users/trishlasinha/Desktop/CS_Baby_Biome/supplementary/figures/")
pdf('beta_diversit_explained_variance_tcam_cs_baby_biome_correction_feeding_mode.pdf', width=8, height=5)



p <- ggplot(data = adonis_results, aes(x = reorder(row.names(adonis_results), R2), y = R2)) +
  geom_bar(stat = "identity", aes(fill = ifelse(`p-value` < 0.05, "Significant", "Non-Significant"))) +
  coord_flip() +
  xlab("Phenotype") +
  theme_bw() +
  scale_fill_manual(values = c("Significant" = "blue", "Non-Significant" = "grey")) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.title = element_blank())


p
dev.off()



selection_2<-selection[!is.na(selection$pre_preg_bmi_mother),]
tcam_2<-tcam[match(row.names(selection_2),row.names(tcam)),]

# Making distance matrix 
Distance=vegdist(tcam_2,method = "euclidean" ) 
Distance=as.matrix(Distance)

# Running adonis for all other factors correcting for feeding mode 
adonis_results<- data.frame(matrix(ncol =4, nrow= ncol(phenotypes)))  
rownames(adonis_results) <- colnames(phenotypes)
colnames(adonis_results) <- c("Df", "F", "R2", "p-value")


adon<-foreach(i=1:ncol(phenotypes),.combine=rbind)%do%{
  #k<-which(!is.na(phenotypes[,i]))
  r<-which(is.na(phenotypes[,i]))
  r2 <- row.names(phenotypes)[r]
  distmat_cleaned <- Distance[!(rownames(Distance) %in% r2),
                              !(colnames(Distance)) %in% r2]  
  phenos2 <- phenotypes[rownames(phenotypes) %in% rownames(distmat_cleaned),]
  
  ad1<-adonis2(distmat_cleaned ~ phenos2[["feeding_mode_pragmatic"]]+ phenos2[["pre_preg_bmi_mother"]] + phenos2[[i]],
               permutations=10000,parallel=8, na.action = na.fail,
               by="margin") 
  print(ad1)
  adonis_results[i,] <- c(ad1$Df[3], ad1$F[3], ad1$R2[3],ad1$'Pr(>F)'[3])
}

adonis_results<-adonis_results %>% drop_na()
setwd("/Users/trishlasinha/Desktop/CS_Baby_Biome/submission/2024_submission/analysis/RESULTS/")
write.table (adonis_results, "adonis_results_tcam_correction_feeding_mode_pre_pre_BMI_01_04_2024.txt", row.names = T, quote = F, sep = "\t")

setwd("/Users/trishlasinha/Desktop/CS_Baby_Biome/supplementary/figures/")
pdf('beta_diversit_explained_variance_tcam_cs_baby_biome_correction_feeding_mode_pre_preg_BMI.pdf', width=8, height=5)



p <- ggplot(data = adonis_results, aes(x = reorder(row.names(adonis_results), R2), y = R2)) +
  geom_bar(stat = "identity", aes(fill = ifelse(`p-value` < 0.05, "Significant", "Non-Significant"))) +
  coord_flip() +
  xlab("Phenotype") +
  theme_bw() +
  scale_fill_manual(values = c("Significant" = "blue", "Non-Significant" = "grey")) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.title = element_blank())


p
dev.off()


######### FIGURE  BETA DIVERSIY FEEDING MODE ########################
setwd("/Users/trishlasinha/Desktop/CS_Baby_Biome/submission/2024_submission/figures/figures_individual/figure_1/")

all<-merge(selection, tcam, by="row.names")
all$feeding_mode_pragmatic<-factor(all$feeding_mode, levels = c("breast_feeding", "mixed_feeding", "formula_feeding"))
tcam_plot<-ggplot(all ,aes(X0,X1, color = feeding_mode_pragmatic))+
  geom_point(size = 2,alpha = 0.8)+
  stat_ellipse(aes(group = feeding_mode_pragmatic, fill = feeding_mode_pragmatic, color = feeding_mode_pragmatic) ,type = "norm",linetype = 2,geom = "polygon",alpha = 0.05,show.legend = F)+
  xlab("Axis 1 (7.94%)")+
  ylab("Axis 2 (7.68%)")+
  scale_fill_manual(values = wes_palette("Moonrise3", n = 5))+
  scale_color_manual(name=NULL, 
                     breaks = c("breast_feeding", "mixed_feeding" ,"formula_feeding"),
                     labels = c("Breast"          , "Mixed", "Formula"),
                     values = wes_palette("Moonrise3", n = 5))+
  theme(plot.subtitle = element_text(vjust = 1),
        plot.caption = element_text(vjust = 1),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(size = 18),
        legend.key = element_rect(fill = NA),
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        panel.background = element_rect(fill = NA))
merged_plot<-ggExtra::ggMarginal(tcam_plot, type = "histogram", groupColour = F, groupFill = TRUE,
                          xparams = list(bins = 60, alpha = 0.8,position = 'identity', color = 'white'),
                          yparams = list(bins = 60, alpha = 0.8,position = 'identity', color = 'white'))


merged_BETA_TCAM_ALL_feeding<-ggarrange(overall, merged_plot)
merged_BETA_TCAM_ALL_feeding


dev.off()


############ Figures BETA DIVERSITY AB ################### 
phenos_imputed_div_pc_CLR$AB<-as.factor(phenos_imputed_div_pc_CLR$AB)

overall_AB <- ggplot(phenos_imputed_div_pc_CLR , aes(PC1, PC2, color = AB)) +
  geom_point(size = 2, alpha = 0.8) +
  stat_ellipse(aes(group = AB, fill = AB, color = AB),
               type = "norm", linetype = 2, geom = "polygon", alpha = 0.05, show.legend = F) +
  xlab(paste("PCo1=", round(my_var_CLR[1], digits = 2), "%", sep = "")) +
  ylab(paste("PCo2=", round(my_var_CLR[2], digits = 2), "%", sep = "")) +
  scale_fill_manual(values = c("#54b7a6", "#fa2b13"))+
  scale_color_manual(name=NULL,
                     breaks = c("No", "Yes"),
                     labels = c("AB-", "AB+"),
                     values = c("#54b7a6", "#fa2b13")) +
  theme(plot.subtitle = element_text(vjust = 1),
        plot.caption = element_text(vjust = 1),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(size = 18),
        legend.key = element_rect(fill = NA),
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        panel.background = element_rect(fill = NA))


overall_AB_extra <- ggExtra::ggMarginal(overall_AB , type = "histogram", groupColour = F, groupFill = TRUE,
                               xparams = list(bins = 60, alpha = 0.8, position = 'identity', color = 'white'),
                               yparams = list(bins = 60, alpha = 0.8, position = 'identity', color = 'white'))

overall_AB_extra



tcam_plot_AB<-ggplot(all ,aes(X0,X1, color = rand_AB))+
  geom_point(size = 2,alpha = 0.8)+
  stat_ellipse(aes(group = rand_AB, fill = rand_AB, color = rand_AB) ,type = "norm",linetype = 2,geom = "polygon",alpha = 0.05,show.legend = F)+
  xlab("Axis 1 (7.94%)")+
  ylab("Axis 2 (7.68%)")+
  scale_fill_manual(values = c("#54b7a6", "#fa2b13"))+
  scale_color_manual(name=NULL,
                     breaks = c("No", "Yes"),
                     labels = c("AB-", "AB+"),
                     values = c("#54b7a6", "#fa2b13")) +
  theme(plot.subtitle = element_text(vjust = 1),
        plot.caption = element_text(vjust = 1),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(size = 18),
        legend.key = element_rect(fill = NA),
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        panel.background = element_rect(fill = NA))
tcam_extra_plot_AB<-ggExtra::ggMarginal(tcam_plot_AB, type = "histogram", groupColour = F, groupFill = TRUE,
                                 xparams = list(bins = 60, alpha = 0.8,position = 'identity', color = 'white'),
                                 yparams = list(bins = 60, alpha = 0.8,position = 'identity', color = 'white'))
tcam_extra_plot_AB

merged_BETA_TCAM_ALL_AB<-ggarrange(overall_AB_extra, tcam_extra_plot_AB)
merged_BETA_TCAM_ALL_AB


merged_BETA_TCAM_ALL_AB <- ggarrange(overall, merged_plot, overall_AB_extra, tcam_extra_plot_AB,
                                     labels = c("D", "E", "F", "G"), 
                                     ncol = 2, nrow = 2)  

merged_BETA_TCAM_ALL_AB
setwd("/Users/trishlasinha/Desktop/CS_Baby_Biome/submission/2024_submission/figures/figures_individual/figure_1/")
save(overall, merged_plot, overall_AB_extra, tcam_extra_plot_AB, file="figure_1_D_E_F_G.RData")


