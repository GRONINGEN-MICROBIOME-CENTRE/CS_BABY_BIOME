#############################################################
# Here we perform a combined  analysis of CS Baby Biome 
# at week 4 together with Lifelines NEXT and MAMI trial  

# Author: Trishla Sinha 
# Date: 20-10-2023
# Last update 01-04-2024

##################################################
rm (list = ls())
setwd("~/Desktop/CS_Baby_Biome/submission/2024_submission/analysis/TAXA_MP4_2022/")

##############################
# Functions
##############################

linear_model_taxa_simple <- function(metadata, ID, CLR_transformed_data, pheno_list) {
  
  df <- metadata
  row.names(df) <- df[,ID]
  df<-merge(df, CLR_transformed_data, by='row.names')
  row.names(df) <- df$Row.names
  df$Row.names <- NULL
  
  Prevalent= c(colnames(CLR_transformed_data))
  #pheno_list= phenotypes
  
  Overall_result_phenos =tibble() 
  
  for (Bug in Prevalent){
    if (! Bug %in% colnames(df)){ next }
    Bug2 = paste(c("`",Bug, "`"), collapse="")
    for ( pheno in pheno_list){
      pheno2 = paste(c("`",pheno, "`"), collapse="")
      df[is.na(df[colnames(df) == pheno]) == F, ID] -> To_keep
      df_pheno = filter(df, !!sym(ID) %in% To_keep )
      
      Model2 = as.formula(paste( c(Bug2,  " ~ clean_reads+ cohort + ",pheno2), collapse="" ))
      lm(Model2, df_pheno) -> resultmodel2
      
      as.data.frame(summary(resultmodel2)$coefficients)[grep(pheno, row.names(as.data.frame(summary(resultmodel2)$coefficients))),] -> Summ_simple
      Summ_simple %>% rownames_to_column("Feature") %>% as_tibble() %>% mutate(P = Summ_simple$`Pr(>|t|)`, Bug =Bug, Pheno=pheno, Model="simple") -> temp_output
      rbind(Overall_result_phenos, temp_output) -> Overall_result_phenos
    }
  }
  
  p=as.data.frame(Overall_result_phenos)
  p$FDR<-p.adjust(p$P, method = "BH")
  
  return(p)
}


linear_model_taxa_cor_feeding <- function(metadata, ID, CLR_transformed_data, pheno_list) {
  
  df <- metadata
  row.names(df) <- df[,ID]
  df<-merge(df, CLR_transformed_data, by='row.names')
  row.names(df) <- df$Row.names
  df$Row.names <- NULL
  
  Prevalent= c(colnames(CLR_transformed_data))
  #pheno_list= phenotypes
  
  Overall_result_phenos =tibble() 
  
  for (Bug in Prevalent){
    if (! Bug %in% colnames(df)){ next }
    Bug2 = paste(c("`",Bug, "`"), collapse="")
    for ( pheno in pheno_list){
      pheno2 = paste(c("`",pheno, "`"), collapse="")
      df[is.na(df[colnames(df) == pheno]) == F, ID] -> To_keep
      df_pheno = filter(df, !!sym(ID) %in% To_keep )
      
      Model2 = as.formula(paste( c(Bug2,  " ~ clean_reads+ cohort + feeding_mode + ",pheno2), collapse="" ))
      lm(Model2, df_pheno) -> resultmodel2
      
      as.data.frame(summary(resultmodel2)$coefficients)[grep(pheno, row.names(as.data.frame(summary(resultmodel2)$coefficients))),] -> Summ_simple
      Summ_simple %>% rownames_to_column("Feature") %>% as_tibble() %>% mutate(P = Summ_simple$`Pr(>|t|)`, Bug =Bug, Pheno=pheno, Model="simple") -> temp_output
      rbind(Overall_result_phenos, temp_output) -> Overall_result_phenos
    }
  }
  
  p=as.data.frame(Overall_result_phenos)
  p$FDR<-p.adjust(p$P, method = "BH")
  
  return(p)
}


##############################
# Loading libraries
##############################
library(tidyverse)
library(stringr)
library(vegan)
library(RColorBrewer)
library(reshape2)
library(pheatmap)
library(vtable)
library(dplyr)
library(plyr)
library(ggpubr)
library(wesanderson)
library(ggpubr)
library(foreach)
library(EnvStats)


##############################
# Input data
##############################

metadata<-read.delim("~/Desktop/CS_Baby_Biome/CS_BABY_BIOME_AMSTERDAM_COMBINED/metdata_combined.txt")
row.names(metadata)<-metadata$ID

metadata$feeding_mode<-factor(metadata$feeding_mode, levels = c("breast_feeding", "mixed_feeding", "formula_feeding"))
metadata$antibiotics<-as.factor(metadata$antibiotics)
summary (metadata$antibiotics)
taxa <-read.delim("~/Desktop/CS_Baby_Biome/submission/2024_submission/analysis/TAXA_MP4_2022/taxa_clean_mp4_oct_2022_01_04_2024.txt")
overlap<-intersect(row.names(taxa), row.names(metadata))
taxa<-taxa[overlap,]
metadata<-metadata[overlap,]

#Derive only species level data 
taxa_sp=taxa
taxa_sp=taxa_sp[,grep("s__",colnames(taxa_sp))]
taxa_sp=taxa_sp[,-grep("t__",colnames(taxa_sp))]
colnames(taxa_sp)=gsub(".*s__","",colnames(taxa_sp))

##############################
# ANALYSIS
##############################

##############################
# Alpha diversity & richness 
##############################

my_alpha_diversity=data.frame(Shannon=diversity(taxa_sp, index = "shannon"), Simpson=diversity(taxa_sp, index = "simpson"))
richness <- ddply(taxa_sp,~rownames(taxa_sp),function(x) {
  data.frame(RICHNESS=sum(x[-1]>0, na.rm = T))
})
row.names(richness)<-richness$`rownames(taxa_sp)`
richness$`rownames(taxa_sp)`=NULL
my_alpha_diversity<-merge(richness, my_alpha_diversity, by="row.names")
row.names(my_alpha_diversity)<-my_alpha_diversity$Row.names
my_alpha_diversity$Row.names=NULL
my_alpha_diversity$Simpson=NULL

# Association of alpha diversity with all phenotypes 
shannon_diversity_simple_linear <- linear_model_taxa_simple(metadata, "ID", my_alpha_diversity, c("infant_sex", "antibiotics", "antibiotic_concentration", "pre_preg_BMI", "gestational_age",  "gravida", "para", "infant_birthweight", "infant_sex","feeding_mode"))

setwd("/Users/trishlasinha/Desktop/CS_Baby_Biome/submission/2024_submission/supplementary/tables/individual/")
write.table(shannon_diversity_simple_linear, "all_cohorts_shannon_diversity_all_phenotypes.txt", sep="\t", row.names=F, quote = F)

# Association of alpha diversity with AB correcting for feeding mode as this has a influence on alpha diversity
shannon_diversity_linear_cor_feeding <- linear_model_taxa_cor_feeding(metadata, "ID", my_alpha_diversity, c("infant_sex", "antibiotics", "antibiotic_concentration", "pre_preg_BMI", "gestational_age",  "gravida", "para", "infant_birthweight", "infant_sex"))
write.table(shannon_diversity_linear_cor_feeding, "all_cohorts_shannon_diversity_all_phenotypes_cor_feeding.txt", sep="\t", row.names=F, quote = F)

merged<-merge(my_alpha_diversity, metadata, by="row.names")
row.names(merged)<-merged$Row.names

setwd("/Users/trishlasinha/Desktop/CS_Baby_Biome/submission/2024_submission/figures/figures_individual/figure_3/")
#pdf('all_cohorts_shannon_diversity_combined_feeding_mode.pdf', width=5, height=5)
shannon_feeding<-ggplot(merged, aes(feeding_mode, y = Shannon, fill = feeding_mode, color = feeding_mode)) +
  scale_fill_manual(values = wes_palette("Moonrise3", n = 5)) +
  scale_color_manual(values = wes_palette("Moonrise3", n = 5)) +
  geom_boxplot(alpha = 0.4, outlier.colour = NA) +
  geom_point(alpha = 0.6, position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0)) +
  ggtitle("") +
  theme_bw() +
  labs(x = "", y = "Shannon Diversity Index") +
  theme(
    axis.title.y = element_text(color = "black", size = 18, ),
    axis.text.y = element_text( size = 10),
    axis.text.x = element_text(size = 18, angle = 60, hjust = 1)
  ) +
  scale_x_discrete(labels = c("Breast feeding", "Mixed feeding", "Formula feeding")) +
  guides(fill = FALSE, color = FALSE)


shannon_feeding

# Specify the groups for comparison (breastfeeding and formula feeding)

#pdf('all_cohorts_shannon_diversity_combined_feeding_mode_antibiotics.pdf', width=5, height=5)


shannon_ab<-ggplot(merged, aes(antibiotics, y = Shannon, fill = antibiotics, color = antibiotics)) +
  scale_fill_manual(values = wes_palette("GrandBudapest1", n = 2)) +
  scale_color_manual(values = wes_palette("GrandBudapest1", n = 2)) +
  geom_boxplot(alpha = 0.4, outlier.colour = NA) +
  geom_point(alpha = 0.6, position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0)) +
  ggtitle("") +
  theme_bw() +
  labs(x = "", y = "Shannon Diversity Index") +
  theme(
    axis.title.y = element_text(color = "black", size = 18),
    axis.text.y = element_text( size = 10),
    axis.text.x = element_text(size = 18, angle = 60, hjust = 1),
  ) +
  guides(fill = FALSE, color = FALSE) +
  scale_x_discrete(labels = c("no" = "AB-", "yes" = "AB+")) 
shannon_ab
#dev.off()


##############################
# Beta diversity 
##############################
taxa_sp_filt = taxa_sp[,(colSums(taxa_sp>0.05)/nrow(taxa_sp))>0.1]
my_pseudocount_normal=min(taxa_sp_filt[taxa_sp_filt!=0])/2
taxa_sp_filt_CLR <- decostand(taxa_sp_filt, "clr", pseudocount=my_pseudocount_normal)
my_distance_CLR=vegdist(taxa_sp_filt, method = "aitchison", pseudocount=my_pseudocount_normal) 
mypcoa_CLR=cmdscale(my_distance_CLR, k = 20, eig = T)
my_var_CLR=round(mypcoa_CLR$eig*100/sum(mypcoa_CLR$eig),2)[1:20]
barplot (my_var_CLR)
mypcoa_df_CLR=as.data.frame(mypcoa_CLR$points)
names(mypcoa_df_CLR) <- c('PC1','PC2','PC3','PC4','PC5', 'PC6','PC7','PC8','PC9','PC10',
                          'PC11','PC12','PC13','PC14','PC15', 'PC16','PC17','PC18','PC19','PC20')
phenos_imputed_div_pc_CLR=merge(mypcoa_df_CLR,merged, by="row.names")
row.names (phenos_imputed_div_pc_CLR)<-phenos_imputed_div_pc_CLR$Row.names

pdf('all_cohorts_beta_diversity_feeding_mode.pdf', width=6, height=5)

overall<-ggplot(phenos_imputed_div_pc_CLR ,aes(PC1,PC2, color = feeding_mode))+
  geom_point(size = 2,alpha = 0.8)+
  stat_ellipse(aes(group = feeding_mode, fill = feeding_mode, color = feeding_mode) ,type = "norm",linetype = 2,geom = "polygon",alpha = 0.05,show.legend = F)+
  xlab(paste("PCo1=",round(my_var_CLR[1],digits = 2),"%",sep = ""))+
  ylab(paste("PCo2=",round(my_var_CLR[2],digits = 2),"%",sep = ""))+
  scale_fill_manual(values = wes_palette("Moonrise3", n = 5))+
  scale_color_manual(name=NULL, 
                     breaks = c("breast_feeding", "mixed_feeding" ,"formula_feeding"),
                     labels = c("Breast"          , "Mixed", "Formula"),
                     values = wes_palette("Moonrise3", n = 5))+
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.key = element_rect(fill = NA), 
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        panel.background = element_rect(fill = NA))
overall<-ggExtra::ggMarginal(overall, type = "histogram", groupColour = F, groupFill = TRUE,
                             xparams = list(bins = 60, alpha = 0.8,position = 'identity', color = 'white'),
                             yparams = list(bins = 60, alpha = 0.8,position = 'identity', color = 'white'))


overall

dev.off()

setwd("/Users/trishlasinha/Desktop/CS_Baby_Biome/submission/2024_submission/figures/figures_individual/figure_3/")

pdf('all_cohorts_beta_diversity.pdf', width=6, height=5)


cohort<-ggplot(phenos_imputed_div_pc_CLR ,aes(PC1,PC2, color = cohort))+
  geom_point(size = 2,alpha = 0.8)+
  stat_ellipse(aes(group = cohort, fill = cohort, color = cohort) ,type = "norm",linetype = 2,geom = "polygon",alpha = 0.05,show.legend = F)+
  xlab(paste("PCo1=",round(my_var_CLR[1],digits = 2),"%",sep = ""))+
  ylab(paste("PCo2=",round(my_var_CLR[2],digits = 2),"%",sep = ""))+
  scale_fill_manual(values=c("#e60049", "#0bb4ff", "#50e991"))+
  scale_color_manual(name=NULL, 
                     breaks = c("amsterdam", "cs_baby_biome" ,"llnext"),
                     labels = c("MAMI trial"  , "CS Baby Biome", "LLNEXT"),
                     values = c("#e60049", "#0bb4ff", "#50e991"))+
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.key = element_rect(fill = NA), 
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        panel.background = element_rect(fill = NA))
cohort<-ggExtra::ggMarginal(cohort, type = "histogram", groupColour = F, groupFill = TRUE,
                             xparams = list(bins = 60, alpha = 0.8,position = 'identity', color = 'white'),
                             yparams = list(bins = 60, alpha = 0.8,position = 'identity', color = 'white'))


cohort

dev.off()



setwd("/Users/trishlasinha/Desktop/CS_Baby_Biome/figures/figure_3/")

pdf('all_ab_beta_diversity.pdf', width=6, height=5)

ab<-ggplot(phenos_imputed_div_pc_CLR ,aes(PC1,PC2, color = antibiotics))+
  geom_point(size = 2,alpha = 0.8)+
  stat_ellipse(aes(group = antibiotics, fill = antibiotics, color = antibiotics) ,type = "norm",linetype = 2,geom = "polygon",alpha = 0.05,show.legend = F)+
  xlab(paste("PCo1=",round(my_var_CLR[1],digits = 2),"%",sep = ""))+
  ylab(paste("PCo2=",round(my_var_CLR[2],digits = 2),"%",sep = ""))+
  scale_fill_manual(values = c("#54b7a6", "#fa2b13"))+
  scale_color_manual(name=NULL, 
                     breaks = c("no", "yes"),
                     labels = c("AB-"          , "AB+"),
                     values = c("#54b7a6", "#fa2b13"))+
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.key = element_rect(fill = NA), 
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        panel.background = element_rect(fill = NA))
ab_extra<-ggExtra::ggMarginal(ab, type = "histogram", groupColour = F, groupFill = TRUE,
                            xparams = list(bins = 60, alpha = 0.8,position = 'identity', color = 'white'),
                            yparams = list(bins = 60, alpha = 0.8,position = 'identity', color = 'white'))


ab_extra
dev.off()



Distance<-as.matrix(my_distance_CLR)
phenotypes<-metadata[,c("infant_sex", "antibiotics", "antibiotic_concentration", "pre_preg_BMI", "gestational_age",  "gravida", "para", "infant_birthweight", "infant_sex","feeding_mode", "cohort")]

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
  
  ad1<-adonis2(distmat_cleaned ~ phenos2[[i]],
               permutations=10000,parallel=8, na.action = na.fail,
               by="margin") 
  print(ad1)
  adonis_results[i,] <- c(ad1$Df[1], ad1$F[1], ad1$R2[1],ad1$'Pr(>F)'[1])
}


setwd("~/Desktop/CS_Baby_Biome/submission/2024_submission/supplementary/tables/individual/")
write.table(adonis_results, "overall_composition_adonis_results_all_cohorts.txt", sep="\t", row.names=T, quote = F)

# Running adonis for all other factors correcting for feeding mode 
adonis_results_cor_feeding<- data.frame(matrix(ncol =4, nrow= ncol(phenotypes)))  
rownames(adonis_results_cor_feeding) <- colnames(phenotypes)
colnames(adonis_results_cor_feeding) <- c("Df", "F", "R2", "p-value")


adon<-foreach(i=1:ncol(phenotypes),.combine=rbind)%do%{
  #k<-which(!is.na(phenotypes[,i]))
  r<-which(is.na(phenotypes[,i]))
  r2 <- row.names(phenotypes)[r]
  distmat_cleaned <- Distance[!(rownames(Distance) %in% r2),
                              !(colnames(Distance)) %in% r2]  
  phenos2 <- phenotypes[rownames(phenotypes) %in% rownames(distmat_cleaned),]
  
  ad1<-adonis2(distmat_cleaned ~ phenos2$feeding_mode + phenos2[[i]],
               permutations=10000,parallel=8, na.action = na.fail,
               by="margin") 
  print(ad1)
  adonis_results_cor_feeding[i,] <- c(ad1$Df[2], ad1$F[2], ad1$R2[2],ad1$'Pr(>F)'[2])
}

setwd("~/Desktop/CS_Baby_Biome/submission/2024_submission/supplementary/tables/individual/")
write.table(adonis_results_cor_feeding, "overall_composition_adonis_results_all_cohorts_cor_feeding.txt", sep="\t", row.names=T, quote = F)


##############################
# Taxanomic analysis
##############################

taxa_simple_linear <- linear_model_taxa_simple(metadata, "ID", taxa_sp_filt_CLR, c("cohort" ,"infant_sex", "antibiotics", "antibiotic_concentration", "pre_preg_BMI", "gestational_age",  "gravida", "para", "infant_birthweight", "infant_sex","feeding_mode"))
setwd("~/Desktop/CS_Baby_Biome/submission/2024_submission/supplementary/tables/individual/")
write.table(taxa_simple_linear, "all_cohorts_taxa_simple_linear.txt", sep = "\t", quote = F, row.names = F)

taxa_cor_feeding_linear<- linear_model_taxa_cor_feeding(metadata, "ID", taxa_sp_filt_CLR, c("cohort" ,"infant_sex", "antibiotics", "antibiotic_concentration", "pre_preg_BMI", "gestational_age",  "gravida", "para", "infant_birthweight", "infant_sex"))
setwd("~/Desktop/CS_Baby_Biome/submission/2024_submission/supplementary/tables/individual/")
write.table(taxa_cor_feeding_linear, "all_cohorts_taxa_cor_feeding_linear.txt", sep = "\t", quote = F, row.names = F)



all<-merge(taxa_sp_filt_CLR, metadata, by="row.names")

pdf('./taxa_associations_AB/top_taxa_sig_AB_cor_cohort_feeding_B.dentium.pdf', width=6, height=5)
ggplot(all, aes(antibiotics, y=Bifidobacterium_dentium, fill=antibiotics, color=antibiotics)) +
  scale_fill_manual(values = wes_palette("GrandBudapest1", n = 2)) +
  scale_color_manual(values = wes_palette("GrandBudapest1", n = 2)) +
  geom_boxplot(alpha=0.4, outlier.colour = NA) +
  geom_point(alpha=0.6, position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0)) +
  ggtitle("") +
  theme_bw() +
  labs(x="Antibiotics", y = "CLR trans Bifidobacterium dentium") +
  theme(
    plot.title = element_text(color="black", size=18, face="bold"),
    axis.title.x = element_text(color="black", size=18, face="bold"),
    axis.title.y = element_text(color="black", size=18, face="bold"),
    axis.text.y = element_text(face="bold", size=10),
    axis.text.x = element_text(size=18, angle = 60, hjust = 1),
    legend.position = "none" # remove the legend
  ) +
  stat_n_text(
    y.pos = -5, # specify where in the y-axis the sample size should be denoted
    color = "black", # choose any color
    text.box = TRUE
  ) +
  annotate("text", x = 1.7, y = max(all$Bifidobacterium_dentium), label = "p = 0.01", hjust = 1, vjust = 1, size = 6) # add the p-value
dev.off()


##############################
# AR resistance genes 
##############################

AR_GENES<-read.delim("~/Desktop/CS_Baby_Biome/submission/2024_submission/analysis/shortbred_CARD_merged_01_04_2024.csv")
AR_GENES$NG_ID<-substr(AR_GENES$Sample, 0, 13) 
AR_GENES$NG_ID<-gsub("_kneadda", "", AR_GENES$NG_ID)
AR_GENES$NG_ID<-gsub("_", "", AR_GENES$NG_ID)
AR_GENES<-AR_GENES[!duplicated(AR_GENES$NG_ID),]
row.names(AR_GENES)<-AR_GENES$NG_ID
AR_GENES$Sample=NULL
AR_GENES$NG_ID=NULL

metadata_infants<-metadata # Only infant metadata
row.names(metadata_infants)<-metadata_infants$ID_seq_AR
AR_GENES=AR_GENES[row.names(AR_GENES)%in% rownames(metadata_infants),] # Only taxa present in infants 
AR_GENES$total_AB_load<-rowSums(AR_GENES)
AR_GENES <- AR_GENES[match(rownames(metadata_infants), rownames(AR_GENES)),]

AR_load<-AR_GENES %>% select(total_AB_load)

#AR_load_mixed_all <- mixed_models(metadata_infants, "ID", AR_load, c("cohort"))
AR_load_all <- linear_model_taxa_simple(metadata_infants, "ID", AR_load, c("infant_sex", "antibiotics", "antibiotic_concentration", "pre_preg_BMI", "gestational_age",  "gravida", "para", "infant_birthweight", "infant_sex","feeding_mode"))
AR_load_all_cor_feeding <- linear_model_taxa_cor_feeding(metadata_infants, "ID", AR_load, c("infant_sex", "antibiotics", "antibiotic_concentration", "pre_preg_BMI", "gestational_age",  "gravida", "para", "infant_birthweight", "infant_sex","feeding_mode"))
setwd("~/Desktop/CS_Baby_Biome/submission/2024_submission/supplementary/tables/individual/")
write.table(AR_load_all, "all_cohorts_AB_load_linear.txt", sep = "\t", quote = F, row.names = F)

AR_load_PRC <- linear_model_taxa_simple(metadata_infants, "ID", PRC_load, c("infant_sex", "antibiotics", "antibiotic_concentration", "pre_preg_BMI", "gestational_age",  "gravida", "para", "infant_birthweight", "infant_sex","feeding_mode"))
AR_load_LAQ <- linear_model_taxa_simple(metadata_infants, "ID", LAQ_1_load, c("infant_sex", "antibiotics", "antibiotic_concentration", "pre_preg_BMI", "gestational_age",  "gravida", "para", "infant_birthweight", "infant_sex","feeding_mode"))
AR_load_CTX_M_3 <- linear_model_taxa_simple(metadata_infants, "ID", CTX_M_3, c("infant_sex", "antibiotics", "antibiotic_concentration", "pre_preg_BMI", "gestational_age",  "gravida", "para", "infant_birthweight", "infant_sex","feeding_mode"))
AR_load_CMY <- linear_model_taxa_simple(metadata_infants, "ID", CMY, c("infant_sex", "antibiotics", "antibiotic_concentration", "pre_preg_BMI", "gestational_age",  "gravida", "para", "infant_birthweight", "infant_sex","feeding_mode"))
AR_load_TEM_1 <- linear_model_taxa_simple(metadata_infants, "ID", TEM_1, c("infant_sex", "antibiotics", "antibiotic_concentration", "pre_preg_BMI", "gestational_age",  "gravida", "para", "infant_birthweight", "infant_sex","feeding_mode"))

df<-merge(metadata_infants, AR_load, by="row.names")
df$antibiotics<-as.factor(df$antibiotics)



pdf('all_cohorts_AB_load.pdf', width=5, height=5)
AR_load_all<-ggplot(df, aes(antibiotics, y = total_AB_load, fill = antibiotics, color = antibiotics)) +
  scale_fill_manual(values = wes_palette("GrandBudapest1", n = 2)) +
  scale_color_manual(values = wes_palette("GrandBudapest1", n = 2)) +
  geom_boxplot(alpha = 0.4, outlier.colour = NA) +
  geom_point(alpha = 0.6, position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0)) +
  ggtitle("") +
  theme_bw() +
  labs(x = "", y = "Total AR load") +
  scale_x_discrete(labels = c("no" = "AB-", "yes" = "AB+")) +
  theme(
    axis.title.y = element_text(color = "black", size = 18, ),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 18, angle = 60, hjust = 1),
    legend.position = "none"
  )
AR_load_all

dev.off()

#### Figures 3 A, B, C, D, E, F 

setwd("/Users/trishlasinha/Desktop/CS_Baby_Biome/submission/2024_submission/figures/figures_individual/figure_3/")

pdf('figure_3_C&D.pdf', width=12, height=8)
Figure_3_C_D<-ggarrange(overall,ab_extra,
                      labels = c("C", "D"), 
                      ncol = 2, nrow = 1)
Figure_3_C_D
dev.off()

pdf('figure_3_E_F_G.pdf', width=12, height=8)
Figure_3_E_F_G<-ggarrange(shannon_feeding, shannon_ab, AR_load_all,
          labels = c("E", "F", "G"), 
         ncol = 3, nrow = 1)
Figure_3_E_F_G
dev.off()
