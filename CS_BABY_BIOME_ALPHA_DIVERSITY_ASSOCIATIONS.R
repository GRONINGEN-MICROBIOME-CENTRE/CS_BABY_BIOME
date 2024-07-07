# Title: "CS_BABY_BIOME_ALPHA_DIVERSITY_ANALYSIS"
# Author: "Trishla Sinha"
# Date: "17/02/2023"
# Last update: "01/04/2024"

rm (list = ls())
setwd("/Users/trishlasinha/Desktop/CS_Baby_Biome/submission/2024_submission/analysis/RESULTS/")


#load packages 
library(tidyverse)
library(stringr)
library(vegan)
library(RColorBrewer)
library(wesanderson)
library(reshape2)
library(dplyr)
library(ggpubr)
library(wesanderson)
library(ggpubr)
library(lmerTest)
library(plyr)


# Load functions
mixed_models <- function(metadata, ID, CLR_transformed_data, pheno_list) {
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
    #Prevalence = sum(as.numeric(as_vector(select(df, Bug)) > 0)) / dim(df)[1]
    # print (c(Bug, Prevalence))
    Bug2 = paste(c("`",Bug, "`"), collapse="")
    for ( pheno in pheno_list){
      pheno2 = paste(c("`",pheno, "`"), collapse="")
      df[is.na(df[colnames(df) == pheno]) == F, ID] -> To_keep
      df_pheno = filter(df, !!sym(ID) %in% To_keep )
      Model0 = as.formula(paste( c(Bug2,  " ~ read_depth + DNA_concentration_ng_ul + Timepoint_numeric + (1|CS_BABY_BIOME_ID)"), collapse="" )) 
      lmer(Model0, df_pheno) -> resultmodel0
      base_model=resultmodel0
      Model2 = as.formula(paste( c(Bug2,  " ~ read_depth + DNA_concentration_ng_ul + Timepoint_numeric + ",pheno2, "+ (1|CS_BABY_BIOME_ID)"), collapse="" ))
      lmer(Model2, df_pheno, REML = F) -> resultmodel2
      M = "Mixed"
      as.data.frame(anova(resultmodel2, base_model))['resultmodel2','Pr(>Chisq)']->p_simp
      as.data.frame(summary(resultmodel2)$coefficients)[grep(pheno, row.names(as.data.frame(summary(resultmodel2)$coefficients))),] -> Summ_simple
      Summ_simple %>% rownames_to_column("Feature") %>% as_tibble() %>% mutate(P = p_simp, Model_choice = M, Bug =Bug, Pheno=pheno, Model="simple") -> temp_output
      rbind(Overall_result_phenos, temp_output) -> Overall_result_phenos
    }
  }
  
  p=as.data.frame(Overall_result_phenos)
  p$FDR<-p.adjust(p$P, method = "BH")
  
  return(p)
  
}

mixed_models_cor_feeding <- function(metadata, ID, CLR_transformed_data, pheno_list) {
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
    #Prevalence = sum(as.numeric(as_vector(select(df, Bug)) > 0)) / dim(df)[1]
    # print (c(Bug, Prevalence))
    Bug2 = paste(c("`",Bug, "`"), collapse="")
    for ( pheno in pheno_list){
      pheno2 = paste(c("`",pheno, "`"), collapse="")
      df[is.na(df[colnames(df) == pheno]) == F, ID] -> To_keep
      df_pheno = filter(df, !!sym(ID) %in% To_keep )
      Model0 = as.formula(paste( c(Bug2,  " ~ read_depth + DNA_concentration_ng_ul + Timepoint_numeric + feeding_mode + (1|CS_BABY_BIOME_ID)"), collapse="" )) 
      lmer(Model0, df_pheno) -> resultmodel0
      base_model=resultmodel0
      Model2 = as.formula(paste( c(Bug2,  " ~ read_depth + DNA_concentration_ng_ul + Timepoint_numeric + feeding_mode+",pheno2, "+ (1|CS_BABY_BIOME_ID)"), collapse="" ))
      lmer(Model2, df_pheno, REML = F) -> resultmodel2
      M = "Mixed"
      as.data.frame(anova(resultmodel2, base_model))['resultmodel2','Pr(>Chisq)']->p_simp
      as.data.frame(summary(resultmodel2)$coefficients)[grep(pheno, row.names(as.data.frame(summary(resultmodel2)$coefficients))),] -> Summ_simple
      Summ_simple %>% rownames_to_column("Feature") %>% as_tibble() %>% mutate(P = p_simp, Model_choice = M, Bug =Bug, Pheno=pheno, Model="simple") -> temp_output
      rbind(Overall_result_phenos, temp_output) -> Overall_result_phenos
    }
  }
  
  p=as.data.frame(Overall_result_phenos)
  p$FDR<-p.adjust(p$P, method = "BH")
  
  return(p)
  
}

#MERGING ALL PHENOTYPIC DATA TOGETHER 

metadata<-read.delim("~/Desktop/CS_Baby_Biome/submission/2024_submission/analysis/Metadata_EGA_CS_BABY_BIOME.txt")
metadata_late<-read.delim("~/Desktop/CS_Baby_Biome/OLD_ANALYSIS/metadata_cs_baby_biome_later_tp.txt")

metadata[sapply(metadata, is.character)] <- lapply(metadata[sapply(metadata, is.character)],  #convert character columns to factors
                                                   as.factor)
metadata$Timepoint_categorical<-factor(metadata$Timepoint_categorical, levels = c("W01", "W02", "W03", "W04", "W05", "W06", "M06", "M12"))
metadata$growth_p_limited<-factor(metadata$growth_p_limited, levels = c("P<10", "P10-P50", "P51-90", ">P90"))
metadata$feeding_mode<-factor(metadata$feeding_mode, levels = c("breast_feeding", "mixed_feeding", "formula_feeding"))
taxa <-read.delim("~/Desktop/CS_Baby_Biome/submission/2024_submission/analysis/TAXA_MP4_2022/taxa_clean_mp4_oct_2022_01_04_2024.txt")


# FILTERING TAXA AT SUB SPECIES, SPECIES & GENUS LEVEL 
metadata_infants<-metadata %>% drop_na(Timepoint_numeric) # Only infant metadata
metadata_infants <- metadata_infants[!(metadata_infants$Timepoint_categorical %in% c("M06", "M12")), ]

row.names(metadata_infants)<-metadata_infants$bioSampleId
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


taxa_sp_filt = taxa_sp[,(colSums(taxa_sp!=0)/nrow(taxa_sp))>0.000]
taxa_genus_filt = taxa_genus[,(colSums(taxa_genus!=0)/nrow(taxa_genus))>0.00]
taxa_sub_species_filt = taxa_sub_species[,(colSums(taxa_sub_species!=0)/nrow(taxa_sub_species))>0.00]

table (metadata_infants$rand_AB)


#### ASSOCIATIONS WITH ALPHA DIVERSITY & RICHNESS ####

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
metadata_infants$ID<-row.names(metadata_infants)
my_alpha_diversity$RICHNESS=as.numeric(my_alpha_diversity$RICHNESS)

# Association of alpha diversity with AB (and other phenotypes)
shannon_diversity_mixed_all <- mixed_models(metadata_infants, "ID", my_alpha_diversity, c("cefazoline_measurement_mg_L", "rand_AB", "pre_preg_bmi_mother",   "preg_gest_age", "preg_weight_gain", "gravida", "para", "infant_birthweight", "infant_sex", "mother_age_at_delivery","APGAR_1", "APGAR_5", "growth_p_limited", "feeding_mode","living_situation", "cats_dogs"))

# Association of alpha diversity with AB (and other phenotypes) correcting for feeding mode (as this is an important factor influencing alpha diversity in previous studies) 
shannon_diversity_mixed_cor_feeding <- mixed_models_cor_feeding(metadata_infants, "ID", my_alpha_diversity, c("cefazoline_measurement_mg_L", "rand_AB", "pre_preg_bmi_mother",   "preg_gest_age", "preg_weight_gain", "gravida", "para", "infant_birthweight", "infant_sex", "mother_age_at_delivery","APGAR_1", "APGAR_5", "growth_p_limited", "feeding_mode","living_situation", "cats_dogs"))

setwd("/Users/trishlasinha/Desktop/CS_Baby_Biome/submission/2024_submission/supplementary/tables/individual/")
write.table(shannon_diversity_mixed_all, "all_phenotypes_shannon_cor_rd_dna_conc_tp_diversity_richness_infants_2024.txt", sep="\t", row.names=F, quote = F)
write.table(shannon_diversity_mixed_cor_feeding, "all_phenotypes_cor_rd_dna_conc_tp_feeding_shannon_diversity_richness_infants_2024.txt", sep="\t", row.names=F, quote = F)

all<-merge(my_alpha_diversity, metadata_infants, by="row.names")

shannon_AB_CS_Baby_Biome <- ggplot(all, aes(Timepoint_categorical, y=Shannon, fill=rand_AB, color=rand_AB)) +
  scale_fill_manual(values = c("#54b7a6", "#fa2b13"), labels = c("AB-", "AB+")) +
  scale_color_manual(values = c("#54b7a6", "#fa2b13"), labels = c("AB-", "AB+")) +
  geom_boxplot(alpha=0.4, outlier.colour = NA) +
  geom_point(alpha=0.6,
             position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0)) +
  ggtitle("") +
  theme_bw() + labs(x="", y = "Shannon Diversity Index", fill="", color="") +
  theme(
    axis.title.y = element_text(color="black", size=18),
    axis.text.y = element_text(size=12),
    axis.text.x = element_text(size=12, angle = 60, hjust = 1))


richness_AB_CS_Baby_Biome <- ggplot(all, aes(Timepoint_categorical, y=RICHNESS, fill=rand_AB, color=rand_AB)) +
  scale_fill_manual(values = c("#54b7a6", "#fa2b13"), labels = c("AB-", "AB+")) +
  scale_color_manual(values = c("#54b7a6", "#fa2b13"), labels = c("AB-", "AB+")) +
  geom_boxplot(alpha=0.4, outlier.colour = NA) +
  geom_point(alpha=0.6,
             position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0)) +
  ggtitle("") +
  theme_bw() + labs(x="", y = "Species Richness", fill="", color="") +
  theme(
    axis.title.y = element_text(color="black", size=18),
    axis.text.y = element_text(size=12),
    axis.text.x = element_text(size=12, angle = 60, hjust = 1))

# Saving part 1, and B of figure 1 
setwd("/Users/trishlasinha/Desktop/CS_Baby_Biome/submission/2024_submission/figures/figures_individual/figure_1/")
save(shannon_AB_CS_Baby_Biome, richness_AB_CS_Baby_Biome, file="figure_1_B_C.RData")
shannon_AB_CS_Baby_Biome

### LATER TIMEPOINTS & ALPHA DIVERSITY & RICHNESS ###
rm (list = ls())
metadata<-read.delim("~/Desktop/CS_Baby_Biome/submission/2024_submission/analysis/Metadata_EGA_CS_BABY_BIOME.txt")
#metadata_late<-read.delim("~/Desktop/CS_Baby_Biome/OLD_ANALYSIS/metadata_cs_baby_biome_later_tp.txt")

metadata[sapply(metadata, is.character)] <- lapply(metadata[sapply(metadata, is.character)],  #convert character columns to factors
                                                   as.factor)

metadata$Timepoint_categorical<-factor(metadata$Timepoint_categorical, levels = c("W01", "W02", "W03", "W04", "W05", "W06", "M06", "M12"))
metadata$growth_p_limited<-factor(metadata$growth_p_limited, levels = c("P<10", "P10-P50", "P51-90", ">P90"))
metadata$feeding_mode<-factor(metadata$feeding_mode, levels = c("breast_feeding", "mixed_feeding", "formula_feeding"))

# Filtering only later timepoints 
metadata_infants <- metadata[(metadata$Timepoint_categorical %in% c("M06", "M12")), ]


# FILTERING TAXA AT SUB SPECIES, SPECIES & GENUS LEVEL 
taxa <-read.delim("~/Desktop/CS_Baby_Biome/submission/2024_submission/analysis/TAXA_MP4_2022/taxa_clean_mp4_oct_2022_01_04_2024.txt")


# Only infant metadata
row.names(metadata_infants)<-metadata_infants$bioSampleId
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


taxa_sp_filt = taxa_sp[,(colSums(taxa_sp!=0)/nrow(taxa_sp))>0.000]
taxa_genus_filt = taxa_genus[,(colSums(taxa_genus!=0)/nrow(taxa_genus))>0.00]
taxa_sub_species_filt = taxa_sub_species[,(colSums(taxa_sub_species!=0)/nrow(taxa_sub_species))>0.00]

table (metadata_infants$rand_AB)


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

all<-merge(my_alpha_diversity, metadata_infants, by="row.names")
all_m6<-all[all$Timepoint_categorical=="M06",]

result_m6_shannon<-wilcox.test(all_m6$Shannon ~ all_m6$rand_AB)
p_value_m6_shannon <-result_m6_shannon$p.value

result_m6_richness <- wilcox.test(all_m6$RICHNESS ~ all_m6$rand_AB)
p_value_m6_richness<-result_m6_richness$p.value

shannon_AB_CS_Baby_Biome_m6<-ggplot(all_m6, aes(rand_AB, y=Shannon, fill=rand_AB, color=rand_AB)) +
  scale_fill_manual(values = c("#54b7a6", "#fa2b13"))+
  scale_color_manual(values = c("#54b7a6", "#fa2b13"))+
  geom_boxplot(alpha=0.4, outlier.colour = NA)+
  geom_point(alpha=0.6,
             position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0))+
  # geom_jitter() +
  ggtitle("")+
  theme_bw()+labs(x="", y = "Shannon Diversity Index Month 6", fill="Antibiotic", color="Antibiotic")+
  theme(
    plot.title = element_text(color="black", size=12, face="bold"),
    axis.title.x = element_text(color="black", size=12, face="bold"),
    axis.title.y = element_text(color="black", size=12, face="bold"),
    axis.text.y = element_text(face="bold", size=10),
    axis.text.x = element_text(size=12, angle = 60, hjust = 1))+
  annotate("text", x = 1.5, y = max(all_m6$Shannon) - 0.0, 
           label = paste("p =", round(p_value_m6_shannon, 3)), size = 4) 
#strip.text.x = element_text(size = 10))
shannon_AB_CS_Baby_Biome_m6

richness_AB_CS_Baby_Biome_m6<-ggplot(all_m6, aes(rand_AB, y=RICHNESS, fill=rand_AB, color=rand_AB)) +
  scale_fill_manual(values = c("#54b7a6", "#fa2b13"))+
  scale_color_manual(values = c("#54b7a6", "#fa2b13"))+
  geom_boxplot(alpha=0.4, outlier.colour = NA)+
  geom_point(alpha=0.6,
             position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0))+
  # geom_jitter() +
  ggtitle("")+
  theme_bw()+labs(x="", y = "Richness Month 6", fill="Antibiotic", color="Antibiotic")+
  theme(
    plot.title = element_text(color="black", size=12, face="bold"),
    axis.title.x = element_text(color="black", size=12, face="bold"),
    axis.title.y = element_text(color="black", size=12, face="bold"),
    axis.text.y = element_text(face="bold", size=10),
    axis.text.x = element_text(size=12, angle = 60, hjust = 1))+
  annotate("text", x = 1.5, y = max(all_m6$RICHNESS) - 0.2, 
           label = paste("p =", round(p_value_m6_richness, 3)), size = 4)  
#strip.text.x = element_text(size = 10))
richness_AB_CS_Baby_Biome_m6




all_m12<-all[all$Timepoint_categorical=="M12",]

result_m12_shannon<-wilcox.test(all_m12$Shannon ~all_m12$rand_AB)
p_value_shannon_m12<-result_m12_shannon$p.value

result_m12_richness<-wilcox.test(all_m12$RICHNESS ~ all_m12$rand_AB)
p_value_richness_m12<-result_m12_richness$p.value

shannon_AB_CS_Baby_Biome_m12<-ggplot(all_m12, aes(rand_AB, y=Shannon, fill=rand_AB, color=rand_AB)) +
  scale_fill_manual(values = c("#54b7a6", "#fa2b13"))+
  scale_color_manual(values = c("#54b7a6", "#fa2b13"))+
  geom_boxplot(alpha=0.4, outlier.colour = NA)+
  geom_point(alpha=0.6,
             position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0))+
  # geom_jitter() +
  ggtitle("")+
  theme_bw()+labs(x="", y = "Shannon Diversity Index Month 12", fill="Antibiotic", color="Antibiotic")+
  theme(
    plot.title = element_text(color="black", size=12, face="bold"),
    axis.title.x = element_text(color="black", size=12, face="bold"),
    axis.title.y = element_text(color="black", size=12, face="bold"),
    axis.text.y = element_text(face="bold", size=10),
    axis.text.x = element_text(size=12, angle = 60, hjust = 1))+
  annotate("text", x = 1.5, y = max(all_m12$Shannon) - 0.1, 
           label = paste("p =", round(p_value_shannon_m12, 3)), size = 4) 
#strip.text.x = element_text(size = 10))
shannon_AB_CS_Baby_Biome_m12


richness_AB_CS_Baby_Biome_m12<-ggplot(all_m12, aes(rand_AB, y=RICHNESS, fill=rand_AB, color=rand_AB)) +
  scale_fill_manual(values = c("#54b7a6", "#fa2b13"))+
  scale_color_manual(values = c("#54b7a6", "#fa2b13"))+
  geom_boxplot(alpha=0.4, outlier.colour = NA)+
  geom_point(alpha=0.6,
             position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0))+
  # geom_jitter() +
  ggtitle("")+
  theme_bw()+labs(x="", y = "Richness Month 12", fill="Antibiotic", color="Antibiotic")+
  theme(
    plot.title = element_text(color="black", size=12, face="bold"),
    axis.title.x = element_text(color="black", size=12, face="bold"),
    axis.title.y = element_text(color="black", size=12, face="bold"),
    axis.text.y = element_text(face="bold", size=10),
    axis.text.x = element_text(size=12, angle = 60, hjust = 1))+
  annotate("text", x = 1.5, y = max(all_m12$RICHNESS) - 0.2, 
           label = paste("p =", round(p_value_richness_m12, 3)), size = 4)  
#strip.text.x = element_text(size = 10))
richness_AB_CS_Baby_Biome_m12



setwd("/Users/trishlasinha/Desktop/CS_Baby_Biome/submission/2024_submission/supplementary/figures/")
pdf('sup_shannon_diversity_richness_later_timepoints_m6_M12_cs_baby_biome.pdf', width=8, height=8)

late_sp<-ggarrange(shannon_AB_CS_Baby_Biome_m6, richness_AB_CS_Baby_Biome_m6, shannon_AB_CS_Baby_Biome_m12, richness_AB_CS_Baby_Biome_m12, common.legend = T)

late_sp

dev.off()
