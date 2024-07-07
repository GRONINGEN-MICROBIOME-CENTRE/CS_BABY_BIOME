# Title: "Taxa_associations_phenotypes"
# Author: "Trishla Sinha"
# Date: "14/03/2023"
# Last update: "01/04/2023"

  
#Loading packages
rm (list = ls())

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
library(pheatmap)


#Loading functions

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

# Loading phenotypic data 

metadata<-read.delim("~/Desktop/CS_Baby_Biome/submission/2024_submission/analysis/Metadata_EGA_CS_BABY_BIOME.txt")
metadata[sapply(metadata, is.character)] <- lapply(metadata[sapply(metadata, is.character)],  #convert character columns to factors
                                                   as.factor)
metadata$Timepoint_categorical<-factor(metadata$Timepoint_categorical, levels = c("W01", "W02", "W03", "W04", "W05", "W06", "M06", "M12"))
metadata$growth_p_limited<-factor(metadata$growth_p_limited, levels = c("P<10", "P10-P50", "P51-90", ">P90"))
metadata$feeding_mode<-factor(metadata$feeding_mode, levels = c("breast_feeding", "mixed_feeding", "formula_feeding"))


metadata_infants<-metadata %>% drop_na(Timepoint_numeric) # Only infant metadata
metadata_infants <- metadata_infants[!(metadata_infants$Timepoint_categorical %in% c("M06", "M12")), ]
row.names(metadata_infants)<-metadata_infants$bioSampleId

#FILTERING TAXA AT SUB SPECIES, SPECIES, GENUS LEVELS
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


#FILTERING TAXA
taxa_sp$CS_Baby_Biome_ID<-substr(row.names (taxa_sp), 1, 5)
unique_counts <- sapply(taxa_sp, function(x) length(unique(taxa_sp$CS_Baby_Biome_ID[x >0.05]))) # RA of 0.05%
taxa_sp_filt <- taxa_sp[, unique_counts >= 2] # In atleast 2 infants 
taxa_sp_filt$CS_Baby_Biome_ID=NULL


# Transform data 
my_pseudocount_normal=min(taxa_sp_filt[taxa_sp_filt!=0])/2
taxa_sp_filt<-taxa_sp_filt[match(row.names(taxa_sp_filt),row.names(metadata_infants)),]
taxa_sp_filt_CLR<-decostand(taxa_sp_filt, "clr", pseudocount=my_pseudocount_normal)

metadata_infants$ID<-row.names(metadata_infants)

#Associations of phenotype in species 
taxa_mixed_all <- mixed_models(metadata_infants, "ID", taxa_sp_filt_CLR, c("cefazoline_measurement_mg_L", "rand_AB", "pre_preg_bmi_mother",   "preg_gest_age", "preg_weight_gain", "gravida", "para", "infant_birthweight", "infant_sex", "mother_age_at_delivery","APGAR_1", "APGAR_5", "growth_p_limited", "feeding_mode","living_situation", "cats_dogs"))

setwd("/Users/trishlasinha/Desktop/CS_Baby_Biome/submission/2024_submission/supplementary/tables/individual/")
write.table(taxa_mixed_all, "species_associations_mixed_all_raw_01_04_2024.txt", sep="\t", row.names=F, quote = F)

# Association of alpha diversity with AB correcting for feeding mode as this has an influence on multiple taxa
taxa_mixed_cor_feeding <- mixed_models_cor_feeding(metadata_infants, "ID", taxa_sp_filt_CLR, c("cefazoline_measurement_mg_L", "rand_AB", "pre_preg_bmi_mother",   "preg_gest_age", "preg_weight_gain", "gravida", "para", "infant_birthweight", "infant_sex", "mother_age_at_delivery","APGAR_1", "APGAR_5", "growth_p_limited", "feeding_mode","living_situation", "cats_dogs"))

setwd("/Users/trishlasinha/Desktop/CS_Baby_Biome/submission/2024_submission/supplementary/tables/individual/")
taxa_mixed_cor_feeding <- mixed_models_cor_feeding(metadata_infants, "ID", taxa_sp_filt_CLR, c("cefazoline_measurement_mg_L", "rand_AB", "pre_preg_bmi_mother",   "preg_gest_age", "preg_weight_gain", "gravida", "para", "infant_birthweight", "infant_sex", "mother_age_at_delivery","APGAR_1", "APGAR_5", "growth_p_limited", "feeding_mode","living_situation", "cats_dogs"))

write.table(taxa_mixed_cor_feeding, "species_associations_mixed_cor_feeding_01_04_2024.txt", sep="\t", row.names=F, quote = F)
write.table(taxa_mixed_all, "species_associations_mixed_all_raw_01_04_2024.txt", sep="\t", row.names=F, quote = F)


########################### PHEATMAP ###############################
# Subset rows with AB info 
subset_df <- taxa_mixed_cor_feeding[grepl("rand_ABYes", taxa_mixed_cor_feeding$Feature), ]
subset_df$FDR<-p.adjust(subset_df$P, method = "BH")
signicant<-subset_df [subset_df$P<0.05,]
common_cols <- intersect(names(taxa_sp_filt_CLR), signicant$Bug)
taxa_signicant<- subset(taxa_sp_filt_CLR, select = common_cols)
mat=as.matrix(taxa_signicant)
met<-metadata_infants[,c("rand_AB","CS_BABY_BIOME_ID", "feeding_mode_pragmatic" )]
met$rand_AB=as.factor(met$rand_AB)
met$CS_BABY_BIOME_ID=as.factor(met$CS_BABY_BIOME_ID)
met$feeding_mode_pragmatic=as.factor(met$feeding_mode_pragmatic)

write.table(mat, "taxa.txt", sep = "\t", quote = F, row.names = T)
write.table(met, "metadata.txt", sep = "\t", quote = F, row.names = T)

all=merge(met,mat, by="row.names")
row.names(all)<-all$Row.names
all$Row.names=NULL
all_sorted<-all[order(all[,'rand_AB']), ]

annotation_2<- select(filter(all_sorted),c( "rand_AB", "feeding_mode_pragmatic"))
names (annotation_2)[1]<-"AB"
names (annotation_2)[2]<-"feeding_mode"
levels(annotation_2$feeding_mode) <- c("breast_feeding", "mixed_feeding", "formula_feeding") 
#levels(annotation_2$AB) <- c("0", "1") 
mat=all_sorted[,c(4:15)]
my_colors = c("#54b7a6", "#fa2b13")
my_colors_2 = c("#85d4e2", "#f4b5bc", "#9c964a")

#mycolors_all <- list(AB = my_colors, feeding_mode=my_colors_2)

mycolors_all <- list(AB = c("No" = my_colors[1], "Yes" = my_colors[2]),
                     feeding_mode = setNames(my_colors_2, levels(annotation_2$feeding_mode)))

#names(my_colors) <- unique(annotation_2$AB)

setwd("/Users/trishlasinha/Desktop/CS_Baby_Biome/submission/2024_submission/figures/figures_individual/figure_2/")
pdf('species_AB_cs_baby_biome.pdf', width=8, height=10)
pheatmap(mat,
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         annotation_row = annotation_2, 
         show_rownames = FALSE,
         annotation_colors = mycolors_all,
         angle_col = 90,   # Rotate x-axis labels by degrees
         main = "Species significantly associated with AB group")


#write.table(taxa_sp_filt_CLR, "Species_names_for_associations_CS_BABY_BIOME_07_05_2023.txt", sep="\t", quote = F)

dev.off()


