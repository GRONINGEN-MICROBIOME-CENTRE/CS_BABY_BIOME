# Title: "CS_BABY_BIOME_ANTIBIOTIC_RESISTANCE_GENES_ANALYSIS"
# Author: "Trishla Sinha"
# Date: "13/02/2023"
# Last update: "01/04/2024"

library(tidyverse)
library(lmerTest)

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


setwd("/Users/trishlasinha/Desktop/CS_Baby_Biome/submission/2024_submission/analysis/RESULTS/")
metadata<-read.delim("~/Desktop/CS_Baby_Biome/submission/2024_submission/analysis/Metadata_EGA_CS_BABY_BIOME.txt")
metadata[sapply(metadata, is.character)] <- lapply(metadata[sapply(metadata, is.character)],  #convert character columns to factors
                                                   as.factor)
metadata$Timepoint_categorical<-factor(metadata$Timepoint_categorical, levels = c("W01", "W02", "W03", "W04", "W05", "W06", "M06", "M12"))
metadata$growth_p_limited<-factor(metadata$growth_p_limited, levels = c("P<10", "P10-P50", "P51-90", ">P90"))
metadata$feeding_mode<-factor(metadata$feeding_mode, levels = c("breast_feeding", "mixed_feeding", "formula_feeding")) 
  
  
AR_GENES<-read.delim("/Users/trishlasinha/Desktop/CS_Baby_Biome/submission/2024_submission/analysis/shortbred_CARD_merged_01_04_2024.csv")

AR_GENES$NG_ID<-substr(AR_GENES$Sample, 0, 13) 
AR_GENES$NG_ID<-gsub("_", "", AR_GENES$NG_ID)
AR_GENES<-AR_GENES[!duplicated(AR_GENES$NG_ID),]
row.names(AR_GENES)<-AR_GENES$NG_ID
AR_GENES$Sample=NULL
AR_GENES$NG_ID=NULL

metadata_infants<-metadata %>% drop_na(Timepoint_numeric) # Only infant metadata
metadata_infants <- metadata_infants[!(metadata_infants$Timepoint_categorical %in% c("M06", "M12")), ]

row.names(metadata_infants)<-metadata_infants$bioSampleId
metadata_infants$ID<-row.names(metadata_infants)

missing_rows <- setdiff(rownames(metadata_infants), rownames(AR_GENES))

AR_GENES=AR_GENES[row.names(AR_GENES)%in% rownames(metadata_infants),] # Only AR genes present in infants 
AR_GENES$total_AB_load<-rowSums(AR_GENES)
AR_GENES <- AR_GENES[match(rownames(metadata_infants), rownames(AR_GENES)),]
str(AR_GENES$total_AB_load)

AR_load<-AR_GENES %>% select(total_AB_load)
filtered_AR_GENES <- AR_GENES[, !(colSums(AR_GENES) == 0)]
filtered_AR_GENES <- filtered_AR_GENES[, colSums(filtered_AR_GENES != 0) >= 10]


# Run mixed models 
AR_load_mixed_all <- mixed_models(metadata_infants, "ID", AR_load, c("cefazoline_measurement_mg_L", "rand_AB", "pre_preg_bmi_mother",   "preg_gest_age", "preg_weight_gain", "gravida", "para", "infant_birthweight", "infant_sex", "mother_age_at_delivery","APGAR_1", "APGAR_5", "growth_p_limited", "feeding_mode","living_situation", "cats_dogs"))

AR_genes_mixed_all <- mixed_models(metadata_infants, "ID", filtered_AR_GENES, c("cefazoline_measurement_mg_L", "rand_AB", "pre_preg_bmi_mother",   "preg_gest_age", "preg_weight_gain", "gravida", "para", "infant_birthweight", "infant_sex", "mother_age_at_delivery","APGAR_1", "APGAR_5", "growth_p_limited", "feeding_mode","living_situation", "cats_dogs"))

setwd("/Users/trishlasinha/Desktop/CS_Baby_Biome/supplementary/tables/individual")

write.table(AR_load_mixed_all, "Total_AR_load_all_phenotypes_CS_BABY_BIOME_01_04_2024.txt", sep="\t", row.names=F, quote = F)
write.table(AR_genes_mixed_all, "AR_GENES_ALL_all_phenotypes_CS_BABY_BIOME_01_04_2024.txt", sep="\t", row.names=F, quote = F)

all<-merge(AR_load, metadata_infants, by="row.names")
all$AB <-all$rand_AB

setwd("/Users/trishlasinha/Desktop/CS_Baby_Biome/submission/2024_submission/figures/figures_individual/figure_2/")


pdf('AR_load_AB_cs_baby_biome.pdf', width=6, height=6)

ggplot(all, aes(Timepoint_categorical, y=total_AB_load, fill=AB, color=AB)) +
  scale_fill_manual(values = c("#54b7a6", "#fa2b13"))+
  scale_color_manual(values = c("#54b7a6", "#fa2b13"))+
  geom_boxplot(alpha=0.4, outlier.colour = NA)+
  geom_point(alpha=0.6,
             position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0))+
  # geom_jitter() +
  ggtitle("")+
  theme_bw()+labs(x="", y = "Total AR load")+
  theme(
    plot.title = element_text(color="black", size=22, face="bold"),
    axis.title.x = element_text(color="black", size=22, face="bold"),
    axis.title.y = element_text(color="black", size=22, face="bold"),
    axis.text.y = element_text(face="bold", size=10),
    axis.text.x = element_text(size=22, angle = 60, hjust = 1))
#strip.text.x = element_text(size = 10))
dev.off()

