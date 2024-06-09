#################### CS Baby Biome SCFA Analysis ########################

# Date: 29th January, 2024
# Last update: 15th March 2024
# Author: Trishla Sinha 

library(tidyverse)
library(vegan)
library(lmerTest)
library(wesanderson)
library(patchwork)


# Load functions 


mixed_models_timepoint <- function(metadata, ID, CLR_transformed_data, pheno_list) {
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
      Model0 = as.formula(paste( c(Bug2,  " ~ (1|CS_BABY_BIOME_ID)"), collapse="" )) 
      lmer(Model0, df_pheno) -> resultmodel0
      base_model=resultmodel0
      Model2 = as.formula(paste( c(Bug2,  " ~ ",pheno2, "+ (1|CS_BABY_BIOME_ID)"), collapse="" ))
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
      Model0 = as.formula(paste( c(Bug2,  " ~ Timepoint + (1|CS_BABY_BIOME_ID)"), collapse="" )) 
      lmer(Model0, df_pheno) -> resultmodel0
      base_model=resultmodel0
      Model2 = as.formula(paste( c(Bug2,  " ~ Timepoint + ",pheno2, "+ (1|CS_BABY_BIOME_ID)"), collapse="" ))
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
      Model0 = as.formula(paste( c(Bug2,  " ~ Timepoint + feeding_mode + (1|CS_BABY_BIOME_ID)"), collapse="" )) 
      lmer(Model0, df_pheno) -> resultmodel0
      base_model=resultmodel0
      Model2 = as.formula(paste( c(Bug2,  " ~ Timepoint+ feeding_mode+",pheno2, "+ (1|CS_BABY_BIOME_ID)"), collapse="" ))
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


# Load variables 
SCFA<-read.delim("~/Desktop/CS_Baby_Biome/submission/2024_submission/analysis/lab/SCFA_TS.txt")
SCFA<-na.omit(SCFA)
names (SCFA) [2] <-"CS_BABY_BIOME_ID"
AB<-read.delim("~/Desktop/CS_Baby_Biome/submission/2024_submission/analysis/lab/AB_samples.txt")
feeding<-read.delim("~/Desktop/CS_Baby_Biome/submission/2024_submission/analysis/lab/feeding_samples.txt")

analyses_df<-left_join(SCFA, AB)
analyses_df<-left_join(analyses_df, feeding)
row.names(analyses_df)<-analyses_df$SAMPLE_ID

analyses_df$Timepoint<-factor(analyses_df$Timepoint, levels = c("W01", "W04"))
analyses_df$feeding_mode_pragmatic<-factor(analyses_df$feeding_mode_pragmatic, levels = c("breast_feeding", "mixed_feeding", "formula_feeding"))

phenotypes<-analyses_df[, c("rand_AB", "feeding_mode_pragmatic", "CS_BABY_BIOME_ID", "Timepoint")]
names (phenotypes)[1]<-"AB"
names (phenotypes)[2]<-"feeding_mode"
phenotypes$ID<-row.names(phenotypes)
SCFA_to_analyze<-analyses_df[, c("Acetate_mM", "Acetate_umolperg", "Propionate_mM", "Propionate_umolperg", "Butyrate_mM",  "Butyrate_umolperg" )]


# Association of bile acids with AB and feeding 
SCFA_mixed_all <- mixed_models(phenotypes, "ID", SCFA_to_analyze, c("AB", "feeding_mode"))

# now correcting for feeding 
SCFA_mixed_all_cor_feeding<-mixed_models_cor_feeding(phenotypes, "ID", SCFA_to_analyze, c("AB"))

# Looking at only timepoint 
SCFA_mixed_all_timepoint<-mixed_models_timepoint(phenotypes, "ID", SCFA_to_analyze, c("Timepoint"))


# Timepoint 
long_data <-  analyses_df%>%
  pivot_longer(cols = c("Acetate_mM", "Propionate_mM",  "Butyrate_mM" ), 
               names_to = "Variable", 
               values_to = "Value")


ggplot(long_data, aes(x = Timepoint, y = Value, fill = Variable, color=Variable)) +
  geom_boxplot(alpha=0.4, outlier.colour = NA)  +
  geom_point(alpha=0.6,
             position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0)) +
  facet_wrap(~ Variable, scales = "free_y") +
  theme_bw() +
  labs(x = "Timepoint", y = "SCFA") +
  theme(legend.position = "none")


# Timepoint 
long_data <-  analyses_df%>%
  pivot_longer(cols = c("Acetate_umolperg", "Propionate_umolperg",  "Butyrate_umolperg" ), 
               names_to = "Variable", 
               values_to = "Value")


p1<-ggplot(long_data, aes(x = Timepoint, y = Value, fill = Timepoint, color=Timepoint)) +
  geom_boxplot(alpha=0.4, outlier.colour = NA)  +
  geom_point(alpha=0.6,
             position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0)) +
  facet_wrap(~ Variable, scales = "free_y") +
  theme_bw() +
  scale_color_manual(values = c("indianred2", "mediumspringgreen", "purple2" ))+
  scale_fill_manual(values = c("indianred2", "mediumspringgreen","purple2" ))+
  labs(x = "Timepoint", y = "umol/g") +
  theme(legend.position = "none")


p2<-ggplot(long_data, aes(x = Timepoint, y = Value, fill = rand_AB, color=rand_AB)) +
  geom_boxplot(alpha = 0.4, outlier.colour = NA, position = position_dodge(width = 0.8)) +
  geom_point(alpha = 0.6, position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.8)) +
  scale_fill_manual(values = c("#54b7a6", "#fa2b13")) +
  scale_color_manual(values = c("#54b7a6", "#fa2b13")) +
  theme_bw() +
  labs(x = "Timepoint", y = "umol/g") +
  theme(legend.position = "bottom") +
  facet_wrap(~ Variable, scales = "free_y", ncol = 3)



p3<-ggplot(long_data, aes(x = Timepoint, y = Value, fill = feeding_mode_pragmatic, color=feeding_mode_pragmatic)) +
  geom_boxplot(alpha = 0.4, outlier.colour = NA, position = position_dodge(width = 0.8)) +
  geom_point(alpha = 0.6, position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.8)) +
  scale_fill_manual(values = wes_palette("Moonrise3", n = 5))+
  scale_color_manual(values = wes_palette("Moonrise3", n = 5))+
  theme_bw() +
  labs(x = "Timepoint", y = "umol/g") +
  theme(legend.position = "bottom") +
  facet_wrap(~ Variable, scales = "free_y", ncol = 3)

write.table(SCFA_mixed_all, "CS_BABY_BIOME_results_mixed_models_SCFA_with_AB_feeding.txt", sep="\t", row.names=F, quote = F)


combined_plot <- p1/p2 / p3

combined_plot





