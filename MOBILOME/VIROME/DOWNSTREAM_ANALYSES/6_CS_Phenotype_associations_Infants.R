################################################################################
##### CS Baby Biome: Associations with  phenotypes in CS Baby Biome infants
### Author(s):Asier Fernández-Pato
### Last updated: 20th June, 2024
################################################################################

#****************
# Load modules
#****************
library(lme4)
library(lmerTest)
library(tidyverse)
library(dplyr)
library(ggplot2)

#****************
# Define functions
#****************

# Functions to perform the association between phenotypes and abundances of vOTUs (mixed-model)
# Different functions are defined:
#--mixed_models: Abundances vs phenotypes (correcting for DNA_conc and read_depth)
#--mixed_models_cor_age: Abundances vs phenotypes (correcting for DNA_conc, read_depth and infant age)
#--mixed_models_cor_age_feeding: Abundances vs phenotypes (correcting for DNA_conc, read_depth, infant age and feeding mode)
#--mixed_models_cor_age_bacteria: Abundances vs phenotypes (correcting for DNA_conc, read_depth, infant age and bacterial host abundance)

#*Note that the functions allow to use vOTU abundances (individual) or vOTU abundances grouped by predicted host

mixed_models <- function(metadata, ID, CLR_transformed_data, pheno_list) {
  df <- metadata
  row.names(df) <- df[,ID]
  df<-merge(df, CLR_transformed_data, by='row.names')
  row.names(df) <- df$Row.names
  df$Row.names <- NULL
  
  Prevalent= c(colnames(CLR_transformed_data))
  Overall_result_phenos =tibble() 
  for (bug_index in 1:length(Prevalent)) {
    Bug <- Prevalent[bug_index]
    cat("Bug", bug_index, "/", length(Prevalent), "\n")
    if (! Bug %in% colnames(df)){ next }
    Bug2 = paste(c("`",Bug, "`"), collapse="")
    for ( pheno in pheno_list){
      pheno2 = paste(c("`",pheno, "`"), collapse="")
      df[is.na(df[colnames(df) == pheno]) == F, ID] -> To_keep
      df_pheno = filter(df, !!sym(ID) %in% To_keep )
      Model0 = as.formula(paste( c(Bug2,  " ~ read_depth + DNA_concentration_ng_ul  + (1|CS_BABY_BIOME_ID)"), collapse="" )) 
      lmer(Model0, df_pheno, REML = F) -> resultmodel0
      base_model=resultmodel0
      Model2 = as.formula(paste( c(Bug2,  " ~ read_depth + DNA_concentration_ng_ul + ",pheno2, "+ (1|CS_BABY_BIOME_ID)"), collapse="" ))
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

mixed_models_cor_age <- function(metadata, ID, CLR_transformed_data, pheno_list) {
  df <- metadata
  row.names(df) <- df[,ID]
  df<-merge(df, CLR_transformed_data, by='row.names')
  row.names(df) <- df$Row.names
  df$Row.names <- NULL
  
  Prevalent= c(colnames(CLR_transformed_data))
  Overall_result_phenos =tibble() 
  for (bug_index in 1:length(Prevalent)) {
    Bug <- Prevalent[bug_index]
    cat("Bug", bug_index, "/", length(Prevalent), "\n")
    if (! Bug %in% colnames(df)){ next }
    Bug2 = paste(c("`",Bug, "`"), collapse="")
    for ( pheno in pheno_list){
      pheno2 = paste(c("`",pheno, "`"), collapse="")
      df[is.na(df[colnames(df) == pheno]) == F, ID] -> To_keep
      df_pheno = filter(df, !!sym(ID) %in% To_keep )
      Model0 = as.formula(paste( c(Bug2,  " ~ read_depth + DNA_concentration_ng_ul + Timepoint_numeric + (1|CS_BABY_BIOME_ID)"), collapse="" )) 
      lmer(Model0, df_pheno, REML = F) -> resultmodel0
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


mixed_models_cor_age_feeding <- function(metadata, ID, CLR_transformed_data, pheno_list) {
  df <- metadata
  row.names(df) <- df[,ID]
  df<-merge(df, CLR_transformed_data, by='row.names')
  row.names(df) <- df$Row.names
  df$Row.names <- NULL
  
  Prevalent= c(colnames(CLR_transformed_data))
  Overall_result_phenos =tibble() 
  
  for (bug_index in 1:length(Prevalent)) {
    Bug <- Prevalent[bug_index]
    cat("Bug", bug_index, "/", length(Prevalent), "\n")
    if (! Bug %in% colnames(df)){ next }
    Bug2 = paste(c("`",Bug, "`"), collapse="")
    for ( pheno in pheno_list){
      pheno2 = paste(c("`",pheno, "`"), collapse="")
      df[is.na(df[colnames(df) == pheno]) == F, ID] -> To_keep
      df_pheno = filter(df, !!sym(ID) %in% To_keep )
      Model0 = as.formula(paste( c(Bug2,  " ~ read_depth + DNA_concentration_ng_ul + Timepoint_numeric + feeding_mode_pragmatic + (1|CS_BABY_BIOME_ID)"), collapse="" )) 
      lmer(Model0, df_pheno, REML = F) -> resultmodel0
      base_model=resultmodel0
      Model2 = as.formula(paste( c(Bug2,  " ~ read_depth + DNA_concentration_ng_ul + Timepoint_numeric + feeding_mode_pragmatic +",pheno2, "+ (1|CS_BABY_BIOME_ID)"), collapse="" ))
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

mixed_models_cor_age_bacteria <- function(metadata, viral_metadata, ID, CLR_transformed_data, pheno_list) {
  df <- metadata
  row.names(df) <- df[,ID]
  df<-merge(df, CLR_transformed_data, by='row.names')
  row.names(df) <- df$Row.names
  df$Row.names <- NULL
  
  Prevalent= c(colnames(CLR_transformed_data))
  Overall_result_phenos =tibble()
  
  for (bug_index in 1:length(Prevalent)) {
    Bug <- Prevalent[bug_index]
    # If working with vOTU abundances -> Extract vOTu host from Viral metadata
    # If working with vOTU abundances grouped by host -> Extract host from from "Bug" name
    Host <- ifelse(Bug %in% viral_metadata$Virus_ID, 
                   viral_metadata[viral_metadata$Virus_ID == Bug, "Bacterial_genus_host"], 
                   sub("_vOTUs$", "", Bug))
    
    cat("Bug", bug_index, "/", length(Prevalent), "\n")
    if (! Bug %in% colnames(df)){ next }
    if (! Host %in% colnames(df)){ next } #exclude associations with vOTUs without host assignment
    Bug2 = paste(c("",Bug, ""), collapse="")
    for ( pheno in pheno_list){
      pheno2 = paste(c("",pheno, ""), collapse="")
      df[is.na(df[colnames(df) == pheno]) == F, ID] -> To_keep
      df_pheno = filter(df, !!sym(ID) %in% To_keep )
      Model0 = as.formula(paste( c(Bug2, " ~ read_depth + DNA_concentration_ng_ul + Timepoint_numeric +", Host, " +(1|CS_BABY_BIOME_ID)"), collapse="" ))
      lmer(Model0, df_pheno, REML = F) -> resultmodel0
      base_model=resultmodel0
      Model2 = as.formula(paste( c(Bug2, " ~ read_depth + DNA_concentration_ng_ul + Timepoint_numeric +", Host, "+", pheno2, "+ (1|CS_BABY_BIOME_ID)"), collapse="" ))
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


# Set working directory 
setwd("~/Desktop/PhD/Projects/CS Baby Biome/Final_analysis/VIRUSES/")

# Read abundance table and metadata tables
Abundance_table_infants <- read.delim("Abundance_table/CS_Abundance_Table_INFANTS_17052024.txt")
Sample_metadata_infants <- read.delim("Metadata_CS/CS_Baby_Biome_Sample_Metadata_Infants_17052024.txt")
Virus_metadata <- read.delim("Metadata_CS/CS_Baby_Biome_Viral_Metadata_17052024.txt")

##************************************************************************
# 1. Phenotype selection and processing: Infants
##************************************************************************

# Select only phenotypes of interest from metadata
colnames(Sample_metadata_infants)
Sample_metadata_infants <- Sample_metadata_infants[,c("DNA_concentration_ng_ul","read_depth", "bioSampleId",
                                                      "CS_BABY_BIOME_ID", "Timepoint_categorical", "Timepoint_numeric", 
                                                      "preg_gest_age","pre_preg_bmi_mother", "infant_birthweight","infant_sex",
                                                      "feeding_mode_pragmatic", "living_situation", "cats_dogs", "rand_AB","richness",
                                                      "shannon", "Virome.composition.infant", "Viral_relab_temperates"), ]
# Convert character variables to factors
Sample_metadata_infants[sapply(Sample_metadata_infants, is.character)] <- lapply(Sample_metadata_infants[sapply(Sample_metadata_infants,
                                                                                                                is.character)],  #convert character columns to factors
                                                   as.factor)


##################################
# Process vOTU abundance table: LOG transformation and filter by prevalence
##################################
Abundance_table_infants <- t(Abundance_table_infants)

# Calculate the prevalence of each vOTU
vOTU_nonZero  <- colSums(Abundance_table_infants>0)/nrow(Abundance_table_infants) 

# Get the list of vOTUs present in more than 5% samples
vOTU_keep <- colnames(Abundance_table_infants)[vOTU_nonZero > 0.05]


# Generate LOG transformed abundance tables
## CLR transformation not suitable here: data hardly compositional after read mapping
## Log transformation can help with data normality necessary for linear models
pseudocount_vOTU <- min(Abundance_table_infants[Abundance_table_infants > 0]) /2
Abundance_table_infants_log <- Abundance_table_infants
Abundance_table_infants_log[Abundance_table_infants_log  == 0] <- pseudocount_vOTU
Abundance_table_infants_log <- log(Abundance_table_infants_log)


##************************************************************************
# 2. Association analysis - Infants: vOTU abundances vs phenotypes
##************************************************************************

##################################
# vOTU abundances vs phenotypes 
##################################
phenotypes_to_test <- colnames(Sample_metadata_infants)[c(6:14)]
vOTU_associations_mixed <- mixed_models(Sample_metadata_infants, "bioSampleId", Abundance_table_infants_log_filtered, 
                                            phenotypes_to_test)

# Remove from results multiple levels of the same phenotypes (factors) and check the number of significant phenotypes
vOTU_associations_mixed <- vOTU_associations_mixed %>%
  distinct(Pheno, Bug, .keep_all = TRUE) %>%
  arrange(FDR)

# Check significant results
vOTU_associations_mixed_sig <- vOTU_associations_mixed [vOTU_associations_mixed$FDR < 0.05,]

# Save results
write.table(vOTU_associations_mixed, "7_ASSOCIATION_PHENOTYPES/vOTU_associations_mixed_20062024.txt",
            sep="\t", row.names=F, quote = F)


##################################
# vOTU abundances vs phenotypes correcting for age
##################################
phenotypes_to_test <-  phenotypes_to_test[!phenotypes_to_test =="Timepoint_numeric"] # exclude age
vOTU_associations_mixed_age <- mixed_models_cor_age(Sample_metadata_infants, "bioSampleId", Abundance_table_infants_log_filtered, 
                                            phenotypes_to_test)

# Remove from results multiple levels of the same phenotypes (factors) and check the number of significant phenotypes
vOTU_associations_mixed_age <- vOTU_associations_mixed_age %>%
  distinct(Pheno, Bug, .keep_all = TRUE) %>%
  arrange(FDR)

# Check significant results
vOTU_associations_mixed_age_sig <- vOTU_associations_mixed_age [vOTU_associations_mixed_age$FDR < 0.05,]

# Save results
write.table(vOTU_associations_mixed_age, "7_ASSOCIATION_PHENOTYPES/vOTU_associations_mixed_age_20062024.txt",
            sep="\t", row.names=F, quote = F)

##################################
# vOTU abundances vs phenotypes correcting for feeding mode + age
##################################
phenotypes_to_test <-  phenotypes_to_test[!phenotypes_to_test =="feeding_mode_pragmatic"] 
vOTU_associations_mixed_age_feeding <- mixed_models_cor_age_feeding(Sample_metadata_infants, "bioSampleId", Abundance_table_infants_log_filtered, 
                                                                   phenotypes_to_test)

# Remove from results multiple levels of the same phenotypes (factors) and check the number of significant phenotypes
vOTU_associations_mixed_age_feeding <- vOTU_associations_mixed_age_feeding %>%
  distinct(Pheno, Bug, .keep_all = TRUE) %>%
  arrange(FDR)

# Check significant results
vOTU_associations_mixed_age_feeding_sig <- vOTU_associations_mixed_age_feeding [vOTU_associations_mixed_age_feeding$FDR < 0.05,]

# Save results
write.table(vOTU_associations_mixed_age_feeding, "7_ASSOCIATION_PHENOTYPES/vOTU_associations_mixed_age_feeding_20062024.txt",
            sep="\t", row.names=F, quote = F)


##################################
# Add bacterial abundance information
##################################
# Load bacterial abundance data
bacterial_abundances <- read.delim("5_BACTERIAL_ANALYSIS/Metaphlan_4_genus_21062024.txt")
bacterial_abundances_infants <- bacterial_abundances[rownames(bacterial_abundances) %in% Sample_metadata_infants$bioSampleId,]

# Perform CLR transformation
my_pseudocount_normal_bact <- min(bacterial_abundances_infants[bacterial_abundances_infants!=0])/2
bacterial_abundances_infants_CLR <- decostand(bacterial_abundances_infants, "clr", pseudocount=my_pseudocount_normal_bact)

# Merge with sample metadata
bacterial_abundances_infants_CLR$bioSampleId <- rownames(bacterial_abundances_infants_CLR)
Sample_metadata_infants <-  dplyr::left_join(Sample_metadata_infants, bacterial_abundances_infants_CLR, by = "bioSampleId")


##################################
# vOTU abundances vs phenotypes correcting for age + bacterial host
##################################
phenotypes_to_test <- colnames(Sample_metadata_infants)[c(6:14)]
phenotypes_to_test <-  phenotypes_to_test[!phenotypes_to_test =="Timepoint_numeric"] # exclude age

# Run the association
vOTU_associations_mixed_age_bacteria <- mixed_models_cor_age_bacteria(Sample_metadata_infants, Virus_metadata,
                                                                      "bioSampleId", Abundance_table_infants_log_filtered, 
                                                                                      phenotypes_to_test)

# Remove from results multiple levels of the same phenotypes (factors) and check the number of significant phenotypes
vOTU_associations_mixed_age_bacteria <- vOTU_associations_mixed_age_bacteria %>%
  distinct(Pheno, Bug, .keep_all = TRUE) %>%
  arrange(FDR)

# Check significant results
vOTU_associations_mixed_age_bacteria_sig <- vOTU_associations_mixed_age_bacteria [vOTU_associations_mixed_age_bacteria$FDR < 0.05,]

# Save results
write.table(vOTU_associations_mixed_age_bacteria, "7_ASSOCIATION_PHENOTYPES/vOTU_associations_mixed_age_bacteria_20062024.txt",
            sep="\t", row.names=F, quote = F)


##************************************************************************
# 3. Association analysis - Infants: vOTU abundances (by predicted host) vs phenotypes
##************************************************************************

# Load vOTU abundances grouped by bacterial host at genus level
vOTU_bacterial_genus_abundance <- read.delim("6_HOST_ASSIGNMENT/Bacterial_genus_host_vOTU_abundance_infants_17052024.txt")

# Calculate the prevalence of each predicted host
genus_nonZero  <- colSums(vOTU_bacterial_genus_abundance>0)/nrow(vOTU_bacterial_genus_abundance) 

# Get the list of genus hosts present in more than 5% samples
genus_keep <- colnames(vOTU_bacterial_genus_abundance)[genus_nonZero > 0.05]

# Generate LOG transformed abundance tables
## CLR transformation not suitable here: data hardly compositional after read mapping
## Log transformation can help with data normality necessary for linear models
pseudocount_vOTU_host <- min(vOTU_bacterial_genus_abundance[vOTU_bacterial_genus_abundance > 0]) /2
vOTU_bacterial_genus_abundance_log <- vOTU_bacterial_genus_abundance
vOTU_bacterial_genus_abundance_log[vOTU_bacterial_genus_abundance_log  == 0] <- pseudocount_vOTU_host
vOTU_bacterial_genus_abundance_log <- log(vOTU_bacterial_genus_abundance_log)

# Remove the vOTUs that are not prevalent
vOTU_bacterial_genus_abundance_log_filtered <- as.data.frame(vOTU_bacterial_genus_abundance_log[,genus_keep])

#Order abundance table according to Sample metadata
vOTU_bacterial_genus_abundance_log_filtered <- vOTU_bacterial_genus_abundance_log_filtered[match(Sample_metadata_infants$bioSampleId,
                                                                                                 rownames(vOTU_bacterial_genus_abundance_log_filtered)), ]

# Add  suffix to colnames (to know that these are vOTU abundances grouped by host)
colnames(vOTU_bacterial_genus_abundance_log_filtered) <- paste0(colnames(vOTU_bacterial_genus_abundance_log_filtered), "_vOTUs")

##################################
# vOTU abundances (by predicted host) vs phenotypes 
##################################
phenotypes_to_test <- colnames(Sample_metadata_infants)[c(6:14)]
vOTU_genus_associations_mixed  <- mixed_models(Sample_metadata_infants, "bioSampleId", vOTU_bacterial_genus_abundance_log_filtered, 
                                        phenotypes_to_test)

# Remove from results multiple levels of the same phenotypes (factors) and check the number of significant phenotypes
vOTU_genus_associations_mixed <- vOTU_genus_associations_mixed %>%
  distinct(Pheno, Bug, .keep_all = TRUE) %>%
  arrange(FDR)

# Check significant results
vOTU_genus_associations_mixed_sig <- vOTU_genus_associations_mixed [vOTU_genus_associations_mixed$FDR < 0.05,]

# Save results
write.table(vOTU_genus_associations_mixed, "7_ASSOCIATION_PHENOTYPES/vOTU_genus_associations_mixed_20062024.txt",
            sep="\t", row.names=F, quote = F)

##################################
# vOTU abundances (by predicted host) vs phenotypes correcting for age
##################################
phenotypes_to_test <-  phenotypes_to_test[!phenotypes_to_test =="Timepoint_numeric"] # exclude age
vOTU_genus_associations_mixed_age <- mixed_models_cor_age(Sample_metadata_infants, "bioSampleId", vOTU_bacterial_genus_abundance_log_filtered, 
                                                    phenotypes_to_test)

# Remove from results multiple levels of the same phenotypes (factors) and check the number of significant phenotypes
vOTU_genus_associations_mixed_age <- vOTU_genus_associations_mixed_age %>%
  distinct(Pheno, Bug, .keep_all = TRUE) %>%
  arrange(FDR)

# Check significant results
vOTU_genus_associations_mixed_age_sig <- vOTU_genus_associations_mixed_age [vOTU_genus_associations_mixed_age$FDR < 0.05,]

# Save results
write.table(vOTU_genus_associations_mixed_age, "7_ASSOCIATION_PHENOTYPES/vOTU_genus_associations_mixed_age_20062024.txt",
            sep="\t", row.names=F, quote = F)

##################################
# vOTU abundances (by predicted host) vs phenotypes correcting for feeding mode + age
##################################

phenotypes_to_test <-  phenotypes_to_test[!phenotypes_to_test =="feeding_mode_pragmatic"] 
vOTU_genus_associations_mixed_age_feeding <- mixed_models_cor_age_feeding(Sample_metadata_infants, "bioSampleId", vOTU_bacterial_genus_abundance_log_filtered, 
                                                                    phenotypes_to_test)

# Remove from results multiple levels of the same phenotypes (factors) and check the number of significant phenotypes
vOTU_genus_associations_mixed_age_feeding <- vOTU_genus_associations_mixed_age_feeding %>%
  distinct(Pheno, Bug, .keep_all = TRUE) %>%
  arrange(FDR)

# Check significant results
vOTU_genus_associations_mixed_age_feeding_sig <- vOTU_genus_associations_mixed_age_feeding [vOTU_genus_associations_mixed_age_feeding$FDR < 0.05,]

# Save results
write.table(vOTU_genus_associations_mixed_age_feeding, "7_ASSOCIATION_PHENOTYPES/vOTU_genus_associations_mixed_age_feeding_20062024.txt",
            sep="\t", row.names=F, quote = F)

##################################
# vOTU abundances (by predicted host) vs phenotypes correcting for age + bacterial host
##################################
phenotypes_to_test <- colnames(Sample_metadata_infants)[c(6:14)]
phenotypes_to_test <-  phenotypes_to_test[!phenotypes_to_test =="Timepoint_numeric"] # exclude age

# Run the association
vOTU_genus_associations_mixed_age_bacteria <- mixed_models_cor_age_bacteria(Sample_metadata_infants, Virus_metadata, "bioSampleId",
                                                                            vOTU_bacterial_genus_abundance_log_filtered, 
                                                                                      phenotypes_to_test)

# Remove from results multiple levels of the same phenotypes (factors) and check the number of significant phenotypes
vOTU_genus_associations_mixed_age_bacteria <- vOTU_genus_associations_mixed_age_bacteria %>%
  distinct(Pheno, Bug, .keep_all = TRUE) %>%
  arrange(FDR)

# Check significant results
vOTU_genus_associations_mixed_age_bacteria_sig <- vOTU_genus_associations_mixed_age_bacteria [vOTU_genus_associations_mixed_age_bacteria$FDR < 0.05,]

# Save results
write.table(vOTU_genus_associations_mixed_age_bacteria, "7_ASSOCIATION_PHENOTYPES/vOTU_genus_associations_mixed_age_bacteria_20062024.txt",
            sep="\t", row.names=F, quote = F)


##************************************************************************
# 4. Association analysis - Infants: Diversity vs phenotypes
##************************************************************************

# Add diversity and composition variables to the abundance table to be tested against the phenotypes (remove from Sample_metadata)
Diversity_infants <- Sample_metadata_infants[, c("richness", "shannon", "Virome.composition.infant", "Viral_relab_temperates")]
colnames(Diversity_infants) <- c("Richness", "Shannon", "Virome_composition", "Viral_relab_temperates")
rownames(Diversity_infants) <- Sample_metadata_infants$bioSampleId

Sample_metadata_infants <- Sample_metadata_infants[, !(colnames(Sample_metadata_infants) 
                                                       %in% c( "richness", 
                                                              "Virome.composition.infant",
                                                              "Viral_relab_temperates"))]


##################################
# Diversity and composition vs phenotypes 
##################################
phenotypes_to_test <- colnames(Sample_metadata_infants)[c(6:14)]
diversity_comp_associations_mixed <- mixed_models(Sample_metadata_infants, "bioSampleId", Diversity_infants, 
                                        phenotypes_to_test)

# Remove from results multiple levels of the same phenotypes (factors) and check the number of significant phenotypes
diversity_comp_associations_mixed <- diversity_comp_associations_mixed %>%
  distinct(Pheno, Bug, .keep_all = TRUE) %>%
  arrange(FDR)

# Check significant results
diversity_comp_associations_mixed_sig <- diversity_comp_associations_mixed [diversity_comp_associations_mixed$FDR < 0.05,]

# Save results
write.table(diversity_comp_associations_mixed, "7_ASSOCIATION_PHENOTYPES/diversity_comp_associations_mixed_20062024.txt",
            sep="\t", row.names=F, quote = F)

##################################
# Diversity and composition vs phenotypes correcting for age
##################################
phenotypes_to_test <-  phenotypes_to_test[!phenotypes_to_test =="Timepoint_numeric"] # exclude age

diversity_comp_associations_mixed_age <- mixed_models_cor_age(Sample_metadata_infants, "bioSampleId", Diversity_infants, 
                                                  phenotypes_to_test)

# Remove from results multiple levels of the same phenotypes (factors) and check the number of significant phenotypes
diversity_comp_associations_mixed_age <- diversity_comp_associations_mixed_age %>%
  distinct(Pheno, Bug, .keep_all = TRUE) %>%
  arrange(FDR)

# Check significant results
diversity_comp_associations_mixed_age_sig <- diversity_comp_associations_mixed_age [diversity_comp_associations_mixed_age$FDR < 0.05,]

# Save results
write.table(diversity_comp_associations_mixed_age, "7_ASSOCIATION_PHENOTYPES/diversity_comp_associations_mixed_age_20062024.txt",
            sep="\t", row.names=F, quote = F)


##################################
# Diversity and composition vs phenotypes correcting for feeding mode + age
##################################
phenotypes_to_test <-  phenotypes_to_test[!phenotypes_to_test =="feeding_mode_pragmatic"] 
diversity_comp_associations_mixed_age_feeding <- mixed_models_cor_age_feeding(Sample_metadata_infants, "bioSampleId", Diversity_infants, 
                                                      phenotypes_to_test)

# Remove from results multiple levels of the same phenotypes (factors) and check the number of significant phenotypes
diversity_comp_associations_mixed_age_feeding <- diversity_comp_associations_mixed_age_feeding %>%
  distinct(Pheno, Bug, .keep_all = TRUE) %>%
  arrange(FDR)

# Check significant results
diversity_comp_associations_mixed_age_feeding_sig <- diversity_comp_associations_mixed_age_feeding [diversity_comp_associations_mixed_age_feeding$FDR < 0.05,]

# Save results
write.table(diversity_comp_associations_mixed_age_feeding, "7_ASSOCIATION_PHENOTYPES/diversity_comp_associations_mixed_age_feeding_20062024.txt",
            sep="\t", row.names=F, quote = F)


##************************************************************************
# 5. Plot significant results
##************************************************************************

# Reload Sample_metadata
Sample_metadata_infants <- read.delim("Metadata_CS/CS_Baby_Biome_Sample_Metadata_Infants_17052024.txt")

Sample_metadata_infants$feeding_mode_pragmatic <- factor(Sample_metadata_infants$feeding_mode_pragmatic,
                                                         levels = c("breast_feeding", "mixed_feeding", "formula_feeding"), 
                                                         labels = c("BM", "Mixed", "Formula"))

Sample_metadata_infants$living_situation <- factor(Sample_metadata_infants$living_situation,
                                                         levels = c("city", "village", "farm"), 
                                                         labels = c("City", "Village", "Farm"))

##################################
# Associations with diversity/richness
##################################

pdf('7_ASSOCIATION_PHENOTYPES/PLOTS/Richness_feeding_mode.pdf', width=3,2, height=3.2)
ggplot(Sample_metadata_infants, aes(x=feeding_mode_pragmatic, y=richness)) +
  geom_violin(aes(fill=feeding_mode_pragmatic), alpha=0.7, width=0.8, trim = F) +
  geom_jitter(aes(color=feeding_mode_pragmatic), alpha=0.7) + 
  geom_boxplot(width=0.2, fill="white", color="black", outliers = FALSE) + 
  labs(x = 'Feeding mode', y = 'Viral Richness') + 
  scale_fill_manual(values = c("#FFBE81", "#DFF0D8","#C7C7FF")) + 
  scale_color_manual(values = c("#FFBE81", "#DFF0D8","#C7C7FF")) + 
  theme_classic() +
  theme(plot.title = element_text(size=18, hjust = 0.5, face="bold"), 
        axis.text = element_text(size=14),
        axis.title = element_text(size=15), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.position = "none") 
dev.off()


pdf('7_ASSOCIATION_PHENOTYPES/PLOTS/Richness_cats_dogs.pdf', width=2.8, height=3.2)
ggplot(Sample_metadata_infants[!is.na(Sample_metadata_infants$cats_dogs),], aes(x=cats_dogs, y=richness)) +
  geom_violin(aes(fill=cats_dogs), alpha=0.7, width=0.8, trim = FALSE) + 
  geom_jitter(aes(color=cats_dogs), alpha=0.7) + # Jittered points
  geom_boxplot(width=0.2, fill="white", color="black", outliers = FALSE) + 
  labs(x = 'Feeding mode', y = 'Viral Richness') + # Labels
  scale_fill_manual(values=c("#8B4513", "#228B22")) + 
  scale_color_manual(values=c("#8B4513", "#228B22")) + 
  theme_classic() + # Classic theme
  theme(plot.title = element_text(size=18, hjust = 0.5, face="bold"), 
        axis.text = element_text(size=14),
        axis.title = element_text(size=15), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.position = "none") 
dev.off()

##################################
# Associations with vOTU abundances
##################################

Sample_metadata_infants <- cbind(Sample_metadata_infants, Abundance_table_infants_log_filtered[,c("MGV_120587", "CS_297",
                                                                                              "ELGV_7255", "ELGV_11626")])


pdf('7_ASSOCIATION_PHENOTYPES/PLOTS/CS_297_living_situation.pdf', width=2.8, height=3.2)
ggplot(Sample_metadata_infants[!is.na(Sample_metadata_infants$living_situation),], aes(x=living_situation, y=CS_297)) +
  geom_violin(aes(fill=living_situation), alpha=0.7, width=0.8, trim = F) +
  geom_jitter(aes(color=living_situation), alpha=0.7) + 
  geom_boxplot(width=0.2, fill="white", color="black", outliers = FALSE) + 
  labs(x = 'Living situation', y = 'CS_297 Abundance') + 
  #scale_fill_manual(values=c("#008080", "#4B0082")) + 
  #scale_color_manual(values=c("#008080", "#4B0082")) +
  theme_classic() +
  theme(plot.title = element_text(size=18, hjust = 0.5, face="bold"), 
        axis.text = element_text(size=14),
        axis.title = element_text(size=15), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.position = "none") 
dev.off()


##################################
# Associations with vOTU abundances grouped by bacterial host
##################################

Sample_metadata_infants <- cbind(Sample_metadata_infants, vOTU_bacterial_genus_abundance_log_filtered[,c("Clostridium_vOTUs", "Lacticaseibacillus_vOTUs",
                                                                                                         "Enterococcus_D_vOTUs", "Bifidobacterium_vOTUs",
                                                                                                         "Limosilactobacillus_vOTUs")])
Sample_metadata_infants <- cbind(Sample_metadata_infants, bacterial_abundances_infants_CLR[,c("Clostridium", "Lacticaseibacillus",
                                                                                              "Enterococcus", "Bifidobacterium")])

pdf('7_ASSOCIATION_PHENOTYPES/PLOTS/Clostridium_vOTUs_by_host_feeding.pdf', width=3.2, height=3.2)
ggplot(Sample_metadata_infants[!is.na(Sample_metadata_infants$feeding_mode_pragmatic), ], aes(x = feeding_mode_pragmatic, y = Clostridium_vOTUs)) +
  geom_violin(aes(fill=feeding_mode_pragmatic), alpha=0.7, width=0.8, trim = F) +
  geom_jitter(aes(color=feeding_mode_pragmatic), alpha=0.7) + 
  geom_boxplot(width=0.2, fill="white", color="black", outliers = FALSE) + 
  labs(x = 'Feeding mode', y = 'Clostridium vOTUs Abundance') + 
  scale_fill_manual(values = c("#FFBE81", "#DFF0D8","#C7C7FF")) + 
  scale_color_manual(values = c("#FFBE81", "#DFF0D8","#C7C7FF")) + 
  theme_classic() +
  theme(plot.title = element_text(size = 18, hjust = 0.5, face = "bold"), 
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 15), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.position = "none")
dev.off()

pdf('7_ASSOCIATION_PHENOTYPES/PLOTS/Clostridium_bacteria_feeding.pdf', width=3.2, height=3.2)
ggplot(Sample_metadata_infants[!is.na(Sample_metadata_infants$feeding_mode_pragmatic), ], aes(x = feeding_mode_pragmatic, y = Clostridium)) +
  geom_violin(aes(fill=feeding_mode_pragmatic), alpha=0.7, width=0.8, trim = F) +
  geom_jitter(aes(color=feeding_mode_pragmatic), alpha=0.7) + 
  geom_boxplot(width=0.2, fill="white", color="black", outliers = FALSE) + 
  labs(x = 'Feeding mode', y = 'Clostridium Abundance') + 
  scale_fill_manual(values = c("#FFBE81", "#DFF0D8","#C7C7FF")) + 
  scale_color_manual(values = c("#FFBE81", "#DFF0D8","#C7C7FF")) + 
  theme_classic() +
  theme(plot.title = element_text(size = 18, hjust = 0.5, face = "bold"), 
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 15), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.position = "none")
dev.off()

pdf('7_ASSOCIATION_PHENOTYPES/PLOTS/Lacticaseibacillus_vOTUs_by_host_feeding.pdf', width=3.2, height=3.2)
ggplot(Sample_metadata_infants[!is.na(Sample_metadata_infants$feeding_mode_pragmatic), ], aes(x = feeding_mode_pragmatic, y = Lacticaseibacillus_vOTUs)) +
  geom_violin(aes(fill=feeding_mode_pragmatic), alpha=0.7, width=0.8, trim = F) +
  geom_jitter(aes(color=feeding_mode_pragmatic), alpha=0.7) + 
  geom_boxplot(width=0.2, fill="white", color="black", outliers = FALSE) + 
  labs(x = 'Feeding mode', y = 'Lacticaseibacillus vOTUs Abundance') + 
  scale_fill_manual(values = c("#FFBE81", "#DFF0D8","#C7C7FF")) + 
  scale_color_manual(values = c("#FFBE81", "#DFF0D8","#C7C7FF")) + 
  theme_classic() +
  theme(plot.title = element_text(size = 18, hjust = 0.5, face = "bold"), 
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 15), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.position = "none")
dev.off()

pdf('7_ASSOCIATION_PHENOTYPES/PLOTS/Lacticaseibacillus_bacteria_feeding.pdf', width=3.2, height=3.2)
ggplot(Sample_metadata_infants[!is.na(Sample_metadata_infants$feeding_mode_pragmatic), ], aes(x = feeding_mode_pragmatic, y = Lacticaseibacillus)) +
  geom_violin(aes(fill=feeding_mode_pragmatic), alpha=0.7, width=0.8, trim = F) +
  geom_jitter(aes(color=feeding_mode_pragmatic), alpha=0.7) + 
  geom_boxplot(width=0.2, fill="white", color="black", outliers = FALSE) + 
  labs(x = 'Feeding mode', y = 'Lacticaseibacillus Abundance') + 
  scale_fill_manual(values = c("#FFBE81", "#DFF0D8","#C7C7FF")) + 
  scale_color_manual(values = c("#FFBE81", "#DFF0D8","#C7C7FF")) + 
  theme_classic() +
  theme(plot.title = element_text(size = 18, hjust = 0.5, face = "bold"), 
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 15), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.position = "none")
dev.off()


pdf('7_ASSOCIATION_PHENOTYPES/PLOTS/Limosilactobacillus_vOTUs_by_pets.pdf', width=3.2, height=3.2)
ggplot(Sample_metadata_infants[!is.na(Sample_metadata_infants$cats_dogs), ], aes(x = cats_dogs, y = Limosilactobacillus_vOTUs)) +
  geom_violin(aes(fill=cats_dogs), alpha=0.7, width=0.8, trim = F) +
  geom_jitter(aes(color=cats_dogs), alpha=0.7) + 
  geom_boxplot(width=0.2, fill="white", color="black", outliers = FALSE) + 
  labs(x = 'Pets', y = 'Limosilactobacillus vOTUs Abundance') + 
  scale_fill_manual(values = c("#8B4513", "#228B22")) + 
  scale_color_manual(values = c("#8B4513", "#228B22")) + 
  theme_classic() +
  theme(plot.title = element_text(size = 18, hjust = 0.5, face = "bold"), 
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 15), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.position = "none")
dev.off()

pdf('7_ASSOCIATION_PHENOTYPES/PLOTS/Bifidobacterium_vOTUs_pre_preg_BMI.pdf', width=2.8, height=3.2)
ggplot(Sample_metadata_infants[!is.na(Sample_metadata_infants$pre_preg_bmi_mother), ], aes(x = Bifidobacterium_vOTUs, y = pre_preg_bmi_mother)) +
  geom_point(alpha = 0.7) + 
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue", alpha = 0.3) +
  labs(x ='Bifidobacterium vOTU Abundance' , y = 'Pre-pregnancy BMI') +
  theme_classic() +
  theme(plot.title = element_text(size = 18, hjust = 0.5, face = "bold"), 
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 15), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.position = "none")
dev.off()

pdf('7_ASSOCIATION_PHENOTYPES/PLOTS/Limosilactobacillus_vOTUs_pre_preg_BMI.pdf', width=2.8, height=3.2)
ggplot(Sample_metadata_infants[!is.na(Sample_metadata_infants$pre_preg_bmi_mother), ], aes(x = Limosilactobacillus_vOTUs, y = pre_preg_bmi_mother)) +
  geom_point(alpha = 0.7) + 
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue", alpha = 0.3) +
  labs(x ='Limosilactobacillus vOTU Abundance' , y = 'Pre-pregnancy BMI') +
  theme_classic() +
  theme(plot.title = element_text(size = 18, hjust = 0.5, face = "bold"), 
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 15), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.position = "none")
dev.off()

pdf('7_ASSOCIATION_PHENOTYPES/PLOTS/Enterococcus_D_vOTUs_pre_preg_BMI.pdf', width=2.8, height=3.2)
ggplot(Sample_metadata_infants[!is.na(Sample_metadata_infants$pre_preg_bmi_mother), ], aes(x = Enterococcus_D_vOTUs, y = pre_preg_bmi_mother)) +
  geom_point(alpha = 0.7) + 
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue", alpha = 0.3) +
  labs(x ='Enterococcus D vOTU Abundance' , y = 'Pre-pregnancy BMI') +
  theme_classic() +
  theme(plot.title = element_text(size = 18, hjust = 0.5, face = "bold"), 
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 15), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.position = "none")
dev.off()


