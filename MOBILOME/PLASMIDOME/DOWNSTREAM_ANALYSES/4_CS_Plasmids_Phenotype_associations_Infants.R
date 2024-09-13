################################################################################
##### CS Baby Biome: Associations with phenotypes in CS Baby Biome infants (plasmids)
### Author(s):Asier Fernández-Pato
### Last updated: 26th June, 2023
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

# Functions to perform the association between Phenotypes and abundances of plasmids (mixed-model)
# Different functions are defined:
#--mixed_models: Abundances vs Phenotypes (correcting for DNA_conc and read_depth)
#--mixed_models_cor_feeding: Abundances vs Phenotypes (correcting for DNA_conc, read_depth, feeding mode)

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

mixed_models_cor_feeding <- function(metadata, ID, CLR_transformed_data, pheno_list) {
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
      Model0 = as.formula(paste( c(Bug2,  " ~ read_depth + DNA_concentration_ng_ul + feeding_mode_pragmatic + (1|CS_BABY_BIOME_ID)"), collapse="" )) 
      lmer(Model0, df_pheno, REML = F) -> resultmodel0
      base_model=resultmodel0
      Model2 = as.formula(paste( c(Bug2,  " ~ read_depth + DNA_concentration_ng_ul + feeding_mode_pragmatic + ",pheno2, "+ (1|CS_BABY_BIOME_ID)"), collapse="" ))
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
setwd("~/Desktop/PhD/Projects/CS Baby Biome/Final_analysis/PLASMIDS/")

# Read abundance table and metadata tables
Plasmid_abundance_infants <- read.delim("ABUNDANCE_TABLES/CS_Abundance_Table_INFANTS_25062024.txt")
Sample_metadata_infants <- read.delim("METADATA_TABLES/CS_Baby_Biome_Sample_Metadata_Infants_25062024.txt")
Plasmid_metadata <- read.delim("METADATA_TABLES/CS_Baby_Biome_Plasmid_Metadata_Infants_25062024.txt")

#Select only samples from W01-W06
Plasmid_abundance_infants <-Plasmid_abundance_infants[, colnames(Plasmid_abundance_infants) 
                                                      %in% Sample_metadata_infants$NG_ID]

##************************************************************************
# 1. Phenotype selection and processing: Infants
##************************************************************************
  
# Select only phenotypes of interest from metadata
colnames(Sample_metadata_infants)
Sample_metadata_infants <- Sample_metadata_infants[,c("DNA_concentration_ng_ul","read_depth", "bioSampleId",
                                                      "CS_BABY_BIOME_ID", "Timepoint_categorical", "Timepoint_numeric", 
                                                      "preg_gest_age","pre_preg_bmi_mother", "infant_birthweight","infant_sex",
                                                      "feeding_mode_pragmatic", "living_situation", "cats_dogs", "rand_AB","plasmid_richness",
                                                      "plasmid_shannon", "Plasmid.composition.infant", "Plasmid_ab"), ]
# Convert character variables to factors
Sample_metadata_infants[sapply(Sample_metadata_infants, is.character)] <- lapply(Sample_metadata_infants[sapply(Sample_metadata_infants,
                                                                                                                is.character)],  #convert character columns to factors
                                                   as.factor)


##################################
# Process plasmid abundance table: LOG transformation and filter by prevalence
##################################
Plasmid_abundance_infants <- t(Plasmid_abundance_infants)

# Calculate the prevalence of each PTU
plasmid_nonZero  <- colSums(Plasmid_abundance_infants>0)/nrow(Plasmid_abundance_infants) 

# Get the list of PTUs present in more than 10% samples
plasmid_keep <- colnames(Plasmid_abundance_infants)[plasmid_nonZero > 0.1] #118

# Generate LOG transformed abundance tables
## CLR transformation not suitable here: data hardly compositional after read mapping
## Log transformation can help with data normality necessary for linear models
pseudocount_PTU <- min(Plasmid_abundance_infants[Plasmid_abundance_infants!=0])/2
Plasmid_abundance_infants_log <- Plasmid_abundance_infants
Plasmid_abundance_infants_log[Plasmid_abundance_infants_log  == 0] <- pseudocount_PTU
Plasmid_abundance_infants_log <- log(Plasmid_abundance_infants_log)

# Remove the PTUs that are not prevalent
Plasmid_abundance_infants_log_filtered <- as.data.frame(Plasmid_abundance_infants_log[,plasmid_keep])

#Order abundance table according to Sample metadata
Plasmid_abundance_infants_log_filtered <- Plasmid_abundance_infants_log_filtered[match(Sample_metadata_infants$bioSampleId,
                                                                                       rownames(Plasmid_abundance_infants_log_filtered)), ]


##************************************************************************
# 2. Association analysis - Infants: PTU abundances vs phenotypes
##************************************************************************

##################################
# Plasmid abundances vs phenotypes 
##################################
phenotypes_to_test <- colnames(Sample_metadata_infants)[c(6:14)]
Plasmid_associations_mixed <- mixed_models(Sample_metadata_infants, "bioSampleId", Plasmid_abundance_infants_log_filtered, 
                                            phenotypes_to_test)

# Remove from results multiple levels of the same phenotypes (factors) and check the number of significant phenotypes
Plasmid_associations_mixed <- Plasmid_associations_mixed %>%
  distinct(Pheno, Bug, .keep_all = TRUE) %>%
  arrange(FDR)

# Check significant results
Plasmid_associations_mixed_sig <- Plasmid_associations_mixed [Plasmid_associations_mixed$FDR < 0.05,]

# Save results
write.table(Plasmid_associations_mixed, "4_ASSOCIATION_PHENOTYPES/Plasmid_associations_mixed_25062024.txt",
            sep="\t", row.names=F, quote = F)


##************************************************************************
# 2. Association analysis - Infants: Diversity and Composition vs phenotypes
##************************************************************************

# Add diversity and composition variables to the abundance table to be tested against the phenotypes (remove from Sample_metadata)
Diversity_infants <- Sample_metadata_infants[, c("plasmid_richness","plasmid_shannon", "Plasmid.composition.infant", "Plasmid_ab")]
colnames(Diversity_infants) <- c("Richness", "Shannon", "Plasmid_composition", "Plasmid_abundance")
rownames(Diversity_infants) <- Sample_metadata_infants$bioSampleId

Sample_metadata_infants <- Sample_metadata_infants[, !(colnames(Sample_metadata_infants) 
                                                       %in% c("plasmid_richness",
                                                              "plasmid_shannon",
                                                              "Plasmid.composition.infant",
                                                              "Plasmid_ab"))]


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
write.table(diversity_comp_associations_mixed, "4_ASSOCIATION_PHENOTYPES/diversity_comp_associations_mixed_25062024.txt",
            sep="\t", row.names=F, quote = F)


##################################
# Diversity and composition vs phenotypes correcting for feeding
##################################
phenotypes_to_test <- colnames(Sample_metadata_infants)[c(6:14)]
phenotypes_to_test <-  phenotypes_to_test[!phenotypes_to_test =="feeding_mode_pragmatic"] 


diversity_comp_associations_mixed_feeding <- mixed_models_cor_feeding(Sample_metadata_infants, "bioSampleId", Diversity_infants, 
                                                              phenotypes_to_test)

# Remove from results multiple levels of the same phenotypes (factors) and check the number of significant phenotypes
diversity_comp_associations_mixed_feeding <- diversity_comp_associations_mixed_feeding %>%
  distinct(Pheno, Bug, .keep_all = TRUE) %>%
  arrange(FDR)

# Check significant results
diversity_comp_associations_mixed_feeding_sig <- diversity_comp_associations_mixed_feeding [diversity_comp_associations_mixed_feeding$FDR < 0.05,]

# Save results
write.table(diversity_comp_associations_mixed_feeding, "4_ASSOCIATION_PHENOTYPES/diversity_comp_associations_mixed_feeding_25062024.txt",
            sep="\t", row.names=F, quote = F)


##************************************************************************
# 3. Plot significant results
##************************************************************************

# Reload Sample_metadata
Sample_metadata_infants <- read.delim("METADATA_TABLES/CS_Baby_Biome_Sample_Metadata_Infants_25062024.txt")

Sample_metadata_infants$feeding_mode_pragmatic <- factor(Sample_metadata_infants$feeding_mode_pragmatic,
                                                         levels = c("breast_feeding", "mixed_feeding", "formula_feeding"), 
                                                         labels = c("Breastmilk", "Mixed", "Formula"))

Sample_metadata_infants$living_situation <- factor(Sample_metadata_infants$living_situation,
                                                         levels = c("city", "village", "farm"), 
                                                         labels = c("City", "Village", "Farm"))

##################################
# Associations with diversity/richness
##################################

pdf('4_ASSOCIATION_PHENOTYPES/PLOTS/Plasmid_Richness_Pets_Infants.pdf', width=2.9, height=3.2) 
ggplot(Sample_metadata_infants[!is.na(Sample_metadata_infants$cats_dogs), ], aes(x = cats_dogs, y = plasmid_richness)) +
  geom_violin(aes(fill=cats_dogs), alpha=0.7, width=0.8, trim = F) +
  geom_jitter(aes(color=cats_dogs), alpha=0.7) + 
  geom_boxplot(width=0.2, fill="white", color="black", outliers = FALSE) + 
  labs(x = 'Pets', y = "Plasmid richness") + 
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

pdf('4_ASSOCIATION_PHENOTYPES//PLOTS/Plasmid_Richness_Feeding_Infants.pdf', width=3.4, height=3.1)
ggplot(Sample_metadata_infants[!is.na(Sample_metadata_infants$feeding_mode_pragmatic), ], aes(x = feeding_mode_pragmatic, y = plasmid_richness)) +
  geom_violin(aes(fill=feeding_mode_pragmatic), alpha=0.7, width=0.8, trim = F) +
  geom_jitter(aes(color=feeding_mode_pragmatic), alpha=0.7) + 
  geom_boxplot(width=0.2, fill="white", color="black", outliers = FALSE) + 
  labs(x = 'Feeding mode', y = 'Plasmid richness') + 
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

