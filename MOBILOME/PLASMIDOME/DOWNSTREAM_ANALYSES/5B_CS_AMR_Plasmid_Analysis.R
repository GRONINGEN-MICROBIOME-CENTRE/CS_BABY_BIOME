################################################################################
##### CS Baby Biome: Plasmid AMR analysis - associations
### Author(s): Asier Fernández-Pato
### Last updated: 27th June, 2024
################################################################################

#****************
# Load libraries
#****************
library(tidyverse)
library(lme4)
library(lmerTest)
library(ggplot2)
library(ComplexHeatmap)


#****************
# Define functions
#****************

# Functions to perform the association between Phenotypes and abundances of plasmids (mixed-model)
# Different functions are defined:
#--mixed_models: Abundances vs Phenotypes (correcting for DNA_conc and read_depth)

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


# Set working directory 
setwd("~/Desktop/PhD/Projects/CS Baby Biome/Final_analysis/PLASMIDS/")

# Read abundance table and metadata tables
Plasmid_abundance_infants <- read.delim("ABUNDANCE_TABLES/CS_Abundance_Table_INFANTS_25062024.txt")
Sample_metadata_infants <- read.delim("METADATA_TABLES/CS_Baby_Biome_Sample_Metadata_Infants_25062024.txt")
Plasmid_metadata_infants <- read.delim("METADATA_TABLES/CS_Baby_Biome_Plasmid_Metadata_Infants_25062024.txt")

#Select only samples from W01-W06
Plasmid_abundance_infants <-Plasmid_abundance_infants[rownames(Plasmid_abundance_infants) %in% 
                                                        Plasmid_metadata_infants$Plasmid_ID,
                                                      colnames(Plasmid_abundance_infants) 
                                                      %in% Sample_metadata_infants$NG_ID]


Plasmid_abundance_infants <- Plasmid_abundance_infants[rowSums(Plasmid_abundance_infants) != 0, ]
Plasmid_metadata_infants <- Plasmid_metadata_infants[Plasmid_metadata_infants$Plasmid_ID %in% rownames(Plasmid_abundance_infants), ]

  
##************************************************************************
# 0. Phenotype selection and AMR data processing: Infants
##************************************************************************

# Read non-transformed AMR abundance tables 
AM_abundance <- read.delim("ABUNDANCE_TABLES/CS_Baby_Biome_AM_Abundance_Table.txt", check.names = F)
AMR_Gene_abundance <- read.delim("ABUNDANCE_TABLES/CS_Baby_Biome_AMR_Gene_Abundance_Table.txt", check.names = F)
AMR_Class_abundance <- read.delim("ABUNDANCE_TABLES/CS_Baby_Biome_AMR_Class_Abundance_Table.txt", check.names = F)
AMR_Mechanism_abundance <- read.delim("ABUNDANCE_TABLES/CS_Baby_Biome_AMR_Mechanism_Abundance_Table.txt", check.names = F)

# Read AMR count data
AM_presence <- read.delim("ABUNDANCE_TABLES/CS_Baby_Biome_AM_Count_Table.txt", check.names = F)
AMR_Gene_presence <- read.delim("ABUNDANCE_TABLES/CS_Baby_Biome_AMR_Gene_Count_Table.txt", check.names = F)
AMR_Class_presence <- read.delim("ABUNDANCE_TABLES/CS_Baby_Biome_AMR_Class_Count_Table.txt", check.names = F)
AMR_Mechanism_presence <- read.delim("ABUNDANCE_TABLES/CS_Baby_Biome_AMR_Mechanism_Count_Table.txt", check.names = F)

# Subset only samples in metadata
AM_abundance <- AM_abundance[rownames(AM_abundance) %in% Sample_metadata_infants$bioSampleId,]
AMR_Gene_abundance <- AMR_Gene_abundance[rownames(AMR_Gene_abundance) %in% Sample_metadata_infants$bioSampleId,]
AMR_Class_abundance <- AMR_Class_abundance[rownames(AMR_Class_abundance) %in% Sample_metadata_infants$bioSampleId,]
AMR_Mechanism_abundance <- AMR_Mechanism_abundance[rownames(AMR_Mechanism_abundance) %in% Sample_metadata_infants$bioSampleId,]

AM_presence <- AM_presence[rownames(AM_presence) %in% Sample_metadata_infants$bioSampleId,]
AMR_Gene_presence <- AMR_Gene_presence[rownames(AMR_Gene_presence) %in% Sample_metadata_infants$bioSampleId,]
AMR_Class_presence <- AMR_Class_presence[rownames(AMR_Class_presence) %in% Sample_metadata_infants$bioSampleId,]
AMR_Mechanism_presence <- AMR_Mechanism_presence[rownames(AMR_Mechanism_presence) %in% Sample_metadata_infants$bioSampleId,]


# Select only phenotypes of interest from metadata
Sample_metadata_infants <- Sample_metadata_infants[,c("DNA_concentration_ng_ul","read_depth", "bioSampleId",
                                                      "CS_BABY_BIOME_ID", "Timepoint_categorical", "Timepoint_numeric", 
                                                      "preg_gest_age","pre_preg_bmi_mother", "infant_birthweight","infant_sex",
                                                      "feeding_mode_pragmatic", "living_situation", "cats_dogs", "rand_AB","plasmid_richness",
                                                      "plasmid_shannon", "Plasmid.composition.infant", "Plasmid_ab", "AMR_gene_counts"), ]
# Convert character variables to factors
Sample_metadata_infants[sapply(Sample_metadata_infants, is.character)] <- lapply(Sample_metadata_infants[sapply(Sample_metadata_infants,
                                                                                                                is.character)],  #convert character columns to factors
                                                                                 as.factor)

##################################
# Process AMR tables: Log-transform abundances and filter by prevalence
##################################

# Generate log transformed abundance tables
## CLR transformation not suitable here: data hardly compositional after read mapping
## Log transformation can help with data normality necessary for linear models
pseudocount_AM_abundance <- min(AM_abundance[AM_abundance > 0]) /2
pseudocount_AMR_Gene_abundance <- min(AMR_Gene_abundance[AMR_Gene_abundance > 0]) /2
pseudocount_AMR_Class_abundance <- min(AMR_Class_abundance[AMR_Class_abundance > 0]) /2
pseudocount_AMR_Mechanism_abundance <- min(AMR_Mechanism_abundance[AMR_Mechanism_abundance > 0]) /2

AM_abundance_log <- AM_abundance
AM_abundance_log[AM_abundance_log  == 0] <- pseudocount_AM_abundance
AMR_Gene_abundance_log <- AMR_Gene_abundance
AMR_Gene_abundance_log[AMR_Gene_abundance_log  == 0] <- pseudocount_AMR_Gene_abundance
AMR_Class_abundance_log <- AMR_Class_abundance
AMR_Class_abundance_log[AMR_Class_abundance_log  == 0] <- pseudocount_AMR_Class_abundance
AMR_Mechanism_abundance_log <- AMR_Mechanism_abundance
AMR_Mechanism_abundance_log[AMR_Mechanism_abundance_log  == 0] <- pseudocount_AMR_Mechanism_abundance

AM_abundance_log <- log(AM_abundance_log)
AMR_Gene_abundance_log <- log(AMR_Gene_abundance_log)
AMR_Class_abundance_log <- log(AMR_Class_abundance_log)
AMR_Mechanism_abundance_log <- log(AMR_Mechanism_abundance_log)

# Select only features with a prevalence > 10%
AM_abundance_prev_list <- colnames(AM_abundance)[colSums(AM_abundance != 0) > 159*0.1] 
AMR_Gene_abundance_prev_list <- colnames(AMR_Gene_abundance)[colSums(AMR_Gene_abundance != 0) > 159*0.1] 
AMR_Gene_presence_prev_list <- colnames(AMR_Gene_presence)[colSums(AMR_Gene_presence != 0) > 159*0.1] # to include "AMR_gene_counts" 
AMR_Class_abundance_prev_list <- colnames(AMR_Class_abundance)[colSums(AMR_Class_abundance != 0) > 159*0.1] 
AMR_Mechanism_abundance_prev_list <- colnames(AMR_Mechanism_abundance)[colSums(AMR_Mechanism_abundance != 0) > 159*0.1] 

AM_abundance_prev <- AM_abundance[,AM_abundance_prev_list]
AMR_Gene_abundance_prev <- AMR_Gene_abundance[,AMR_Gene_abundance_prev_list]
AMR_Class_abundance_prev  <- AMR_Class_abundance[,AMR_Class_abundance_prev_list]
AMR_Mechanism_abundance_prev <- AMR_Mechanism_abundance[,AMR_Mechanism_abundance_prev_list]

AM_abundance_log_prev <- AM_abundance_log[,AM_abundance_prev_list]
AMR_Gene_abundance_log_prev <- AMR_Gene_abundance_log[,AMR_Gene_abundance_prev_list]
AMR_Class_abundance_log_prev  <- AMR_Class_abundance_log[,AMR_Class_abundance_prev_list]
AMR_Mechanism_abundance_log_prev <- AMR_Mechanism_abundance_log[,AMR_Mechanism_abundance_prev_list]

AM_presence_prev <- AM_presence[,AM_abundance_prev_list]
AMR_Gene_presence_prev <- AMR_Gene_presence[,AMR_Gene_presence_prev_list]
AMR_Class_presence_prev  <- AMR_Class_presence[,AMR_Class_abundance_prev_list]
AMR_Mechanism_presence_prev <- AMR_Mechanism_presence[,AMR_Mechanism_abundance_prev_list]

# Merge AMR data 
# In this case, we will only use the abundance of antimicrobials for the association analyses
# We add add also total AMR gene abundance and AMR gene count
AMR_abundance_data <- data.frame(cbind(AM_abundance_log_prev, "AMR_gene_abundance"= AMR_Gene_abundance_log_prev$AMR_gene_abundance),
                                 check.names = F)
AMR_count_data <- data.frame(cbind(AM_presence_prev,"AMR_Gene_counts" = AMR_Gene_presence_prev$AMR_gene_counts),
                             check.names = F)

##************************************************************************
# 1. General statistics
##************************************************************************

# Estimate number and proportion of plasmids with AMRs in infants (W1-W6)
length(which(Plasmid_metadata_infants$AMR_genes>0))
100*(length(which(Plasmid_metadata_infants$AMR_genes>0)) / nrow(Plasmid_metadata_infants))

# Estimate total number of AMRs in infant plasmids (W1-W6)
sum(Plasmid_metadata_infants$AMR_genes)

# Estimate number of different AMR genes found in infant plasmids
length(colSums(AMR_Gene_abundance)[colSums(AMR_Gene_abundance) != 0])

# Estimate the proportion of plasmid genes that are AMRs
100*(sum(Plasmid_metadata_infants$AMR_genes) / sum(Plasmid_metadata_infants$n_genes))

# Association between AMR gene presence and plasmid mobility
table <- table(Plasmid_metadata_infants$AMR_genes_binary, Plasmid_metadata_infants$Mobility)
chi_square_test <- chisq.test(table)
proportions <- prop.table(table, margin = 2) 
print(chi_square_test)
print(proportions)

proportions_df <- as.data.frame(as.table(proportions))
colnames(proportions_df) <- c("AMR_genes_binary", "Mobility", "Proportion")

# Estimate which AM/Am classes has resistance genes more frequently among infants
names(sort(colSums(AM_presence), decreasing = TRUE))
names(sort(colSums(AMR_Class_presence), decreasing = TRUE))

# Generate barplot: AMR gene  presence vs mobility
pdf('5_AMR_ANALYSIS/PLOTS/AMR_Presence_Mobility.pdf', width=3.5, height=3.8) 
ggplot(proportions_df, aes(x = Mobility, y = Proportion, fill = AMR_genes_binary)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
  geom_text(aes(label = scales::percent(Proportion)), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5, 
            size = 4, 
            color = "black") +
  labs(x = "Plasmid Mobility", y = "Proportion", fill = "ARG Presence") +
  scale_fill_manual(values = c("Yes" = alpha("#C40003", 0.7), "No" = alpha("#0055AA", 0.7))) +
  theme_classic() +
  theme(plot.title = element_text(size = 18, hjust = 0.5, face = "bold"), 
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 15), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13),
        legend.position = "bottom") +
  scale_y_continuous(labels = scales::percent)
dev.off()

# Estimate the ARG density per PTU (ARG per kilobase of plasmid genome)
Plasmid_metadata_infants$AMR_density <- Plasmid_metadata_infants$AMR_genes / Plasmid_metadata_infants$length
Plasmid_metadata_infants$Mobility <- factor(Plasmid_metadata_infants$Mobility)

# Association between AMR density and plasmid mobility
kruskal_test <- kruskal.test(AMR_density ~ Mobility, data = Plasmid_metadata_infants)
kruskal_test
dunn_result <- dunn.test(Plasmid_metadata_infants$AMR_density, Plasmid_metadata_infants$Mobility)
dunn_result$P.adjusted <- p.adjust(dunn_result$P, method = "fdr")

# Generate violin plot: vOTU vs PTU ARM gene content
Virus_metadata_infants <- read.delim("../VIRUSES/Metadata_CS/CS_Baby_Biome_Viral_Metadata_17052024.txt")
Virus_metadata_infants$AMR_gene_proportion <- 100*(Virus_metadata_infants$AMR_genes / Virus_metadata_infants$n_genes)
Virus_abundance_infants <- read.delim("../VIRUSES/Abundance_table/CS_Abundance_Table_INFANTS_17052024.txt")
Virus_metadata_infants <- Virus_metadata_infants[Virus_metadata_infants$Virus_ID %in% rownames(Virus_abundance_infants),]

Virus_metadata_infants$Source <- "Virus"
Plasmid_metadata_infants$Source <- "Plasmid"
combined_data <- rbind(
  data.frame(Source = "Virus", AMR_gene_proportion = Virus_metadata_infants$AMR_gene_proportion),
  data.frame(Source = "Plasmid", AMR_gene_proportion = Plasmid_metadata_infants$AMR_gene_proportion)
)

pdf('5_AMR_ANALYSIS/PLOTS/AMR_vOTU_PTU.pdf', width=3.2, height=3.2)
ggplot(combined_data, aes(x = Source, y = AMR_gene_proportion, fill = Source)) +
  geom_violin(aes(fill=Source), alpha=0.7, width=0.8, trim = FALSE) +
  geom_jitter(aes(color=Source), alpha=0.9, width = 0.2) + 
  geom_boxplot(width=0.2, fill="white", color="black", outlier.shape = NA) + 
  labs(x = 'MGE', y = 'ARG Proportion (log scale)') + 
  scale_fill_manual(values = c("Virus" ="#404080", "Plasmid" = "#69b3a2")) + 
  scale_color_manual(values = c("Virus" ="#404080", "Plasmid" = "#69b3a2")) + 
  scale_y_log10() +
  theme_classic() +
  theme(plot.title = element_text(size = 18, hjust = 0.5, face = "bold"), 
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 15), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.position = "none")
dev.off()

##************************************************************************
# 2. Association analysis - Infants: AM abundances vs phenotypes
##************************************************************************

##################################
# Plasmid abundances vs phenotypes 
##################################
phenotypes_to_test <- colnames(Sample_metadata_infants)[c(6:14)]
AM_associations_mixed <- mixed_models(Sample_metadata_infants, "bioSampleId", AMR_abundance_data, 
                                           phenotypes_to_test)

# Remove from results multiple levels of the same phenotypes (factors) and check the number of significant phenotypes
AM_associations_mixed <- AM_associations_mixed %>%
  distinct(Pheno, Bug, .keep_all = TRUE) %>%
  arrange(FDR)

# Check significant results
AM_associations_mixed_sig <- AM_associations_mixed [AM_associations_mixed$FDR < 0.05,]

# Save results
write.table(AM_associations_mixed, "5_AMR_ANALYSIS/AM_associations_mixed_27062024.txt", sep="\t", row.names=F, quote = F)


##************************************************************************
# 3. Association analysis - Infants: AM counts vs phenotypes
##************************************************************************

##################################
# AM counts vs phenotypes 
##################################

phenotypes_to_test <- colnames(Sample_metadata_infants)[c(6:14)]
AM_count_associations_mixed <- mixed_models(Sample_metadata_infants, "bioSampleId", AMR_count_data, 
                                      phenotypes_to_test)

# Remove from results multiple levels of the same phenotypes (factors) and check the number of significant phenotypes
AM_count_associations_mixed <- AM_count_associations_mixed %>%
  distinct(Pheno, Bug, .keep_all = TRUE) %>%
  arrange(FDR)

# Check significant results
AM_count_associations_mixed_sig <- AM_count_associations_mixed [AM_count_associations_mixed$FDR < 0.05,]

# Save results
write.table(AM_count_associations_mixed, "5_AMR_ANALYSIS/AM_count_associations_mixed_27062024.txt", sep="\t", row.names=F, quote = F)


##************************************************************************
# 4. Generation of plots
##************************************************************************

# Plots exploring total AMR gene abundance and presence 
# Load AMR abundance without log transformation
AMR_abundance_data <- data.frame(cbind(AM_abundance_prev, "AMR_gene_abundance"= AMR_Gene_abundance_prev$AMR_gene_abundance),
                                 check.names = F)

Sample_metadata_infants$AMR_gene_abundance <- AMR_abundance_data$AMR_gene_abundance
Sample_metadata_infants$AMR_gene_count <- AMR_count_data$AMR_gene_counts
Sample_metadata_infants$Type <- rep("Infant", nrow(Sample_metadata_infants))

pdf('5_AMR_ANALYSIS/PLOTS/Total_AMR_abundance.pdf', width=3.2, height=3.2)
ggplot(Sample_metadata_infants, 
       aes(x = rand_AB, y = AMR_gene_abundance, fill = rand_AB)) +
  geom_violin(aes(fill=rand_AB), alpha=0.7, width=0.8, trim = F) +
  geom_jitter(aes(color=rand_AB), alpha=0.7) + 
  geom_boxplot(width=0.2, fill="white", color="black", outlier.shape = NA) + 
  scale_fill_manual(values=c("#0055AA", "#C40003")) + 
  scale_color_manual(values=c("#0055AA", "#C40003")) +
  labs(x = "AB Group", y = "AMR Gene Abundance", fill = "AB Group") +
  theme_classic() +
  theme(plot.title = element_text(size=18, hjust = 0.5, face="bold"), 
        axis.text = element_text(size=14),
        axis.title = element_text(size=15), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.position = "none") 
dev.off()

pdf('5_AMR_ANALYSIS/PLOTS/Total_AMR_count.pdf', width=3.2, height=3.2)
ggplot(Sample_metadata_infants, 
       aes(x = rand_AB, y = AMR_gene_counts, fill = rand_AB)) +
  geom_violin(aes(fill=rand_AB), alpha=0.7, width=0.8, trim = F) +
  geom_jitter(aes(color=rand_AB), alpha=0.7) + 
  geom_boxplot(width=0.2, fill="white", color="black", outlier.shape = NA) + 
  scale_fill_manual(values=c("#0055AA", "#C40003")) + 
  scale_color_manual(values=c("#0055AA", "#C40003")) +
  labs(x = "AB Group", y = "AMR Gene Count", fill = "AB Group") +
  theme_classic() +
  theme(plot.title = element_text(size=18, hjust = 0.5, face="bold"), 
        axis.text = element_text(size=14),
        axis.title = element_text(size=15), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.position = "none") 
dev.off()

pdf('5_AMR_ANALYSIS/PLOTS/Total_AMR_abundance_timepoint.pdf', width=3.9, height=3.5)
ggplot(Sample_metadata_infants, 
       aes(x=Timepoint_categorical, y=AMR_gene_abundance, fill = rand_AB)) +
  geom_boxplot(aes(fill=rand_AB), alpha=0.7, width=0.8) +
  geom_point(aes(color=rand_AB), size=1.2, position=position_jitterdodge(), alpha=0.7) +
  labs(x = "Timepoint", y = "AMR Gene Abundance", fill = "AB Group") +
  scale_fill_manual(values=c("#0055AA", "#C40003")) + 
  scale_color_manual(values=c("#0055AA", "#C40003")) +
  theme_classic()  +
  theme(plot.title = element_text(size=18, hjust = 0.5, face="bold"), 
        axis.text = element_text(size=14),
        axis.title = element_text(size=16), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.position = "bottom",
        legend.text = element_text(size=12),
        legend.title = element_text(size=13)) +
  guides(color = FALSE)
dev.off()

pdf('5_AMR_ANALYSIS/PLOTS/Total_AMR_count_timepoint.pdf', width=3.9, height=3.5)
ggplot(Sample_metadata_infants, 
       aes(x=Timepoint_categorical, y=AMR_gene_counts, fill = rand_AB)) +
  geom_boxplot(aes(fill=rand_AB), alpha=0.7, width=0.8) +
  geom_point(aes(color=rand_AB), size=1.2, position=position_jitterdodge(), alpha=0.7) +
  labs(x = "Timepoint", y = "AMR Gene Count", fill = "AB Group") +
  scale_fill_manual(values=c("#0055AA", "#C40003")) + 
  scale_color_manual(values=c("#0055AA", "#C40003")) +
  theme_classic()  +
  theme(plot.title = element_text(size=18, hjust = 0.5, face="bold"), 
        axis.text = element_text(size=14),
        axis.title = element_text(size=16), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.position = "bottom",
        legend.text = element_text(size=12),
        legend.title = element_text(size=13)) +
  guides(color = FALSE)
dev.off()

pdf('5_AMR_ANALYSIS/PLOTS/Total_AMR_count_timepoint_ungrouped.pdf', width=3.7, height=3.1)
ggplot(Sample_metadata_infants, 
       aes(x=Timepoint_categorical, y=AMR_gene_counts)) +
  geom_violin(aes(fill=Type), alpha=0.5, width=0.8, trim = F) +
  geom_jitter(alpha = 0.7, color ="#008080") +
  geom_boxplot(width=0.2, fill="white", color="black") + 
  labs(x = "Timepoint", y = "AMR Gene Count") +
  scale_fill_manual(values=c("#008080")) +
  theme_classic()  +
  theme(plot.title = element_text(size=18, hjust = 0.5, face="bold"), 
        axis.text = element_text(size=14),
        axis.title = element_text(size=16), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.position = "none")
dev.off()

pdf('2_DIVERSITY/PLOTS/Plasmid_Richness_Infants.pdf', width=3.6, height=3.1)
ggplot(Sample_metadata_infants
       , aes(x=Timepoint_categorical, y=plasmid_richness)) + 
  geom_violin(aes(fill=Type), alpha=0.5, width=0.8, trim = F) +
  geom_jitter(alpha = 0.7, color ="#008080") +
  geom_boxplot(width=0.2, fill="white", color="black") + 
  labs(x = 'Timepoint', y = 'Plasmid richness') + 
  scale_fill_manual(values=c("#008080")) + #"#008080", "#4B0082"
  theme_classic() + 
  #ylim(0,6) +
  theme(plot.title = element_text(size=18, hjust = 0.5, face="bold"), 
        axis.text = element_text(size=14),
        axis.title = element_text(size=16), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.position = "none")
dev.off()


Sample_metadata_infants$feeding_mode_pragmatic <- factor(
  Sample_metadata_infants$feeding_mode_pragmatic,
  levels = c("breast_feeding", "mixed_feeding", "formula_feeding"),
  labels = c("breastmilk", "mixed", "formula")
)

pdf('5_AMR_ANALYSIS/PLOTS/Total_AMR_abundance_timepoint_FEEDING.pdf', width=3.9, height=3.5)
ggplot(Sample_metadata_infants[!is.na(Sample_metadata_infants$feeding_mode_pragmatic), ],
       aes(x=Timepoint_categorical, y=AMR_gene_abundance, fill = feeding_mode_pragmatic)) +
  geom_boxplot(aes(fill=feeding_mode_pragmatic), alpha=0.7, width=0.8) +
  geom_point(aes(color=feeding_mode_pragmatic), size=1.2, position=position_jitterdodge(), alpha=0.7) +
  labs(x = "Timepoint", y = "AMR Gene Abundance", fill = "Feeding mode") +
  scale_fill_manual(values=c("#FFA255", "#BFE2C5", "#9A9AFF")) + 
  scale_color_manual(values=c("#FFA255", "#BFE2C5", "#9A9AFF")) +
  theme_classic()  +
  theme(plot.title = element_text(size=18, hjust = 0.5, face="bold"), 
        axis.text = element_text(size=14),
        axis.title = element_text(size=16), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.position = "bottom") +
  guides(color = FALSE)
dev.off()

pdf('5_AMR_ANALYSIS/PLOTS/Total_AMR_count_timepoint_FEEDING.pdf', width=3.9, height=3.5)
ggplot(Sample_metadata_infants[!is.na(Sample_metadata_infants$feeding_mode_pragmatic), ],
       aes(x=Timepoint_categorical, y=AMR_gene_counts, fill = feeding_mode_pragmatic)) +
  geom_boxplot(aes(fill=feeding_mode_pragmatic), alpha=0.7, width=0.8) +
  geom_point(aes(color=feeding_mode_pragmatic), size=1.2, position=position_jitterdodge(), alpha=0.7) +
  labs(x = "Timepoint", y = "AMR Gene Count", fill = "Feeding mode") +
  scale_fill_manual(values=c("#FFA255", "#BFE2C5", "#9A9AFF")) + 
  scale_color_manual(values=c("#FFA255", "#BFE2C5", "#9A9AFF")) +
  theme_classic()  +
  theme(plot.title = element_text(size=18, hjust = 0.5, face="bold"), 
        axis.text = element_text(size=14),
        axis.title = element_text(size=16), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.position = "bottom") +
  guides(color = FALSE)
dev.off()


# Plots showing association found with phenotypes
Sample_metadata_infants$streptomycin <- AMR_abundance_data$streptomycin
Sample_metadata_infants$gentamicin <- AMR_abundance_data$gentamicin
Sample_metadata_infants$`nalidixic acid` <- AMR_abundance_data$`nalidixic acid`

pdf('5_AMR_ANALYSIS/PLOTS/Streptomycin_Gestational_age.pdf', width=2.8, height=3.2)
ggplot(Sample_metadata_infants[!is.na(Sample_metadata_infants$streptomycin), ], aes(x = streptomycin , y = preg_gest_age)) +
  geom_point(alpha = 0.7) + 
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue", alpha = 0.3) +
  labs(x ='Streptomycin Resistance', y = 'Gestational Age') +
  theme_classic() +
  theme(plot.title = element_text(size = 18, hjust = 0.5, face = "bold"), 
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 15), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.position = "none")
dev.off()


pdf('5_AMR_ANALYSIS//PLOTS/Gentamicin_birthweight.pdf', width=2.8, height=3.2)
ggplot(Sample_metadata_infants[!is.na(Sample_metadata_infants$gentamicin), ], aes(x = gentamicin , y = infant_birthweight)) +
  geom_point(alpha = 0.7) + 
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue", alpha = 0.3) +
  labs(x ='Gentamicin Resistance', y = 'Infant birthweight (g)') +
  theme_classic() +
  theme(plot.title = element_text(size = 18, hjust = 0.5, face = "bold"), 
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 15), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.position = "none")
dev.off()

pdf('5_AMR_ANALYSIS//PLOTS/nalidixic_acid_BMI.pdf', width=2.8, height=3.2)
ggplot(Sample_metadata_infants[!is.na(Sample_metadata_infants$`nalidixic acid`), ], aes(x = `nalidixic acid` , y = pre_preg_bmi_mother)) +
  geom_point(alpha = 0.7) + 
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue", alpha = 0.3) +
  labs(x ='Nalidixic acid Resistance', y = 'Infant birthweight (g)') +
  theme_classic() +
  theme(plot.title = element_text(size = 18, hjust = 0.5, face = "bold"), 
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 15), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.position = "none")
dev.off()

################################################
# Heatmap of AMR counts and abundances
################################################

# Prepare column annotations
Sample_metadata_infants_heatmap <- Sample_metadata_infants[, c("Timepoint_categorical",
                                                               "feeding_mode_pragmatic", "rand_AB")]

rownames(Sample_metadata_infants_heatmap) <- Sample_metadata_infants$bioSampleId


# Prepare AMR class abundances
AMR_abundance_infants_heatmap <- AMR_Class_abundance_prev 
rownames(AMR_abundance_infants_heatmap) <- Sample_metadata_infants$NG_ID
AMR_abundance_infants_heatmap <- data.frame(t(AMR_abundance_infants_heatmap))
colnames(AMR_abundance_infants_heatmap) <- Sample_metadata_infants$bioSampleId

rows_to_modify <- grep(".antibiotic", rownames(AMR_abundance_infants_heatmap))
rownames(AMR_abundance_infants_heatmap)[rows_to_modify] <- gsub(".antibiotic", "", rownames(AMR_abundance_infants_heatmap)[rows_to_modify])
rownames(AMR_abundance_infants_heatmap) <- gsub("\\.", " ", rownames(AMR_abundance_infants_heatmap))
rownames(AMR_abundance_infants_heatmap) [7] <- "antiseptics"

# Prepare AMR class counts
AMR_count_infants_heatmap <- AMR_Class_presence_prev
rownames(AMR_count_infants_heatmap) <- Sample_metadata_infants$NG_ID
AMR_count_infants_heatmap <- data.frame(t(AMR_count_infants_heatmap))
colnames(AMR_count_infants_heatmap) <- Sample_metadata_infants$bioSampleId

rows_to_modify <- grep(".antibiotic", rownames(AMR_count_infants_heatmap))
rownames(AMR_count_infants_heatmap)[rows_to_modify] <- gsub(".antibiotic", "",
                                                            rownames(AMR_count_infants_heatmap)[rows_to_modify])
rownames(AMR_count_infants_heatmap) <- gsub("\\.", " ", rownames(AMR_count_infants_heatmap))
rownames(AMR_count_infants_heatmap) [7] <- "antiseptics"

# Create the heatmap excluding columns from the legend (default euclidean clustering)
annotation_colors <- list(
  rand_AB = c("No" = alpha("#0055AA", 0.7), "Yes" = alpha("#C40003", 0.7)),
  feeding_mode_pragmatic = c("breastmilk" = alpha("#FFBE81", 0.7), "formula" = alpha("#C7C7FF", 0.7),
                             "mixed" = alpha("#DFF0D8", 0.7)),
  Timepoint_categorical = c("W01" = alpha("#007ED3", 0.7), "W02" = alpha("#7FD2FF", 0.7),
                            "W03" = alpha("#96B964", 0.7), "W04" = alpha("#EAC862", 0.7),
                            "W05" = alpha("#FF9D1E", 0.7), "W06" = alpha("#C40003", 0.7))
)

pdf('5_AMR_ANALYSIS/PLOTS/Heatmap_AMR_class_abundances.pdf', width=10, height=7.2)
ComplexHeatmap::pheatmap(AMR_abundance_infants_heatmap,
                         annotation_col = Sample_metadata_infants_heatmap, 
                         color = hcl.colors(50, "Oslo"),
                         show_colnames = FALSE, 
                         annotation_legend = TRUE, 
                         fontsize_row = 12, 
                         cluster_rows = T, cluster_cols = T,
                         annotation_colors = annotation_colors)
dev.off()

pdf('5_AMR_ANALYSIS/PLOTS/Heatmap_AMR_class_counts.pdf', width=10, height=7.2)
ComplexHeatmap::pheatmap(AMR_count_infants_heatmap, 
                         annotation_col = Sample_metadata_infants_heatmap, 
                         color = rev(hcl.colors(50, "Oslo")),
                         show_colnames = FALSE, 
                         annotation_legend = TRUE, 
                         fontsize_row = 12, 
                         cluster_rows = T, cluster_cols = T,
                         annotation_colors = annotation_colors)
dev.off()

