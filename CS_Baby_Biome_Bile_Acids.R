#################### CS Baby Biome Bile Acids Analysis ########################

# Date: 24th of January, 2024
# Last update: 24th March, 2024
# Author: Trishla Sinha 

library(tidyverse)
library(vegan)
library(lmerTest)
library(wesanderson)


# Load functions 

urso.group    <- c('UDCA', 'GUDCA', 'TUDCA')
ca.group      <- c('CA', "GCA", 'TCA', 'DCA','GDCA', 'TDCA')
cdca.group    <- c('CDCA','GCDCA', 'TCDCA',"LCA", 'GLCA', 'TLCA') # 'TLCA_3S','GLCA_3S'
udca.group    <- c('UDCA', 'GUDCA',"TUDCA")
unconjugated  <- c('CA','CDCA','DCA','LCA')
conjugated    <- c('TCA','GCA','TCDCA','GCDCA','TDCA', 'GDCA', 'TLCA','GLCA') # 'TLCA_3S','GLCA_3S'
ca.deconjugated    <- c('CA','TCA', 'GCA')
ca.dehydroxylation <- c('DCA','TDCA', 'GDCA')
ca.deconjugated    <- c('DCA','TDCA', 'GDCA')
ca.dehydroxylation <- c('CA','TCA', 'GCA')
glycine.group <- c('GCA', 'GCDCA','GDCA', 'GLCA') # 'GLCA_3S'
taurine.group <- c('TCA','TCDCA','TDCA','TLCA') # 'TLCA_3S'

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
bile_acids<-read.delim("~/Desktop/CS_Baby_Biome/submission/2024_submission/analysis/lab/bile_acid_measurements_TS.txt")
link<-read.delim("~/Desktop/CS_Baby_Biome/submission/2024_submission/analysis/lab/lab_log_TS.txt")
link<-link[,c(1,4,11)]
bile_acids_link<-left_join(link, bile_acids)
names (bile_acids_link)[1]<-"CS_BABY_BIOME_ID"
AB<-read.delim("~/Desktop/CS_Baby_Biome/submission/2024_submission/analysis/lab/AB_samples.txt")
feeding<-read.delim("~/Desktop/CS_Baby_Biome/submission/2024_submission/analysis/lab/feeding_samples.txt")

analyses_df<-left_join(bile_acids_link, AB)
analyses_df <- analyses_df[!is.na(analyses_df$rand_AB), ]
analyses_df$SAMPLE_ID<-paste0(analyses_df$CS_BABY_BIOME_ID, "_", analyses_df$Timepoint)
row.names(analyses_df)<-analyses_df$SAMPLE_ID
analyses_df<-left_join(analyses_df, feeding)
row.names(analyses_df)<-analyses_df$SAMPLE_ID

#write.table(analyses_df, "Bile_acids_CS_BABY_BIOME_02_04_2024.txt", sep = "\t", quote = F)


bile_acids_only <- grep("^(?!.*birthcard).*CA.*$", names(analyses_df), ignore.case = TRUE, perl = TRUE, value = TRUE)
only_phenos <- analyses_df[, !names(analyses_df) %in% bile_acids_only]
bile_acids_selection<-analyses_df[,bile_acids_only]
bile_acids_conc<-bile_acids_selection
bile_acids_prop <- bile_acids_conc / rowSums(bile_acids_conc) * 100

#write.table(only_phenos, "metadata_cs_baby_biome.txt", sep="\t", row.names=F, quote = F )

all_basic<- analyses_df[, !names(analyses_df) %in% bile_acids_only]

all_basic_bile_acids_prop <- merge(bile_acids_prop, all_basic, by="row.names")
row.names(all_basic_bile_acids_prop)<-all_basic_bile_acids_prop$Row.names
all_basic_bile_acids_prop$Row.names=NULL
all_basic_bile_acids_prop<- all_basic_bile_acids_prop[match(rownames(bile_acids_prop), rownames(all_basic_bile_acids_prop)), ]
all_basic_bile_acids_prop$Timepoint<-factor(all_basic_bile_acids_prop$Timepoint, levels = c("W01", "W04",  "W06"))
all_basic_bile_acids_prop$feeding_mode_pragmatic<-factor(all_basic_bile_acids_prop$feeding_mode_pragmatic, levels = c("breast_feeding", "mixed_feeding", "formula_feeding"))

all_basic_bile_acids_prop$total_primary<-all_basic_bile_acids_prop$CA+all_basic_bile_acids_prop$CDCA+all_basic_bile_acids_prop$TCDCA+all_basic_bile_acids_prop$GCDCA+all_basic_bile_acids_prop$TCA+all_basic_bile_acids_prop$GCA

all_basic_bile_acids_prop$total_secondary<-all_basic_bile_acids_prop$DCA+all_basic_bile_acids_prop$TDCA+all_basic_bile_acids_prop$GDCA+all_basic_bile_acids_prop$UDCA+all_basic_bile_acids_prop$TUDCA+all_basic_bile_acids_prop$GUDCA+
                                           all_basic_bile_acids_prop$LCA+all_basic_bile_acids_prop$TLCA+all_basic_bile_acids_prop$GLCA+all_basic_bile_acids_prop$LCA.3S+all_basic_bile_acids_prop$TLCA.3S.quant+all_basic_bile_acids_prop$GLCA.3S.quant



phenotypes<-all_basic_bile_acids_prop[, c("rand_AB", "feeding_mode_pragmatic", "CS_BABY_BIOME_ID", "Timepoint")]
names (phenotypes)[1]<-"AB"
names (phenotypes)[2]<-"feeding_mode"
phenotypes$ID<-row.names(phenotypes)
bile_acids_to_analyze<-all_basic_bile_acids_prop[, c(1:17, 26:27)]

# Add ratio's to calculations 
bile_acids_to_analyze$ratio_pri_secon <- bile_acids_to_analyze$total_primary / (bile_acids_to_analyze$total_primary + bile_acids_to_analyze$total_secondary)
bile_acids_to_analyze$CA_CDCA_ratio <- rowSums(bile_acids_to_analyze[,ca.group])/rowSums(bile_acids_to_analyze[,cdca.group])
#bile_acids_to_analyze$Unconjugated_conjugated_ratio <- rowSums(bile_acids_to_analyze[,unconjugated])/rowSums(bile_acids_to_analyze[,conjugated])
#bile_acids_to_analyze$CA_dehydro_deconju_ratio      <- rowSums(bile_acids_to_analyze[,ca.dehydroxylation])/rowSums(bile_acids_to_analyze[,ca.deconjugated])
#bile_acids_to_analyze$Taurine_glycine_ratio         <- rowSums(bile_acids_to_analyze[,taurine.group])/rowSums(bile_acids_to_analyze[,glycine.group])


# Association of bile acids with AB and feeding 
BA_mixed_all <- mixed_models(phenotypes, "ID", bile_acids_to_analyze, c("AB", "feeding_mode"))

# now correcting for feeding 
BA_mixed_all_cor_feeding<-mixed_models_cor_feeding(phenotypes, "ID", bile_acids_to_analyze, c("AB"))

# Looking at only timepoint 
BA_mixed_all_timepoint<-mixed_models_timepoint(phenotypes, "ID", bile_acids_to_analyze, c("Timepoint"))

write.table(BA_mixed_all, "BA_mixed_all.txt", sep = "\t", quote = F)

write.table(BA_mixed_all_timepoint, "BA_mixed_all_timepoint.txt", sep = "\t", quote = F)

# Timepoint 
long_data_tp <- all_basic_bile_acids_prop %>%
  pivot_longer(cols = c("TLCA.3S.quant", "GCA", "GCDCA"), 
               names_to = "Variable", 
               values_to = "Value")


A=ggplot(long_data_tp, aes(x = Timepoint, y = Value, fill = Timepoint, color=Timepoint)) +
  geom_boxplot(alpha=0.4, outlier.colour = NA)  +
  geom_point(alpha=0.6,
             position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0)) +
  facet_wrap(~ Variable, scales = "free_y") +
  theme_bw() +
  scale_color_manual(values = c("indianred2", "mediumspringgreen", "purple2" ))+
  scale_fill_manual(values = c("indianred2", "mediumspringgreen","purple2" ))+
  labs(x = "Timepoint", y = "Relative abundance of Bile Acid") +
  theme(legend.position = "none")


long_data_ab <- all_basic_bile_acids_prop %>%
  pivot_longer(cols = c("GCDCA","GUDCA","GLCA"), 
               names_to = "Variable", 
               values_to = "Value") %>%
  select(Timepoint, rand_AB, Variable, Value)



C=ggplot(long_data_ab, aes(x = Timepoint, y = Value, fill = rand_AB, color=rand_AB)) +
  geom_boxplot(alpha = 0.4, outlier.colour = NA, position = position_dodge(width = 0.8)) +
  geom_point(alpha = 0.6, position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.8)) +
  scale_fill_manual(values = c("#54b7a6", "#fa2b13")) +
  scale_color_manual(values = c("#54b7a6", "#fa2b13")) +
  theme_bw() +
  labs(x = "Timepoint", y = "Relative abundance of Bile Acid") +
  theme(legend.position = "bottom") +
  facet_wrap(~ Variable, scales = "free_y", ncol = 3)




long_data_feed <- all_basic_bile_acids_prop %>%
  pivot_longer(cols = c("TCDCA","TCA", "CA"), 
               names_to = "Variable", 
               values_to = "Value") %>%
  select(Timepoint, feeding_mode_pragmatic, Variable, Value)



B=ggplot(long_data_feed, aes(x = Timepoint, y = Value, fill = feeding_mode_pragmatic, color=feeding_mode_pragmatic)) +
  geom_boxplot(alpha = 0.4, outlier.colour = NA, position = position_dodge(width = 0.8)) +
  geom_point(alpha = 0.6, position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.8)) +
  scale_fill_manual(values = wes_palette("Moonrise3", n = 5, name=""))+
  scale_color_manual(values = wes_palette("Moonrise3", n = 5, name=""))+
  theme_bw() +
  labs(x = "Timepoint", y = "Relative abundance of Bile Acid") +
  theme(legend.position = "none") +
  facet_wrap(~ Variable, scales = "free_y", ncol = 3)


B

all_ba_p_dis<-vegdist(bile_acids_prop, method = 'canberra', na.rm = T)
all_ba_p_dis_mds<-cmdscale(all_ba_p_dis, k=5, eig = T)
all_ba_p_pcoa <- data.frame(all_ba_p_dis_mds$points)

p_ba_p_pcoa_timepoint<-ggplot(all_ba_p_pcoa ,aes(X1,X2, color = all_basic_bile_acids_prop$Timepoint))+
  geom_point(size = 2,alpha = 0.5)+
  stat_ellipse(aes(group = all_basic_bile_acids_prop$Timepoint, fill = all_basic_bile_acids_prop$Timepoint, color = all_basic_bile_acids_prop$Timepoint) ,type = "norm",linetype = 2,geom = "polygon",alpha = 0.05,show.legend = F)+
  xlab(paste("PCo1=",round(all_ba_p_dis_mds$eig[1],digits = 2),"%",sep = ""))+
  ylab(paste("PCo2=",round(all_ba_p_dis_mds$eig[2],digits = 2),"%",sep = ""))+
  scale_color_manual(name=NULL, 
                     breaks = c("W01",  "W04", "W06"),
                     labels = c("W01          ", "W04       ", "W06"),
                     values = c("indianred2", "mediumspringgreen","purple2" ))+
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
p_ba_p_pcoa_timepoint<-ggExtra::ggMarginal(p_ba_p_pcoa_timepoint, type = "histogram", groupColour = F, groupFill = TRUE,
                                           xparams = list(bins = 60, alpha = 0.5,position = 'identity', color = 'white'),
                                           yparams = list(bins = 60, alpha = 0.5,position = 'identity', color = 'white'))

p_ba_p_pcoa_timepoint

p_ba_p_pcoa_feeding<-ggplot(all_ba_p_pcoa ,aes(X1,X2, color = all_basic_bile_acids_prop$feeding_mode_pragmatic))+
  geom_point(size = 2,alpha = 0.8)+
  stat_ellipse(aes(group =all_basic_bile_acids_prop$feeding_mode_pragmatic, fill = all_basic_bile_acids_prop$feeding_mode_pragmatic, color = all_basic_bile_acids_prop$feeding_mode_pragmatic) ,type = "norm",linetype = 2,geom = "polygon",alpha = 0.05,show.legend = F)+
  xlab(paste("PCo1=",round(all_ba_p_dis_mds$eig[1],digits = 2),"%",sep = ""))+
  ylab(paste("PCo2=",round(all_ba_p_dis_mds$eig[2],digits = 2),"%",sep = ""))+
  scale_fill_manual(values = wes_palette("Moonrise3", n = 5))+
  scale_color_manual(name=NULL, 
                     breaks = c("breast_feeding", "mixed_feeding" ,"formula_feeding"),
                     labels = c("Breast", "Mixed", "Formula"),
                     values = wes_palette("Moonrise3", n = 5))+
  theme(plot.subtitle = element_text(vjust = 1),
        plot.caption = element_text(vjust = 1),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.key = element_rect(fill = NA),
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        panel.background = element_rect(fill = NA))
p_ba_p_pcoa_feeding<-ggExtra::ggMarginal(p_ba_p_pcoa_feeding, type = "histogram", groupColour = F, groupFill = TRUE,
                                 xparams = list(bins = 60, alpha = 0.8,position = 'identity', color = 'white'),
                                 yparams = list(bins = 60, alpha = 0.8,position = 'identity', color = 'white'))
p_ba_p_pcoa_feeding




p_ba_p_pcoa_rand_AB<-ggplot(all_ba_p_pcoa ,aes(X1,X2, color = all_basic_bile_acids_prop$rand_AB))+
  geom_point(size = 2,alpha = 0.8)+
  stat_ellipse(aes(group =all_basic_bile_acids_prop$rand_AB, fill = all_basic_bile_acids_prop$rand_AB, color = all_basic_bile_acids_prop$rand_AB) ,type = "norm",linetype = 2,geom = "polygon",alpha = 0.05,show.legend = F)+
  xlab(paste("PCo1=",round(all_ba_p_dis_mds$eig[1],digits = 2),"%",sep = ""))+
  ylab(paste("PCo2=",round(all_ba_p_dis_mds$eig[2],digits = 2),"%",sep = ""))+
  scale_fill_manual(values = c("#54b7a6", "#fa2b13"))+
  scale_color_manual(name=NULL, 
                     breaks = c("No", "Yes"),
                     labels = c("AB-"          , "AB+"),
                     values = c("#54b7a6", "#fa2b13"))+
  theme(plot.subtitle = element_text(vjust = 1),
        plot.caption = element_text(vjust = 1),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.key = element_rect(fill = NA),
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        panel.background = element_rect(fill = NA))

p_ba_p_pcoa_rand_AB<-ggExtra::ggMarginal(p_ba_p_pcoa_rand_AB, type = "histogram", groupColour = F, groupFill = TRUE,
                                         xparams = list(bins = 60, alpha = 0.8,position = 'identity', color = 'white'),
                                         yparams = list(bins = 60, alpha = 0.8,position = 'identity', color = 'white'))
p_ba_p_pcoa_rand_AB

library(gridExtra)


grid.arrange(p_ba_p_pcoa_timepoint,A,  p_ba_p_pcoa_feeding, B, p_ba_p_pcoa_rand_AB, C, nrow = 3)





adonis2(bile_acids_prop~all_basic_bile_acids_prop$Timepoint)
adonis2(bile_acids_prop~all_basic_bile_acids_prop$feeding_mode_pragmatic)
adonis2(bile_acids_prop~all_basic_bile_acids_prop$rand_AB)
adonis2(bile_acids_prop~all_basic_bile_acids_prop$feeding_mode_pragmatic+ all_basic_bile_acids_prop$rand_AB)


setwd("/Users/trishlasinha/Desktop/CS_Baby_Biome/supplementary/tables/individual")
write.table(BA_mixed_all_timepoint, "BA_associations_timepoint_31_03_2024.txt", sep="\t", row.names=F, quote = F)
write.table(BA_mixed_all, "BA_associations_AB_feeding_31_03_2024.txt", sep="\t", row.names=F, quote = F)
write.table(BA_mixed_all_cor_feeding, "BA_associations_cor_feeding_AB_31_03_2024.txt", sep="\t", row.names=F, quote = F)
