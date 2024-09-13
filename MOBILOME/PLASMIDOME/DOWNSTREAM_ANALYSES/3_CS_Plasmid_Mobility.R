################################################################################
##### CS Baby Biome: Plasmid mobility analysis: CONJscan and oriT
### Author(s): Asier Fernández-Pato
### Last updated: 22th June, 2024
################################################################################

# Set working directory 
setwd("~/Desktop/PhD/Projects/CS Baby Biome/Final_analysis/PLASMIDS/")

#****************
# Load libraries
#****************
library(dplyr)
library(plyranges)
library(data.table)


##************************************************************************
# 1. Load results from CONJscan and oriT blastN and process them
#*************************************************************************
CONJcan_result <- read.delim("3_MOBILITY/Plasmid_CONJ_MOB_elements.tsv") 
oriT_result <- read.delim("3_MOBILITY/oriT_blast_with_coverage.tsv") 

# Select only oriT matches with target coverage of 50% 
oriT_result <- oriT_result[oriT_result$qcov > 50, ]

# Add a new ID combining orti ID and sequence ID
oriT_result <- oriT_result %>%
  group_by(sseqid) %>%
  mutate(unique_sseqid = paste0(sseqid, "_", row_number())) %>%
  ungroup()

# Filter overlapping matches selecting only the one with lowest e-value
oriT_result <- oriT_result[,c("qseqid", "sseqid", "pident", "length", "mismatch", 
                              "gapopen", "qstart", "qend", "sstart", "send",
                              "evalue", "bitscore")]

test <- oriT_result %>% dplyr::select(qseqid, sseqid, sstart, send, evalue) 

test <- test %>% mutate(i_start = case_when(
  sstart < send ~ sstart,
  send < sstart ~ send
))

test <- test %>% mutate(i_end = case_when(
  sstart < send ~ send,
  send < sstart ~ sstart
))

test <- test %>% mutate(strand = case_when(
  sstart < send ~ "+",
  send < sstart ~ "-"
))

test <- test %>% select(qseqid, sseqid, i_start, i_end, evalue, strand)
colnames(test) <- c("querynames", "seqnames", "start", "end", "evalue", "strand")
test_irange <- test %>% as_granges()
test_disjoin <- reduce(test_irange,with.revmap=TRUE)
list_revmap <- as.data.frame(mcols(test_disjoin))

filtered_data <- c()
for(i in 1:nrow(list_revmap)){
  filtered_data <- c(filtered_data, 
                     (slice(list_revmap, i) %>% 
                        unlist(use.names=FALSE))[which.min(slice(test,  slice(list_revmap, i) %>%
                                                                   unlist(use.names=FALSE))$evalue)])
}

oriT_result_processed <- slice(test, filtered_data)
colnames(oriT_result_processed) <- c("oriT_ID", "Plasmid_ID", "start", "end", "evalue", "strand")

# Merge CONJscan and oriT results
Mobility_results <- left_join(CONJcan_result, oriT_result_processed, by="Plasmid_ID")
Mobility_results <- Mobility_results[, c("Plasmid_ID","CONJ_MOB_elements","oriT_ID" )]

Mobility_results <- Mobility_results %>%
  group_by(Plasmid_ID, CONJ_MOB_elements) %>%
  summarize(oriT_ID = paste(na.omit(oriT_ID), collapse = ","), .groups = 'drop') %>%
  ungroup()

# Replace empty strings with NA
Mobility_results <- Mobility_results %>%
  mutate_all(~na_if(., ""))

# Classify plasmids as: conjugative, mobilizable or non-mobilizable
##- Plasmids encoding complete conjugation systems —> Conjugative
##- Plasmids encoding a relaxase or an oriT —> mobilizable
##- Remaining plasmids —> Non-mobilizable
Mobility_results <- Mobility_results %>%
  mutate(
    Plasmid_mobility = case_when(
      grepl("T4SS", CONJ_MOB_elements) ~ "Conjugative",
      grepl("MOB", CONJ_MOB_elements) ~ "Mobilizable",
      !is.na(oriT_ID) ~ "Mobilizable",
      TRUE ~ "Non-mobilizable"
    )
  )

#****************
# Write results
#****************
write.table(Mobility_results,"3_MOBILITY/CS_Plasmid_Mobility.txt", sep = "\t", row.names = F, quote = FALSE) # Plasmid metadata


