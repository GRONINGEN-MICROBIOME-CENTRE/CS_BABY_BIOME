# Load R packages
library(stats)
library(Biostrings)
library(data.table)

# Get FASTA file of plasmid genomes from command line arguments
args <- commandArgs(trailingOnly = TRUE)
plasmid_sequences <- readDNAStringSet(args[1]) # FASTA file with the plasmid sequences 

# Get sequence IDs
plasmid_origin <- data.frame(names(plasmid_sequences))

# Add DB information
colnames(plasmid_origin)[1] <- "plasmid_seq"
plasmid_origin$Origin <- NA
plasmid_origin[,"Origin"] [grep("^IMGPR", plasmid_origin$plasmid_seq)] <- "IMG_circular_plasmid" # 3472
plasmid_origin[,"Origin"] [grep("RNODE|type_circular", plasmid_origin$plasmid_seq)] <- "Circular_plasmid" #750
plasmid_origin[,"Origin"][which(is.na(plasmid_origin$Origin))] <- "Plasmid_fragment" #2640


# Add new simplified plasmid IDs 
plasmid_origin$
  new_ID <- ave(plasmid_origin$Origin, plasmid_origin$Origin, 
                                                 FUN = function(x) paste0(x, "_", seq_along(x)))


# Replace the IDs in the FASTA file with the new simplified IDs
names(plasmid_sequences) <- plasmid_origin$new_ID

# Save final table with the X sequences used as input for dereplication and their DB of origin 
write.table(plasmid_origin,"Renamed_plasmids.txt", sep = "\t", 
            row.names = F, col.names = F, quote = FALSE)

# Save the updated FASTA file
writeXStringSet(plasmid_sequences, "rep_seqs_all_plasmids_renamed.fna")
