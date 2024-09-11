# Get file with clusters that have plasmid fragments
args <- commandArgs(trailingOnly = TRUE)
clusters <- read.delim(args[1], sep = ",", header = F)

# Get a vector of plasmid fragment names
fragments <- as.vector(grep("fragment", unlist(clusters), value = T))

# Save list of plasmid fragments to remove
write.table(fragments,"Plasmid_fragments_to_remove.txt", sep = "\t", 
            row.names = F, col.names = F, quote = FALSE)
