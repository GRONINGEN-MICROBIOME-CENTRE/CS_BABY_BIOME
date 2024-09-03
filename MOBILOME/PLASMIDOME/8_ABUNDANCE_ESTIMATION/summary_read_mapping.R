# Get names of plasmid genomes and samples from command line arguments
args <- commandArgs(trailingOnly = TRUE)
plasmid_names <- readLines(args[1]) # vector with the names of plasmid genomes
sample_names <- readLines(args[2]) # vector of sample names
BED_output <- args[3] # path to BED coverage output
reads_table <- read.delim(args[4], header=T) # file with number of reads per sample

#Set sample names as row names in read_table (set colnames[2] to "clean_reads" if not done before)
rownames(reads_table) <- reads_table[, 1]

# Create an empty matrix to store the results
results_table <- matrix(0, nrow = length(plasmid_names), ncol = length(sample_names), dimnames = list(plasmid_names, sample_names))

# Loop through each sample in sample_names
for (sample in sample_names) {
    # Set the location of the result files from Bedtools coverage output
    file_path <- paste0(BED_output, sample, '.coverage.txt')
    # Read the table 
    cov_table <- read.table(file_path, sep = '\t', row.names = 1, header = F, stringsAsFactors = F)
    # Subset the table to only include rows corresponding to the genomes in plasmid_names
    cov_table <- cov_table[plasmid_names, ]
    # Find the indices of plasmid genomes where the coverage is greater than or equal to 75%
    idx <- which(cov_table[, 6] >= 0.75)
    # Calculate the abundance values for results_table (column 3 = number of reads mapped; column 5 = plasmid sequence lenght)
    results_table[idx, sample] <- cov_table[idx, 3] / cov_table[idx, 5] / reads_table[sample, 'clean_reads'] * 10^9
}

# Write abundance table
write.table(results_table, file="CS_Baby_Biome_Plasmid_Abundance_Table.txt", sep = "\t")
