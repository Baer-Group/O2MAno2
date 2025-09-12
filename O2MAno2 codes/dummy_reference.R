# Set current folder of the R script as the working directory
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
getwd()
library(GenomicRanges)
library(IRanges)
library(Biostrings)
library(dplyr)

# Step 1: Load the reference genome (multi-chromosome FASTA)
genome <- readDNAStringSet("c_elegans.PRJNA13758.WS292.genomic.fa")

# Initialize empty data frames to log base substitutions and indels
mutation_log <- data.frame(Chromosome = character(), Position = numeric(), 
                           Mutation_Type = character(), Original_Sequence = character(), 
                           New_Sequence = character(), stringsAsFactors = FALSE)

vcf_data <- data.frame(CHROM = character(), POS = numeric(), ID = character(),
                       REF = character(), ALT = character(), QUAL = character(), 
                       FILTER = character(), INFO = character(), stringsAsFactors = FALSE)

# Random base generator
random_base <- function(current_base) {
  bases <- c("A", "T", "C", "G")
  bases <- bases[bases != current_base]  # Exclude the current base
  return(sample(bases, 1))  # Randomly pick another base
}

# Placeholder function to check if a region is coding or non-coding
is_coding_region <- function(chr_name, position) {
  return(position <= 50000)  # For demonstration purposes only
}

# Indel length generator based on coding or non-coding region
indel_length <- function(region_type) {
  if (region_type == "coding") {
    return(sample(c(3, 6, 9, 12, 15, 18), 1))  # Only multiples of 3 for coding regions
  } else {
    return(sample(1:20, 1))  # Any length for non-coding regions
  }
}

# Step 2: Iterate over all chromosomes in the genome
set.seed(123)  # For reproducibility
for (chr_name in names(genome)) {
  genome_seq <- as.character(genome[[chr_name]])
  
  # Genome length for the current chromosome
  genome_length <- nchar(genome_seq)
  
  # Step 3: Generate 200 random base substitutions
  substitution_positions <- sample(1:genome_length, 200, replace = FALSE)
  
  # Log substitutions without modifying the genome
  for (pos in substitution_positions) {
    original_base <- substr(genome_seq, pos, pos)
    new_base <- random_base(original_base)
    
    # Log the substitution in the mutation log
    mutation_log <- rbind(mutation_log, data.frame(Chromosome = chr_name, 
                                                   Position = pos, 
                                                   Mutation_Type = "Substitution", 
                                                   Original_Sequence = original_base, 
                                                   New_Sequence = new_base))
    
    # Log substitution in VCF format
    vcf_data <- rbind(vcf_data, data.frame(CHROM = chr_name, 
                                           POS = pos, 
                                           ID = ".", 
                                           REF = original_base, 
                                           ALT = new_base, 
                                           QUAL = ".", FILTER = ".", INFO = "."))
  }
  
  # Step 4: Generate 1800 random indels (with random size between 1 and 20 bp, and frame-preserving indels in coding regions)
  indel_positions <- sample(1:genome_length, 1800, replace = FALSE)
  
  for (pos in indel_positions) {
    region_type <- ifelse(is_coding_region(chr_name, pos), "coding", "non-coding")
    indel_type <- sample(c("Insertion", "Deletion"), 1)
    indel_len <- indel_length(region_type)
    
    if (indel_type == "Insertion") {
      # Generate a random sequence to insert
      new_sequence <- paste0(sample(c("A", "T", "C", "G"), indel_len, replace = TRUE), collapse = "")
      
      # Use the reference base at the position for context (we won't modify it)
      ref_base <- substr(genome_seq, pos - 1, pos - 1)  # Get the base before the insertion site
      
      # Log insertion in the mutation log
      mutation_log <- rbind(mutation_log, data.frame(Chromosome = chr_name, 
                                                     Position = pos, 
                                                     Mutation_Type = "Insertion", 
                                                     Original_Sequence = "-", 
                                                     New_Sequence = new_sequence))
      
      # Log insertion in VCF format (REF is the adjacent base, ALT is the adjacent base + inserted sequence)
      vcf_data <- rbind(vcf_data, data.frame(CHROM = chr_name, 
                                             POS = pos - 1,  # Adjacent base position
                                             ID = ".", 
                                             REF = ref_base, 
                                             ALT = paste0(ref_base, new_sequence), 
                                             QUAL = ".", FILTER = ".", INFO = "."))
    } else {
      # Use the reference base for the deletion context
      ref_base <- substr(genome_seq, pos - 1, pos - 1)
      
      # Get the sequence to delete (including the adjacent base)
      deleted_sequence <- substr(genome_seq, pos - 1, pos + indel_len - 1)
      
      # Log deletion in the mutation log
      mutation_log <- rbind(mutation_log, data.frame(Chromosome = chr_name, 
                                                     Position = pos, 
                                                     Mutation_Type = "Deletion", 
                                                     Original_Sequence = deleted_sequence, 
                                                     New_Sequence = "-"))
      
      # Log deletion in VCF format (REF is the deleted sequence, ALT is the adjacent base)
      vcf_data <- rbind(vcf_data, data.frame(CHROM = chr_name, 
                                             POS = pos - 1,  # Adjacent base position
                                             ID = ".", 
                                             REF = deleted_sequence, 
                                             ALT = ref_base, 
                                             QUAL = ".", FILTER = ".", INFO = "."))
    }
  }
}

# Step 5: Save the log and VCF data to CSV and VCF files
write.csv(mutation_log, "mutation_log_it2.csv", row.names = FALSE)

# Write VCF header
vcf_header <- c("##fileformat=VCFv4.2", 
                paste0("##source=Dummy_Variants"),
                paste0("##reference=c_elegans.PRJNA13758.WS292.genomic.fa"),
                paste0("##contig=<ID=", names(genome), ",length=", width(genome), ">"),
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")

# Write the VCF file with header
writeLines(vcf_header, "dummy_variants_it2.vcf")
write.table(vcf_data, "dummy_variants_it2.vcf", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)

# Output a message when the process is complete
cat("Updated mutation log and VCF have been saved to CSV and VCF files.")


