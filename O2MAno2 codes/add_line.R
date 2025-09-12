# Set current folder of the R script as the working directory
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
getwd()

# Load necessary libraries
library(dplyr)

# Define a function to transform the dataset
generate_mutation_line_file <- function(input_file, gene_file, output_file) {
  # Read the input CSV file
  data <- read.csv(input_file, header = TRUE)
  
  # Read the gene annotation file
  gene_data <- read.csv(gene_file, header = TRUE)
  
  # Extract constant columns (Chromosome, Position, REF, ALT)
  constant_cols <- data[, 1:4]
  
  # Get the line names from the columns
  line_names <- colnames(data)[5:ncol(data)]
  
  # Create an empty data frame to store the results
  mutation_results <- data.frame(Chromosome = character(),
                                 Position = numeric(),
                                 REF = character(),
                                 ALT = character(),
                                 Line = character(),
                                 Gene_name = character(),
                                 Putative_impact = character(),
                                 Annotation = character(),
                                 stringsAsFactors = FALSE)
  
  # Iterate over each row of the dataset
  for (i in 1:nrow(data)) {
    # Extract row data (excluding the constant columns)
    row_data <- data[i, 5:ncol(data)]
    
    # Find columns with "1/1" or "1|1" mutation
    mutated_lines <- which(row_data == "1/1" | row_data == "1|1")
    
    # If there are any mutations, create new rows for each mutation found
    if (length(mutated_lines) > 0) {
      for (line in mutated_lines) {
        # Extract gene information for the current row
        gene_info <- gene_data[i, c("Gene_name", "Putative_impact","Annotation")]
        
        # Create a new row with gene information
        new_row <- c(as.character(constant_cols[i, ]), line_names[line], as.character(gene_info))
        mutation_results <- rbind(mutation_results, new_row)
      }
    }
  }
  
  # Rename columns
  colnames(mutation_results) <- c("Chromosome", "Position", "REF", "ALT", "Line", "Gene_name", "Putative_impact", "Annotation")
  
  # Write the output CSV file
  write.csv(mutation_results, output_file, row.names = FALSE)
}

# Apply the function to the four datasets
input_files <- c("E.bp.GVCFs.snp.10x.mutation.rows.csv", "G.bp.GVCFs.snp.10x.mutation.rows.csv", 
                 "E.bp.GVCFs.indel.10x.mutation.rows.csv", "G.bp.GVCFs.indel.10x.mutation.rows.csv")
gene_files <- c("genes_from_snp_10xE.csv", "genes_from_snp_10xG.csv", 
                "genes_from_indel_10xE.csv", "genes_from_indel_10xG.csv")
output_files <- c("E.snp.10x.mutation.lines.csv", "G.snp.10x.mutation.lines.csv", 
                  "E.indel.10x.mutation.lines.csv", "G.indel.10x.mutation.lines.csv")

for (i in 1:length(input_files)) {
  generate_mutation_line_file(input_files[i], gene_files[i], output_files[i])
}

