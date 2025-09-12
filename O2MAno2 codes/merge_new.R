## Load necessary libraries
library(dplyr)

# Function to process SNP and INDEL files and merge them
process_and_merge_files <- function(snp_file, indel_file, output_file) {
  
  # Read the SNP file
  snp_df <- read.csv(snp_file)
  
  
  # Read the INDEL file
  indel_df <- read.csv(indel_file)
  
  # Create "Deletions" by summing neg__inf_16, neg_15_2, neg_1
  indel_df <- indel_df %>%
    mutate(Deletions = neg_inf_16 + neg_15_2 + neg_1)
  
  # Create "Insertions" by summing pos_1, pos_2_15, pos_16_inf
  indel_df <- indel_df %>%
    mutate(Insertions = pos_1 + pos_2_15 + pos_16_inf)
  
  # Select relevant columns from INDEL file
  indel_df <- indel_df %>%
    select(Sample, Deletions, Insertions, Total.indels)
  
  # Merge the SNP and INDEL files based on the "Sample" column
  combined_df <- merge(snp_df, indel_df, by = "Sample")
  
  # Create a new column "Total_mutations" by summing "Total.snp" and "Total.indels"
  combined_df <- combined_df %>%
    mutate(Total_mutations = Total.snp + Total.indels)
  
  # Write the combined dataframe to the new CSV file
  write.csv(combined_df, output_file, row.names = FALSE)
}

# Function to process all pairs of SNP and INDEL files in a directory
process_all_pairs <- function(directory) {
  # Get all the files in the directory
  files <- list.files(directory, pattern = "\\.csv$", full.names = TRUE)
  
  # Separate SNP and INDEL files based on the presence of 'snp' and 'indel' in filenames
  snp_files <- files[grepl("snp", files)]
  indel_files <- files[grepl("indel", files)]
  
  # Loop over SNP files to find corresponding INDEL files
  for (snp_file in snp_files) {
    # Create the base filename by removing the 'snp' portion from the snp_file
    base_name <- sub("snp", "combined", snp_file)
    
    # Find the corresponding indel file by replacing 'snp' with 'indel'
    indel_file <- sub("snp", "indel", snp_file)
    
    # Check if the indel file exists in the list
    if (indel_file %in% indel_files) {
      # Construct output filename
      output_file <- sub("snp", "combined", snp_file)
      
      # Process and merge the SNP and INDEL files
      process_and_merge_files(snp_file, indel_file, output_file)
      
      # Print a message indicating successful merge
      cat("Processed and merged:", basename(snp_file), "and", basename(indel_file), "\n")
    }
  }
}

# Example usage:
# Set the directory where your CSV files are located
directory <- getwd()

# Process all pairs in the directory
process_all_pairs(directory)
