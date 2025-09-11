#Set current folder of the R script as the working directory
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
getwd()

# Install required packages if not already installed
if(!require(readr)) install.packages('readr')
if(!require(dplyr)) install.packages('dplyr')

# Load libraries
library(readr)
library(dplyr)

# Function to convert CSV to BED
convert_csv_to_bed <- function(input_csv, output_bed) {
  # Read the CSV file
  data <- read_csv(input_csv)
  
  # Create the BED columns: Start (0-based) and End (exclusive)
  bed_data <- data %>%
    mutate(Start = Position - 1,  # Convert to 0-based start
           End = Position) %>%    # End is exclusive, same as Position for single-point mutations
    select(Chromosome, Start, End)
  
  # Write the BED file (tab-separated, no header)
  write.table(bed_data, file = output_bed, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  print(paste("Converted", input_csv, "to BED format:", output_bed))
}

# Get the list of all CSV files in the current directory
csv_files <- list.files(pattern = "*.csv")  # Adjust this pattern if necessary

# Loop through each file and convert to BED
for (csv_file in csv_files) {
  # Define the output BED filename (replacing .csv with .bed)
  bed_file <- sub(".csv$", ".bed", csv_file)
  
  # Call the function to convert CSV to BED
  convert_csv_to_bed(csv_file, bed_file)
}
# Function to convert reference CSV (with Start and End) to BED
convert_range_csv_to_bed <- function(input_csv, output_bed) {
  # Read the CSV file
  data <- read_csv(input_csv)
  
  # Ensure it has Chromosome, Start, and End columns
  bed_data <- data %>%
    select(Chromosome, Start, End)
  
  # Write to BED file (no changes to Start and End since they already define ranges)
  write.table(bed_data, file = output_bed, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  print(paste("Converted", input_csv, "to BED format:", output_bed))
}

# Example usage for reference file
convert_range_csv_to_bed("briggsae_mono_repeats.csv", "briggsae_mono_repeats.bed")




# Function to convert BED back to CSV
convert_bed_to_csv <- function(input_bed, output_csv) {
  # Read the BED file
  bed_data <- read.table(input_bed, sep = "\t", header = FALSE)
  colnames(bed_data) <- c("Chromosome", "Start", "End")
  
  # Convert Start and End to Position (if necessary, depending on your use case)
  bed_data$Position <- bed_data$End
  
  # Save as CSV
  write.csv(bed_data, output_csv, row.names = FALSE)
}

# Example usage to convert back to CSV
convert_bed_to_csv("briggsae_mono.bed", "briggsae_mono.csv")
convert_bed_to_csv("briggsae_non_mono.bed", "briggsae_non_mono.csv")
