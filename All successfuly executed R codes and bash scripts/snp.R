# Load necessary library

library(dplyr)



# Step 1: Get all CSV files in the directory (you can modify the pattern if needed)

csv_files <- list.files(pattern = "*.csv")



# Step 2: Process each CSV file

for (file in csv_files) {

  

  # Read the current SNP data CSV file

  df = read.csv(file, stringsAsFactors = FALSE)
snp_data = as.data.frame(lapply(df, function(x) gsub("\\|", "/", x)), stringsAsFactors = FALSE)
  

  # Split the data based on the combination of REF and ALT

  GC_to_AT = filter(snp_data, (REF == "G" & ALT == "A") | (REF == "C" & ALT == "T"))

  GC_to_CG = filter(snp_data, (REF == "G" & ALT == "C") | (REF == "C" & ALT == "G"))

  GC_to_TA = filter(snp_data, (REF == "G" & ALT == "T") | (REF == "C" & ALT == "A"))

  AT_to_GC = filter(snp_data, (REF == "A" & ALT == "G") | (REF == "T" & ALT == "C"))

  AT_to_CG = filter(snp_data, (REF == "A" & ALT == "C") | (REF == "T" & ALT == "G"))

  AT_to_TA = filter(snp_data, (REF == "A" & ALT == "T") | (REF == "T" & ALT == "A"))

  

  # List of datasets of different base substitution types

  dataset_names = c("GC_to_AT", "GC_to_CG", "GC_to_TA", "AT_to_GC", "AT_to_CG", "AT_to_TA")

  

  # Create an empty list to store the results

  output_dataframes = list()

  

  # Loop over each dataset

  for (dataset_name in dataset_names) {

    

    # Get the actual dataset using get() to reference the dataframe by name

    snp_data_split = get(dataset_name)

    

    # Count the total number of 1/1 mutations for each sample

    mt_counts_11 = colSums(snp_data_split[-c(1, 2, 3, 4)] == "1/1", na.rm = TRUE)

    

    # Create a summary table with sample names and the dataset name as the column header

    mt_summary_11 = data.frame(

      Sample = names(mt_counts_11),

      stringsAsFactors = FALSE

    )

    

    # Dynamically set the column name to the dataset name and assign mutation counts

    mt_summary_11[[dataset_name]] = mt_counts_11

    

    # Store the resulting dataframe in the list using the dataset name as the key

    output_dataframes[[dataset_name]] = mt_summary_11

  }

  

  # Combine the six dataframes into a single dataframe

  bs_data = Reduce(function(x, y) merge(x, y, by = "Sample", all = TRUE), output_dataframes)

  

  # Add the 8th column ("Total snp") as the sum of the values from the 2nd to 7th columns

  bs_data$`Total snp` <- rowSums(bs_data[, 2:7], na.rm = TRUE)

  

  # Generate output file name based on the input file

  output_file <- sub(".csv", "_total.csv", file)

  

  # Save the summary table to a CSV file

  write.csv(bs_data, output_file, row.names = FALSE)

  

  # Print progress

  cat("Processed:", file, "->", output_file, "\n")

}


