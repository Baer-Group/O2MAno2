# Load the necessary library

library(dplyr)



# Step 1: Get all CSV files in the directory (you can modify the pattern if needed)

csv_files <- list.files(pattern = "*.csv")



# Step 2: Loop through each CSV file and process the indel data

for (file in csv_files) {

  

  # Read the current indel data CSV file

  df <- read.csv(file, stringsAsFactors = FALSE)

indel_data = as.data.frame(lapply(df, function(x) gsub("\\|", "/", x)), stringsAsFactors = FALSE)  

  # List of indel categories based on length ranges

  indel_categories <- list(

    "neg_20_11" = filter(indel_data, ALT != REF & (nchar(REF) - nchar(ALT)) >= 11 & (nchar(REF) - nchar(ALT)) <= 20),

    "neg_10_6" = filter(indel_data, ALT != REF & (nchar(REF) - nchar(ALT)) >= 6 & (nchar(REF) - nchar(ALT)) <= 10),

    "neg_5_2" = filter(indel_data, ALT != REF & (nchar(REF) - nchar(ALT)) >= 2 & (nchar(REF) - nchar(ALT)) <= 5),

    "neg_1" = filter(indel_data, ALT != REF & (nchar(REF) - nchar(ALT)) == 1),

    "pos_1" = filter(indel_data, ALT != REF & (nchar(ALT) - nchar(REF)) == 1),

    "pos_2_5" = filter(indel_data, ALT != REF & (nchar(ALT) - nchar(REF)) >= 2 & (nchar(ALT) - nchar(REF)) <= 5),

    "pos_6_10" = filter(indel_data, ALT != REF & (nchar(ALT) - nchar(REF)) >= 6 & (nchar(ALT) - nchar(REF)) <= 10),

    "pos_11_20" = filter(indel_data, ALT != REF & (nchar(ALT) - nchar(REF)) >= 11 & (nchar(ALT) - nchar(REF)) <= 20)

  )

  

  # Create an empty list to store the results for each category

  output_indel_dataframes <- list()



  # Step 3: Loop over each category in the list

  for (indel_name in names(indel_categories)) {

    

    # Get the actual dataset for the current category

    indel_data_filtered <- indel_categories[[indel_name]]

    

    # Step 4: Count the total number of 1/1 mutations for each sample

    indel_counts_11 <- colSums(indel_data_filtered[-c(1, 2, 3, 4)] == "1/1", na.rm = TRUE)

    

    # Step 5: Create a summary table with sample names and the indel category as the column header

    indel_summary_11 <- data.frame(

      Sample = names(indel_counts_11),

      stringsAsFactors = FALSE

    )

    

    # Step 6: Dynamically set the column name to the indel category and assign mutation counts

    indel_summary_11[[indel_name]] <- indel_counts_11

    

    # Step 7: Store the resulting dataframe in the list using the indel category as the key

    output_indel_dataframes[[indel_name]] <- indel_summary_11

  }



  # Step 8: Combine the indel dataframes into a single dataframe

  indel_data_combined <- Reduce(function(x, y) merge(x, y, by = "Sample", all = TRUE), output_indel_dataframes)



  # Step 9: Create a total count column for all indel categories

  indel_data_combined$`Total indels` <- rowSums(indel_data_combined[, 2:9], na.rm = TRUE)



  # Generate output file name based on the input file

  output_file <- sub(".csv", "_total.csv", file)



  # Step 10: Save the indel summary table to a CSV file

  write.csv(indel_data_combined, output_file, row.names = FALSE)

  

  # Print progress

  cat("Processed:", file, "->", output_file, "\n")

}


