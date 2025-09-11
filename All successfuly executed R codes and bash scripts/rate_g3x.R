# Load necessary libraries

library(dplyr)

library(tidyr)



# Step 1: Load the coverage data

coverage_data <- read.csv("G3x_coverage_results.csv")



# Step 2: Define the function to process each *total.csv file

process_file <- function(file_name, coverage_data) {

  

  # Step 3: Load the *total.csv data

  data <- read.csv(file_name)

  

  # Step 4: Merge the coverage data with the current *total.csv file based on the Sample column

  # Assuming that both files have the same "Sample" column to match on

  merged_data <- data %>% 

    left_join(coverage_data, by = "Sample")

  

  # Step 5: Replace zeros with the pseudocount (0.00000000000001)

  merged_data[merged_data == 0] <- 0.00000000000001

  

  # Step 6: Apply the formula for each data point: value / (Coverage * 1.5 * 100286401)

  genome_length <- 100286401

  merged_data <- merged_data %>% 

    mutate(across(-c(Sample, Coverage), ~ . / (Coverage * 1.5 * genome_length)))

  

  # Step 7: Remove specific rows (H7GWFDSX7_EG8072_BaerAncestor_POOLRET74_a/b/c)

  rows_to_delete <- c("H7GWFDSX7_EG8072_BaerAncestor_POOLRET74_a", 

                      "H7GWFDSX7_EG8072_BaerAncestor_POOLRET74_b", 

                      "H7GWFDSX7_EG8072_BaerAncestor_POOLRET74_c")

  

  # Delete rows where any column matches the specified rows

  merged_data <- merged_data[!apply(merged_data, 1, function(row) any(row %in% rows_to_delete)), ]

  

  # Step 8: Write the new file

  # Generate new file name by replacing "total" with "rate"

  new_file_name <- gsub("total", "rate", file_name)

  

  # Save the resulting data as a new CSV file

  write.csv(merged_data, new_file_name, row.names = FALSE)

  

  # Return the new file name for confirmation

  return(new_file_name)

}



# Step 9: Get a list of all files that match the pattern "*total.csv"

file_list <- list.files(pattern = "*total.csv")



# Step 10: Apply the process_file function to all files, using the coverage data

processed_files <- lapply(file_list, process_file, coverage_data = coverage_data)



# Print the names of the generated files

print(processed_files)


