# Load necessary libraries
library(dplyr)
library(tidyr)

# Step 1: Load the coverage data
coverage_data <- read.csv("E3X_coverage_results.csv")


# Step 2: Define the function to process each *total.csv file

process_file <- function(file_name, coverage_data) {

  
  # Step 3: Load the *total.csv data

  data <- read.csv(file_name)

  

  # Step 4: Merge the coverage data with the current *total.csv file based on the Sample column

  # Assuming that both files have the same "Sample" column to match on

  merged_data <- data %>% 

    left_join(coverage_data, by = "Sample")

 

  # Step 5: Remove specific rows (H7GWFDSX7_EG8072_BaerAncestor_POOLRET74_a/b/c)

  rows_to_delete <- c("H7GWFDSX7_EG8072_BaerAncestor_POOLRET74_a", 

                      "H7GWFDSX7_EG8072_BaerAncestor_POOLRET74_b", 

                      "H7GWFDSX7_EG8072_BaerAncestor_POOLRET74_c")

  

  # Delete rows where any column matches the specified rows

  merged_data <- merged_data[!apply(merged_data, 1, function(row) any(row %in% rows_to_delete)), ]

  

  # Step 6: Write the new file

  # Generate new file name by replacing "total" with "rate"

  new_file_name <- gsub("total", "cov", file_name)

  

  # Save the resulting data as a new CSV file

  write.csv(merged_data, new_file_name, row.names = FALSE)

  

  # Return the new file name for confirmation

  return(new_file_name)

}



# Step 7: Get a list of all files that match the pattern "*total.csv"

file_list <- list.files(pattern = "*total.csv")



# Step 8: Apply the process_file function to all files, using the coverage data

processed_files <- lapply(file_list, process_file, coverage_data = coverage_data)



# Print the names of the generated files

print(processed_files)


