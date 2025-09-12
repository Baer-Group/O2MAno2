# Set current folder of the R script as the working directory
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
getwd()

# Load necessary library
library(dplyr)

# Step 1: Get all CSV files in the directory
csv_files <- list.files(pattern = "*.csv")

# Step 2: Process each CSV file
for (file in csv_files) {
  
  # Read the current SNP data CSV file
  df = read.csv(file, stringsAsFactors = FALSE)
  df2 = as.data.frame(lapply(df, function(x) gsub("\\|", "/", x)), stringsAsFactors = FALSE)
  
  # Define heterozygote genotypes
  ht = c("0/1", "1/0")
  
  # Identify the three ancestors
  anc = c("H7GWFDSX7_EG8072_BaerAncestor_POOLRET74_a",
          "H7GWFDSX7_EG8072_BaerAncestor_POOLRET74_b",
          "H7GWFDSX7_EG8072_BaerAncestor_POOLRET74_c")
  
  # Remove rows where the ancestors contain heterozygous values or NA
  df3 = df2 %>%
    filter(!apply(.[, anc], 1, function(x) any(x %in% ht) || any(is.na(x))))
  
  # Step 1: Keep rows where there's one 0/0 or 0/1, with rest being 1/1 or NA
  keep_rows = function(row) {
    gt = as.character(row[-c(1, 2, 3, 4)])  # Exclude first four columns
    count_0_0 = sum(gt == "0/0", na.rm = TRUE)
    count_0_1 = sum(gt == "0/1", na.rm = TRUE)
    count_1_1 = sum(gt == "1/1", na.rm = TRUE)
    total_na = sum(is.na(gt))
    
    # Keep if there's exactly one 0/0 or 0/1 and the rest are 1/1 or NA
    if ((count_0_0 == 1 && count_1_1 + total_na == (length(gt) - 1)) ||
        (count_0_1 == 1 && count_1_1 + total_na == (length(gt) - 1))) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
  
  # Apply the function to each row and filter the rows with unique genotypes
  no_recall_rows = df3[apply(df3, 1, keep_rows), ]
  
  # Save the filtered SNP data to *no_recall.csv files
  output_no_recall <- sub(".csv", "_no_recall.csv", file)  # Generate output file name
  write.csv(no_recall_rows, output_no_recall, row.names = FALSE)
  
  # Step 2: Keep rows with less than 10 NA and save to *strict.csv
  strict_rows <- no_recall_rows %>%
    filter(rowSums(is.na(.)) < 10)
  
  output_strict <- sub("_no_recall.csv", "_strict.csv", output_no_recall)  # Generate output file name
  write.csv(strict_rows, output_strict, row.names = FALSE)
  
  # Step 3: Keep rows with no NA and only one 0/0 or 0/1, saving to *zero_na.csv
  zero_na_rows = strict_rows %>%
    filter(rowSums(is.na(.)) == 0) %>%
    filter(apply(., 1, keep_rows))  # Reuse the keep_rows function
  
  output_zero_na <- sub("_strict.csv", "_zero_na.csv", output_strict)  # Generate output file name
  write.csv(zero_na_rows, output_zero_na, row.names = FALSE)
  
  # Print progress
  cat("Processed:", file, "->", output_no_recall, "\n")
  cat("Processed:", file, "->", output_strict, "\n")
  cat("Processed:", file, "->", output_zero_na, "\n")
}
