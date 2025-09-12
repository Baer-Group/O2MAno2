# Load necessary library
library(dplyr)

# Define the folder path containing the CSV files
folder_path <- getwd()

# Get a list of all files with E and G prefixes
e_files <- list.files(path = folder_path, pattern = "^E.*\\.csv$", full.names = TRUE)
g_files <- list.files(path = folder_path, pattern = "^G.*\\.csv$", full.names = TRUE)

# Function to extract the part of the filename without the 'E' or 'G' prefix and without the redundant part
strip_filename <- function(filename) {
  sub("^[EG]\\.combined\\.", "", basename(filename))
}

# Loop through each E file
for (e_file in e_files) {
  # Strip the prefix to find the corresponding G file
  base_name <- strip_filename(e_file)
  g_file <- file.path(folder_path, paste0("G.combined.", base_name))
  
  if (file.exists(g_file)) {
    # Read both E and G files
    df_e <- read.csv(e_file)
    df_g <- read.csv(g_file)
    
    # Combine the rows from E and G
    combined_df <- bind_rows(df_e, df_g)
    
    # Define the new filename by stripping redundant parts
    new_file_name <- paste0(folder_path, "/", base_name)
    
    # Write the combined file
    write.csv(combined_df, new_file_name, row.names = FALSE)
    
    print(paste("Combined", basename(e_file), "and", basename(g_file), "into", basename(new_file_name)))
  } else {
    print(paste("No matching G file found for", basename(e_file)))
  }
}

print("All files processed and combined.")


# Ensure the folder path ends with a slash to correctly form the file path
if (!grepl("/$", folder_path)) {
  folder_path <- paste0(folder_path, "/")
}

# Get a list of the 12 combined files
file_list <- list.files(path = folder_path, pattern = ".*\\.csv$", full.names = TRUE)

# Print the folder path to verify
print(paste("Folder path:", folder_path))

# Function to extract the line number from the sample name
extract_line <- function(sample_name) {
  line <- NA
  
  # Search for E or G in the sample name and determine the line number
  if (grepl("E[0-9]{3}", sample_name)) {
    line_number <- as.numeric(sub(".*E([0-9]{3}).*", "\\1", sample_name))
    line <- line_number
  } else if (grepl("G[0-9]{3}", sample_name)) {
    line_number <- as.numeric(sub(".*G([0-9]{3}).*", "\\1", sample_name))
    line <- line_number + 200
  }
  
  return(line)
}

# Loop through each file
for (file in file_list) {
  # Print the file being processed
  print(paste("Processing file:", file))
  
  # Read the CSV file
  df <- read.csv(file)
  
  # Extract Line based on the Sample name
  df <- df %>%
    rowwise() %>%
    mutate(Line = extract_line(Sample)) %>%
    ungroup()
  
  # Assign Replicate values based on row numbers using row_number()
  df <- df %>%
    mutate(Replicate = ifelse(row_number() <= 100, 1, 2))
  
  # Relocate Line and Replicate to be the 2nd and 3rd columns
  df <- df %>%
    relocate(Line, Replicate, .after = Sample)
  
  # Define the new file name with "_new.csv" suffix
  new_file <- sub("\\_cov.csv$", ".csv", file)
  
  # Print the new file path to confirm it's being written
  print(paste("Writing new file:", new_file))
  
  # Write the updated CSV file
  write.csv(df, new_file, row.names = FALSE)
}

print("All files processed with Line and Replicate columns added.")
