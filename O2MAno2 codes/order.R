library(dplyr)



# Step 1: Get all the mutation CSV files

mutation_files <- list.files(pattern = "*.csv")



# Step 2: Separate files into E and G groups

E_files <- grep("^E", mutation_files, value = TRUE)

G_files <- grep("^G", mutation_files, value = TRUE)



# Step 3: Process each pair of twins

for (e_file in E_files) {

  

  # Find the corresponding G twin (replace the first character with 'G')

  g_file <- sub("^E", "G", e_file)

  

  # Ensure the G twin exists

  if (g_file %in% G_files) {

    

    # Read both E and G files

    e_mt_data <- read.csv(e_file, stringsAsFactors = FALSE)

    g_mt_data <- read.csv(g_file, stringsAsFactors = FALSE)

    

    # Step 4: Extract positions (assuming the second column contains positions)

    e_positions <- e_mt_data[, 2]

    g_positions <- g_mt_data[, 2]

    

    # Step 5: Find shared positions between E and G datasets

    shared_positions <- intersect(e_positions, g_positions)

    

    # Step 6: Filter matching data (O1) for both E and G datasets

    filtered_e_mt_data <- e_mt_data %>% filter(e_mt_data[, 2] %in% shared_positions)

    filtered_g_mt_data <- g_mt_data %>% filter(g_mt_data[, 2] %in% shared_positions)

    

    # Step 7: Find non-matching data (O2) for both E and G datasets

    non_matching_e_mt_data <- e_mt_data %>% filter(!(e_mt_data[, 2] %in% shared_positions))

    non_matching_g_mt_data <- g_mt_data %>% filter(!(g_mt_data[, 2] %in% shared_positions))

    

    # Step 8: Define file names for saving

    e_o1_file <- sub(".csv", "_O1.csv", e_file)

    e_o2_file <- sub(".csv", "_O2.csv", e_file)

    g_o1_file <- sub(".csv", "_O1.csv", g_file)

    g_o2_file <- sub(".csv", "_O2.csv", g_file)

    

    # Step 9: Write the filtered data to CSV files

    write.csv(filtered_e_mt_data, e_o1_file, row.names = FALSE)

    write.csv(non_matching_e_mt_data, e_o2_file, row.names = FALSE)

    write.csv(filtered_g_mt_data, g_o1_file, row.names = FALSE)

    write.csv(non_matching_g_mt_data, g_o2_file, row.names = FALSE)

    

    cat("Processed twins:", e_file, "and", g_file, "\n")

    cat("Generated:", e_o1_file, e_o2_file, g_o1_file, g_o2_file, "\n")

  }

}


