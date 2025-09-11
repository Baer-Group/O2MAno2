library(dplyr)



# Step 1: Load the tropicalis_mono_repeats.csv file

tropicalis_mono_repeats <- read.csv("tropicalis_mono_repeats.csv", stringsAsFactors = FALSE)



# Ensure correct data types for tropicalis_mono_repeats

tropicalis_mono_repeats$Chromosome <- as.character(tropicalis_mono_repeats$Chromosome)

tropicalis_mono_repeats$Start <- as.numeric(tropicalis_mono_repeats$Start)

tropicalis_mono_repeats$End <- as.numeric(tropicalis_mono_repeats$End)



# Step 2: Define a function to check if a mutation is in a mononucleotide repeat region

is_in_repeat <- function(chromosome, position, repeats_df) {

  # Debugging: Print the current chromosome and position being checked

  print(paste("Checking Chromosome:", chromosome, "Position:", position))

  

  # Filter rows where Chromosome matches and Position is within Start and End range

  match <- repeats_df %>%

    filter(Chromosome == chromosome & Start <= position & position <= End)

  

  # Return TRUE if any match is found, FALSE otherwise

  return(nrow(match) > 0)

}



# Step 3: Get all the mutation CSV files

mutation_files <- list.files(pattern = "*.csv") # Adjust if needed



# Step 4: Loop through each mutation file

for (file in mutation_files) {

  

  # Read the mutation file

  mutations <- read.csv(file, stringsAsFactors = FALSE)

  

  # Debugging: Print column names to check if "Position" and "Chromosome" columns exist

  print(paste("Processing file:", file))

  print(colnames(mutations))

  

  # Ensure correct data types for mutations

  if ("Position" %in% colnames(mutations) & "Chromosome" %in% colnames(mutations)) {

    mutations$Chromosome <- as.character(mutations$Chromosome)

    mutations$Position <- as.numeric(mutations$Position)

    

    # Initialize empty data frames for storing mononucleotide and non-mono mutations

    mono_mutations <- data.frame()

    non_mono_mutations <- data.frame()

    

    # Step 5: Iterate through each mutation

    for (i in 1:nrow(mutations)) {

      chromosome <- mutations$Chromosome[i]

      position <- mutations$Position[i]

      

      # Ensure chromosome and position are not missing

      if (!is.na(chromosome) & !is.na(position)) {

        

        # Debugging: Print the current row's chromosome and position values

        print(paste("Processing Chromosome:", chromosome, "Position:", position))

        

        # Check if the mutation is in a repeat region

        if (is_in_repeat(chromosome, position, tropicalis_mono_repeats)) {

          # Debugging: Print if the mutation is classified as mono

          print(paste("Position", position, "is classified as mono"))

          mono_mutations <- rbind(mono_mutations, mutations[i, ])

        } else {

          # Debugging: Print if the mutation is classified as non-mono

          print(paste("Position", position, "is classified as non-mono"))

          non_mono_mutations <- rbind(non_mono_mutations, mutations[i, ])

        }

        

      } else {

        # Debugging: Print a message if chromosome or position is missing

        print(paste("Skipping row", i, "due to missing Chromosome or Position"))

      }

    }

    

    # Step 6: Write the mono and non-mono mutations to separate files

    mono_file <- sub(".csv", "_mono.csv", file)

    non_mono_file <- sub(".csv", "_non_mono.csv", file)

    

    write.csv(mono_mutations, mono_file, row.names = FALSE)

    write.csv(non_mono_mutations, non_mono_file, row.names = FALSE)

    

    cat("Processed:", file, "->", mono_file, "and", non_mono_file, "\n")

  } else {

    cat("Error: 'Position' or 'Chromosome' column not found in file:", file, "\n")

  }

}


