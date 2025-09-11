library(dplyr)



# Step 1: Load the mono_repeats.csv file

mono_repeats <- read.csv("mono_repeats.csv", stringsAsFactors = FALSE)



# Step 2: Define a function to check if a mutation is in a mononucleotide repeat region

is_in_repeat <- function(chromosome, position, repeats_df) {

  # Filter rows where Chromosome matches and Position is within Start and End range

  match <- repeats_df %>%

    filter(Chromosome == chromosome & Start <= position & position <= End)

  

  # Return TRUE if any match is found, FALSE otherwise

  return(nrow(match) > 0)

}



# Step 3: Get all the mutation CSV files

mutation_files <- list.files(pattern = "*x.csv") # Adjust if eeeded



# Step 4: Loop through each mutation file

for (file in mutation_files) {

  

  # Read the mutation file

  mutations <- read.csv(file, stringsAsFactors = FALSE)

  

  # Initialize empty data frames for storing mononucleotide and non-mono mutations

  mono_mutations <- data.frame()

  non_mono_mutations <- data.frame()

  

  # Step 5: Iterate through each mutation

  for (i in 1:nrow(mutations)) {

    chromosome <- mutations$Chromosome[i]

    position <- mutations$Position[i]

    

    # Check if the mutation is in a repeat region

    if (is_in_repeat(chromosome, position, mono_repeats)) {

      mono_mutations <- rbind(mono_mutations, mutations[i, ])

    } else {

      non_mono_mutations <- rbind(non_mono_mutations, mutations[i, ])

    }

  }

  

  # Step 6: Write the mono and non-mono mutations to separate files

  mono_file <- sub(".csv", ".mono.csv", file)

  non_mono_file <- sub(".csv", ".non_mono.csv", file)

  

  write.csv(mono_mutations, mono_file, row.names = FALSE)

  write.csv(non_mono_mutations, non_mono_file, row.names = FALSE)

  

  cat("Processed:", file, "->", mono_file, "and", non_mono_file, "\n")

}


