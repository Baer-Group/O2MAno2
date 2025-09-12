library(dplyr)

# Define the path to your TRF .dat file
dat_file <- "c_elegans.PRJNA13758.WS292.genomic.fa.2.5.5.75.10.0.2.dat"

# Read the file into a list of lines
lines <- readLines(dat_file)

# Initialize an empty list to store repeat information
repeat_data_list <- list()
current_chromosome <- NULL  # To store the current chromosome or sequence

# Loop through each line in the file and extract relevant data
for (i in seq_along(lines)) {
  line <- lines[i]
  
  # Check if the line starts with 'Sequence': this typically contains chromosome information
  if (grepl("^Sequence", line)) {
    # Extract the chromosome name (everything after "Sequence: ")
    current_chromosome <- sub("^Sequence: ", "", line)
  }
  
  # Check if the line starts with a digit (indicating a repeat data line)
  if (grepl("^[0-9]", line)) {
    # Split the line into fields
    fields <- unlist(strsplit(line, "\\s+"))
    
    # Extract relevant fields from the line
    start <- as.integer(fields[1])
    end <- as.integer(fields[2])
    period_size <- as.integer(fields[3])
    copy_number <- as.numeric(fields[4])
    percent_matches <- as.numeric(fields[6])
    repeat_seq <- fields[length(fields)]
    
    # Append the data as a list (faster than rbind)
    repeat_data_list[[length(repeat_data_list) + 1]] <- list(
      Chromosome = current_chromosome,  # Associate chromosome name with repeat
      Start = start,
      End = end,
      Period_Size = period_size,
      Copy_Number = copy_number,
      Percent_Matches = percent_matches,
      Repeat_Sequence = repeat_seq
    )
  }
}

# Convert the list to a dataframe
repeat_data <- do.call(rbind, lapply(repeat_data_list, as.data.frame, stringsAsFactors = FALSE))

# Display the dataframe structure
str(repeat_data)

# Extend the Start and End positions by 2 base pairs
repeat_data <- repeat_data %>%
  mutate(Start = as.numeric(Start) - 2,
         End = as.numeric(End) + 2)

# Remove rows where Copy_Number is less than 5
repeats_filtered <- repeat_data[repeat_data$Copy_Number >= 5, ]

# Save the extended repeat regions to a CSV file
write.csv(repeats_filtered, "mono_di_repeats.csv", row.names = FALSE)

### Filter for mononucleotide repeats (Period_Size == 1)
mono_repeats <- repeats_filtered %>%
  filter(Period_Size == 1)

# Output filtered mono repeat data
write.csv(mono_repeats, "mono_repeats.csv", row.names = FALSE)

### Filter for mononucleotide repeats (Period_Size == 1)
di_repeats <- repeats_filtered %>%
  filter(Period_Size == 2)

# Output filtered mono repeat data
write.csv(di_repeats, "di_repeats.csv", row.names = FALSE)

# Summing the values in the Copy_Number column
total_sum <- sum(mono_repeats$Copy_Number)
