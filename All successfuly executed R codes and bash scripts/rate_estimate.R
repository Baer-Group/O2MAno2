#Set current folder of the R script as the working directory
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
getwd()
# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)
library(gridExtra) # For arranging multiple plots

# Step 1: Define the function to load and process each file for GLM
process_rate_file <- function(file_name, start_col, end_col, custom_order, rename_labels) {
  
  # Load the renamed CSV file
  data <- read.csv(file_name)
  
  # Select only the relevant columns (2nd to 7th for SNP, 2nd to 9th for Indel)
  mutation_data <- data[, start_col:end_col]
  
  # Rename the columns based on rename_labels
  colnames(mutation_data) <- rename_labels
  
  # Convert the data into long format for GLM
  long_data <- pivot_longer(mutation_data, cols = everything(), names_to = "Category", values_to = "MutationRate")
  
  # Ensure that the factor levels are in the correct order based on the custom_order
  long_data$Category <- factor(long_data$Category, levels = custom_order)
  
  # Step 2: Fit a General Linear Model (GLM) for Mutation Rate estimation
  glm_model <- glm(MutationRate ~ Category, data = long_data, family = gaussian())
  
  # Step 3: Calculate mean and standard error for each category
  glm_summary <- long_data %>%
    group_by(Category) %>%
    summarise(
      mean_rate = mean(MutationRate),
      se_rate = sd(MutationRate) / sqrt(n())
    )
  
  # Return the summarized data and the long data (for plotting)
  return(glm_summary)
}

# Step 4: Define the function to categorize Dataset into Shared O1 Line, E Line, and G Line
categorize_dataset <- function(file_name) {
  if (grepl("O1", basename(file_name))) {
    return("Shared O1 Line")
  } else if (grepl("^E", basename(file_name))) {
    return("E Line")
  } else if (grepl("^G", basename(file_name))) {
    return("G Line")
  } else {
    return("Unknown")
  }
}


# Step 5: Define the function to compare three datasets and generate the plot
compare_three_datasets <- function(file1, file2, file3, start_col, end_col, y_label = "Mutation rate", file_type = "snp", custom_order, rename_labels, plot_title) {
  
  # Load and process each file
  data1 <- process_rate_file(file1, start_col, end_col, custom_order, rename_labels) %>% mutate(Dataset = categorize_dataset(file1))
  data2 <- process_rate_file(file2, start_col, end_col, custom_order, rename_labels) %>% mutate(Dataset = categorize_dataset(file2))
  data3 <- process_rate_file(file3, start_col, end_col, custom_order, rename_labels) %>% mutate(Dataset = categorize_dataset(file3))
  
  # Combine the processed data from all three files
  combined_data <- bind_rows(data1, data2, data3)
  
  # Plot the results with error bars and grouped by Dataset
  plot <- ggplot(combined_data, aes(x = Category, y = mean_rate, fill = Dataset)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(aes(ymin = mean_rate - se_rate, ymax = mean_rate + se_rate), width = 0.2, position = position_dodge(0.9)) +
    labs(x = "Category", y = y_label, title = plot_title, fill = "Strains Compared") +  # Change the legend title to "Strains Compared"
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Tilt the x-axis labels
  
  return(plot)
}

# Custom order and renaming for SNPs
snp_custom_order <- c("G:C→A:T","G:C→C:G","G:C→T:A", "A:T→G:C",   "A:T→C:G", "A:T→T:A")
snp_rename_labels <- c("G:C→A:T","G:C→C:G","G:C→T:A", "A:T→G:C",   "A:T→C:G", "A:T→T:A")

# Custom order and renaming for Indels
indel_custom_order <- c("-20 to -11", "-10 to -6", "-5 to -2", "-1", "1", "2 to 5", "6 to 10", "11 to 20")
indel_rename_labels <- c("-20 to -11", "-10 to -6", "-5 to -2", "-1", "1", "2 to 5", "6 to 10", "11 to 20")

# # For 3x coverage
plot1 <- compare_three_datasets("E.bp.GVCFs.snp.3x.mutation.rows_O1_rate.csv", "E.bp.GVCFs.snp.3x.mutation.rows_O2_rate.csv", "G.bp.GVCFs.snp.3x.mutation.rows_O2_rate.csv",
                                start_col = 2, end_col = 7, y_label = "SNP rate", 
                                file_type = "snp", custom_order = snp_custom_order, rename_labels = snp_rename_labels, "Whole Genome")


plot2 <- compare_three_datasets("E.bp.GVCFs.indel.3x.mutation.rows_O1_rate.csv", "E.bp.GVCFs.indel.3x.mutation.rows_O2_rate.csv", "G.bp.GVCFs.indel.3x.mutation.rows_O2_rate.csv",
                                start_col = 2, end_col = 9, y_label = "Indel rate", 
                                file_type = "indel", custom_order = indel_custom_order, rename_labels = indel_rename_labels, "Whole Genome")

plot3 <- compare_three_datasets("E.bp.GVCFs.snp.3x.mutation.rows_non_mono_O1_rate.csv", "E.bp.GVCFs.snp.3x.mutation.rows_non_mono_O2_rate.csv", "G.bp.GVCFs.snp.3x.mutation.rows_non_mono_O2_rate.csv",
                                start_col = 2, end_col = 7, y_label = "SNP rate", 
                                file_type = "snp", custom_order = snp_custom_order, rename_labels = snp_rename_labels, "Non Mono")

plot4 <- compare_three_datasets("E.bp.GVCFs.indel.3x.mutation.rows_non_mono_O1_rate.csv", "E.bp.GVCFs.indel.3x.mutation.rows_non_mono_O2_rate.csv", "G.bp.GVCFs.indel.3x.mutation.rows_non_mono_O2_rate.csv", 
                                start_col = 2, end_col = 9, y_label = "Indel rate", 
                                file_type = "indel", custom_order = indel_custom_order, rename_labels = indel_rename_labels, "Non Mono")

plot5 <- compare_three_datasets("E.bp.GVCFs.snp.3x.mutation.rows_mono_O1_rate.csv", "E.bp.GVCFs.snp.3x.mutation.rows_mono_O2_rate.csv", "G.bp.GVCFs.snp.3x.mutation.rows_mono_O2_rate.csv",
                                start_col = 2, end_col = 7, y_label = "SNP rate", 
                                file_type = "snp", custom_order = snp_custom_order, rename_labels = snp_rename_labels, "Mono")

plot6 <- compare_three_datasets("E.bp.GVCFs.indel.3x.mutation.rows_mono_O1_rate.csv", "E.bp.GVCFs.indel.3x.mutation.rows_mono_O2_rate.csv", "G.bp.GVCFs.indel.3x.mutation.rows_mono_O2_rate.csv", 
                                start_col = 2, end_col = 9, y_label = "Indel rate", 
                                file_type = "indel", custom_order = indel_custom_order, rename_labels = indel_rename_labels, "Mono")



# Arrange the plots into a 3x2 grid
compound_plot <- grid.arrange(plot1, plot2, plot3, plot4, plot5, plot6, nrow = 3, ncol = 2, top = "Mutation rate and Spectrum at 10X coverage")



# For 10x coverage
plot1 <- compare_three_datasets("E.bp.GVCFs.snp.10x.mutation.rows_O1_rate.csv", "E.bp.GVCFs.snp.10x.mutation.rows_O2_rate.csv", "G.bp.GVCFs.snp.10x.mutation.rows_O2_rate.csv",
                                start_col = 2, end_col = 7, y_label = "SNP rate", 
                                file_type = "snp", custom_order = snp_custom_order, rename_labels = snp_rename_labels, "Whole Genome")


plot2 <- compare_three_datasets("E.bp.GVCFs.indel.10x.mutation.rows_O1_rate.csv", "E.bp.GVCFs.indel.10x.mutation.rows_O2_rate.csv", "G.bp.GVCFs.indel.10x.mutation.rows_O2_rate.csv",
                                start_col = 2, end_col = 9, y_label = "Indel rate", 
                                file_type = "indel", custom_order = indel_custom_order, rename_labels = indel_rename_labels, "Whole Genome")

plot3 <- compare_three_datasets("E.bp.GVCFs.snp.10x.mutation.rows_non_mono_O1_rate.csv", "E.bp.GVCFs.snp.10x.mutation.rows_non_mono_O2_rate.csv", "G.bp.GVCFs.snp.10x.mutation.rows_non_mono_O2_rate.csv",
                                start_col = 2, end_col = 7, y_label = "SNP rate", 
                                file_type = "snp", custom_order = snp_custom_order, rename_labels = snp_rename_labels, "Non Mono")

plot4 <- compare_three_datasets("E.bp.GVCFs.indel.10x.mutation.rows_non_mono_O1_rate.csv", "E.bp.GVCFs.indel.10x.mutation.rows_non_mono_O2_rate.csv", "G.bp.GVCFs.indel.10x.mutation.rows_non_mono_O2_rate.csv", 
                                start_col = 2, end_col = 9, y_label = "Indel rate", 
                                file_type = "indel", custom_order = indel_custom_order, rename_labels = indel_rename_labels, "Non Mono")

plot5 <- compare_three_datasets("E.bp.GVCFs.snp.10x.mutation.rows_mono_O1_rate.csv", "E.bp.GVCFs.snp.10x.mutation.rows_mono_O2_rate.csv", "G.bp.GVCFs.snp.10x.mutation.rows_mono_O2_rate.csv",
                                start_col = 2, end_col = 7, y_label = "SNP rate", 
                                file_type = "snp", custom_order = snp_custom_order, rename_labels = snp_rename_labels, "Mono")

plot6 <- compare_three_datasets("E.bp.GVCFs.indel.10x.mutation.rows_mono_O1_rate.csv", "E.bp.GVCFs.indel.10x.mutation.rows_mono_O2_rate.csv", "G.bp.GVCFs.indel.10x.mutation.rows_mono_O2_rate.csv", 
                                start_col = 2, end_col = 9, y_label = "Indel rate", 
                                file_type = "indel", custom_order = indel_custom_order, rename_labels = indel_rename_labels, "Mono")

# Define the file names for each plot
# E.bp.GVCFs.indel.3x.mutation.rows_mono_O1_rate.csv
