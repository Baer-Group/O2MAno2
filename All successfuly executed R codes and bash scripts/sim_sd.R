# Set current folder of the R script as the working directory
library(rstudioapi)
library(ggplot2)
setwd(dirname(getActiveDocumentContext()$path))
getwd()

# Details of 'variance in mutation rate' simulation
## This script simulates the second order MA experiment data assuming
## the null hypothesis of one mutation rate is true
## It incorporates sampling variance and unequal coverage

# Input files and constants
dataset_o1 <- read.csv("3x.cleaned_O1.csv")
dataset_o2 <- read.csv("3x.cleaned_O2.csv")
ref <- 100286401

# Define column indices
column_generations <- 12
column_coverage  <- 11
column_snps      <- 10
column_indels    <- 16
column_mutations <- 17

# Number of lines (data points)
N <- 181

# Function to simulate mean, variance, and SD for a given dataset and column index
simulate_stats <- function(dataset, column_index) {
  # maximum likelihood estimate of rate per genome per generation
  rate_mu <- mean(
    dataset[1:N, column_index] /
      (dataset[1:N, column_coverage] * dataset[1:N, column_generations])
  )
  
  Mean_sim <- numeric(10000)
  Var_sim  <- numeric(10000)
  SD_sim   <- numeric(10000)
  
  for (k in seq_len(10000)) {
    rates <- numeric(N)
    for (i in seq_len(N)) {
      # simulate raw counts over generations
      raw <- sum(rpois(dataset[i, column_generations], rate_mu))
      # apply coverage sampling: keep fraction that falls in covered regions
      kept <- sum(runif(raw) < dataset[i, column_coverage])
      # convert to rate per genome per generation
      rates[i] <- kept / (dataset[i, column_coverage] * dataset[i, column_generations] * ref)
    }
    Mean_sim[k] <- mean(rates)
    Var_sim[k]  <- var(rates)
    SD_sim[k]   <- sd(rates)
  }
  
  list(mean = Mean_sim, variance = Var_sim, sd = SD_sim)
}

# Run simulations for O1 and O2 datasets
o1_mutations <- simulate_stats(dataset_o1, column_mutations)
o1_indels    <- simulate_stats(dataset_o1, column_indels)
o1_snps      <- simulate_stats(dataset_o1, column_snps)

o2_mutations <- simulate_stats(dataset_o2, column_mutations)
o2_indels    <- simulate_stats(dataset_o2, column_indels)
o2_snps      <- simulate_stats(dataset_o2, column_snps)

# Save all 10,000 simulation results to CSV
i_df <- data.frame(
  o1_mut_mean = o1_mutations$mean,
  o1_mut_var  = o1_mutations$variance,
  o1_mut_sd   = o1_mutations$sd,
  o1_ind_mean = o1_indels$mean,
  o1_ind_var  = o1_indels$variance,
  o1_ind_sd   = o1_indels$sd,
  o1_snp_mean = o1_snps$mean,
  o1_snp_var  = o1_snps$variance,
  o1_snp_sd   = o1_snps$sd,
  o2_mut_mean = o2_mutations$mean,
  o2_mut_var  = o2_mutations$variance,
  o2_mut_sd   = o2_mutations$sd,
  o2_ind_mean = o2_indels$mean,
  o2_ind_var  = o2_indels$variance,
  o2_ind_sd   = o2_indels$sd,
  o2_snp_mean = o2_snps$mean,
  o2_snp_var  = o2_snps$variance,
  o2_snp_sd   = o2_snps$sd
)
write.csv(i_df, "mutation_simulation_results.csv", row.names = FALSE)

options(scipen = 0)
####################SNP#######################

# Calculate observed SD of rates from original data
obs_o1_snp <- sd(
  dataset_o1[1:N, column_snps] /
    (dataset_o1[1:N, column_coverage] * dataset_o1[1:N, column_generations] * ref)
)
obs_o2_snp <- sd(
  dataset_o2[1:N, column_snps] /
    (dataset_o2[1:N, column_coverage] * dataset_o2[1:N, column_generations] * ref)
)

# Compute p-values: fraction of simulated SDs greater than observed
p_o1 <- mean(o1_snps$sd > obs_o1_snp)
p_o2 <- mean(o2_snps$sd > obs_o2_snp)

# Plot distribution of simulated SD for SNP rate comparison
plot_data <- data.frame(
  sd  = c(o1_snps$sd, o2_snps$sd),
  set = rep(c("O1", "O2"), each = length(o1_snps$sd))
)

p <- ggplot(plot_data, aes(x = sd, fill = set)) +
  geom_histogram(position = "identity", alpha = 0.6, bins = 50) +
  # observed lines and labels
  geom_vline(xintercept = obs_o1_snp, color = "red",   linetype = "dashed", size = 1) +
  annotate("text", x = obs_o1_snp, y = Inf, label = "Observed O1", color = "red", angle = 90, vjust = -0.5,hjust = 1.1) +
  annotate("text", x = obs_o1_snp, y = 0,   label = paste0("pO1=", signif(p_o1,3)), color = "black", angle = 0, vjust = 1.2,hjust = -0.1) +
  geom_vline(xintercept = obs_o2_snp, color = "blue",  linetype = "dashed", size = 1) +
  annotate("text", x = obs_o2_snp, y = Inf, label = "Observed O2", color = "blue", angle = 90, vjust = -0.5,hjust = 1.1) +
  annotate("text", x = obs_o2_snp, y = 0,   label = paste0("pO2=", signif(p_o2,3)), color = "black", angle = 0, vjust = 1.2,hjust = -0.1) +
  # legend without p-values
  scale_fill_manual(name = "Dataset",
                    values = c("O1" = "#F8766D", "O2" = "#00BFC4")) +
  labs(x = "Simulated SD in SNP rate",
       y = "Frequency",
       title = "SNP Rate SD: O1 vs O2") +
  theme_minimal()

ggsave("snp_sd_comparison.png", plot = p, width = 8, height = 5)

####################INDEL#######################

# Calculate observed SD of rates from original data
obs_o1_indel <- sd(
  dataset_o1[1:N, column_indels] /
    (dataset_o1[1:N, column_coverage] * dataset_o1[1:N, column_generations] * ref)
)
obs_o2_indel <- sd(
  dataset_o2[1:N, column_indels] /
    (dataset_o2[1:N, column_coverage] * dataset_o2[1:N, column_generations] * ref)
)

# Compute p-values: fraction of simulated SDs greater than observed
p_o1 <- mean(o1_indels$sd > obs_o1_indel)
p_o2 <- mean(o2_indels$sd > obs_o2_indel)

# Plot distribution of simulated SD for INDEL rate comparison
plot_data <- data.frame(
  sd  = c(o1_indels$sd, o2_indels$sd),
  set = rep(c("O1", "O2"), each = length(o1_indels$sd))
)

p <- ggplot(plot_data, aes(x = sd, fill = set)) +
  geom_histogram(position = "identity", alpha = 0.6, bins = 50) +
  # observed lines and labels
  geom_vline(xintercept = obs_o1_indel, color = "red",   linetype = "dashed", size = 1) +
  annotate("text", x = obs_o1_indel, y = Inf, label = "Observed O1", color = "red", angle = 90, vjust = -0.5,hjust = 1.1) +
  annotate("text", x = obs_o1_indel, y = 0,   label = paste0("pO1=", signif(p_o1,5)), color = "black", angle = 0, vjust = 1.2,hjust = -0.1) +
  geom_vline(xintercept = obs_o2_indel, color = "blue",  linetype = "dashed", size = 1) +
  annotate("text", x = obs_o2_indel, y = Inf, label = "Observed O2", color = "blue", angle = 90, vjust = -0.5,hjust = 1.1) +
  annotate("text", x = obs_o2_indel, y = 0,   label = sprintf("pO2 = %.6f", p_o2), color = "black", angle = 0, vjust = 1.2,hjust = -0.1) +
  # legend without p-values
  scale_fill_manual(name = "Dataset",
                    values = c("O1" = "#F8766D", "O2" = "#00BFC4")) +
  labs(x = "Simulated SD in INDEL rate",
       y = "Frequency",
       title = "INDEL Rate SD: O1 vs O2") +
  theme_minimal()

ggsave("indel_sd_comparison.png", plot = p, width = 8, height = 5)

####################Total Mutations#######################

# Calculate observed SD of rates from original data
obs_o1_mutation <- sd(
  dataset_o1[1:N, column_mutations] /
    (dataset_o1[1:N, column_coverage] * dataset_o1[1:N, column_generations] * ref)
)
obs_o2_mutation <- sd(
  dataset_o2[1:N, column_mutations] /
    (dataset_o2[1:N, column_coverage] * dataset_o2[1:N, column_generations] * ref)
)

# Compute p-values: fraction of simulated SDs greater than observed
p_o1 <- mean(o1_mutations$sd > obs_o1_mutation)
p_o2 <- mean(o2_mutations$sd > obs_o2_mutation)

# Plot distribution of simulated SD for INDEL rate comparison
plot_data <- data.frame(
  sd  = c(o1_mutations$sd, o2_mutations$sd),
  set = rep(c("O1", "O2"), each = length(o1_mutations$sd))
)

p <- ggplot(plot_data, aes(x = sd, fill = set)) +
  geom_histogram(position = "identity", alpha = 0.6, bins = 50) +
  # observed lines and labels
  geom_vline(xintercept = obs_o1_mutation, color = "red",   linetype = "dashed", size = 1) +
  annotate("text", x = obs_o1_mutation, y = Inf, label = "Observed O1", color = "red", angle = 90, vjust = -0.5,hjust = 1.1) +
  annotate("text", x = obs_o1_mutation, y = 0,   label = paste0("pO1=", signif(p_o1,5)), color = "black", angle = 0, vjust = 1.2,hjust = 0.1) +
  geom_vline(xintercept = obs_o2_mutation, color = "blue",  linetype = "dashed", size = 1) +
  annotate("text", x = obs_o2_mutation, y = Inf, label = "Observed O2", color = "blue", angle = 90, vjust = -0.5,hjust = 1.1) +
  annotate("text", x = obs_o2_mutation, y = 0,   label = sprintf("pO2 = %.6f", p_o2), color = "black", angle = 0, vjust = 1.2,hjust = 0.1) +
  # legend without p-values
  scale_fill_manual(name = "Dataset",
                    values = c("O1" = "#F8766D", "O2" = "#00BFC4")) +
  labs(x = "Simulated SD in Total Mutation rate",
       y = "Frequency",
       title = "Total Mutation Rate SD: O1 vs O2") +
  theme_minimal()

ggsave("mutation_sd_comparison.png", plot = p, width = 8, height = 5)