# Set current folder of the R script as the working directory
library(rstudioapi)
library(ggplot2)
library(scales)
library(matrixStats)
setwd(dirname(getActiveDocumentContext()$path))
getwd()

# Details of 'variance in mutation rate' simulation
## This script simulates the second-order MA experiment data assuming
## that the null hypothesis of one mutation rate is true.
## We simulate across 300 total generations (split into O1 = first 150 and O2 = next 150),
## drop mutations according to “no-recall” probabilities per 300-generation total,
## but still compare O1 vs. O2 differences.

# 1) Load data and set constants
dataset_o1 <- read.csv("3x.cleaned_O1.csv")
dataset_o2 <- read.csv("3x.cleaned_O2.csv")

# ref: genome length in bp
ref <- 100286401

# Column indices (must match CSV structure)
column_generations <- 12   # Number of generations (150 in O1 and 150 in O2)
column_coverage    <- 11   # Coverage fraction (e.g., 0.98)
column_snps        <- 10   # Observed SNP counts
column_indels      <- 16   # Observed indel counts
column_mutations   <- 17   # Observed total‐mutation counts

# 2) Merge O1 and O2 into a 300-generation dataset
dataset_all <- dataset_o1

# Sum mutation counts across O1 and O2
cols_to_sum <- c(column_snps, column_indels, column_mutations)
for (col in cols_to_sum) {
  dataset_all[[col]] <- dataset_o1[[col]] + dataset_o2[[col]]
}

# Set the total generations to 300
dataset_all[[column_generations]] <- dataset_o1[[column_generations]] + dataset_o2[[column_generations]]
# (i.e., 150 + 150 = 300 for each line)

# Coverage stays as in O1/O2 (same column holds coverage fraction)

# Preview
head(dataset_all)

# 3) Function to compute per-line observed rates (for summary)
calc_rates <- function(df) {
  snp_rate   <- df[[column_snps]] / (df[[column_coverage]] * df[[column_generations]] * ref)
  indel_rate <- df[[column_indels]] / (df[[column_coverage]] * df[[column_generations]] * ref)
  mut_rate   <- df[[column_mutations]] / (df[[column_coverage]] * df[[column_generations]] * ref)
  data.frame(snp_rate, indel_rate, mut_rate)
}

# 4) Compute and report means & SDs on the real (summed) data
rates_all <- calc_rates(dataset_all)
mean_all  <- colMeans(rates_all)
sd_all    <- colSds(as.matrix(rates_all))

mean_all
sd_all

# 5) Compute observed O2 - O1 differences for comparison
mean_o1 <- colMeans(calc_rates(dataset_o1))
mean_o2 <- colMeans(calc_rates(dataset_o2))
observed_diff_snp   <- mean_o2["snp_rate"]   - mean_o1["snp_rate"]
observed_diff_indel <- mean_o2["indel_rate"] - mean_o1["indel_rate"]
observed_diff_mut   <- mean_o2["mut_rate"]   - mean_o1["mut_rate"]

# 6) Simulation function: split into O1 (150 gen) and O2 (150 gen),
#    but estimate the “true” per‐bp per‐gen rate from the combined 300‐generation data,
#    then simulate true counts in each half, apply no-recall drop, and compute rates.
simulate_split <- function(dataset, column_index, iterations = 10000) {
  # CHANGED: define “no-recall” probabilities (per 300-generation total)
  no_recall_snp   <- 1.7798 / 100   # 0.017798 chance a true SNP is dropped
  no_recall_indel <- 2.2574 / 100   # 0.022574 chance a true indel is dropped
  
  # pick correct no-recall rate based on column being simulated
  no_recall <- if (column_index == column_snps) {
    no_recall_snp
  } else if (column_index == column_indels) {
    no_recall_indel
  } else {
    0  # no extra drop for total‐mutation simulation
  }
  
  # 6a) Estimate “true” μ_bp_per_gen from the real combined (300-gen) data:
  #     Each observed count_i = true_mu_bp_per_gen * cov_i * 300 * ref 
  true_mu_bp_per_gen <- mean(
    dataset[[column_index]] /
      (dataset[[column_coverage]] * dataset[[column_generations]] * ref )
  )
  
  n <- nrow(dataset)
  
  # storage for simulated O1 and O2 means and variances
  mean_o1_sim <- numeric(iterations)
  mean_o2_sim <- numeric(iterations)
  var_o1_sim  <- numeric(iterations)
  var_o2_sim  <- numeric(iterations)
  
  # 6b) Monte Carlo loop
  for (k in seq_len(iterations)) {
    rates_o1 <- numeric(n)
    rates_o2 <- numeric(n)
    
    for (i in seq_len(n)) {
      cov_i <- dataset[i, column_coverage]      # e.g. 0.98
      gen1  <- 150                              # first half
      gen2  <- 150                              # second half
      
      # 6b-i) Simulate TRUE mutations in O1 (150 gen)
      lambda1 <- true_mu_bp_per_gen * gen1 * ref
      true_count1 <- rpois(1, lambda1)
      # apply no-recall drop
      observed1 <- rbinom(1, true_count1, 1 - no_recall)
      # compute observed rate for O1
      rates_o1[i] <- observed1 / (cov_i * gen1 * ref)
      
      # 6b-ii) Simulate TRUE mutations in O2 (150 gen)
      lambda2 <- true_mu_bp_per_gen  * gen2 * ref
      true_count2 <- rpois(1, lambda2)
      # apply no-recall drop
      observed2 <- rbinom(1, true_count2, 1 - no_recall)
      # compute observed rate for O2
      rates_o2[i] <- observed2 / (cov_i * gen2 * ref)
    }
    
    mean_o1_sim[k] <- mean(rates_o1)
    mean_o2_sim[k] <- mean(rates_o2)
    var_o1_sim[k]  <- var(rates_o1)
    var_o2_sim[k]  <- var(rates_o2)
  }
  
  list(
    mean_o1 = mean_o1_sim,
    mean_o2 = mean_o2_sim,
    var_o1  = var_o1_sim,
    var_o2  = var_o2_sim
  )
}

# 7) Run simulations for SNPs, indels, and total mutations
set.seed(2025)  # for reproducibility

sim_snp <- simulate_split(dataset_all, column_snps,   iterations = 10000)
sim_indel <- simulate_split(dataset_all, column_indels, iterations = 10000)
sim_mut <- simulate_split(dataset_all, column_mutations, iterations = 10000)

# 8) Compute simulated differences (O2 - O1) for each iteration
diff_snp_sim   <- sim_snp$mean_o2 - sim_snp$mean_o1
diff_indel_sim <- sim_indel$mean_o2 - sim_indel$mean_o1
diff_mut_sim   <- sim_mut$mean_o2 - sim_mut$mean_o1

# 9) Plot histograms of simulated differences with observed lines

## SNPs
diff_df_snp <- data.frame(diff_snp = diff_snp_sim)
ggplot(diff_df_snp, aes(x = diff_snp)) +
  geom_histogram(bins = 50, fill = "cyan", color = "white") +
  geom_vline(xintercept = observed_diff_snp, color = "firebrick", linetype = "dashed", linewidth = 1) +
  annotate("text",
           x = observed_diff_snp,
           y = Inf,
           label = paste0("Obs Δ = ", format(observed_diff_snp, scientific = TRUE, digits = 2)),
           color = "firebrick",
           vjust = 1.5,
           hjust = 1.1) +
  scale_x_continuous(labels = scientific_format(digits = 2)) +
  labs(title = "Simulated vs Observed Δ SNP Rate (O2 - O1)",
       x     = "Δ Mean SNP Rate",
       y     = "Frequency") +
  theme_bw()

## Indels
diff_df_indel <- data.frame(diff_indel = diff_indel_sim)
ggplot(diff_df_indel, aes(x = diff_indel)) +
  geom_histogram(bins = 50, fill = "aquamarine", color = "white") +
  geom_vline(xintercept = observed_diff_indel, color = "firebrick", linetype = "dashed", linewidth = 1) +
  annotate("text",
           x = observed_diff_indel,
           y = Inf,
           label = paste0("Obs Δ = ", format(observed_diff_indel, scientific = TRUE, digits = 2)),
           color = "firebrick",
           vjust = 1.5,
           hjust = 1.1) +
  scale_x_continuous(labels = scientific_format(digits = 2)) +
  labs(title = "Simulated vs Observed Δ Indel Rate (O2 - O1)",
       x     = "Δ Mean Indel Rate",
       y     = "Frequency") +
  theme_bw()

## Total Mutations
diff_df_mut <- data.frame(diff_mut = diff_mut_sim)
ggplot(diff_df_mut, aes(x = diff_mut)) +
  geom_histogram(bins = 50, fill = "blue", color = "white") +
  geom_vline(xintercept = observed_diff_mut, color = "firebrick", linetype = "dashed", linewidth = 1) +
  annotate("text",
           x = observed_diff_mut,
           y = Inf,
           label = paste0("Obs Δ = ", format(observed_diff_mut, scientific = TRUE, digits = 2)),
           color = "firebrick",
           vjust = 1.5,
           hjust = 1.1) +
  scale_x_continuous(labels = scientific_format(digits = 2)) +
  labs(title = "Simulated vs Observed Δ Total Mutation Rate (O2 - O1)",
       x     = "Δ Mean Total Mutation Rate",
       y     = "Frequency") +
  theme_bw()
