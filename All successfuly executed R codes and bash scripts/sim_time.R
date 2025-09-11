# Set current folder of the R script as the working directory
library(rstudioapi)
library(ggplot2)
library(scales)
library(dplyr)
library(matrixStats)
setwd(dirname(getActiveDocumentContext()$path))
getwd()

# Details of 'variance in mutation rate' simulation
## This script simulates the second-order MA experiment data assuming
## that the null hypothesis of one mutation rate is true.
## We simulate across **1200 days total** (split into O1 = first 600 days and O2 = next 600 days),
## drop mutations according to “no-recall” probabilities per 1200-day total,
## but still compare O1 vs. O2 differences.

# 1) Load data and set constants
dataset_o1 <- read.csv("3x.cleaned_O1.csv")
dataset_o2 <- read.csv("3x.cleaned_O2.csv")

# ref: genome length in bp
ref <- 100286401

# CHANGED: simulate over 1200 days total, 600 days per half
t_total <- 1200    # total days
half_t  <- t_total / 2   # 600 days for O1 and 600 for O2

# Column indices (must match CSV structure)
column_coverage  <- 11   # Coverage fraction (e.g., 0.98)
column_snps      <- 10   # Observed SNP counts
column_indels    <- 16   # Observed indel counts
column_mutations <- 17   # Observed total‐mutation counts

# 2) Merge O1 and O2 into a single dataset covering 1200 days
dataset_all <- dataset_o1

# Sum mutation counts across O1 and O2
cols_to_sum <- c(column_snps, column_indels, column_mutations)
for (col in cols_to_sum) {
  dataset_all[[col]] <- dataset_o1[[col]] + dataset_o2[[col]]
}

# 3) Function to compute per-line observed rates (for summary)
calc_rates <- function(df, days) {
  snp_rate   <- df[[column_snps]]   / (df[[column_coverage]] * days * ref)
  indel_rate <- df[[column_indels]] / (df[[column_coverage]] * days * ref)
  mut_rate   <- df[[column_mutations]] / (df[[column_coverage]] * days * ref)
  data.frame(snp_rate, indel_rate, mut_rate)
}

# 4) Compute and report means & SDs on the real (combined) data over 1200 days
rates_all <- calc_rates(dataset_all, t_total)
mean_all  <- colMeans(rates_all)
sd_all    <- colSds(as.matrix(rates_all))

mean_all
sd_all

# 5) Compute observed O2 - O1 differences for comparison (each 600-day half)
rates_o1 <- calc_rates(dataset_o1, half_t)
rates_o2 <- calc_rates(dataset_o2, half_t)
mean_o1  <- colMeans(rates_o1)
mean_o2  <- colMeans(rates_o2)

observed_diff_snp   <- mean_o2["snp_rate"]   - mean_o1["snp_rate"]
observed_diff_indel <- mean_o2["indel_rate"] - mean_o1["indel_rate"]
observed_diff_mut   <- mean_o2["mut_rate"]   - mean_o1["mut_rate"]

# 6) Simulation function with single-step “no-recall” drop,
#    using fixed 1200 days to estimate true μ, then splitting into 600-day halves.
simulate_split <- function(dataset, column_index, iterations = 10000) {
  # define “no-recall” probabilities over 1200 days
  no_recall_snp   <- 1.7798 / 100   # 0.017798 fraction of true SNPs dropped
  no_recall_indel <- 2.2574 / 100   # 0.022574 fraction of true indels dropped
  
  no_recall <- if (column_index == column_snps) {
    no_recall_snp
  } else if (column_index == column_indels) {
    no_recall_indel
  } else {
    0
  }
  
  # CHANGED: estimate true μ_bp_per_day using t_total = 1200 days
  true_mu_bp_per_day <- mean(
    dataset[[column_index]] /
      (dataset[[column_coverage]]  * t_total * ref)
  )
  
  n <- nrow(dataset)
  mean_o1_sim <- numeric(iterations)
  mean_o2_sim <- numeric(iterations)
  var_o1_sim  <- numeric(iterations)
  var_o2_sim  <- numeric(iterations)
  
  for (k in seq_len(iterations)) {
    rates_o1 <- numeric(n)
    rates_o2 <- numeric(n)
    
    for (i in seq_len(n)) {
      cov_i <- dataset[i, column_coverage]
      # CHANGED: use half_t = 600 days for each half
      gen1 <- half_t
      gen2 <- half_t
      
      p_keep <- (1 - no_recall)
      
      # O1 half
      lambda1     <- true_mu_bp_per_day * ref * gen1
      true_count1 <- rpois(1, lambda1)
      observed1   <- rbinom(1, true_count1, cov_i)
      rates_o1[i] <- observed1 / (cov_i * gen1 * ref)
      
      # O2 half
      lambda2     <- true_mu_bp_per_day * ref * gen2
      true_count2 <- rpois(1, lambda2)
      observed2   <- rbinom(1, true_count2, cov_i)
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

# 7) Run simulations
set.seed(4)
sim_snp   <- simulate_split(dataset_all, column_snps,   iterations = 10000)
sim_indel <- simulate_split(dataset_all, column_indels, iterations = 10000)
sim_mut   <- simulate_split(dataset_all, column_mutations, iterations = 10000)

# 8) Compute simulated Δ = O2 – O1
diff_snp_sim   <- sim_snp$mean_o2   - sim_snp$mean_o1
diff_indel_sim <- sim_indel$mean_o2 - sim_indel$mean_o1
diff_mut_sim   <- sim_mut$mean_o2   - sim_mut$mean_o1

# 9) Plot histograms for each category

## SNPs
data.frame(diff_snp = diff_snp_sim) %>%
  ggplot(aes(x = diff_snp)) +
  geom_histogram(bins = 50, fill = "cyan", color = "white") +
  geom_vline(xintercept = observed_diff_snp, color = "firebrick",
             linetype = "dashed", linewidth = 1) +
  annotate("text", x = observed_diff_snp, y = Inf,
           label = paste0("Obs Δ = ", format(observed_diff_snp,
                                             scientific = TRUE, digits = 2)),
           color = "firebrick", vjust = 1.5, hjust = 1.1) +
  labs(title = "Simulated vs Observed Δ SNP Rate (O2–O1)",
       x = "Δ Mean SNP Rate", y = "Frequency") +
  scale_x_continuous(labels = scientific_format(digits = 2)) +
  theme_bw()

## Indels
data.frame(diff_indel = diff_indel_sim) %>%
  ggplot(aes(x = diff_indel)) +
  geom_histogram(bins = 50, fill = "aquamarine", color = "white") +
  geom_vline(xintercept = observed_diff_indel, color = "firebrick",
             linetype = "dashed", linewidth = 1) +
  annotate("text", x = observed_diff_indel, y = Inf,
           label = paste0("Obs Δ = ", format(observed_diff_indel,
                                             scientific = TRUE, digits = 2)),
           color = "firebrick", vjust = 1.5, hjust = 1.1) +
  labs(title = "Simulated vs Observed Δ Indel Rate (02–01)",
       x = "Δ Mean Indel Rate", y = "Frequency") +
  scale_x_continuous(labels = scientific_format(digits = 2)) +
  theme_bw()

## Total mutations
data.frame(diff_mut = diff_mut_sim) %>%
  ggplot(aes(x = diff_mut)) +
  geom_histogram(bins = 50, fill = "blue", color = "white") +
  geom_vline(xintercept = observed_diff_mut, color = "firebrick",
             linetype = "dashed", linewidth = 1) +
  annotate("text", x = observed_diff_mut, y = Inf,
           label = paste0("Obs Δ = ", format(observed_diff_mut,
                                             scientific = TRUE, digits = 2)),
           color = "firebrick", vjust = 1.5, hjust = 1.1) +
  labs(title = "Simulated vs Observed Δ Total Mutation Rate (600d–600d)",
       x = "Δ Mean Total Mutation Rate", y = "Frequency") +
  scale_x_continuous(labels = scientific_format(digits = 2)) +
  theme_bw()

