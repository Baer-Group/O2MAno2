# =======================================
# O1 vs O2 variance sims (O2 = 150 gens)
# =======================================

suppressPackageStartupMessages({
  library(rstudioapi)
  library(ggplot2)
})

# Set working directory to this script's folder (if in RStudio)
try({
  setwd(dirname(getActiveDocumentContext()$path))
}, silent = TRUE)

set.seed(42)
options(scipen = 0)

# -------------------------
# 1) Inputs and parameters
# -------------------------
dataset_o1 <- read.csv("3x_O1.csv", stringsAsFactors = FALSE)
dataset_o2 <- read.csv("3x_O2.csv", stringsAsFactors = FALSE)

# Column indices (as in your previous scripts)
column_generations <- 12
column_coverage    <- 11
column_snps        <- 10
column_indels      <- 16
column_mutations   <- 17

# Callable genome size (bp)
ref <- 100286401

# Simulation replicates
NSIM <- 10000

# FP means per line for each class (used only for O2 simulation and for the 300-gen expectation printout)
fp_lambda_for_column <- function(column_index, column_snps, column_indels) {
  if (column_index == column_snps) {
    2L   # SNP
  } else if (column_index == column_indels) {
    0.33L   # INDEL
  } else {
    2.33L   # TOTAL (SNP + INDEL)
  }
}

# -------------------------------------------------------
# 2) Per-dataset pooled MLE (ratio of sums) for the rate
#    If subtract_fp = TRUE (O2), subtract FP per line
#    before summing (clamped at 0). Returns per-genome
#    per-generation rate (not per-bp).
# -------------------------------------------------------
compute_uniform_rate_dataset <- function(dataset,
                                         column_index,
                                         column_coverage, column_generations,
                                         column_snps, column_indels,
                                         subtract_fp = FALSE) {
  fp_per_line <- if (subtract_fp) fp_lambda_for_column(column_index, column_snps, column_indels) else 0L
  counts_vec  <- pmax(dataset[, column_index] - fp_per_line, 0L)
  expo_vec    <- dataset[, column_coverage] * dataset[, column_generations]
  rate_mu     <- sum(counts_vec) / sum(expo_vec)
  rate_mu
}

# -------------------------------------------------------
# 3) Observed SDs (per-bp per-gen)
#    O1: from O1 alone; O2: from O2 alone (true observed O2 SD)
# -------------------------------------------------------
observed_sd_dataset <- function(dataset, column_index, col_cov, col_gen, ref) {
  sd(dataset[, column_index] / (dataset[, col_cov] * dataset[, col_gen] * ref))
}

# -------------------------------------------------------
# 4) Simulators
#    - O1: simulate using mu_o1, thin by O1 coverage, normalize by
#      (cov1 * gens1 * ref). No FP.
#    - O2: simulate *only* the O2 150 generations:
#      Poisson(gens2 * mu_o2) -> thin by O2 coverage -> add FP -> normalize by
#      (cov2 * gens2 * ref).
# -------------------------------------------------------
simulate_o1 <- function(dataset, mu_o1,
                        column_index, col_cov, col_gen, ref,
                        nsim = NSIM) {
  N <- nrow(dataset)
  Mean_sim <- numeric(nsim)
  Var_sim  <- numeric(nsim)
  SD_sim   <- numeric(nsim)
  
  for (k in seq_len(nsim)) {
    rates <- numeric(N)
    for (i in seq_len(N)) {
      gens1 <- dataset[i, col_gen]
      cov1  <- dataset[i, col_cov]
      # True mutations across gens1 (per-genome)
      raw_true  <- rpois(1, gens1 * mu_o1)
      kept_true <- rbinom(1, size = raw_true, prob = cov1)    # thin by coverage
      rates[i]  <- kept_true / (cov1 * gens1 * ref)           # per-bp per-gen
    }
    Mean_sim[k] <- mean(rates)
    Var_sim[k]  <- var(rates)
    SD_sim[k]   <- sd(rates)
  }
  list(mean = Mean_sim, variance = Var_sim, sd = SD_sim)
}

simulate_o2_150 <- function(dataset_o2, mu_o2,
                            column_index, col_cov_o2, col_gen_o2, ref,
                            nsim = NSIM) {
  N <- nrow(dataset_o2)
  lambda_fp <- fp_lambda_for_column(column_index, column_snps, column_indels)
  
  Mean_sim <- numeric(nsim)
  Var_sim  <- numeric(nsim)
  SD_sim   <- numeric(nsim)
  
  for (k in seq_len(nsim)) {
    rates <- numeric(N)
    for (i in seq_len(N)) {
      gens2 <- dataset_o2[i, col_gen_o2]    # should be ~150
      cov2  <- dataset_o2[i, col_cov_o2]
      
      # True mutations across O2 gens only (per-genome)
      raw_true  <- rpois(1, gens2 * mu_o2)
      kept_true <- rbinom(1, size = raw_true, prob = cov2)
      
      # Add FP at O2 observation
      fp_i <- rpois(1, lambda_fp)
      
      kept_total <- kept_true + fp_i
      rates[i]   <- kept_total / (cov2 * gens2 * ref)   # per-bp per-gen (O2 segment)
    }
    Mean_sim[k] <- mean(rates)
    Var_sim[k]  <- var(rates)
    SD_sim[k]   <- sd(rates)
  }
  list(mean = Mean_sim, variance = Var_sim, sd = SD_sim)
}

# -------------------------------------------------------
# 5) Estimate rates (per-genome per generation)
#    O1: from O1 data (no FP subtraction)
#    O2: from O2 data with FP subtracted per line
# -------------------------------------------------------
mu_o1_snp <- compute_uniform_rate_dataset(dataset_o1, column_snps,
                                          column_coverage, column_generations,
                                          column_snps, column_indels,
                                          subtract_fp = FALSE)
mu_o1_ind <- compute_uniform_rate_dataset(dataset_o1, column_indels,
                                          column_coverage, column_generations,
                                          column_snps, column_indels,
                                          subtract_fp = FALSE)
mu_o1_tot <- compute_uniform_rate_dataset(dataset_o1, column_mutations,
                                          column_coverage, column_generations,
                                          column_snps, column_indels,
                                          subtract_fp = FALSE)

mu_o2_snp <- compute_uniform_rate_dataset(dataset_o2, column_snps,
                                          column_coverage, column_generations,
                                          column_snps, column_indels,
                                          subtract_fp = TRUE)
mu_o2_ind <- compute_uniform_rate_dataset(dataset_o2, column_indels,
                                          column_coverage, column_generations,
                                          column_snps, column_indels,
                                          subtract_fp = TRUE)
mu_o2_tot <- compute_uniform_rate_dataset(dataset_o2, column_mutations,
                                          column_coverage, column_generations,
                                          column_snps, column_indels,
                                          subtract_fp = TRUE)

cat("Per-genome per-generation MLEs:\n")
cat(sprintf("  O1 SNP   = %.6f   (per-bp: %.3e)\n", mu_o1_snp, mu_o1_snp / ref))
cat(sprintf("  O1 INDEL = %.6f   (per-bp: %.3e)\n", mu_o1_ind, mu_o1_ind / ref))
cat(sprintf("  O1 TOTAL = %.6f   (per-bp: %.3e)\n", mu_o1_tot, mu_o1_tot / ref))
cat(sprintf("  O2 SNP   = %.6f   (per-bp: %.3e)  [FP subtracted]\n", mu_o2_snp, mu_o2_snp / ref))
cat(sprintf("  O2 INDEL = %.6f   (per-bp: %.3e)  [FP subtracted]\n", mu_o2_ind, mu_o2_ind / ref))
cat(sprintf("  O2 TOTAL = %.6f   (per-bp: %.3e)  [FP subtracted]\n\n", mu_o2_tot, mu_o2_tot / ref))

# -------------------------------------------------------
# 6) Report the 300-gen uniform expectations (for reference)
#    E[count_300] = 150*mu_o1 + 150*mu_o2 + E[FP]
#    (per line, per-genome). We just print the scalar expectation.
# -------------------------------------------------------
report_300_expectation <- function(mu_o1, mu_o2, column_index) {
  lambda_fp <- fp_lambda_for_column(column_index, column_snps, column_indels)
  150 * mu_o1 + 150 * mu_o2 + lambda_fp
}
cat("Expected observed counts over 300 gens (per line, per-genome):\n")
cat(sprintf("  SNP   : %.3f\n", report_300_expectation(mu_o1_snp, mu_o2_snp, column_snps)))
cat(sprintf("  INDEL : %.3f\n", report_300_expectation(mu_o1_ind, mu_o2_ind, column_indels)))
cat(sprintf("  TOTAL : %.3f\n\n", report_300_expectation(mu_o1_tot, mu_o2_tot, column_mutations)))

# -------------------------------------------------------
# 7) Simulate SDs (O1 and O2 separately, per the new rules)
# -------------------------------------------------------
o1_snps  <- simulate_o1(dataset_o1, mu_o1_snp, column_snps,
                        column_coverage, column_generations, ref, nsim = NSIM)
o1_ind   <- simulate_o1(dataset_o1, mu_o1_ind, column_indels,
                        column_coverage, column_generations, ref, nsim = NSIM)
o1_tot   <- simulate_o1(dataset_o1, mu_o1_tot, column_mutations,
                        column_coverage, column_generations, ref, nsim = NSIM)

o2_snps  <- simulate_o2_150(dataset_o2, mu_o2_snp, column_snps,
                            column_coverage, column_generations, ref, nsim = NSIM)
o2_ind   <- simulate_o2_150(dataset_o2, mu_o2_ind, column_indels,
                            column_coverage, column_generations, ref, nsim = NSIM)
o2_tot   <- simulate_o2_150(dataset_o2, mu_o2_tot, column_mutations,
                            column_coverage, column_generations, ref, nsim = NSIM)

# -------------------------------------------------------
# 8) Observed SDs (true O1, true O2)
# -------------------------------------------------------
obs_o1_snp   <- observed_sd_dataset(dataset_o1, column_snps,     column_coverage, column_generations, ref)
obs_o1_ind   <- observed_sd_dataset(dataset_o1, column_indels,   column_coverage, column_generations, ref)
obs_o1_total <- observed_sd_dataset(dataset_o1, column_mutations, column_coverage, column_generations, ref)

obs_o2_snp   <- observed_sd_dataset(dataset_o2, column_snps,     column_coverage, column_generations, ref)
obs_o2_ind   <- observed_sd_dataset(dataset_o2, column_indels,   column_coverage, column_generations, ref)
obs_o2_total <- observed_sd_dataset(dataset_o2, column_mutations, column_coverage, column_generations, ref)

# -------------------------------------------------------
# 9) p-values (right-tail: simulated SD > observed SD)
# -------------------------------------------------------
p_o1_snp <- mean(o1_snps$sd > obs_o1_snp)
p_o1_ind <- mean(o1_ind$sd  > obs_o1_ind)
p_o1_tot <- mean(o1_tot$sd  > obs_o1_total)

p_o2_snp <- mean(o2_snps$sd > obs_o2_snp)
p_o2_ind <- mean(o2_ind$sd  > obs_o2_ind)
p_o2_tot <- mean(o2_tot$sd  > obs_o2_total)

# -------------------------------------------------------
# 10) Overlay plots (O1 & O2 together, like before)
#     Show first, then save.
# -------------------------------------------------------
plot_sd_overlay <- function(sim_o1_sd, sim_o2_sd, obs_o1_sd, obs_o2_sd, p_o1, p_o2,
                            title, file_out, bins = 50) {
  plot_data <- data.frame(
    sd  = c(sim_o1_sd, sim_o2_sd),
    set = rep(c("O1", "O2"), each = length(sim_o1_sd))
  )
  p <- ggplot(plot_data, aes(x = sd, fill = set)) +
    geom_histogram(position = "identity", alpha = 0.6, bins = bins) +
    # observed lines and labels
    geom_vline(xintercept = obs_o1_sd, color = "red",   linetype = "dashed", size = 1) +
    annotate("text", x = obs_o1_sd, y = Inf, label = "Observed O1",
             color = "red", angle = 90, vjust = -0.5, hjust = 1.1) +
    annotate("text", x = obs_o1_sd, y = 0,   label = sprintf("pO1 = %.3f", p_o1),
             color = "black", angle = 0, vjust = 1.2, hjust = 0.5) +
    geom_vline(xintercept = obs_o2_sd, color = "blue",  linetype = "dashed", size = 1) +
    annotate("text", x = obs_o2_sd, y = Inf, label = "Observed O2",
             color = "blue", angle = 90, vjust = -0.5, hjust = 1.1) +
    annotate("text", x = obs_o2_sd, y = 0,   label = sprintf("pO2 = %.3f", p_o2),
             color = "black", angle = 0, vjust = 1.2, hjust = 0.5) +
    # legend and titles
    scale_fill_manual(name = "Dataset",
                      values = c("O1" = "#F8766D", "O2" = "#00BFC4")) +
    labs(x = "Simulated SD (per-bp per-gen)",
         y = "Frequency",
         title = title) +
    theme_bw(base_size = 18) +        # <— raises the default font for everything
    theme(
      text        = element_text(size = 18),
      axis.title  = element_text(size = 18),
      axis.text   = element_text(size = 16),
      legend.text = element_text(size = 16),
      legend.title= element_text(size = 16),
      plot.title  = element_text(size = 22, face = "bold")
    )
  
  print(p)  # show before saving
  ggsave(file_out, plot = p, width = 8, height = 5, dpi = 500)
}

# SNP
plot_sd_overlay(o1_snps$sd, o2_snps$sd, obs_o1_snp, obs_o2_snp, p_o1_snp, p_o2_snp,
                "SNV Rate SD: O1 vs O2",
                "snp_sd_comparison.png")

# INDEL
plot_sd_overlay(o1_ind$sd,  o2_ind$sd,  obs_o1_ind,  obs_o2_ind,  p_o1_ind, p_o2_ind,
                "INDEL Rate SD: O1 vs O2",
                "indel_sd_comparison.png")

# TOTAL
plot_sd_overlay(o1_tot$sd,  o2_tot$sd,  obs_o1_total, obs_o2_total, p_o1_tot, p_o2_tot,
                "Total Mutation Rate SD: O1 vs O2 (O2 = 150 gens, FP added)",
                "mutation_sd_comparison.png")

# -------------------------------------------------------
# 11) Save simulation draws to CSV (optional)
# -------------------------------------------------------
sim_df <- data.frame(
  o1_snp_mean = o1_snps$mean,  o1_snp_var = o1_snps$variance,  o1_snp_sd = o1_snps$sd,
  o1_ind_mean = o1_ind$mean,   o1_ind_var = o1_ind$variance,   o1_ind_sd = o1_ind$sd,
  o1_tot_mean = o1_tot$mean,   o1_tot_var = o1_tot$variance,   o1_tot_sd = o1_tot$sd,
  o2_snp_mean = o2_snps$mean,  o2_snp_var = o2_snps$variance,  o2_snp_sd = o2_snps$sd,
  o2_ind_mean = o2_ind$mean,   o2_ind_var = o2_ind$variance,   o2_ind_sd = o2_ind$sd,
  o2_tot_mean = o2_tot$mean,   o2_tot_var = o2_tot$variance,   o2_tot_sd = o2_tot$sd
)
write.csv(sim_df, "o1_o2_simulations_150gen_o2.csv", row.names = FALSE)

cat("Done. Plots were shown and saved. CSV: o1_o2_simulations_150gen_o2.csv\n")
