# Set current folder of the R script as the working directory
library(rstudioapi)
library(ggplot2)
library(scales)
library(matrixStats)
setwd(dirname(getActiveDocumentContext()$path))
getwd()

# --- Inputs ----------------------------------------------------
ref <- 100286401  # size of callable/reference genome in bp

# Adjust these if your CSV layout differs
col_gen   <- 12  # generations
col_cov   <- 11  # coverage (fraction)
col_snp   <- 10  # SNP count
col_indel <- 14  # indel count
col_mut   <- 17  # total mutation count

o1 <- read.csv("3x_O1.csv", stringsAsFactors = FALSE)
o2 <- read.csv("3x_O2.csv", stringsAsFactors = FALSE)

# --- Helpers ---------------------------------------------------
rate_vec <- function(df, count_col) {
  df[[count_col]] / (df[[col_cov]] * df[[col_gen]] * ref)
}
sem <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

summarize_rates <- function(df, label) {
  r_snp   <- rate_vec(df, col_snp)
  r_indel <- rate_vec(df, col_indel)
  r_mut   <- rate_vec(df, col_mut)
  data.frame(
    dataset = label,
    class   = c("SNP", "Indel", "Total"),
    mean    = c(mean(r_snp),   mean(r_indel),   mean(r_mut)),
    sem     = c(sem(r_snp),    sem(r_indel),    sem(r_mut)),
    stringsAsFactors = FALSE
  )
}

# --- Compute ---------------------------------------------------
out <- rbind(
  summarize_rates(o1, "O1"),
  summarize_rates(o2, "O2")
)

# View results
print(out)


#######days######
# --- Inputs ---
ref <- 100286401              # callable/reference bp
t_o1 <- 600                   # days (O1)
t_o2 <- 624                   # days (O2)

# Adjust if your CSV layout differs
col_cov   <- 11               # coverage (fraction)
col_snp   <- 10               # SNP count
col_indel <- 14               # indel count
col_mut   <- 17               # total mutation count

o1 <- read.csv("3x_O1.csv", stringsAsFactors = FALSE)
o2 <- read.csv("3x_O2.csv", stringsAsFactors = FALSE)

# --- Helpers ---
rate_vec_time <- function(df, count_col, t_days) {
  df[[count_col]] / (df[[col_cov]] * ref * t_days)
}
sem <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

summarize_rates_time <- function(df, label, t_days) {
  r_snp   <- rate_vec_time(df, col_snp,   t_days)
  r_indel <- rate_vec_time(df, col_indel, t_days)
  r_mut   <- rate_vec_time(df, col_mut,   t_days)
  data.frame(
    dataset = label,
    time_days = t_days,
    class   = c("SNP", "Indel", "Total"),
    mean    = c(mean(r_snp),   mean(r_indel),   mean(r_mut)),
    sem     = c(sem(r_snp),    sem(r_indel),    sem(r_mut)),
    stringsAsFactors = FALSE
  )
}

# --- Compute (per bp per day) ---
out <- rbind(
  summarize_rates_time(o1, "O1", t_o1),
  summarize_rates_time(o2, "O2", t_o2)
)

print(out)


# Details of 'variance in mutation rate' simulation
## This script simulates the second-order MA experiment data assuming
## that the null hypothesis of one mutation rate is true.
## We simulate across 300 total generations (split into O1 = first 150 and O2 = next 150),
## drop mutations according to “no-recall” probabilities per 300-generation total,
## but still compare O1 vs. O2 differences.

# 1) Load data and set constants
dataset_o1 <- read.csv("3x_O1.csv")
dataset_o2 <- read.csv("3x_O2.csv")

# ref: genome length in bp
ref <- 100286401

# Column indices (must match CSV structure)
column_generations <- 12   # Number of generations (150 in O1 and 150 in O2)
column_coverage    <- 11   # Coverage fraction (e.g., 0.98)
column_snps        <- 10   # Observed SNP counts
column_indels      <- 14   # Observed indel counts
column_mutations   <- 17   # Observed total‐mutation counts


# Assumptions for corrections (3X):
FP_snp   <- 2
FP_indel <- 0.333
pFTR_snp   <- 0.022   # 2.2%
pFTR_indel <- 0.025   # 2.5%

point_estimates_epoch <- function(df, epoch = c("O1","O2")) {
  epoch <- match.arg(epoch)
  t  <- df[[column_generations]]
  nS <- df[[column_snps]]
  nI <- df[[column_indels]]
  
  # O2-only FP subtraction; O1 untouched
  adjS <- if (epoch == "O2") pmax(nS - FP_snp,   0) else nS
  adjI <- if (epoch == "O2") pmax(nI - FP_indel, 0) else nI
  
  # FtR added as counts (proportion � observed n)
  FtR_S <- pFTR_snp   * nS
  FtR_I <- pFTR_indel * nI
  
  # Class-specific adjusted counts
  adjS <- adjS + FtR_S
  adjI <- adjI + FtR_I
  adjT <- adjS + adjI  # total is sum of classes, consistent with the formula
  
  # Per-generation point estimates (NO coverage term)
  mu_snp   <- adjS / (ref * t)
  mu_indel <- adjI / (ref * t)
  mu_total <- adjT / (t)
  
  data.frame(mu_snp, mu_indel, mu_total)
}

# --- Compute per-line point estimates for each epoch ---
mu_o1 <- point_estimates_epoch(dataset_o1, "O1")
mu_o2 <- point_estimates_epoch(dataset_o2, "O2")

# --- Report epoch means (point estimates) ---
mu_o1_mean <- colMeans(mu_o1, na.rm = TRUE)
mu_o2_mean <- colMeans(mu_o2, na.rm = TRUE)

mu_o1_mean
mu_o2_mean

# --- Helper: SEM ---
sem <- function(v) {
  v <- v[!is.na(v)]
  if (length(v) < 2) return(NA_real_)
  sd(v) / sqrt(length(v))
}

# --- Epoch means and SEMs ---
o1_mean <- sapply(mu_o1, mean, na.rm = TRUE)
o1_sem  <- sapply(mu_o1, sem)

o2_mean <- sapply(mu_o2, mean, na.rm = TRUE)
o2_sem  <- sapply(mu_o2, sem)

# Tidy summary table (per-generation point estimates; no coverage in denominator)
summary_sem <- data.frame(
  Epoch = rep(c("O1","O2"), each = 3),
  Class = rep(c("SNP","Indel","Total"), times = 2),
  Mean  = c(o1_mean, o2_mean),
  SEM   = c(o1_sem,  o2_sem)
)

summary_sem

###### DAYS #####

# --- Day constants ---
days_o1 <- 600L
days_o2 <- 624L

# Assumptions for corrections (3X):
FP_snp     <- 2
FP_indel   <- 0.333
pFTR_snp   <- 0.022   # 2.2%
pFTR_indel <- 0.025   # 2.5%

# Point estimates per epoch (PER DAY; no coverage in denominator)
point_estimates_epoch <- function(df, epoch = c("O1","O2"),
                                  days_o1 = 600L, days_o2 = 624L) {
  epoch <- match.arg(epoch)
  n <- nrow(df)
  t  <- if (epoch == "O1") rep(days_o1, n) else rep(days_o2, n)
  
  nS <- df[[column_snps]]
  nI <- df[[column_indels]]
  
  # O2-only FP subtraction; O1 untouched
  adjS <- if (epoch == "O2") pmax(nS - FP_snp,   0) else nS
  adjI <- if (epoch == "O2") pmax(nI - FP_indel, 0) else nI
  
  # FtR added as counts (proportion � observed n)
  FtR_S <- pFTR_snp   * nS
  FtR_I <- pFTR_indel * nI
  
  # Class-specific adjusted counts
  adjS <- adjS + FtR_S
  adjI <- adjI + FtR_I
  adjT <- adjS + adjI
  
  # Per-day point estimates (NO coverage term)
  mu_snp_day   <- adjS / (ref * t)
  mu_indel_day <- adjI / (ref * t)
  mu_total_day <- adjT / (ref * t)
  
  data.frame(mu_snp_day, mu_indel_day, mu_total_day)
}

# --- Compute per-line point estimates using DAYS ---
mu_o1 <- point_estimates_epoch(dataset_o1, "O1", days_o1, days_o2)
mu_o2 <- point_estimates_epoch(dataset_o2, "O2", days_o1, days_o2)

# --- Means and SEMs (unchanged logic) ---
sem <- function(v) { v <- v[!is.na(v)]; if (length(v) < 2) return(NA_real_); sd(v) / sqrt(length(v)) }

o1_mean <- sapply(mu_o1, mean, na.rm = TRUE)
o1_sem  <- sapply(mu_o1, sem)

o2_mean <- sapply(mu_o2, mean, na.rm = TRUE)
o2_sem  <- sapply(mu_o2, sem)

summary_sem_days <- data.frame(
  Epoch = rep(c("O1","O2"), each = 3),
  Class = rep(c("SNP","Indel","Total"), times = 2),
  Mean_per_day = c(o1_mean, o2_mean),
  SEM_per_day  = c(o1_sem,  o2_sem)
)

summary_sem_days





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


# 4) Compute and report means & SEMs on the real (summed) data
rates_all <- calc_rates(dataset_all)
mean_all  <- colMeans(rates_all, na.rm = TRUE)

# SEM = SD / sqrt(n)  (n = non-NA count per column)
n_all   <- colSums(!is.na(rates_all))
sem_all <- colSds(as.matrix(rates_all), na.rm = TRUE) / sqrt(n_all)

mean_all
sem_all
Overall_mean <- mean_all
Overall_SEM  <- sem_all


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
  
  
  # 6a) Estimate “true” μ_bp_per_gen from the real combined (300-gen) data:
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
      lambda1 <- true_mu_bp_per_gen * cov_i * gen1 * ref
      observed1 <- rpois(1, lambda1)
      
      # compute observed rate for O1
      rates_o1[i] <- observed1 / (cov_i * gen1 * ref)
      
      # 6b-ii) Simulate TRUE mutations in O2 (150 gen)
      lambda2 <- true_mu_bp_per_gen * cov_i * gen2 * ref
      observed2 <- rpois(1, lambda2)
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
set.seed(2)  # for reproducibility

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

