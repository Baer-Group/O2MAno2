# Set current folder of the R script as the working directory
library(rstudioapi)
library(ggplot2)
library(scales)
library(matrixStats)
setwd(dirname(getActiveDocumentContext()$path))
getwd()

# Details of 'variance in mutation rate' simulation
## This script simulates the second order MA experiment data assuming
## that the null hypothesis of one mutation rate is true
## It tries to accommodate for sampling variance and unequal sampling

# 1) Load data and set constants
dataset_o1 <- read.csv("3x.cleaned_O1.csv")
dataset_o2 <- read.csv("3x.cleaned_O2.csv")

# ref, column_coverage, column_generations, etc. all defined
# Define column indices
ref                <- 100286401
column_generations <- 12
column_coverage    <- 11
column_snps        <- 10
column_indels      <- 16
column_mutations   <- 17

# 2) Function to compute per-line rates
calc_rates <- function(df) {
  # per-line rates for SNPs, indels, total muts
  snp_rate   <- (df[[column_snps]]   * 100) / (df[[column_coverage]] * df[[column_generations]] * ref)
  indel_rate <- (df[[column_indels]] * 100) / (df[[column_coverage]] * df[[column_generations]] * ref)
  mut_rate   <- (df[[column_mutations]] * 100) / (df[[column_coverage]] * df[[column_generations]] * ref)
  
  data.frame(snp_rate, indel_rate, mut_rate)
}

# 3) Compute and report means
rates_o1 <- calc_rates(dataset_o1)
rates_o2 <- calc_rates(dataset_o2)

mean_o1 <- colMeans(rates_o1)
mean_o2 <- colMeans(rates_o2)

sd_o1 <- colSds(as.matrix( rates_o1[ , sapply(rates_o1, is.numeric) ] ))
sd_o2 <- colSds(as.matrix( rates_o2[ , sapply(rates_o2, is.numeric) ] ))

sd(rates_o1[,1])

# 4) Compute observed difference 
observed_diff_snp <- mean_o2["snp_rate"] - mean_o1["snp_rate"]
observed_diff_indel <- mean_o2["indel_rate"] - mean_o1["indel_rate"]
observed_diff_mut <- mean_o2["mut_rate"] - mean_o1["mut_rate"]

# 5) Make a copy of O1 to hold the aggregated data
dataset_all <- dataset_o1

# 6) Sum the 4 columns that change across O1-O2 (indices 10, 12, 16, 17)
cols_to_sum <- c(10, 12, 16, 17)
for (col in cols_to_sum) {
  dataset_all[[col]] <- dataset_o1[[col]] + dataset_o2[[col]]
}


# Preview the first few rows
head(dataset_all)



# Initialize arrays for means and variances
Var_total <- array()
Mean_total <- array()
Var_indels <- array()
Mean_indels <- array()
Var_snps <- array()
Mean_snps <- array()


simulate_split <- function(dataset, column_index) {
  # CHANGED: define false‐negative rates
  fn_rate_snp   <- 0.223 / 100    ## 0.223% SNPs missed  
  fn_rate_indel <- 0.5228 / 100   ## 0.5228% indels missed
  
  # pick the correct FN rate based on which column we’re simulating
  fn_rate <- if (column_index == column_snps) {
    fn_rate_snp
  } else if (column_index == column_indels) {
    fn_rate_indel
  } else {
    0  # no FN for total‐mutations simulation, or adjust if desired
  }                             ## CHANGED
  
  # 1) Compute the overall mean μ from the real data (per‐gen rate)
  total_mu <- mean(
    (dataset[, column_index] * 100) /
      (dataset[, column_coverage] * dataset[, column_generations])
  )
  
  n <- nrow(dataset)
  iterations <- 10000
  
  # Prepare storage
  mean_o1 <- numeric(iterations)
  mean_o2 <- numeric(iterations)
  var_o1  <- numeric(iterations)
  var_o2  <- numeric(iterations)
  
  # 2) Simulate
  for (k in seq_len(iterations)) {
    rates_o1 <- numeric(n)
    rates_o2 <- numeric(n)
    
    for (i in seq_len(n)) {
      cov_i <- dataset[i, column_coverage]
      
      # Simulate first 150 generations
      muts1 <- sum(rpois(150, total_mu))
      pass1 <- sum(runif(muts1, 0, 100) < cov_i)
      
      # CHANGED: apply false-negative filter to the “passed” mutations
      observed1 <- sum(rbinom(pass1, 1, 1 - fn_rate))
      rates_o1[i] <- (observed1 * 100) / (cov_i * 150 * ref)
      
      # Simulate second 150 generations
      muts2 <- sum(rpois(150, total_mu))
      pass2 <- sum(runif(muts2, 0, 100) < cov_i)
      
      # CHANGED: apply false-negative filter here as well
      observed2 <- sum(rbinom(pass2, 1, 1 - fn_rate))
      rates_o2[i] <- (observed2 * 100) / (cov_i * 150 * ref)
    }
    
    mean_o1[k] <- mean(rates_o1)
    mean_o2[k] <- mean(rates_o2)
    var_o1[k]  <- var(rates_o1)
    var_o2[k]  <- var(rates_o2)
  }
  
  list(
    mean_o1 = mean_o1,
    mean_o2 = mean_o2,
    var_o1  = var_o1,
    var_o2  = var_o2
  )
}


# 7) Run for mutations, indels, and snps
sim_mut   <- simulate_split(dataset_all, column_mutations)
sim_indel <- simulate_split(dataset_all, column_indels)
sim_snp   <- simulate_split(dataset_all, column_snps)


# 8) Compute differences (O2 - O1) for each iteration
diff_mut   <- sim_mut$mean_o2   - sim_mut$mean_o1
diff_indel <- sim_indel$mean_o2 - sim_indel$mean_o1
diff_snp   <- sim_snp$mean_o2   - sim_snp$mean_o1

####SNP
sd
# 9) Build simulation dataframe
diff_df <- as.data.frame(diff_snp)  
# or if already in memory: diff_df <- data.frame(diff_indel = diff_indel)

# 10) Plot histogram of simulated SNP differences with observed line
ggplot(diff_df, aes(x = diff_snp)) +
  geom_histogram(
    bins = 50,
    fill = "cyan",
    color = "white"
  ) +
  geom_vline(
    xintercept = observed_diff_snp,
    color = "firebrick",
    linetype = "dashed",
    linewidth = 2
  ) +
  annotate(
    "text",
    x = observed_diff_snp,
    y = Inf,
    label = paste0("Obs Δ = ", 
                   format(observed_diff_snp, scientific = TRUE, digits = 2)),
    color = "firebrick",
    vjust = 1.5,
    hjust = 1.1
  ) +
  scale_x_continuous(
    labels = scientific_format(digits = 2)
  ) +
  labs(
    title = "Simulated vs Observed Difference in SNP Rate (O2 - O1)",
    x     = "Difference in Mean SNP Rate",
    y     = "Frequency"
  ) +
  theme_bw()



#####Indel
# 9) Build simulation dataframe
diff_df <- as.data.frame(diff_indel)  


# 10) Plot histogram of simulated Indel differences with observed line
ggplot(diff_df, aes(x = diff_indel)) +
  geom_histogram(
    bins = 50,
    fill = "aquamarine",
    color = "white"
  ) +
  geom_vline(
    xintercept = observed_diff_indel,
    color = "firebrick",
    linetype = "dashed",
    linewidth = 2
  ) +
  annotate(
    "text",
    x = observed_diff_indel,
    y = Inf,
    label = paste0("Obs Δ = ", 
                   format(observed_diff_indel, scientific = TRUE, digits = 2)),
    color = "firebrick",
    vjust = 1.5,
    hjust = 1.1
  ) +
  scale_x_continuous(
    labels = scientific_format(digits = 2)
  ) +
  labs(
    title = "Simulated vs Observed Difference in INDEL Rate (O2 - O1)",
    x     = "Difference in Mean INDEL Rate",
    y     = "Frequency"
  ) +
  theme_bw()



##### Total
# 9) Build simulation dataframe
diff_df <- as.data.frame(diff_mut)  


# 10) Plot histogram of simulated total mutation differences with observed line
ggplot(diff_df, aes(x = diff_mut)) +
  geom_histogram(
    bins = 50,
    fill = "blue",
    color = "white"
  ) +
  geom_vline(
    xintercept = observed_diff_mut,
    color = "firebrick",
    linetype = "dashed",
    linewidth = 2
  ) +
  annotate(
    "text",
    x = observed_diff_mut,
    y = Inf,
    label = paste0("Obs Δ = ", 
                   format(observed_diff_mut, scientific = TRUE, digits = 2)),
    color = "firebrick",
    vjust = 1.5,
    hjust = 1.1
  ) +
  scale_x_continuous(
    labels = scientific_format(digits = 2)
  ) +
  labs(
    title = "Simulated vs Observed Difference in Total Mutation Rate (O2 - O1)",
    x     = "Difference in Mean Mutation Rate",
    y     = "Frequency"
  ) +
  theme_bw()
