# ============================================================
#Simulate Δ = mean(O2) - mean(O1) from μ_O1, grow FP until
# observed Δ > 95% of simulated Δ (i.e., > 95th percentile)
# ============================================================

suppressPackageStartupMessages({
  library(rstudioapi)
  library(ggplot2)
  library(scales)
  library(dplyr)
})

# --- Setup ---------------------------------------------------
try({ setwd(dirname(getActiveDocumentContext()$path)) }, silent = TRUE)
set.seed(2025)
options(scipen = 0)

ref  <- 100286401           # callable bp
gcap <- 150                 # cap generations at 150 if needed
NSIM <- 10000               # simulations per FP level
MAX_FP <- 20                # upper bound search for FP (per line)

# --- Load ----------------------------------------------------
dataset_o1 <- read.csv("3x_O1.csv", stringsAsFactors = FALSE)
dataset_o2 <- read.csv("3x_O2.csv", stringsAsFactors = FALSE)

# Column indices
col_gen   <- 12  # generations (per dataset, usually ~150)
col_cov   <- 11  # coverage fraction
col_snp   <- 10  # observed SNP count
col_indel <- 16  # observed indel count
col_mut   <- 17  # observed total mutation count

# ---------- Helpers -----------------------------------------
# pooled MLE (ratio of sums) for μ from O1 only (per-genome per-generation)
mu_from_o1 <- function(dataset, column) {
  counts <- dataset[[column]]
  expo   <- dataset[[col_cov]] * dataset[[col_gen]]
  sum(counts) / sum(expo)
}

# Compute observed mean rates (per-bp per-gen) and Δ = mean(O2)-mean(O1)
observed_delta <- function(col_index) {
  r1 <- dataset_o1[[col_index]] /
    (dataset_o1[[col_cov]] * dataset_o1[[col_gen]] * ref)
  r2 <- dataset_o2[[col_index]] /
    (dataset_o2[[col_cov]] * dataset_o2[[col_gen]] * ref)
  c(mean_o1 = mean(r1), mean_o2 = mean(r2), delta = mean(r2) - mean(r1))
}

# Simulate Δ using μ_O1 for both O1 and O2, with FP added to O2 lines
simulate_delta_once <- function(mu_o1, col_index, fp_lambda) {
  n   <- min(nrow(dataset_o1), nrow(dataset_o2))
  g1  <- pmin(dataset_o1[[col_gen]], gcap)[1:n]
  g2  <- pmin(dataset_o2[[col_gen]], gcap)[1:n]
  c1  <- dataset_o1[[col_cov]][1:n]
  c2  <- dataset_o2[[col_cov]][1:n]
  
  # O1 true -> thin by coverage -> per-bp-per-gen rate
  # (Poisson thinning equivalent; use Binomial for clarity)
  rate_o1 <- numeric(n)
  rate_o2 <- numeric(n)
  
  for (i in seq_len(n)) {
    # O1
    raw1   <- rpois(1, g1[i] * mu_o1)
    kept1  <- rbinom(1, raw1, c1[i])
    rate_o1[i] <- kept1 / (c1[i] * g1[i] * ref)
    
    # O2 (same Î¼_O1), then add FP
    raw2   <- rpois(1, g2[i] * mu_o1)
    kept2  <- rbinom(1, raw2, c2[i])
    fp_i   <- rpois(1, fp_lambda)
    rate_o2[i] <- (kept2 + fp_i) / (c2[i] * g2[i] * ref)
  }
  
  mean(rate_o2) - mean(rate_o1)
}

simulate_delta <- function(mu_o1, col_index, fp_lambda, nsim = NSIM) {
  v <- numeric(nsim)
  for (k in seq_len(nsim)) v[k] <- simulate_delta_once(mu_o1, col_index, fp_lambda)
  v
}


# ================================
# Overlapping FP histograms
# ================================

make_fp_overlap <- function(mu_o1, col_index, label, outfile,
                            fps = 0:3, nsim = NSIM,
                            palette = c("#1b9e77","#d95f02","#7570b3","#e7298a")) {
  # observed Δ
  obs <- observed_delta(col_index)
  obs_delta <- as.numeric(obs["delta"])
  
  # simulate Δ for each FP
  sims <- lapply(fps, function(fp) simulate_delta(mu_o1, col_index, fp, nsim = nsim))
  names(sims) <- as.character(fps)
  
  # long df + stats
  df_long <- dplyr::bind_rows(
    lapply(names(sims), function(nm) {
      data.frame(diff = sims[[nm]], fp = factor(as.integer(nm), levels = fps))
    })
  )
  stats <- dplyr::bind_rows(
    lapply(names(sims), function(nm) {
      v <- sims[[nm]]
      tibble::tibble(fp = as.integer(nm), q95 = as.numeric(quantile(v, 0.95)))
    })
  )
  
  # axis range
  xr <- range(df_long$diff); xr <- xr + c(-1,1)*diff(xr)*0.05
  
  # math title: Δ = μ̄₂ − μ̄₁
  title_expr <- bquote(.(label) * ": " * Delta == bar(mu)[2] - bar(mu)[1])
  
  p <- ggplot(df_long, aes(x = diff)) +
    # overlapping histograms; keep the FP color legend
    geom_histogram(aes(fill = fp), bins = 50, position = "identity",
                   alpha = 0.35, color = NA, show.legend = TRUE) +
    scale_fill_manual(values = palette, name = "fp") +
    # observed Δ (red dashed) — legend entry
    geom_vline(aes(xintercept = obs_delta,
                   linetype = "Observed Δ", color = "Observed Δ"),
               linewidth = 1, show.legend = TRUE) +
    # per-FP 95th percentile lines COLORED BY FP (no legend for these)
    geom_vline(data = stats,
               aes(xintercept = q95, color = factor(fp, levels = fps)),
               linetype = "dotted", linewidth = 0.9, show.legend = FALSE) +
    # add a single gray dotted key to legend (no line drawn, just the key)
    geom_vline(aes(xintercept = -Inf,
                   linetype = "95th percentile", color = "95th percentile"),
               linewidth = 0.9, show.legend = TRUE) +
    scale_color_manual(
      values = c("Observed Δ" = "firebrick", "95th percentile" = "gray30"),
      breaks = c("Observed Δ", "95th percentile"),
      name = NULL
    ) +
    scale_linetype_manual(
      values = c("Observed Δ" = "dashed", "95th percentile" = "dotted"),
      breaks = c("Observed Δ", "95th percentile"),
      name = NULL
    ) +
    scale_x_continuous(limits = xr, labels = scales::scientific_format(digits = 2)) +
    labs(title = title_expr,
         x = "Simulated Δ (per-bp per-gen)",
         y = "Frequency") +
    guides(
      # remove the black border on FP color swatches
      fill = guide_legend(
        order = 1,
        override.aes = list(alpha = 0.6, colour = NA, color = NA, linetype = 0, linewidth = 0)
      ),
      # keep the two-item line legend
      color = guide_legend(order = 2),
      linetype = guide_legend(order = 2)
    ) +
    theme_bw()
  
  print(p)                       # show before saving
  ggsave(outfile, plot = p, width = 10, height = 6.5, dpi = 500)
  
  invisible(list(plot = p, stats = stats, obs_delta = obs_delta))
}






# ---- Make overlapping histograms for FP = 0,1,2,3 ----
ov_snp <- make_fp_overlap(mu_o1_snp, col_snp,   "SNV",   "SNP_FP0-3_overlap.png")
ov_ind <- make_fp_overlap(mu_o1_ind, col_indel, "INDEL", "INDEL_FP0-3_overlap.png")
ov_tot <- make_fp_overlap(mu_o1_tot, col_mut,   "TOTAL", "TOTAL_FP0-3_overlap.png")

##############
# =======================================
# Find FP cutoff where right-tail p = 0.05
# (Observed Δ equals 95th percentile)
# =======================================

# Evaluate p_right and q95 at a given FP lambda
eval_p_right <- function(mu_o1, col_index, fp_lambda, nsim = NSIM) {
  obs <- observed_delta(col_index)
  obs_delta <- as.numeric(obs["delta"])
  sims <- simulate_delta(mu_o1, col_index, fp_lambda, nsim = nsim)
  p_right <- mean(sims > obs_delta)
  q95 <- as.numeric(quantile(sims, 0.95))
  list(p_right = p_right, q95 = q95, obs_delta = obs_delta)
}

# Bracket and bisect to find the smallest FP with p_right >= target (default 0.05)
find_fp_cutoff <- function(mu_o1, col_index, target = 0.05,
                           nsim_bracket = NSIM, nsim_bisect = NSIM,
                           fp_lo = 0, fp_hi_init = 1, fp_hi_max = 50,
                           tol = 0.01, max_iter = 40, label = "") {
  
  # Evaluate at lo
  lo <- eval_p_right(mu_o1, col_index, fp_lo, nsim = nsim_bracket)
  message(sprintf("[%s] FP=%.3f  p_right=%.6f  q95=%s  obsΔ=%s",
                  label, fp_lo, lo$p_right,
                  format(lo$q95, scientific = TRUE, digits = 2),
                  format(lo$obs_delta, scientific = TRUE, digits = 2)))
  
  # If already above target at FP=0, cutoff is ~0
  if (lo$p_right >= target) {
    return(list(fp_cutoff = 0, p_right = lo$p_right, bracket = c(0, 0),
                note = "p_right(0) >= target; cutoff at 0"))
  }
  
  # Expand high bound until p_right >= target or hit fp_hi_max
  fp_hi <- fp_hi_init
  hi <- eval_p_right(mu_o1, col_index, fp_hi, nsim = nsim_bracket)
  while (hi$p_right < target && fp_hi < fp_hi_max) {
    message(sprintf("[%s] Expanding: FP=%.3f  p_right=%.6f", label, fp_hi, hi$p_right))
    fp_hi <- min(fp_hi * 2, fp_hi_max)
    hi <- eval_p_right(mu_o1, col_index, fp_hi, nsim = nsim_bracket)
  }
  
  if (hi$p_right < target) {
    return(list(fp_cutoff = NA_real_, p_right = hi$p_right, bracket = c(fp_lo, fp_hi),
                note = "Did not reach target within fp_hi_max"))
  }
  
  message(sprintf("[%s] Bracketed: [%.3f, %.3f] with p_right(lo)=%.4f, p_right(hi)=%.4f",
                  label, fp_lo, fp_hi, lo$p_right, hi$p_right))
  
  # Bisection
  iter <- 0
  while ((fp_hi - fp_lo) > tol && iter < max_iter) {
    fp_mid <- 0.5 * (fp_lo + fp_hi)
    mid <- eval_p_right(mu_o1, col_index, fp_mid, nsim = nsim_bisect)
    iter <- iter + 1
    message(sprintf("[%s] Iter %02d: FP=%.4f  p_right=%.6f",
                    label, iter, fp_mid, mid$p_right))
    if (mid$p_right >= target) {
      fp_hi <- fp_mid
      hi <- mid
    } else {
      fp_lo <- fp_mid
      lo <- mid
    }
  }
  
  fp_star <- fp_hi
  list(fp_cutoff = fp_star, p_right = hi$p_right,
       bracket = c(fp_lo, fp_hi),
       iters = iter,
       note = sprintf("Converged with tol=%.3f (smallest FP s.t. p_right>=%.2f)", tol, target))
}

# ---------- Run it for SNP / INDEL / TOTAL ----------
# Make sure μ_O1 is available
mu_o1_snp <- mu_from_o1(dataset_o1, col_snp)
mu_o1_ind <- mu_from_o1(dataset_o1, col_indel)
mu_o1_tot <- mu_from_o1(dataset_o1, col_mut)

# You can bump nsim_* for more stable estimates (e.g., 20000–50000)
cut_snp <- find_fp_cutoff(mu_o1_snp, col_snp,   target = 0.05,
                          nsim_bracket = 15000, nsim_bisect = 30000,
                          fp_lo = 0, fp_hi_init = 0.5, fp_hi_max = 10,
                          tol = 0.01, max_iter = 40, label = "SNP")
cut_ind <- find_fp_cutoff(mu_o1_ind, col_indel, target = 0.05,
                          nsim_bracket = 15000, nsim_bisect = 30000,
                          fp_lo = 0, fp_hi_init = 1.0, fp_hi_max = 20,
                          tol = 0.01, max_iter = 40, label = "INDEL")
cut_tot <- find_fp_cutoff(mu_o1_tot, col_mut,   target = 0.05,
                          nsim_bracket = 15000, nsim_bisect = 30000,
                          fp_lo = 0, fp_hi_init = 0.5, fp_hi_max = 20,
                          tol = 0.01, max_iter = 40, label = "TOTAL")

print(list(SNP = cut_snp, INDEL = cut_ind, TOTAL = cut_tot))

