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


# Make sure μ_O1 is available
mu_o1_snp <- mu_from_o1(dataset_o1, col_snp)
mu_o1_ind <- mu_from_o1(dataset_o1, col_indel)
mu_o1_tot <- mu_from_o1(dataset_o1, col_mut)


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

make_fp_overlap <- function(mu_o1,
                            col_index,
                            title_str = "[A]",
                            outfile   = "FP_overlap.png",
                            fps       = 0:3,
                            nsim      = NULL,
                            base_color = "#1f77b4",
                            # new args:
                            fp_cutoff = NA_real_,     # numeric; from find_fp_cutoff(...)
                            fp_q95    = NA_real_,     # numeric; mean FP to compute a single q95 line
                            shade_from = 0.95, shade_to = 0.05) {
  
  # ----- local helpers (single-hue light→dark palette) -----
  mix_col <- function(col, bg = "#FFFFFF", f = 0.5) {
    c1 <- grDevices::col2rgb(col); c2 <- grDevices::col2rgb(bg)
    m  <- (1 - f) * c1 + f * c2
    grDevices::rgb(m[1,]/255, m[2,]/255, m[3,]/255)
  }
  gen_shades <- function(base, n, from = 0.70, to = 0.10) {
    fs <- seq(from, to, length.out = n)
    vapply(fs, function(f) mix_col(base, "#FFFFFF", f), character(1))
  }
  
  # ----- nsim default -----
  if (is.null(nsim)) nsim <- if (!is.null(get0("NSIM", inherits = TRUE))) get0("NSIM", inherits = TRUE) else 10000
  
  # ----- observed Δ -----
  obs <- observed_delta(col_index)
  obs_delta <- as.numeric(obs["delta"])
  
  # ----- simulate Δ for each FP (for the overlay histograms) -----
  sims <- lapply(fps, function(fp) simulate_delta(mu_o1, col_index, fp, nsim = nsim))
  names(sims) <- as.character(fps)
  
  df_long <- dplyr::bind_rows(
    lapply(names(sims), function(nm) {
      data.frame(diff = sims[[nm]], fp = factor(as.integer(nm), levels = fps))
    })
  )
  
  # ----- single q95 line computed at fp_q95 (if provided) -----
  q95_meanFP <- NA_real_
  if (!is.na(fp_q95)) {
    v_mid <- simulate_delta(mu_o1, col_index, fp_q95, nsim = nsim)
    q95_meanFP <- as.numeric(stats::quantile(v_mid, 0.95))
  }
  
  # ----- palette & plot limits -----
  palette <- gen_shades(base_color, length(fps), from = shade_from, to = shade_to)
  xr <- range(df_long$diff); xr <- xr + c(-1, 1) * diff(xr) * 0.05
  
  # rough y max for placing an annotation near the baseline
  y_max <- max(hist(df_long$diff, breaks = 50, plot = FALSE)$counts)
  ann_y <- y_max * 0.06  # a bit above the x-axis
  
  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = diff)) +
    ggplot2::geom_histogram(ggplot2::aes(fill = fp), bins = 50, position = "identity",
                            alpha = 0.90, color = NA, show.legend = TRUE) +
    ggplot2::scale_fill_manual(values = palette, name = "FP") +
    
    # Observed Δ (red dashed) — legend entry
    ggplot2::geom_vline(ggplot2::aes(xintercept = obs_delta,
                                     linetype = "Observed Δ", color = "Observed Δ"),
                        linewidth = 1, show.legend = TRUE) +
    
    # Single q95 line at the mean FP (gray dotted) — legend entry
    {if (!is.na(q95_meanFP))
      ggplot2::geom_vline(ggplot2::aes(xintercept = q95_meanFP,
                                       linetype = "mean FP(95p)",
                                       color = "mean FP(95p)"),
                          linewidth = 0.9, show.legend = TRUE)
    } +
    
    # annotate FP cutoff under the observed Δ line
    {if (!is.na(fp_cutoff))
      ggplot2::annotate("text",
                        x = obs_delta, y = ann_y,
                        label = paste0("FP cutoff = ", format(fp_cutoff, digits = 3)),
                        hjust = 0.5, vjust = 0.2, color = "firebrick", size = 5)
    } +
    
    ggplot2::scale_color_manual(values = c("Observed Δ" = "firebrick",
                                           "mean FP(95p)" = "gray30"),
                                breaks = c("Observed Δ", "mean FP(95p)"),
                                name = NULL) +
    ggplot2::scale_linetype_manual(values = c("Observed Δ" = "dashed",
                                              "mean FP(95p)" = "dotted"),
                                   breaks = c("Observed Δ", "mean FP(95p)"),
                                   name = NULL) +
    
    ggplot2::scale_x_continuous(limits = xr,
                                labels = scales::scientific_format(digits = 2)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.12, 0.05))) +
    ggplot2::labs(title = title_str,
                  x = "Simulated Δ (per-bp per-gen)",
                  y = "Frequency") +
    ggplot2::guides(
      fill     = ggplot2::guide_legend(order = 1, override.aes = list(alpha = 1)),
      color    = ggplot2::guide_legend(order = 2, override.aes = list(fill = NA)),
      linetype = ggplot2::guide_legend(order = 2, override.aes = list(fill = NA))
    ) +
    ggplot2::theme_bw(base_size = 18) +
    ggplot2::theme(
      legend.key        = ggplot2::element_rect(fill = "white", colour = NA),
      legend.background = ggplot2::element_rect(fill = "white", colour = "grey85"),
      text         = ggplot2::element_text(size = 18),
      axis.title   = ggplot2::element_text(size = 20),
      axis.text    = ggplot2::element_text(size = 18),
      legend.text  = ggplot2::element_text(size = 20),
      legend.title = ggplot2::element_text(size = 20),
      plot.title   = ggplot2::element_text(size = 22, face = "bold")
    )
  
  print(p)
  ggplot2::ggsave(outfile, plot = p, width = 10, height = 6.5, dpi = 500)
  invisible(list(plot = p, q95_meanFP = q95_meanFP, obs_delta = obs_delta))
}


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


# cutoffs (as before)
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


# base colors (pick any you like)
col_snp_base <- "turquoise"
col_ind_base <- "slateblue2"
col_tot_base <- "tan"

# mean FP per panel for the single q95 line
fp_q95_snp <- 2.00
fp_q95_ind <- 0.33
fp_q95_tot <- fp_q95_snp + fp_q95_ind  # 2.33

ov_snp <- make_fp_overlap(mu_o1_snp, col_snp,   title_str = "[A]",
                          outfile   = "SNP_FP0-3_overlap.png",
                          fps = 0:3, base_color = col_snp_base,
                          fp_cutoff = cut_snp$fp_cutoff, fp_q95 = fp_q95_snp)

ov_ind <- make_fp_overlap(mu_o1_ind, col_indel, title_str = "[B]",
                          outfile   = "INDEL_FP0-3_overlap.png",
                          fps = 0:3, base_color = col_ind_base,
                          fp_cutoff = cut_ind$fp_cutoff, fp_q95 = fp_q95_ind)

ov_tot <- make_fp_overlap(mu_o1_tot, col_mut,   title_str = "[C]",
                          outfile   = "TOTAL_FP0-3_overlap.png",
                          fps = 0:3, base_color = col_tot_base,
                          fp_cutoff = cut_tot$fp_cutoff, fp_q95 = fp_q95_tot)



