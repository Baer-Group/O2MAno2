suppressPackageStartupMessages({
  library(rstudioapi)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(forcats)
  library(scales)
  library(readr)
  library(zoo)
})

# --- optional: set working directory to this script's folder (RStudio only)
try({ setwd(dirname(getActiveDocumentContext()$path)) }, silent = TRUE)

# =========================
# O1 vs O2 variant locations per chromosome
#   - E and G combined within each generation
#   - 2 panels per chromosome (SNP on top, Indel on bottom)
#   - O1 = #0070C0, O2 = #E46C0A
# =========================

# ---- Inputs (adjust pathnames if needed)
files <- list(
  snp_O1   = c("E.snp.3x_O1.csv",   "G.snp.3x_O1.csv"),
  snp_O2   = c("E.snp.3x_O2.csv",   "G.snp.3x_O2.csv"),
  indel_O1 = c("E.indel.3x_O1.csv", "G.indel.3x_O1.csv"),
  indel_O2 = c("E.indel.3x_O2.csv", "G.indel.3x_O2.csv")
)

# ---------- TUNABLES ----------
bin_size     <- 250000L       # was 100k; try 250k or 500k for smoother lines
smooth_bins  <- 3L            # rolling-mean over this many bins (set 1 to disable)
show_points  <- FALSE         # TRUE=points; FALSE=rug marks
free_y_scales <- TRUE         # TRUE gives each panel its own y-scale (declutters)

# x ticks every N kb (e.g., 5000 = 5 Mb)
x_tick_kb <- 5000

# ---------- rebuild counts with optional smoothing ----------


counts_df <- points_df %>%
  group_by(Chromosome, Type, Gen, Bin) %>%
  summarise(Count = n(), .groups = "drop") %>%
  arrange(Chromosome, Type, Gen, Bin) %>%
  group_by(Chromosome, Type, Gen) %>%
  mutate(Count_smooth = if (smooth_bins > 1) 
    zoo::rollmean(Count, k = smooth_bins, fill = NA, align = "center") 
    else Count) %>%
  ungroup()

# attach smoothed y for dots/rugs
points_df <- points_df %>%
  left_join(counts_df %>% select(Chromosome, Type, Gen, Bin, Count_smooth),
            by = c("Chromosome","Type","Gen","Bin"))


#######################################Bar #######################
# --- rebuild counts (no smoothing) ---
counts_df_bar <- points_df %>%
  dplyr::group_by(Chromosome, Type, Gen, Bin) %>%
  dplyr::summarise(Count = dplyr::n(), .groups = "drop") %>%
  dplyr::mutate(
    Chromosome = forcats::fct_relevel(Chromosome, c("I","II","III","IV","V","X")),
    Type = factor(Type, levels = c("SNP","Indel")),
    Gen  = factor(Gen,  levels = c("O1","O2"))
  )

# --- bar plot params ---
bar_width_frac <- 0.80   # fraction of bin width used by each bar group
x_tick_mb      <- 5      # tick every 5 Mb (adjust to taste)

p_bar <- ggplot(counts_df_bar, aes(x = Bin, y = Count, fill = Gen)) +
  # one column per generation, dodged within each bin
  geom_col(
    width = bar_width_frac * bin_size,
    position = position_dodge2(
      width = bin_size, preserve = "single",
      padding = (1 - bar_width_frac) / 2
    ),
    alpha = 0.95
  ) +
  facet_grid(rows = vars(Type), cols = vars(Chromosome),
             scales = if (free_y_scales) "free_y" else "fixed",
             switch = "y") +
  scale_fill_manual(
    values = c(O1 = "#0070C0", O2 = "#E46C0A"),
    guide = guide_legend(title = NULL, override.aes = list(alpha = 1))
  ) +
  # x in Mb
  scale_x_continuous(
    breaks = scales::breaks_width(x_tick_mb * 1e6),
    labels = scales::label_number(scale = 1e-6, suffix = " Mb", big.mark = ",")
  ) +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(x = "Genomic position (Mb; binned)", y = NULL) +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid.minor = element_blank(),
    strip.placement = "outside",
    strip.text.x = element_text(face = "bold"),
    strip.text.y.left = element_text(angle = 90, face = "bold", size = 18),
    legend.position = "top",
    legend.text = element_text(size = 16)
  )

print(p_bar)
ggsave("chromosome_panels_bars_O1_vs_O2.png", p_bar, width = 16, height = 9, dpi = 500)



################################## bar #########################


# ---------- plot ----------
# x ticks every N Mb
x_tick_mb <- 5  # try 2, 5, or 10

p <- ggplot() +
  geom_line(
    data = counts_df,
    aes(x = Bin, y = Count_smooth, color = Gen, group = interaction(Chromosome, Type, Gen)),
    linewidth = 0.7, alpha = 0.95, na.rm = TRUE
  ) +
  { if (show_points) {
    geom_point(
      data = points_df,
      aes(x = Bin, y = Count_smooth, color = Gen),
      size = 0.9, alpha = 0.25,
      position = position_jitter(width = bin_size * 0.15, height = 0, seed = 42),
      na.rm = TRUE
    )
  } else {
    geom_rug(
      data = points_df,
      aes(x = Bin, color = Gen),
      sides = "b", alpha = 0.15, linewidth = 0.3, na.rm = TRUE
    )
  }
  } +
  # ← back to 2-row grid (SNP top, Indel bottom)
  facet_grid(rows = vars(Type), cols = vars(Chromosome),
             scales = if (free_y_scales) "free_y" else "fixed",
             switch = "y") +
  scale_color_manual(values = c(O1 = "#0070C0", O2 = "#E46C0A"),
                     guide = guide_legend(title = NULL, override.aes = list(alpha = 1))) +
  # x in Mb
  scale_x_continuous(
    breaks = scales::breaks_width(x_tick_mb * 1e6),
    labels = scales::label_number(scale = 1e-6, big.mark = ",")
  ) +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(x = "Genomic position (in Mb; binned)", y = NULL) +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid.minor = element_blank(),
    strip.placement = "outside",
    strip.text.x = element_text(face = "bold", size = 20 ),
    # ↓ vertical strip labels for SNP/Indel
    strip.text.y.left = element_text(angle = 90, face = "bold", size =20),
    legend.position = "top",
    # Legend text ("O1", "O2")
    legend.text  = element_text(size = 20),                                  # <- bump size here
    legend.title = element_text(size = 20)   
  )


# --- Recompute counts with O2 halved ---
library(zoo)

counts_df_adj <- points_df %>%
  dplyr::group_by(Chromosome, Type, Gen, Bin) %>%
  dplyr::summarise(Count = dplyr::n(), .groups = "drop") %>%
  dplyr::mutate(
    Count_adj = dplyr::if_else(Gen == "O2", Count * 0.5, Count)  # <-- halve O2
  ) %>%
  dplyr::arrange(Chromosome, Type, Gen, Bin) %>%
  dplyr::group_by(Chromosome, Type, Gen) %>%
  dplyr::mutate(
    Count_smooth = if (smooth_bins > 1)
      zoo::rollmean(Count_adj, k = smooth_bins, fill = NA, align = "center")
    else Count_adj
  ) %>%
  dplyr::ungroup()

# Attach adjusted/smoothed y to points so dots/rugs sit on the adjusted line
points_df_adj <- points_df %>%
  dplyr::left_join(
    counts_df_adj %>% dplyr::select(Chromosome, Type, Gen, Bin, Count_smooth),
    by = c("Chromosome","Type","Gen","Bin")
  )

# --- Plot (same style as before) ---
p_halfO2 <- ggplot() +
  geom_line(
    data = counts_df_adj,
    aes(x = Bin, y = Count_smooth, color = Gen,
        group = interaction(Chromosome, Type, Gen)),
    linewidth = 0.7, alpha = 0.95, na.rm = TRUE
  ) +
  {
    if (show_points) {
      geom_point(
        data = points_df_adj,
        aes(x = Bin, y = Count_smooth, color = Gen),
        size = 0.9, alpha = 0.25,
        position = position_jitter(width = bin_size * 0.15, height = 0, seed = 42),
        na.rm = TRUE
      )
    } else {
      geom_rug(
        data = points_df_adj,
        aes(x = Bin, color = Gen),
        sides = "b", alpha = 0.15, linewidth = 0.3, na.rm = TRUE
      )
    }
  } +
  facet_grid(rows = vars(Type), cols = vars(Chromosome),
             scales = if (free_y_scales) "free_y" else "fixed", switch = "y") +
  scale_color_manual(name = NULL, values = c(O1 = "#0070C0", O2 = "#E46C0A"),
                     guide = guide_legend(title = NULL, override.aes = list(alpha = 1))) +
  scale_x_continuous(
    breaks = scales::breaks_width(x_tick_mb * 1e6),
    labels = scales::label_number(scale = 1e-6, big.mark = ",")
  ) +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(
    x = "Genomic position (in Mb; binned)",
    y = NULL
  ) +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid.minor = element_blank(),
    strip.placement  = "outside",
    strip.text.x     = element_text(face = "bold"),
    strip.text.y.left= element_text(angle = 90, face = "bold", size = 20), # your larger facet labels
    legend.position  = "top",
    legend.text      = element_text(size = 20)
  ) +
  guides(color = guide_legend(override.aes = list(size = 2.5)))

print(p_halfO2)
ggsave("chromosome_panels_O1_vs_O2_O2halved_Mb.png", p_halfO2, width = 16, height = 9, dpi = 500)




suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
  library(forcats); library(scales)
})
library(gridExtra) # For arranging multiple plots

# ---- order like your screenshot (top→bottom)
chr_order <- c("X","V","IV","III","II","I")

# helper: make per-chromosome bin indices, fill missing bins with 0s
prep_counts_for_heat <- function(points_df, adjust_O2 = FALSE) {
  base <- points_df %>%
    mutate(
      Chromosome = forcats::fct_relevel(Chromosome, c("X","V","IV","III","II","I")),
      Type = factor(Type, levels = c("SNP","Indel")),
      Gen  = factor(Gen,  levels = c("O1","O2"))
    ) %>%
    dplyr::group_by(Chromosome, Type, Gen, Bin) %>%
    dplyr::summarise(Count = dplyr::n(), .groups = "drop")
  
  if (adjust_O2) base <- base %>% dplyr::mutate(Count = dplyr::if_else(Gen == "O2", Count * 0.5, Count))
  
  base %>%
    dplyr::group_by(Chromosome) %>%
    dplyr::mutate(min_bin = min(Bin), max_bin = max(Bin)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      Bin_idx = as.integer((Bin - min_bin) / bin_size),
      Pos_Mb  = (Bin - min_bin) / 1e6
    ) %>%
    dplyr::group_by(Type, Gen, Chromosome) %>%
    tidyr::complete(
      Bin_idx = tidyr::full_seq(min(Bin_idx):max(Bin_idx), 1),
      fill = list(Count = 0)
    ) %>%
    dplyr::mutate(Pos_Mb = Bin_idx * (bin_size / 1e6)) %>%  # keep Mb aligned to bins
    dplyr::ungroup()
}


# plotting function (looks like your screenshot)
make_species_style_heatmap <- function(counts_df, type_pick = "SNP",
                                       title_str = "Distribution of SNP Hotspots",
                                       x_break_mb = 5) {
  df <- counts_df %>% dplyr::filter(Type == type_pick)
  
  max_mb <- ceiling(max(df$Pos_Mb, na.rm = TRUE))
  breaks <- seq(0, max_mb, by = x_break_mb)
  
  ggplot(df, aes(x = Pos_Mb, y = Chromosome, fill = Count)) +
    geom_tile(width = 0.95, height = 0.9, linewidth = 0.15, color = "grey90") +
    facet_grid(rows = vars(Gen), switch = "y") +
    scale_x_continuous(breaks = breaks, labels = scales::label_number(accuracy = 1), expand = c(0.002, 0.02)) +
    labs(x = "Genomic position (Mb)", y = "Chromosome", title = title_str) +
    scale_fill_viridis_c(name = "Count", option = "magma",
                         trans = "sqrt", breaks = scales::pretty_breaks(5)) +
    theme_minimal(base_size = 18) +
    theme(
      panel.grid = element_blank(),
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0, face = "bold", size = 16),
      legend.position = "right"
    )
}


# ---------- Build datasets ----------
counts_exact <- prep_counts_for_heat(points_df, adjust_O2 = FALSE)
counts_adj   <- prep_counts_for_heat(points_df, adjust_O2 = TRUE)

# ---------- FIGURES ----------
# SNP (exact)
p_snp_exact <- make_species_style_heatmap(counts_exact, type_pick = "SNP",
                                          title_str = "SNP Hotspots(Exact)")
ggsave("heatmap_SNP_exact.png", p_snp_exact, width = 16, height = 9, dpi = 500)

# SNP (O2 halved)
p_snp_adj <- make_species_style_heatmap(counts_adj, type_pick = "SNP",
                                        title_str = "SNP Hotspots")
ggsave("heatmap_SNP_adjusted.png", p_snp_adj, width = 16, height = 9, dpi = 500)

# INDEL (exact)
p_indel_exact <- make_species_style_heatmap(counts_exact, type_pick = "Indel",
                                            title_str = "Indel Hotspots(Exact)")
ggsave("heatmap_Indel_exact.png", p_indel_exact, width = 16, height = 9, dpi = 500)

# INDEL (O2 halved)
p_indel_adj <- make_species_style_heatmap(counts_adj, type_pick = "Indel",
                                          title_str = "Indel Hotspots")
ggsave("heatmap_Indel_adjusted.png", p_indel_adj, width = 16, height = 9, dpi = 500)

# Arrange the plots into a 1x4 grid
compound_plot1 <- grid.arrange(p_snp_exact, p_indel_exact, nrow = 1, ncol = 2)
# Arrange the plots into a 1x4 grid
compound_plot2 <- grid.arrange(p_snp_adj, p_indel_adj, nrow = 1, ncol = 2)
