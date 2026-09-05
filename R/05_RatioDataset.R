# =============================================================================
# 05_RatioDataset.R
# County-level ratio of iNaturalist to VertNet records per species.
# Shows the proportion of records from each dataset within each county,
# ranked by iNaturalist proportion.
#
# DEPENDENCIES: run 01_DataPrep.R first
#
# OUTPUT: Figure_ObsRatio_VertNet_iNat_by_county.pdf / .png
#         [Figure 3]
# =============================================================================

library(tidyverse)
library(here)          # project-root-relative paths
library(grid)
library(gridExtra)
library(cowplot)
library(conflicted)
library(rphylopic)    
conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(dplyr::select)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

DATA_DIR <- here::here("data")
PREP_DIR <- here::here("data", "prepared")
FIG_DIR  <- here::here("figures")

SPECIES_LEVELS <- readRDS(file.path(PREP_DIR, "SPECIES_LEVELS.rds"))

COLOURS <- c("iNat_prop"    = "#D55E00",   
             "VertNet_prop" = "#009E73")   

paper_theme <- theme_bw(base_size = 22) +
  theme(
    axis.title       = element_text(size = 26),
    axis.text        = element_text(size = 20, color = "black"),
    strip.text       = element_text(size = 22, face = "italic"),
    legend.text      = element_text(size = 20),
    legend.title     = element_text(size = 22),
    plot.title       = element_text(size = 28, face = "bold"),
    plot.subtitle    = element_text(size = 22),
    panel.grid.minor = element_blank(),
    panel.spacing    = unit(1, "lines")
  )
# =============================================================================
# Reusable phylopic placement
# =============================================================================
PHYLOPIC_CFG <- list(
  x_frac  = 0.20,   # fraction across x-axis (county-rank) range
  y_frac  = 0.85,   # fraction up the fixed 0-1 y range
  h_frac  = 0.16,   # silhouette height as fraction of the fixed 0-1 y range
  x_hard  = NULL,   # NULL -> use x_frac (x range varies per species)
  y_min   = 0,      # fixed y range for this plot (proportion is always 0-1)
  y_max   = 1
)

phylopic_lookup <- tibble(
  Species       = SPECIES_LEVELS,
  phylopic_name = c("Bufo bufo", "Bufo bufo", "Bufo bufo",
                    "Leiocephalus carinatus", "Leiocephalus carinatus", "Leiocephalus carinatus", 
                    "Rhinoclemmys punctularia", "Rhinoclemmys punctularia")
)

mirrored_species <- c("A. punctatus", "A. woodhousii", "I. nebulifer",
                      "Co. texanus", "Cr. collaris", "P. cornutum",
                      "Te. ornata", "Tr. scripta")

# Build placement table from whatever x/y columns a given plot uses.
# `df` is the plotted data; xvar drives the x range per species. The y range
# is taken from cfg$y_min/cfg$y_max when both are set (fixed-scale plots like
# this one), otherwise derived from yvar.
build_phylopic_positions <- function(df, xvar, yvar = NULL, lookup = phylopic_lookup,
                                     cfg = PHYLOPIC_CFG,
                                     mirrored = mirrored_species,
                                     species_levels = SPECIES_LEVELS) {
  fixed_y <- !is.null(cfg$y_min) && !is.null(cfg$y_max)
  
  ranges <- df %>%
    group_by(Species) %>%
    summarise(
      xmin = min(.data[[xvar]], na.rm = TRUE),
      xmax = max(.data[[xvar]], na.rm = TRUE),
      ymin = if (fixed_y) cfg$y_min else min(.data[[yvar]], na.rm = TRUE),
      ymax = if (fixed_y) cfg$y_max else max(.data[[yvar]], na.rm = TRUE),
      .groups = "drop"
    )
  
  lookup %>%
    left_join(ranges, by = "Species") %>%
    mutate(
      Species    = factor(Species, levels = species_levels),
      x = if (!is.null(cfg$x_hard)) cfg$x_hard
      else xmin + cfg$x_frac * (xmax - xmin),
      y = ymin + cfg$y_frac * (ymax - ymin),
      height     = cfg$h_frac * (ymax - ymin),
      horizontal = Species %in% mirrored
    )
}

add_phylopic_auto <- function(p, positions) {
  for (i in seq_len(nrow(positions))) {
    row <- positions[i, ]
    ann <- data.frame(
      x = row$x, y = row$y, name = row$phylopic_name
    )
    p <- p + geom_phylopic(
      data = ann, aes(x = x, y = y, name = name),
      inherit.aes = FALSE, color = "gray20",
      height = row$height, horizontal = row$horizontal
    )
  }
  p
}

# =============================================================================
# 1. Per-species summary: % counties dominated by each dataset
# =============================================================================
urbanization_raw <- readRDS(file.path(PREP_DIR, "urbanization_raw.rds"))

calc_summary <- function(species_name) {
  urbanization_raw %>%
    filter(Species == species_name, !is.na(county)) %>%
    group_by(Species, Dataset, county) %>%
    count() %>%
    pivot_wider(names_from = Dataset, values_from = n,
                values_fill = 0) %>%
    mutate(
      sp.n         = iNaturalist + VertNet,
      iNat.prop    = iNaturalist / sp.n,
      VertNet.prop = VertNet / sp.n
    ) %>%
    ungroup() %>%
    summarise(
      pct_iNat_dominant    = sum(iNat.prop > 0.5)    / n() * 100,
      pct_VertNet_dominant = sum(VertNet.prop > 0.5) / n() * 100,
      pct_equal            = sum(iNat.prop == 0.5)   / n() * 100
    ) %>%
    mutate(Species = species_name)
}

summaries <- map_dfr(SPECIES_LEVELS, calc_summary)

# =============================================================================
# 2. Build one stacked bar plot per species
# =============================================================================

make_ratio_plot <- function(sp, summaries, urbanization_raw) {
  sp_summary <- summaries %>% filter(Species == sp)
  
  # County-rank data for this species; `rank` becomes the x aesthetic driver.
  ranked <- urbanization_raw %>%
    filter(Species == sp, !is.na(county)) %>%
    group_by(county, Dataset, Species) %>%
    count() %>%
    pivot_wider(names_from = Dataset, values_from = n,
                values_fill = 0) %>%
    mutate(
      total        = iNaturalist + VertNet,
      iNat_prop    = iNaturalist / total,
      VertNet_prop = VertNet / total
    ) %>%
    arrange(desc(iNat_prop)) %>%
    tibble::rowid_to_column(var = "rank")
  
  n_counties <- nrow(ranked)
  
  p <- ranked %>%
    pivot_longer(cols = c(iNat_prop, VertNet_prop),
                 names_to  = "Dataset",
                 values_to = "proportion") %>%
    ggplot(aes(fill = Dataset, y = proportion,
               x = reorder(county, -rank))) +
    geom_bar(position = "fill", stat = "identity", width = 1) +
    geom_hline(yintercept = 0.5, color = "gray99",
               linetype = "dashed", linewidth = 0.4) +
    scale_fill_manual(values = COLOURS) +
    scale_y_continuous(expand = c(0, 0)) +
    paper_theme +
    theme(
      aspect.ratio     = 1,
      legend.position  = "none",
      axis.text.x      = element_blank(),
      axis.ticks.x     = element_blank(),
      axis.text.y      = element_text(size = 14),
      plot.title = ggtext::element_markdown(
        hjust = 0.5, size = 18, lineheight = 0.9, , face = "plain"
      ),
      plot.margin      = margin(0, 0, 0, 0),
    ) +
    annotate("text",
             x = n_counties * 0.95, y = 0.62,
             label    = sprintf("%.0f%%", sp_summary$pct_iNat_dominant),
             hjust    = 1, size = 6, fontface = "bold", color = "white") +
    annotate("text",
             x = n_counties * 0.05, y = 0.38,
             label    = sprintf("%.0f%%", sp_summary$pct_VertNet_dominant),
             hjust    = -0.25, size = 6, fontface = "bold", color = "white") +
    labs(
      title = sprintf("*%s*<br>(*n* = %d)", sp, n_counties),
      y = "", x = ""
    )
  
  # Phylopic silhouette via the shared workflow. x is driven by county rank
  # (per-species range); y is the fixed 0-1 proportion set in PHYLOPIC_CFG.
#  pos <- build_phylopic_positions(ranked, xvar = "rank") %>%
#    filter(Species == sp)
#  add_phylopic_auto(p, pos)
# 
  }

plot_list <- map(SPECIES_LEVELS, make_ratio_plot,
                 summaries   = summaries,
                 urbanization_raw = urbanization_raw) %>%
  set_names(SPECIES_LEVELS)

# =============================================================================
# 3. Combine panels and add shared y-axis label [Figure 3]
# =============================================================================
combo <- cowplot::plot_grid(
  plotlist = plot_list,
  ncol     = 4,
  align    = "hv"
)

y_label <- textGrob(
  "Ratio of VertNet \nto iNaturalist records",
  gp  = gpar(fontface = "plain", fontsize = 18),
  rot = 90
)

figure_ratio <- grid.arrange(arrangeGrob(combo, left = y_label))

ggsave(file.path(FIG_DIR, "fig3_ObsRatio_VertNet_iNat_by_county.pdf"), figure_ratio, width = 12, height = 5.5, units = "in", dpi = 300)
ggsave(file.path(FIG_DIR, "fig3_ObsRatio_VertNet_iNat_by_county.png"), figure_ratio, width = 12, height = 5.5, units = "in", dpi = 300)
