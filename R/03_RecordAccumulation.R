# =============================================================================
# 03_RecordAccumulation.R
# Temporal coverage bias: record accumulation curves over time.
#
# DEPENDENCIES: run 01_DataPrep.R first or load temporal_data.
#
# Expected columns in temporal_data:
#   Species          (character, e.g. "A. punctatus")
#   Dataset          (character: "VertNet" | "iNaturalist")
#   year             (integer)
#   RecordCount      (integer — annual record count per species x dataset)
#   UniqueCollectors (integer — annual unique observer count per species x dataset)
# 
# OUTPUT: record_accumulation.pdf / .png
#         [Figure 1]
# =============================================================================

library(tidyverse)     # dplyr, ggplot2, tidyr, purrr, stringr
library(here)          # project-root-relative paths
library(mgcv)          # gam, nb()
library(patchwork)     # wrap_plots (if needed)
library(cowplot)       # plot_grid
library(sf)            # only if spatial plots added later
library(rphylopic)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

DATA_DIR <- here::here("data")
PREP_DIR <- here::here("data", "prepared")
FIG_DIR  <- here::here("figures")

temporal_data  <- readRDS(file.path(PREP_DIR, "temporal_data.rds"))
SPECIES_LEVELS <- readRDS(file.path(PREP_DIR, "SPECIES_LEVELS.rds"))

# Paper theme
paper_theme <- theme_bw(base_size = 16) +
  theme(
    axis.title       = element_text(size = 16),
    axis.text        = element_text(size = 13, color = "black"),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    strip.text       = element_text(size = 14, face = "italic"),
    legend.text      = element_text(size = 14),
    legend.title     = element_text(size = 15),
    legend.position  = "none",
    plot.title       = element_text(size = 18, face = "bold"),
    plot.subtitle    = element_text(size = 15),
    panel.grid.minor = element_blank(),
    panel.spacing    = unit(1, "lines"),
  )

# Plot colours
COLOURS <- c("iNaturalist" = "#D55E00", "VertNet" = "#009E73")

# Italic facet labels in canonical order: returns a factor whose levels are
# plotmath italic('...') strings, so label_parsed renders italics while the
# facet order stays locked to SPECIES_LEVELS.
ITALIC_LEVELS <- paste0("italic('", SPECIES_LEVELS, "')")
italic_species <- function(x) {
  factor(paste0("italic('", as.character(x), "')"), levels = ITALIC_LEVELS)
}

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

# Fit a negative-binomial GAM for a given response column
fit_nb_gam <- function(data, response_col) {
  f <- as.formula(paste0(response_col, " ~ s(year, by = fDataset, bs = 'cr')"))
  gam(f, family = nb(), data = data)
}

# Generate predictions from a fitted GAM, trimming iNat to observed year range
predict_nb_gam <- function(model, data, species_name) {
  newdat <- with(data, expand.grid(
    year     = seq(min(year), max(year), by = 1),
    fDataset = unique(fDataset)
  ))
  pred          <- predict(model, newdata = newdat, type = "link", se.fit = TRUE)
  newdat$fit    <- pred$fit
  newdat$se.fit <- pred$se.fit
  newdat$Species <- species_name
  
  min_inat <- data %>% filter(fDataset == "iNaturalist") %>%
    summarise(min(year)) %>% pull()
  newdat %>% filter(!(fDataset == "iNaturalist" & year < min_inat))
}

# Extract smooth-term edf and p-value for annotation
get_gam_stats <- function(model, species_name) {
  s <- summary(model)$s.table
  data.frame(
    Species = species_name,
    edf     = round(s[1, "edf"], 2),
    pval    = signif(s[1, "p-value"], 2)
  )
}

# Run fit + predict + stats for all species, given a response column
run_all_species <- function(temporal_data, response_col, species_levels) {
  models <- map(species_levels, function(sp) {
    df <- temporal_data %>%
      filter(!is.na(year), year < 2024, Species == sp) %>%
      mutate(fDataset = as.factor(Dataset))
    fit_nb_gam(df, response_col)
  }) %>% set_names(species_levels)
  
  preds <- map_dfr(species_levels, function(sp) {
    df <- temporal_data %>%
      filter(!is.na(year), year < 2024, Species == sp) %>%
      mutate(fDataset = as.factor(Dataset))
    predict_nb_gam(models[[sp]], df, sp)
  }) %>%
    mutate(est     = exp(fit),
           low     = exp(fit - 1.96 * se.fit),
           high    = exp(fit + 1.96 * se.fit),
           Dataset = as.character(fDataset))
  
  gam_stats <- map_dfr(species_levels, ~ get_gam_stats(models[[.x]], .x))
  
  preds <- preds %>%
    left_join(gam_stats, by = "Species") %>%
    mutate(label = paste0("edf=", edf, "\np=", pval))
  
  list(models = models, preds = preds, gam_stats = gam_stats)
}

# -----------------------------------------------------------------------------
# Reusable phylopic placement
# -----------------------------------------------------------------------------
# Fixed placement fractions — consistent across ALL plots regardless of data.
# Tie x to a hard year (x_hard) when every panel shares an x-axis; set x_hard
# to NULL to fall back to x_frac for plots where the x range varies per panel.
PHYLOPIC_CFG <- list(
  x_frac  = 0.20,   # fraction across x-axis range (used when x_hard is NULL)
  y_frac  = 0.9,   # fraction up y-axis range
  h_frac  = 0.14,   # silhouette height as fraction of y range
  x_hard  = 1835    # hard x value; overrides x_frac when non-NULL
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
# `df` is the plotted data; xvar/yvar are the aesthetics driving the panels.
build_phylopic_positions <- function(df, xvar, yvar, lookup = phylopic_lookup,
                                     cfg = PHYLOPIC_CFG,
                                     mirrored = mirrored_species,
                                     species_levels = SPECIES_LEVELS) {
  ranges <- df %>%
    group_by(Species) %>%
    summarise(
      xmin = min(.data[[xvar]], na.rm = TRUE),
      xmax = max(.data[[xvar]], na.rm = TRUE),
      ymin = min(.data[[yvar]], na.rm = TRUE),
      ymax = max(.data[[yvar]], na.rm = TRUE),
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
      Species = factor(row$Species, levels = levels(positions$Species)),
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
# 1. Prep: year-spine and plot counts
# =============================================================================
plot_counts <- temporal_data %>%
  filter(!is.na(year)) %>%
  count(Species, Dataset, year, name = "n") %>%
  mutate(Species = factor(Species, levels = SPECIES_LEVELS))

# First observed year per species x dataset
firstYear <- temporal_data %>%
  filter(!is.na(year)) %>%
  group_by(Species, Dataset) %>%
  summarise(first_year = min(year, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Dataset, values_from = first_year)

# Year spine: one row per species x dataset x year (fills gaps)
syears <- map_dfr(SPECIES_LEVELS, function(sp) {
  fy <- firstYear %>% filter(Species == sp)
  bind_rows(
    tibble(Species = sp, Dataset = "VertNet",
           first.obs = seq(as.integer(fy$VertNet), 2024)),
    tibble(Species = sp, Dataset = "iNaturalist",
           first.obs = seq(as.integer(fy$iNaturalist), 2024))
  )
})

# =============================================================================
# 2. Make record accumulation figure [Figure 1]
# =============================================================================

# Cumulative records per Species x Dataset across years
accum <- temporal_data %>%
  filter(!is.na(year)) %>%
  group_by(Species, Dataset, year) %>%
  summarise(n = sum(RecordCount, na.rm = TRUE), .groups = "drop") %>%
  arrange(Species, Dataset, year) %>%
  group_by(Species, Dataset) %>%
  mutate(cumulative = cumsum(n)) %>%
  ungroup() %>%
  mutate(Species = factor(Species, levels = SPECIES_LEVELS))

# Plot base 
accum_plot <- ggplot(accum, aes(year, cumulative, color = Dataset)) +
  geom_line(linewidth = 1) +
  facet_wrap(~Species, scales = "free_y", nrow = 2) +
  scale_color_manual(values = COLOURS) +
  scale_y_continuous(breaks = scales::breaks_extended(n = 4),
                     expand = expansion(mult = c(0.02, 0.12))) +
  labs(x = "Year", y = "Cumulative records") +
  paper_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "none",
        aspect.ratio = 1,
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 18),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_text(size = 18))

# Add critters using the reusable placement workflow
phylopic_auto <- build_phylopic_positions(accum, xvar = "year", yvar = "cumulative")
accum_plot    <- add_phylopic_auto(accum_plot, phylopic_auto)

print(accum_plot)
ggsave(file.path(FIG_DIR, "fig1_record_accumulation.pdf"), accum_plot, width = 12, height = 6.5, units = "in", dpi = 300)
ggsave(file.path(FIG_DIR, "fig1_record_accumulation.png"), accum_plot, width = 12, height = 6.5, units = "in", dpi = 300)