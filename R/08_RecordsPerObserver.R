# =============================================================================
# 08_RecordsPerObserver.R
# Records per observer over time, by species and dataset.
#
# DEPENDENCIES: run 01_DataPrep.R first
#
# OUTPUT: records_per_observer.pdf / .png
#         [SI Figure 1]
# =============================================================================

library(tidyverse)     
library(here)         
library(conflicted)   
library(broom)         
library(mgcv)          
library(rphylopic)   

conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(dplyr::select)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
DATA_DIR <- here::here("data")
PREP_DIR <- here::here("data", "prepared")
FIG_DIR  <- here::here("figures")

SPECIES         <- readRDS(file.path(PREP_DIR, "SPECIES.rds"))
SPECIES_LEVELS  <- readRDS(file.path(PREP_DIR, "SPECIES_LEVELS.rds"))
temporal_data   <- readRDS(file.path(PREP_DIR, "temporal_data.rds"))
vertnet_spatial <- readRDS(file.path(PREP_DIR, "vertnet_spatial.rds"))
inat_spatial    <- readRDS(file.path(PREP_DIR, "inat_spatial.rds"))

N_SPECIES <- length(SPECIES_LEVELS)

# Plot colours (matches 03)
COLOURS <- c("iNaturalist" = "#D55E00", "VertNet" = "#009E73")

# -----------------------------------------------------------------------------
# Paper theme (matches 03_RecordAccumulation.R)
# -----------------------------------------------------------------------------
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
    panel.spacing    = unit(1, "lines")
  )

# Italic facet labels in canonical order 
ITALIC_LEVELS <- paste0("italic('", SPECIES_LEVELS, "')")
italic_species <- function(x) {
  factor(paste0("italic('", as.character(x), "')"), levels = ITALIC_LEVELS)
}

# =============================================================================
# 1. Calculate records per observer
# =============================================================================
cat("--- Computing records per observer ---\n")

temporal_data <- temporal_data %>%
  mutate(
    RecordsPerObserver = RecordCount / UniqueCollectors,
    RecordsPerObserver = ifelse(is.infinite(RecordsPerObserver) | is.nan(RecordsPerObserver),
                                NA, RecordsPerObserver)
  )

# =============================================================================
# 2. Reusable phylopic placement
# =============================================================================
PHYLOPIC_CFG <- list(
  x_frac  = 0.20,   # fraction across x-axis range (used when x_hard is NULL)
  y_frac  = 0.9,    # fraction up y-axis range
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

# Build placement table 
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

# `parse_species` converts the facet key to the italic plotmath factor so the
# layers land in the right panels when facets use label_parsed.
add_phylopic_auto <- function(p, positions, parse_species = FALSE) {
  for (i in seq_len(nrow(positions))) {
    row <- positions[i, ]
    sp_key <- if (parse_species) italic_species(as.character(row$Species))
    else factor(row$Species, levels = levels(positions$Species))
    ann <- data.frame(
      Species = sp_key, x = row$x, y = row$y, name = row$phylopic_name
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
# 3. Records-per-observer over time by species [SI Figure 1]
# =============================================================================
plot_data <- temporal_data %>%
  filter(!is.na(RecordsPerObserver), is.finite(RecordsPerObserver)) %>%
  mutate(Species = factor(Species, levels = SPECIES_LEVELS))

# -----------------------------------------------------------------------------
# Significance asterisks + GAM stats table
# -----------------------------------------------------------------------------
# Wrapped in tryCatch since a handful of species x dataset combos may have
# too few distinct years for a k = 5 smooth to fit.
get_smooth_stats <- function(data) {
  m <- tryCatch(
    gam(RecordsPerObserver ~ s(year, k = 5), data = data),
    error = function(e) NULL
  )
  if (is.null(m)) return(tibble(edf = NA_real_, pval = NA_real_))
  s <- summary(m)$s.table
  tibble(edf = round(s[1, "edf"], 2), pval = s[1, "p-value"])
}

gam_stats <- plot_data %>%
  group_by(Species, Dataset) %>%
  group_modify(~ get_smooth_stats(.x)) %>%
  ungroup()

cat("\n=== GAM smooth stats: RecordsPerObserver ~ s(year, k = 5), by Species x Dataset ===\n")
gam_stats %>%
  mutate(pval = signif(pval, 3)) %>%
  arrange(Species, Dataset) %>%
  print(n = Inf)

# Significance asterisks, keyed to the italic (parsed) facet so they land in
# the correct panels. Stacked so the two datasets' asterisks don't overlap.
sig_labels <- gam_stats %>%
  mutate(
    sig_label = case_when(
      is.na(pval)  ~ "",
      pval < 0.001 ~ "***",
      pval < 0.01  ~ "**",
      pval < 0.05  ~ "*",
      TRUE         ~ ""
    ),
    label_vjust = ifelse(Dataset == "iNaturalist", 1.2, 2.6),
    Species     = italic_species(Species)
  )

# Plot 
sifig1 <- plot_data %>%
  mutate(Species = italic_species(Species)) %>%
  ggplot(aes(x = year, y = RecordsPerObserver, color = Dataset, fill = Dataset)) +
  geom_point(alpha = 0.45, size = 1) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5), linewidth = 1) +
  facet_wrap(~ Species, scales = "free_y", nrow = 2,
             drop = FALSE, labeller = label_parsed) +
  scale_color_manual(values = COLOURS) +
  scale_fill_manual(values  = COLOURS) +
  labs(x = "Year", y = "Records per observer / collector", color = NULL, fill = NULL) +
  paper_theme +
  theme(
    legend.position  = "none",
    aspect.ratio     = 1,
    strip.background = element_blank(),
    strip.text       = element_text(face = "italic"),
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18))

# Silhouettes via the shared fraction-based workflow (parse_species = TRUE so
# the layers land in the italic label_parsed panels).
phylopic_positions <- build_phylopic_positions(
  plot_data, xvar = "year", yvar = "RecordsPerObserver"
)
sifig1 <- add_phylopic_auto(sifig1, phylopic_positions, parse_species = TRUE)

# Significance asterisks, top-right corner (opposite the top-left silhouettes)
sifig1 <- sifig1 +
  geom_text(
    data = sig_labels,
    aes(x = Inf, y = Inf, label = sig_label, vjust = label_vjust, color = Dataset),
    inherit.aes = FALSE, hjust = 1.1, size = 6, fontface = "bold", show.legend = FALSE
  )

print(sifig1)

ggsave(file.path(FIG_DIR, "sifig1_records_per_observer.pdf"), sifig1,
       width = 12, height = 6.5, units = "in", dpi = 300)
ggsave(file.path(FIG_DIR, "sifig1_records_per_observer.png"), sifig1,
       width = 12, height = 6.5, units = "in", dpi = 300)
