# =============================================================================
# 04_TemporalBiasGAMs.R
# Negative-binomial GAMs for record counts and unique observers over time.
# Produces a two-panel combined figure: (a) records, (b) observers.
#
# DEPENDENCIES: run 01_DataPrep.R first
#
# OUTPUT: temporal_records_observers.pdf / .png
#         [Figure 2]
# =============================================================================

library(tidyverse)
library(here)         
library(mgcv)
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

temporal_data <- readRDS(file.path(PREP_DIR, "temporal_data.rds"))
SPECIES_LEVELS <- readRDS(file.path(PREP_DIR, "SPECIES_LEVELS.rds"))

COLOURS <- c("iNaturalist" = "#D55E00", "VertNet" = "#009E73")

ITALIC_LEVELS <- paste0("italic('", SPECIES_LEVELS, "')")
italic_species <- function(x) {
  factor(paste0("italic('", as.character(x), "')"), levels = ITALIC_LEVELS)
}

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

# -----------------------------------------------------------------------------
# Reusable phylopic placement
# -----------------------------------------------------------------------------
PHYLOPIC_CFG <- list(
  x_frac  = 0.20,   # fraction across x-axis range (used when x_hard is NULL)
  y_frac  = 0.68,   # fraction up y-axis range (May need to change depending on edf sizing)
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

# Adds one geom_phylopic() layer per row. Looping (rather than one layer with
# aes()-mapped columns) is necessary because `horizontal` only behaves
# reliably as a fixed per-layer argument, not as an aes()-mapped column.
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

# Range source for silhouette placement: union of plotted raw points and the
# GAM ribbon, keyed by the plain (un-parsed) Species so build_phylopic_*
# can join on phylopic_lookup. yvar becomes a common name ("y_range").
placement_data <- function(points_df, preds_df, response_col) {
  pts <- points_df %>%
    filter(!is.na(year)) %>%
    transmute(Species = as.character(Species), year, y_range = .data[[response_col]])
  rib <- preds_df %>%
    transmute(
      Species = gsub("italic\\('|'\\)", "", as.character(Species)),
      year, y_range = high
    )
  bind_rows(pts, rib) %>% filter(!is.na(y_range))
}

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

fit_nb_gam <- function(data, response_col) {
  f <- as.formula(paste0(response_col, " ~ s(year, by = fDataset, bs = 'cr')"))
  gam(f, family = nb(), data = data)
}

predict_nb_gam <- function(model, data, species_name) {
  newdat <- with(data, expand.grid(
    year     = seq(min(year), max(year), by = 1),
    fDataset = unique(fDataset)
  ))
  pred           <- predict(model, newdata = newdat, type = "link", se.fit = TRUE)
  newdat$fit     <- pred$fit
  newdat$se.fit  <- pred$se.fit
  newdat$Species <- species_name
  
  min_inat <- data %>% filter(fDataset == "iNaturalist") %>%
    summarise(min(year)) %>% pull()
  newdat %>% filter(!(fDataset == "iNaturalist" & year < min_inat))
}

get_gam_stats <- function(model, species_name) {
  s <- summary(model)$s.table
  dataset_name <- gsub(".*fDataset", "", rownames(s))
  pv <- s[, "p-value"]
  data.frame(
    Species = species_name,
    Dataset = dataset_name,
    edf     = round(s[, "edf"], 2),
    pval    = pv
  )
}

# Fit GAMs for all species, return predictions and stats
run_all_species <- function(temporal_data, response_col) {
  models <- map(SPECIES_LEVELS, function(sp) {
    d <- temporal_data %>%
      filter(!is.na(year), year < 2024, Species == sp) %>%
      mutate(fDataset = as.factor(Dataset))
    fit_nb_gam(d, response_col)
  }) %>% set_names(SPECIES_LEVELS)
  
  preds <- map_dfr(SPECIES_LEVELS, function(sp) {
    d <- temporal_data %>%
      filter(!is.na(year), year < 2024, Species == sp) %>%
      mutate(fDataset = as.factor(Dataset))
    predict_nb_gam(models[[sp]], d, sp)
  }) %>%
    mutate(
      est     = exp(fit),
      low     = exp(fit - 1.96 * se.fit),
      high    = exp(fit + 1.96 * se.fit),
      Dataset = as.character(fDataset)
    )
  
  stats <- map_dfr(SPECIES_LEVELS, ~ get_gam_stats(models[[.x]], .x))
  
  preds <- preds %>%
    left_join(stats, by = c("Species", "Dataset")) %>%
    mutate(
      label   = paste0("edf=", round(edf, 2)),
      Species = italic_species(Species)
    )
  
  list(models = models, preds = preds, stats = stats)
}

# =============================================================================
# 1. Record count GAMs
# =============================================================================

rec  <- run_all_species(temporal_data, "RecordCount")

effortTime.Records.plot <- temporal_data %>%
  filter(!is.na(year)) %>%
  mutate(Species = italic_species(Species)) %>%
  ggplot(aes(group = Dataset, y = RecordCount, x = year)) +
  geom_point(aes(color = Dataset), size = 1, alpha = 0.4) +
  paper_theme +
  geom_ribbon(data = rec$preds,
              aes(x = year, ymin = low, ymax = high, fill = Dataset),
              inherit.aes = FALSE, linetype = 0, alpha = 0.75) +
  geom_line(data = rec$preds,
            aes(x = year, y = est, group = Dataset, color = Dataset),
            inherit.aes = FALSE) +
  geom_text(
    data = rec$preds %>% group_by(Species, Dataset) %>% slice(1) %>%
      mutate(label_vjust = ifelse(Dataset == "iNaturalist", 1.2, 2.6)),
    aes(x = min(year), y = Inf, label = label, vjust = label_vjust, color = Dataset),
    inherit.aes = FALSE, hjust = 0, size = 6, show.legend = FALSE
  ) +
  scale_color_manual(values = COLOURS) +
  scale_fill_manual(values  = COLOURS) +
  facet_wrap(~ Species, scales = "free_y", nrow = 2, labeller = label_parsed) +
  ylab("Number of records") +
  xlab(NULL) +
  scale_y_continuous(labels = function(x) formatC(x, width = 5, format = "d", flag = " ")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 12),
        legend.position = "none",
        aspect.ratio = 1,
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 18),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_text(size = 18))

phylopic_positions_records <- build_phylopic_positions(
  placement_data(temporal_data, rec$preds, "RecordCount"),
  xvar = "year", yvar = "y_range"
)
effortTime.Records.plot <- add_phylopic_auto(
  effortTime.Records.plot, phylopic_positions_records, parse_species = TRUE
)

print(effortTime.Records.plot)

# =============================================================================
# 2. Unique observer/collector GAMs
# =============================================================================

obs <- run_all_species(temporal_data, "UniqueCollectors")

obs_points <- temporal_data %>%
  filter(!is.na(year)) %>%
  select(Species, Dataset, year, UniqueCollectors) %>%
  distinct()

ObsTime.plot <- obs_points %>%
  mutate(Species = italic_species(Species)) %>%
  ggplot(aes(group = Dataset, y = UniqueCollectors, x = year)) +
  geom_point(aes(color = Dataset), size = 1, alpha = 0.4) +
  paper_theme +
  geom_ribbon(data = obs$preds,
              aes(x = year, ymin = low, ymax = high, fill = Dataset),
              inherit.aes = FALSE, linetype = 0, alpha = 0.75) +
  geom_line(data = obs$preds,
            aes(x = year, y = est, group = Dataset, color = Dataset),
            inherit.aes = FALSE) +
  geom_text(
    data = obs$preds %>% group_by(Species, Dataset) %>% slice(1) %>%
      mutate(label_vjust = ifelse(Dataset == "iNaturalist", 1.2, 2.6)),
    aes(x = min(year), y = Inf, label = label, vjust = label_vjust, color = Dataset),
    inherit.aes = FALSE, hjust = 0, size = 6, show.legend = FALSE
  ) +
  scale_color_manual(values = COLOURS) +
  scale_fill_manual(values  = COLOURS) +
  facet_wrap(~ Species, scales = "free_y", nrow = 2, labeller = label_parsed) +
  xlab("Year") +
  ylab("Number of observers") +
  scale_y_continuous(labels = function(x) formatC(x, width = 5, format = "d", flag = " ")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "none",
        aspect.ratio = 1,
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 18),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_text(size = 18))


phylopic_positions_observers <- build_phylopic_positions(
  placement_data(obs_points, obs$preds, "UniqueCollectors"),
  xvar = "year", yvar = "y_range"
)
ObsTime.plot <- add_phylopic_auto(
  ObsTime.plot, phylopic_positions_observers, parse_species = TRUE
)

print(ObsTime.plot)

# =============================================================================
# 3. Combined two-panel figure [Figure 2]
# =============================================================================

combo <- cowplot::plot_grid(
  ggplotGrob(effortTime.Records.plot),
  ggplotGrob(ObsTime.plot),
  ncol           = 1,
  labels         = c("(a)", "(b)"),
  label_fontface = "plain",
  label_size     = 24,
  hjust          = 0,
  label_x        = 0.01,
  align          = "h"
)

combo <- cowplot::plot_grid(
  effortTime.Records.plot,
  ObsTime.plot,
  ncol           = 1,
  labels         = c("(a)", "(b)"),
  label_fontface = "plain",
  label_size     = 24,
  hjust          = 0,
  label_x        = 0.01,
  align          = "v",
  axis           = "lr"
)

print(combo)

ggsave(file.path(FIG_DIR, "fig2_temporal_records_observers.pdf"), combo, width = 12, height = 13, units = "in", dpi = 300)
ggsave(file.path(FIG_DIR, "fig2_temporal_records_observers.png"), combo, width = 12, height = 13, units = "in", dpi = 300)
