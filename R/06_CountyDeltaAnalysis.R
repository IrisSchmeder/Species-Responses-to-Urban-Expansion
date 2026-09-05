# =============================================================================
# 06_CountyDeltaAnalysis.R
# County-level dataset delta analysis: paired t-tests, z-score maps,
# regional summaries, and dataset-shift GLMM.
#
# DEPENDENCIES: run 01_DataPrep.R first
#
# OUTPUT: paired_ttest_boxplot.pdf / .png AND county_log_diff_map.png
#         [Figure 5] AND [Figure 4]
# =============================================================================

library(tidyverse)    
library(here)         
library(conflicted)   
library(sf)            
library(lme4)         
library(emmeans)       
library(viridis)      
library(patchwork)     
library(scales)       
library(rphylopic)     

conflicted::conflicts_prefer(dplyr::select)
conflicted::conflicts_prefer(dplyr::filter)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
DATA_DIR <- here::here("data")
PREP_DIR <- here::here("data", "prepared")
FIG_DIR  <- here::here("figures")

SPECIES               <- readRDS(file.path(PREP_DIR, "SPECIES.rds"))
SPECIES_LEVELS        <- readRDS(file.path(PREP_DIR, "SPECIES_LEVELS.rds"))
texas_counties        <- readRDS(file.path(PREP_DIR, "texas_counties.rds"))
urbanization_raw      <- readRDS(file.path(PREP_DIR, "urbanization_raw.rds"))
county_delta_combined <- readRDS(file.path(PREP_DIR, "county_delta_combined.rds"))
county_data_filtered  <- readRDS(file.path(PREP_DIR, "county_data_filtered.rds"))

species_list <- unlist(SPECIES, use.names = FALSE)

species_display <- tibble::tribble(
  ~full,             ~label,          ~phylopic_name,
  "A. punctatus",    "A. punctatus",  "Rhinella marina",
  "A. woodhousii",   "A. woodhousii", "Bufo bufo",
  "I. nebulifer",    "I. nebulifer",  "Duttaphrynus brevirostris",
  "Co. texanus",     "Co. texanus",   "Uta stansburiana stejnegeri",
  "Cr. collaris",    "Cr. collaris",  "Crotaphytus collaris",
  "P. cornutum",     "P. cornutum",   "Phrynosoma solare",
  "Te. ornata",      "Te. ornata",    "Terrapene carolina",
  "Tr. scripta",     "Tr. scripta",   "Trachemys scripta elegans"
)

SPECIES_LEVELS <- species_display$full                                   # display order
SPECIES_LABELS <- setNames(species_display$label, species_display$full)  # axis text

# Reversed factor: use when Species maps to the y-axis or the plot is flipped,
# so the canonical SPECIES_LEVELS order reads top-to-bottom.
SPECIES_LEVELS_REV <- rev(SPECIES_LEVELS)
species_rev <- function(x) factor(as.character(x), levels = SPECIES_LEVELS_REV)

# -----------------------------------------------------------------------------
# Shared theme (matches 05_RatioDataset.R)
# -----------------------------------------------------------------------------
paper_theme <- theme_bw(base_size = 24) +
  theme(
    axis.title       = element_text(size = 24),
    axis.text        = element_text(size = 20, color = "black"),
    strip.text       = element_text(size = 22, face = "italic"),
    legend.text      = element_text(size = 22),
    legend.title     = element_text(size = 23),
    plot.title       = element_text(size = 26, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.spacing    = unit(1, "lines")
  )

# -----------------------------------------------------------------------------
# Shared phylopic config (matches the other scripts)
# -----------------------------------------------------------------------------
phylopic_lookup <- tibble(
  Species       = SPECIES_LEVELS,
  phylopic_name = c("Bufo bufo", "Bufo bufo", "Bufo bufo",
                    "Leiocephalus carinatus", "Leiocephalus carinatus", "Leiocephalus carinatus", 
                    "Rhinoclemmys punctularia", "Rhinoclemmys punctularia")
)

mirrored_species <- c("A. punctatus", "A. woodhousii", "I. nebulifer",
                      "Co. texanus", "Cr. collaris", "P. cornutum",
                      "Te. ornata", "Tr. scripta")


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
# run a paired t-test on a species subset
run_paired_ttest <- function(df) {
  tt <- t.test(df$inat_count, df$vertnet_count, paired = TRUE)
  tibble(n_counties   = nrow(df),
         mean_inat    = mean(df$inat_count),
         mean_vertnet = mean(df$vertnet_count),
         mean_diff    = tt$estimate,
         ci_low       = tt$conf.int[1],
         ci_high      = tt$conf.int[2],
         t_statistic  = tt$statistic,
         p_value      = tt$p.value)
}

# =============================================================================
# 1. Time-filtered paired t-test (VertNet pre-2008, iNat 2008+)
# =============================================================================
paired_ttest_filtered <- county_data_filtered %>%
  group_by(Species) %>%
  group_modify(~ run_paired_ttest(.x)) %>%
  ungroup() %>%
  mutate(
    p_adj          = p.adjust(p_value, method = "BH"),
    interpretation = case_when(
      p_adj < 0.05 & mean_diff > 0 ~ "Significantly more iNaturalist",
      p_adj < 0.05 & mean_diff < 0 ~ "Significantly more VertNet",
      TRUE ~ "No significant difference"
    )
  )

print(paired_ttest_filtered)

# =============================================================================
# 2. Boxplot of paired t-test [Figure 5]
# =============================================================================

# -----------------------------------------------------------------------------
# Boxplot: county-level record-count difference per species, colored by
# significance (BH-adjusted) instead of an asterisk. Individual counties are
# shown behind the box as jittered, semi-transparent points, so counties that
# land close together read as a darker patch rather than overlapping solid
# dots. outlier.shape = NA on the boxplot avoids drawing the same points
# twice (once as jitter, once as boxplot outliers).
# -----------------------------------------------------------------------------

paired_ttest_plot_df <- county_data_filtered %>%
  mutate(diff = inat_count - vertnet_count) %>%
  left_join(paired_ttest_filtered %>% select(Species, p_value), by = "Species") %>%
  mutate(sig = p_value < 0.05, Species_y = species_rev(Species),
         # Darker/more saturated variants of the fill colors, used only for
         # the jitter dots and the box outline/whiskers/median - so the
         # outline stays visible against the (lighter) fill instead of
         # blending into it.
         dot_color    = ifelse(sig, "#CC79A7", "gray55"),
         border_color = ifelse(sig, "#994C71", "gray35"))

# Placement table for the discrete species y-axis: one silhouette per species
# at a fixed data-x just inside the left plot limit (set via scale_x_continuous
# below) so it isn't dropped by oob = squish. This is the boxplot analogue of
# add_phylopic_auto - same shared lookup and mirrored set, but the x/y/height
# come from the discrete axis rather than a per-species data range.
img_x <- -140

phylopic_box <- paired_ttest_filtered %>%
  distinct(Species) %>%
  left_join(phylopic_lookup, by = "Species") %>%
  mutate(
    Species_y  = species_rev(Species),
    x          = img_x,
    horizontal = Species %in% mirrored_species,
    height     = 0.40
  )

# add one geom_phylopic layer per species (mirrors add_phylopic_auto's loop,
# but keyed to the discrete y-axis)
add_phylopic_discrete <- function(p, positions) {
  for (i in seq_len(nrow(positions))) {
    row <- positions[i, ]
    ann <- data.frame(x = row$x, Species_y = row$Species_y,
                      name = row$phylopic_name)
    p <- p + geom_phylopic(
      data = ann, aes(x = x, y = Species_y, name = name),
      inherit.aes = FALSE, color = "gray20",
      height = row$height, horizontal = row$horizontal
    )
  }
  p
}

pval_labels <- paired_ttest_plot_df %>%
  distinct(Species_y, sig, p_value)
print(pval_labels)

paired_ttest_boxplot_fig <- ggplot(paired_ttest_plot_df, aes(x = diff, y = Species_y)) +
  geom_vline(xintercept = 0, linewidth = 0.5, color = "gray30", linetype = "dashed") +
  geom_jitter(aes(color = dot_color), height = 0.28, alpha = 0.12, size = 1.4, show.legend = FALSE) +
  geom_boxplot(aes(fill = sig, color = border_color), outlier.shape = NA, width = 0.6, linewidth = 0.4) +
  scale_fill_manual(values = c("FALSE" = "gray55", "TRUE" = "#CC79A7"),
                    labels = c("n.s.", "BH p < 0.05"), name = NULL) +
  scale_color_identity() +
  scale_x_continuous(limits = c(-150, 310), oob = squish,
                     breaks = c(-100, -50, 0, 50, 100, 200, 300),
                     labels = c("-100", "-50", "0", "50", "100", "200", "300")) +
  scale_y_discrete(labels = SPECIES_LABELS) +
  geom_text(
    data = pval_labels,
    aes(x = 295, y = Species_y,
        label = ifelse(p_value < 0.001,
                       "italic(p) < 0.001",
                       sprintf("italic(p)==\"%.3f\"", p_value))),
    parse = TRUE, hjust = 1, size = 6, color = "gray10",
    inherit.aes = FALSE
  ) +
  paper_theme +
  theme(axis.text.y        = element_text(face = "italic", size = 20, color = "gray10"),
        axis.text.x        = element_text(size = 18, color = "gray10"),
        axis.title.x       = element_text(size = 20, color = "black"),
        panel.grid.major.y = element_blank(),
        legend.position    = "none", #removed "top"
        legend.text        = element_text(size = 16, color = "gray10"),
        plot.margin        = margin(10, 10, 10, 40),
        aspect.ratio       = 0.6) +
  labs(x = "County-level record count difference \n(iNaturalist 2008-2024 minus VertNet pre-2008)",
       y = "")

# silhouettes added last so they sit on top of the boxes
paired_ttest_boxplot_fig <- add_phylopic_discrete(paired_ttest_boxplot_fig, phylopic_box)

print(paired_ttest_boxplot_fig)
ggsave(file.path(FIG_DIR, "fig5_paired_ttest_boxplot.pdf"), paired_ttest_boxplot_fig, width = 12, height = 8, units = "in", dpi = 300)
ggsave(file.path(FIG_DIR, "fig5_paired_ttest_boxplot.png"), paired_ttest_boxplot_fig, width = 12, height = 8, units = "in", dpi = 300)

# =============================================================================
# 3. County-level log-difference map [Figure 4]
# =============================================================================
county_prop <- urbanization_raw %>%
  filter(!is.na(county)) %>%
  group_by(Species, county, Dataset) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = Dataset, values_from = n, values_fill = 0L) %>%
  mutate(total     = VertNet + iNaturalist,
         prop_inat = iNaturalist / total,
         diff      = iNaturalist - VertNet,
         log_diff  = sign(diff) * log1p(abs(diff))) %>%
  filter(total > 0)

county_map_prop <- texas_counties %>% left_join(county_prop, by = "county")
log_cap         <- quantile(abs(county_prop$log_diff), 0.95, na.rm = TRUE)
map_plots       <- list()

# Phylopic placement for the maps: same top-right-corner idea and the same
# mirrored logic as the boxplot above, just scaled to the map's own
# coordinate range (a Texas county map's bounding box) instead of a data
# range. Silhouette names come from the shared phylopic_lookup so the maps
# use the SAME silhouette set as Figure 5. All icons use one shared height
# (icon_height_base).
map_bbox <- sf::st_bbox(texas_counties)
img_x <- map_bbox["xmax"] - 0.18 * (map_bbox["xmax"] - map_bbox["xmin"])
img_y <- map_bbox["ymax"] - 0.04 * (map_bbox["ymax"] - map_bbox["ymin"])
icon_height_base <- 0.18 * (map_bbox["ymax"] - map_bbox["ymin"])
phylopic_attributions <- list()

for (sp in species_list) {
  sp_data    <- county_map_prop %>% filter(Species == sp)
  n_counties <- sum(!is.na(sp_data$log_diff))
  if (n_counties == 0) { cat("Skipping", sp, "- no data\n"); next }
  
  sp_phylopic_name <- phylopic_lookup$phylopic_name[phylopic_lookup$Species == sp]
  sp_uuid       <- get_uuid(name = sp_phylopic_name)
  sp_attribution <- get_attribution(sp_uuid)
  phylopic_attributions[[sp]] <- sp_attribution
  cat("Attribution for", sp, "(", sp_phylopic_name, "):\n")
  print(sp_attribution)
  sp_mirrored       <- sp %in% mirrored_species
  sp_height         <- icon_height_base
  
  map_plots[[sp]] <- ggplot() +
    geom_sf(data = texas_counties,
            fill = "gray30", color = "white", linewidth = 0.25) +
    geom_sf(data = sp_data %>% filter(!is.na(log_diff)),
            aes(fill = log_diff), color = "white", linewidth = 0.25) +
    geom_phylopic(
      data = data.frame(x = img_x, y = img_y, name = sp_phylopic_name),
      aes(x = x, y = y, name = name),
      inherit.aes = FALSE, color = "gray20",
      height = sp_height, horizontal = sp_mirrored
    ) +
    scale_fill_gradient2(
      low = "#009E73", mid = "#F7F7F7", high = "#D55E00",
      midpoint = 0, limits = c(-log_cap, log_cap), oob = squish,
      na.value = "gray30",
      name   = "Log Difference (iNaturalist − VertNet)",
      breaks = c(-log_cap, 0, log_cap),
      labels = c(
        sprintf("%.1f \nMore Historical Records\n(VertNet)", -log_cap),
        "0 \nEqual",
        sprintf("%.1f \nMore Contemporary Records\n(iNaturalist)", log_cap)
      ),
      guide = guide_colorbar(
        direction = "horizontal",
        title.position = "top",
        title.hjust = 0.5,
        label.position = "bottom",
        barwidth = unit(28, "lines"),
        barheight = unit(0.8, "lines")
      )
    ) +
    paper_theme +
    theme(axis.text  = element_blank(), axis.ticks  = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          # legend.position = "none",
          legend.position = "bottom",
          legend.title = element_text(size = 16, hjust = 0.5),
          legend.text  = element_text(size = 12),
          plot.title   = element_text(face = "italic", hjust = 0.5, size = 18),
          plot.caption = ggtext::element_markdown(hjust = 0.5, size = 14, color = "gray10")
    ) +
    labs(title = sp, x = NULL, y = NULL,
         caption = paste0("*n* = ", n_counties))
}

credit_lines <- map_chr(names(phylopic_attributions), function(sp) {
  img_info     <- phylopic_attributions[[sp]]$images[[1]]
  sp_name      <- phylopic_lookup$phylopic_name[phylopic_lookup$Species == sp]
  license_note <- if (img_info$license_abbr != "CC0 1.0") {
    paste0(" (", img_info$license_abbr, ")")
  } else {
    ""
  }
  paste0(sp_name, " by ", img_info$attribution, license_note)
})

phylopic_credit <- paste0(
  "Silhouettes: ", paste(unique(credit_lines), collapse = ", "),
  ". All via phylopic.org."
)

cat(phylopic_credit, "\n")
writeLines(phylopic_credit, "phylopic_credit_line.txt")

combined_map <- wrap_plots(map_plots, ncol = 4) /
  guide_area() +
  plot_layout(guides = "collect", heights = c(1, 0.15)) &
  theme(legend.position = "bottom")

combined_map <- combined_map +
  plot_annotation(
    theme = theme(plot.margin = margin(2, 2, 2, 2))
  )

print(combined_map)

ggsave(file.path(FIG_DIR, "fig4_county_log_diff_map.png"), combined_map, width = 12, height = 9, units = "in", dpi = 600)
ggsave(file.path(FIG_DIR, "fig4_county_log_diff_map.pdf"), combined_map, width = 12, height = 9, units = "in", dpi = 600)
