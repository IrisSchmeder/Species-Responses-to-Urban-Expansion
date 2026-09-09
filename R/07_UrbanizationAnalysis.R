# =============================================================================
# 07_UrbanizationAnalysis.R
# Urbanization trends, impervious-surface models
#
# DEPENDENCIES: run 01_DataPrep.R first
#
# OUTPUT: imperv_coef_weighted.pdf / .png
#         [Figure 6]
# =============================================================================

library(tidyverse)     
library(here)          
library(conflicted)    
library(sf)           
library(lme4)         
library(mgcv)        
library(emmeans)      
library(scales)       
library(rphylopic)    

conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(dplyr::select)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
DATA_DIR <- here::here("data")
PREP_DIR <- here::here("data", "prepared")
FIG_DIR  <- here::here("figures")

SPECIES            <- readRDS(file.path(PREP_DIR, "SPECIES.rds"))
SPECIES_LEVELS     <- readRDS(file.path(PREP_DIR, "SPECIES_LEVELS.rds"))
texas_counties     <- readRDS(file.path(PREP_DIR, "texas_counties.rds"))
county_imperv      <- readRDS(file.path(PREP_DIR, "county_imperv.rds"))
urbanization_raw   <- readRDS(file.path(PREP_DIR, "urbanization_raw.rds"))
county_counts      <- readRDS(file.path(PREP_DIR, "county_counts.rds"))
master             <- readRDS(file.path(PREP_DIR, "master.rds"))
inat_spatial       <- readRDS(file.path(PREP_DIR, "inat_spatial.rds"))
imperv_scaling     <- readRDS(file.path(PREP_DIR, "imperv_scaling.rds"))

species_list  <- unlist(SPECIES, use.names = FALSE)
imperv_center <- imperv_scaling$center
imperv_scale  <- imperv_scaling$scale

# Reversed factor: use when Species maps to the y-axis or the plot is flipped,
# so the canonical SPECIES_LEVELS order reads top-to-bottom.
SPECIES_LEVELS_REV <- rev(SPECIES_LEVELS)
species_rev <- function(x) factor(as.character(x), levels = SPECIES_LEVELS_REV)

# =============================================================================
# 1. Custom species display order + abbreviated axis labels
# =============================================================================
species_display <- tibble::tribble(
  ~full,             ~label,          ~phylopic_name,
  "A. punctatus",    "A. punctatus",  "Bufo bufo",
  "A. woodhousii",   "A. woodhousii", "Bufo bufo",
  "I. nebulifer",    "I. nebulifer",  "Bufo bufo",
  "Co. texanus",     "Co. texanus",   "Leiocephalus carinatus",
  "Cr. collaris",    "Cr. collaris",  "Leiocephalus carinatus",
  "P. cornutum",     "P. cornutum",   "Leiocephalus carinatus",
  "Te. ornata",      "Te. ornata",    "Rhinoclemmys punctularia",
  "Tr. scripta",     "Tr. scripta",   "Rhinoclemmys punctularia"
)

SPECIES_LEVELS <- species_display$full                                   # display order
SPECIES_LABELS <- setNames(species_display$label, species_display$full)  # axis text

EPS <- 1e-4
squeeze <- function(x) pmin(pmax(x, EPS), 1 - EPS)

stopifnot("pct_impervious_local" %in% names(urbanization_raw))

pt <- urbanization_raw %>%
  filter(!is.na(pct_impervious_local), !is.na(year),
         !is.na(longitude), !is.na(latitude)) %>%
  mutate(
    y_imperv = squeeze(pct_impervious_local / 100),
    Species  = factor(Species, levels = SPECIES_LEVELS),
    Dataset  = factor(Dataset, levels = c("iNaturalist", "VertNet"))
  )

species_list <- levels(droplevels(pt$Species))

SPECIES_LEVELS_REV <- rev(SPECIES_LEVELS)
species_rev <- function(x) factor(as.character(x), levels = SPECIES_LEVELS_REV)

# =============================================================================
# 2. Per-species rebalancing weight
# =============================================================================
# iNaturalist records scaled so the dataset's
# total weight equals VertNet's. VertNet records get weight 1. Balanced species
# (ratio ~ 1) end up ~unweighted; lopsided species get iNat heavily downweighted.

add_balance_weights <- function(d) {
  n_in <- sum(d$Dataset == "iNaturalist")
  n_vn <- sum(d$Dataset == "VertNet")
  if (n_in == 0 || n_vn == 0) { d$w <- 1; return(d) }
  d %>% mutate(w = ifelse(Dataset == "iNaturalist", n_vn / n_in, 1))
}

# =============================================================================
# 3. Fit one pooled, weighted GAM per species
# =============================================================================
# return the shared linear year slope 
# (+ a smooth model for the figure and non-linearity check).
fit_pooled <- function(df, label) {
  map_dfr(species_list, function(sp) {
    d <- df %>% filter(Species == sp) %>% droplevels()
    n_in <- sum(d$Dataset == "iNaturalist"); n_vn <- sum(d$Dataset == "VertNet")
    
    if (n_in < 50 || n_vn < 30 || n_distinct(d$year) < 5) {
      return(tibble(Species = sp, slope = NA_real_, se = NA_real_,
                    p_value = NA_real_, edf_year = NA_real_, p_smooth = NA_real_,
                    n_inat = n_in, n_vertnet = n_vn, ratio = n_in / max(n_vn, 1),
                    note = "insufficient data"))
    }
    
    d <- add_balance_weights(d) %>% mutate(year_c = year - mean(year))
    
    m_lin <- tryCatch(
      gam(y_imperv ~ year_c + Dataset + s(longitude, latitude, by = Dataset),
          family = betar(), weights = w, method = "REML", gamma = 1.4, data = d),
      error = function(e) NULL)
    m_smooth <- tryCatch(
      gam(y_imperv ~ s(year, k = 3) + Dataset + s(longitude, latitude, by = Dataset),
          family = betar(), weights = w, method = "REML", gamma = 1.4, data = d),
      error = function(e) NULL)
    
    if (is.null(m_lin)) {
      return(tibble(Species = sp, slope = NA_real_, se = NA_real_,
                    p_value = NA_real_, edf_year = NA_real_, p_smooth = NA_real_,
                    n_inat = n_in, n_vertnet = n_vn, ratio = n_in / n_vn,
                    note = "model failed"))
    }
    ct <- summary(m_lin)$p.table
    ss <- if (!is.null(m_smooth)) summary(m_smooth)$s.table else NULL
    yrow <- if (!is.null(ss)) grep("s\\(year\\)", rownames(ss)) else integer(0)
    
    tibble(
      Species  = sp,
      slope    = ct["year_c", "Estimate"],
      se       = ct["year_c", "Std. Error"],
      p_value  = ct["year_c", ncol(ct)],
      edf_year = if (length(yrow)) ss[yrow[1], "edf"]     else NA_real_,
      p_smooth = if (length(yrow)) ss[yrow[1], "p-value"] else NA_real_,
      n_inat   = n_in, n_vertnet = n_vn, ratio = n_in / n_vn,
      note     = NA_character_
    )
  }) %>%
    mutate(window = label, Species = factor(Species, levels = SPECIES_LEVELS))
}

res_1990 <- fit_pooled(pt %>% filter(year >= 1990), "1990-2025") #smaller year filter for smoother trends
res_all  <- fit_pooled(pt,                           "all years")

slopes <- bind_rows(res_1990, res_all) %>%
  filter(is.na(note)) %>%
  group_by(window) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(sig = p_adj < 0.05, Species = factor(Species, levels = SPECIES_LEVELS))

cat("\n=== Dataset-balanced pooled urban trend (logit slope per year) ===\n")
cat("slope < 0 = declining in present-day urban areas; > 0 = increasing\n")
cat("ratio = iNat/VertNet n (near 1 = naturally balanced, e.g. P. cornutum)\n\n")
slopes %>%
  select(window, Species, slope, se, p_adj, sig, n_inat, n_vertnet, ratio) %>%
  arrange(window, Species) %>%
  print(n = Inf, width = Inf)


# =============================================================================
# 4. Forest plot [Figure 6]
# =============================================================================

coef_df <- slopes %>%
  filter(window == "all years") %>%
  mutate(ci_low = slope - 1.96 * se, ci_high = slope + 1.96 * se,
         Species_y = species_rev(Species))

# --- Phylopic silhouettes -----------------------------------------------
# Pulled from species_display (defined near the top of the script) so the
# image, the display label, and the plotting order all stay in one place.
phylopic_lookup <- species_display %>%
  select(Species = full, phylopic_name)

coef_df <- coef_df %>% left_join(phylopic_lookup, by = "Species")

mirrored_species <- c("A. punctatus", "A. woodhousii", "I. nebulifer",
                      "Co. texanus", "Cr. collaris", "P. cornutum",
                      "Te. ornata", "Tr. scripta")

phylopic_df <- coef_df %>%
  distinct(Species, Species_y, phylopic_name) %>%
  mutate(height = ifelse(Species == "Tr. scripta", 0.35, 0.18))

phylopic_df_mirrored   <- phylopic_df %>% filter(Species %in% mirrored_species)
phylopic_df_unmirrored <- phylopic_df %>% filter(!Species %in% mirrored_species)

# x position for the silhouettes: just outside the left edge of the CIs, in
# the expanded margin created by scale_x_continuous() below. Because this
# point becomes part of the plotted data (not a fixed `limits=`), ggplot's
# automatic range + expansion will make room for it without clipping.
img_x <- min(coef_df$ci_low, na.rm = TRUE) -
  diff(range(c(coef_df$ci_low, coef_df$ci_high), na.rm = TRUE)) * 0.24

imperv_coef_fig <- ggplot(coef_df, aes(x = slope, y = Species_y, color = sig)) +
  geom_vline(xintercept = 0, color = "gray30", linewidth = 0.5, linetype = "dashed") +
  geom_errorbar(
    aes(xmin = ci_low, xmax = ci_high),
    height    = 0.25,
    linewidth = 0.9
  ) +
  geom_point(size = 3.6) +
  # Small mirrored silhouettes
  geom_phylopic(
    data = phylopic_df_mirrored,
    aes(x = img_x, y = Species_y, name = phylopic_name),
    inherit.aes = FALSE,
    color = "gray20",
    horizontal = TRUE,
    height = 0.40
  ) +
  scale_color_manual(values = c("FALSE" = "gray55", "TRUE" = "#CC79A7"),
                     labels = c("n.s.", "BH p < 0.05"), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.06, 0.08))) +
  scale_y_discrete(labels = SPECIES_LABELS,
                   expand = expansion(add = c(1.5, 0.6))) +
  annotate("text", x = -Inf, y = -0.1, label = "\u2190 declining in urban",
           hjust = -0.1, vjust = 0, size = 5.5, color = "gray10") +
  annotate("text", x =  Inf, y = -0.1, label = "increasing in urban \u2192",
           hjust = 1.1, vjust = 0, size = 5.5, color = "gray10") +
  labs(
    x = "Dataset-balanced change in present-day urban occurrence\n(logit impervious per year, \u03B2 \u00B1 95% CI; all years)",
    y = ""
  ) +
  theme_bw(base_size = 22) +
  theme(axis.text.y        = element_text(face = "italic", size = 20, color = "gray10"),
        axis.text.x        = element_text(size = 18, color = "gray10"),
        axis.title.x       = element_text(size = 20, color = "black"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        legend.position    = "none",
        legend.text        = element_text(size = 16, color = "gray10"),
        plot.margin        = margin(10, 10, 10, 40),
        aspect.ratio       = 0.6)

print(imperv_coef_fig)

ggsave(file.path(FIG_DIR, "fig6_imperv_coef_weighted.png"), imperv_coef_fig, width = 12, height = 8, units = "in", dpi = 600)
ggsave(file.path(FIG_DIR, "fig6_imperv_coef_weighted.pdf"), imperv_coef_fig, width = 12, height = 8, units = "in", dpi = 600)
