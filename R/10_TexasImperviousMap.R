# =============================================================================
# 10_TexasImperviousMap.R
# Produces map of percent impervious cover by county, major cities, 
# and the top four museums represented in VertNet herpetology records in Texas.
#
# DEPENDENCIES: run 01_DataPrep.R first
#
# OUTPUT: texas_impervious_map.pdf / .png
#         [SI Figure 3]
# =============================================================================

library(tidyverse)
library(here)          # project-root-relative paths
library(sf)
library(ggspatial)
library(patchwork)
library(cowplot)
library(paletteer)
library(shadowtext)
library(tigris)
library(units)
library(ggrepel)
conflicted::conflicts_prefer(lubridate::intersect)

options(tigris_use_cache = TRUE)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

# Paths resolve relative to the project root via here::here().
PREP_DIR <- here::here("data", "prepared")
FIG_DIR  <- here::here("figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

texas_counties <- readRDS(file.path(PREP_DIR, "texas_counties.rds"))
county_imperv  <- readRDS(file.path(PREP_DIR, "county_imperv.rds"))
urban_shp <- st_read(here::here("Data", "TxDOT_UrbanizedAreas",
                                "TxDOT_Urbanized_Areas.shp"))

muse_raw <- tibble(
  name = c("Angelo State", "A&M", "UT", "UTEP"),
  lon  = c(-100.4590, -96.3364, -97.7335, -106.5060),
  lat  = c(  31.4382,  30.6187,  30.2850,  31.7726)
) %>% st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

cities_raw <- tibble(
  name = c("Houston", "Dallas", "San Antonio", "Austin", "El Paso", "Lubbock"),
  lon  = c( -95.3701,  -96.8089,    -98.4946, -97.7431,  -106.4850, -101.8552),
  lat  = c(  29.7601,   32.7792,     29.4252,  30.2672,    31.7619,   33.5845)
) %>% st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

# =============================================================================
# 1. Make boundaries and perform transformations
# =============================================================================
tx <- texas_counties %>%
  st_transform(4326) %>%
  st_make_valid() %>%
  st_union() %>%
  st_as_sf()

urban_poly <- urban_shp %>% st_transform(4326) %>% st_make_valid()

# Transform
TARGET_CRS <- 32139

tx            <- st_transform(tx,            TARGET_CRS)
county_imperv <- st_transform(county_imperv, TARGET_CRS)
urban_poly    <- st_transform(urban_poly,    TARGET_CRS)
cities_sf     <- st_transform(cities_raw,    TARGET_CRS)
muse_sf       <- st_transform(muse_raw,      TARGET_CRS)

muse_overlap    <- muse_sf %>% filter(name %in% c("UTEP", "UT")) %>% arrange(desc(name))
muse_standalone <- muse_sf %>% filter(!name %in% c("UTEP", "UT"))

overlap_coords <- st_coordinates(muse_overlap)

offset_x <- c(80000,      0)   # UTEP east, UT no x shift
offset_y <- c(    0,  25000)   # UTEP no y shift, UT north
muse_overlap_offset <- muse_overlap %>%
  st_set_geometry(
    st_sfc(
      mapply(function(x, y, dx, dy) st_point(c(x + dx, y + dy)),
             overlap_coords[, 1], overlap_coords[, 2],
             offset_x, offset_y,
             SIMPLIFY = FALSE),
      crs = TARGET_CRS
    )
  )

city_label <- tibble(
  name    = c("Houston",  "Dallas",  "San Antonio", "Austin",  "El Paso",  "Lubbock"),
  label_x = c( 1179927,    1050252,      868463,       922218,   138527,     578598),
  label_y = c(3058100,   3354333,     3004704,      3059475,  3249385,    3436438),
  hjust   = c(0,          0,            1,             1,        0,           0)
)

muse_label <- tibble(
  name    = c("Angelo State", "A&M",    "UT",      "UTEP"),
  label_x = c(  664054,        992095,   925095,    169606),
  label_y = c( 3199393,       3122421,  3100469,  3215680),
  hjust   = c(1,               0,         1,         0)
)


city_color   <- "white"
museum_color <- "white"

# =============================================================================
# 2. Make map of impervious cover and features [SI Figure 3]
# =============================================================================
map_bbox_plot <- st_bbox(tx)
pad_x <- 0.02 * (map_bbox_plot["xmax"] - map_bbox_plot["xmin"])
pad_y <- 0.02 * (map_bbox_plot["ymax"] - map_bbox_plot["ymin"])

p_imperv <- ggplot() +
  geom_sf(data = county_imperv, aes(fill = pct_impervious), color = "grey72", linewidth = 0.2) +
  geom_sf(data = tx,         fill = NA, color = "grey20", linewidth = 0.8) +
  geom_sf(data = urban_poly, aes(color = "Urbanized area"), fill = NA, linewidth = 0.75) +
  geom_sf(data = cities_sf,
          aes(shape = "City"),
          color = "black", fill = city_color, size = 5.3, stroke = 0.8) +
  geom_sf(data = muse_standalone,
          aes(shape = "Museum"),
          color = "black", fill = museum_color, size = 5.5, stroke = 0.8) +
  geom_sf(data = muse_overlap,
          shape = 22, color = "black", fill = NA, size = 6.7, stroke = 1.2) +
  geom_shadowtext(
    data     = city_label,
    aes(x = label_x, y = label_y, label = name, hjust = hjust),
    size     = 6.8,
    color    = "white",
    bg.colour = "black",
    bg.r     = 0.15,
    fontface = "bold"
  ) +
  geom_shadowtext(
    data     = muse_label,
    aes(x = label_x, y = label_y, label = name, hjust = hjust),
    size     = 6.8,
    color    = "white",
    bg.colour = "black",
    bg.r     = 0.15,
    fontface = "bold"
  ) +
  scale_color_manual(
    values = c("Urbanized area" = "darkorchid1"),
    name   = "",
    guide  = guide_legend(
      order         = 1,
      override.aes  = list(linewidth = 2.5, shape = NA, fill = NA)
    )
  ) +
  scale_shape_manual(
    values = c("City" = 21, "Museum" = 22),
    name   = "",
    guide  = guide_legend(order = 2)
  ) +
  scale_fill_paletteer_c(
    "grDevices::YlOrRd",
    name      = "% Impervious\nSurface Cover",
    direction = -1,  # low values light/yellow, high values dark red - darkest = most impervious
    guide     = guide_colorbar(order = 3, barheight = unit(4, "cm"), barwidth = unit(0.6, "cm"))
  ) +
  coord_sf(
    xlim   = c(map_bbox_plot["xmin"] - pad_x, map_bbox_plot["xmax"] + pad_x),
    ylim   = c(map_bbox_plot["ymin"] - pad_y, map_bbox_plot["ymax"] + pad_y),
    expand = FALSE
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.title        = element_blank(),
    legend.position   = "right",
    legend.title      = element_text(size = 25, face = "bold"),
    legend.text       = element_text(size = 25),
    legend.key.height = unit(1.35, "cm"),
    legend.key.width  = unit(1.0, "cm")
  )

p_imperv


# Save
ggsave(file.path(FIG_DIR, "texas_impervious_map.png"), p_imperv,
       width = 16, height = 12, dpi = 300)
ggsave(file.path(FIG_DIR, "texas_impervious_map.pdf"), p_imperv,
       width = 16, height = 12, dpi = 300)

