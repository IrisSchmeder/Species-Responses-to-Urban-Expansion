# =============================================================================
# 09_TexasFeaturesMap.R
# Produces map of biotic provinces, urbanized areas, major cities, 
# and the top four museums represented in VertNet herpetology records in Texas.
#
# DEPENDENCIES: run 01_DataPrep.R first
#
# OUTPUT: texas_provinces_map.pdf / .png
#         [SI Figure 2]
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
PREP_DIR    <- here::here("data", "prepared")
FIG_DIR     <- here::here("figures")
RESULTS_DIR <- here::here("results")
dir.create(FIG_DIR,     showWarnings = FALSE, recursive = TRUE)
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

texas_counties <- readRDS(file.path(PREP_DIR, "texas_counties.rds"))

# Environmental layers each live in their own folder under data/
# ecoregions_shp <- st_read(here::here("data", "tx_eco_l3", "tx_eco_l3.shp"))
bioticprovinces_shp <- st_read(here::here("data", "TexasBioticProvinces",
                                          "TexasBioticProvinces.shp"))
urban_shp <- st_read(here::here("data", "TxDOT_UrbanizedAreas",
                                "Urbanized_Area.shp"))

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
# State outline is now derived from the already-prepared
# texas_counties (a single st_union of the county polygons), instead of a
# fresh tigris::states() call.
tx <- texas_counties %>%
  st_transform(4326) %>%
  st_make_valid() %>%
  st_union() %>%
  st_as_sf()

# eco <- ecoregions_shp %>% st_transform(4326) %>% st_make_valid()
eco <- bioticprovinces_shp %>% st_transform(4326) %>% st_make_valid()
eco_name_col <- intersect(names(eco), c("ProvName"))[1]
eco <- eco %>% mutate(eco_name = .data[[eco_name_col]])

# st_intersection()/st_union() on lon-lat (4326) data uses the S2 spherical
# engine by default, which is much stricter about micro duplicate-vertex
# edges than the old planar GEOS engine - it can flag a geometry as invalid
# even right after st_make_valid(), because the intersection itself can
# introduce a new degenerate edge. Switching to the planar engine for this
# step (and re-validating right after the intersection, since that's the
# operation producing the bad geometry) avoids it.
sf::sf_use_s2(FALSE)

eco_tx <- st_intersection(eco, tx) %>%
  st_make_valid() %>%
  group_by(eco_name) %>%
  summarize(geometry = st_union(geometry), .groups = "drop") %>%
  st_as_sf(crs = 4326)

sf::sf_use_s2(TRUE)

urban_poly <- urban_shp %>% st_transform(4326) %>% st_make_valid()

# Transform
TARGET_CRS <- 32139

tx         <- st_transform(tx,         TARGET_CRS)
eco_tx     <- st_transform(eco_tx,     TARGET_CRS)
urban_poly <- st_transform(urban_poly, TARGET_CRS)
cities_sf  <- st_transform(cities_raw, TARGET_CRS)
muse_sf    <- st_transform(muse_raw,   TARGET_CRS)

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
# 2. Make map of ecoregions and features [SI Figure 2]
# =============================================================================
map_bbox_plot <- st_bbox(tx)
pad_x <- 0.02 * (map_bbox_plot["xmax"] - map_bbox_plot["xmin"])
pad_y <- 0.02 * (map_bbox_plot["ymax"] - map_bbox_plot["ymin"])

p_provinces <- ggplot() +
  geom_sf(data = eco_tx,     aes(fill = eco_name), color = "grey72", linewidth = 0.2) +
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
    values = c("Urbanized area" = "green"),
    name   = "",
    guide  = guide_legend(
      order         = 1,
      override.aes  = list(linewidth = 2.5, shape = NA, fill = NA)
    )
  ) +
  scale_shape_manual(
    values = c("City" = 21, "Museum" = 22),
    name   = ""
  ) +
  scale_fill_manual(
    name   = "Biotic Province",
    values = c(
      "Austroriparian" = "#633372",  # dark purple, colors from palette Signac
      "Balconian"      = "#9b3441",  # dark maroon
      "Chihuahuan"     = "#FE9B00",  # orange
      "Kansan"         = "#d8443c",  # red-orange
      "Navahonian"     = "#9f5691",  # purple
      "Tamaulipan"     = "#2b9b81",  # teal green
      "Texan"          = "#de597c"   # magenta/pink
    )
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

p_provinces


# Save
ggsave(file.path(FIG_DIR, "texas_provinces_map.png"), p_provinces,
       width = 16, height = 12, dpi = 300)
ggsave(file.path(FIG_DIR, "texas_provinces_map.pdf"), p_provinces,
       width = 16, height = 12, dpi = 300)

summary_tbl <- eco_tx %>% st_drop_geometry() %>% count(eco_name, sort = TRUE)
write_csv(summary_tbl, file.path(RESULTS_DIR, "texas_provinces_summary.csv"))
