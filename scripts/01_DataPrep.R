# =============================================================================
# 01_DataPrep.R
# All data loading, cleaning, spatial joins, and derived analysis tables.
# Run this first; downstream scripts load the saved outputs.
#
# Inputs:  VertNet Final CSVs (per species) + new iNat CSVs (globbed by
#          genus_species stem), NLCD impervious raster
# Outputs: prepared/ directory with .rds files for each analysis table
# =============================================================================

library(sf)            # spatial operations, st_join, st_transform
library(here)          # project-root-relative paths
library(tidyverse)     # dplyr, ggplot2, readr, tidyr, purrr, lubridate, stringr
library(conflicted)    # resolve namespace conflicts
library(tigris)        # Census TIGER/Line county boundaries
library(terra)         # raster operations (NLCD)
library(exactextractr) # fast zonal statistics for rasters
library(scales)        # squish, percent, etc. (used in downstream plots)

conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(dplyr::select)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

# Paths resolve relative to the project root via here::here(), so the repo
# runs unchanged on any machine as long as the folder layout is preserved.
DATA_DIR  <- here::here("Data", "raw")   # occurrence CSVs (VertNet + iNat)
NLCD_PATH <- here::here("Data", "NLCD", "Annual_NLCD_FctImp_2024_CU_C1V1.tif")
OUT_DIR   <- here::here("Data", "prepared")
FIG_DIR   <- here::here("figures")

# Directory holding the new iNat CSVs (same folder as the occurrence CSVs)
INAT_DIR  <- DATA_DIR

# Year range we will include in analysis
YEAR_MIN <- 1800L
YEAR_MAX <- 2024L

# Species lookup: short filename tag -> binomial label
SPECIES <- list(
  punctatus   = "A. punctatus",
  woodhousii  = "A. woodhousii",
  nebulifer   = "I. nebulifer",
  cophosaurus = "Co. texanus",
  crotaphytus = "Cr. collaris",
  phrynosoma  = "P. cornutum",
  terrapene   = "Te. ornata",
  trachemys   = "Tr. scripta"
)

# iNat file lookup: match species tag -> genus_species stem used in new
# filenames (e.g. "Incilius_nebulifer_30856obs_8March2025.csv"). Files are
# globbed and regex-matched on the stem rather than reconstructed, because
# the raw names carry inconsistent obs counts, dates, and spaces.
INAT_STEMS <- list(
  punctatus   = "Anaxyrus_punctatus",
  woodhousii  = "Anaxyrus_woodhousii",   # fowleri merged in separately below
  nebulifer   = "Incilius_nebulifer",
  cophosaurus = "Cophosaurus_texanus",
  crotaphytus = "Crotaphytus_collaris",
  phrynosoma  = "Phrynosoma_cornutum",
  terrapene   = "Terrapene_ornata",
  trachemys   = "Trachemys_scripta"
)

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

# Build a filename from a species tag and source (VertNet only now)
make_fname <- function(tag, source, final = FALSE) {
  Tag <- paste0(toupper(substr(tag, 1, 1)), substr(tag, 2, nchar(tag)))
  if (final) file.path(DATA_DIR, sprintf("%s_%s_Final.csv", Tag, source))
  else       file.path(DATA_DIR, sprintf("%s_%s.csv",       Tag, source))
}

# Locate the single iNat CSV matching a genus_species stem; error if 0 or >1.
# Pattern tolerates the "_NNNNobs_" convention and any trailing date/spaces.
find_inat_file <- function(stem, dir = INAT_DIR) {
  all_csv <- list.files(dir, pattern = "\\.csv$", full.names = TRUE)
  base    <- basename(all_csv)
  hits    <- all_csv[grepl(paste0("^", stem, "_.*obs"), base, ignore.case = TRUE)]
  if (length(hits) == 0)
    stop(sprintf("No iNat file matching stem '%s' in %s", stem, dir))
  if (length(hits) > 1)
    stop(sprintf("Multiple iNat files matching stem '%s':\n  %s",
                 stem, paste(basename(hits), collapse = "\n  ")))
  hits
}

# Clean year column: coerce to integer and filter to valid range
clean_year <- function(df) {
  df %>%
    mutate(year = as.integer(.data$year)) %>%
    filter(
      !is.na(.data$year),
      .data$year >= YEAR_MIN,
      .data$year <= YEAR_MAX
    )
}

# Find date column
find_date_col <- function(df) {
  candidates <- grep("^observed", names(df), value = TRUE, ignore.case = TRUE)
  if (length(candidates) == 0)
    stop("No 'observed*' date column found")
  if ("observed_on" %in% candidates) return("observed_on")
  if ("observed_1" %in% candidates) return("observed_1")
  candidates[1]
}

# Parse year from iNat date column and clean
parse_inat_year <- function(df) {
  date_col <- find_date_col(df)
  df %>%
    mutate(year = as.integer(substr(.data[[date_col]], 1, 4))) %>%
    filter(!is.na(year), year >= YEAR_MIN, year <= YEAR_MAX)
}

# Spatial join: assign county to iNat records via point-in-polygon
assign_counties <- function(inat_data, counties_sf, species_name) {
  inat_sf <- inat_data %>%
    filter(!is.na(longitude), !is.na(latitude)) %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE) %>%
    st_transform(st_crs(counties_sf))
  
  result <- st_join(inat_sf, counties_sf, join = st_within) %>%
    st_drop_geometry() %>%
    mutate(Species = species_name)
  
  cat(sprintf("[%s] Assigned %d / %d iNat records to counties\n",
              species_name, sum(!is.na(result$county)), nrow(result)))
  result
}

cat("=== 01_DataPrep.R ===\n")
cat("Working directory:", getwd(), "\n\n")

# =============================================================================
# 1. Load VertNet Final CSVs + new iNat CSVs
# =============================================================================
cat("--- Loading VertNet Final (spatial) CSVs ---\n")

vertnet_spatial <- map(names(SPECIES), function(tag) {
  read.csv(make_fname(tag, "VertNet", final = TRUE)) %>% clean_year()
}) %>% set_names(names(SPECIES))

vertnet_spatial <- map(vertnet_spatial, function(df) {
  df %>%
    rename(
      latitude = decimallatitude,
      longitude = decimallongitude
    )
}) %>%
  set_names(names(SPECIES))

# --- New iNat CSVs, globbed by genus_species stem ----------------------------
# don't worry if this throws NA warning on year parse, the column still builds
cat("--- Loading new iNat CSVs (globbed by stem) ---\n")

inat_spatial_raw <- map(names(SPECIES), function(tag) {
  f <- find_inat_file(INAT_STEMS[[tag]])
  cat(sprintf("  [%s] -> %s\n", tag, basename(f)))
  read.csv(f) %>% parse_inat_year()
}) %>% set_names(names(SPECIES))

# --- Merge A. fowleri records into A. woodhousii -----------------------------
cat("--- Merging A. fowleri into A. woodhousii ---\n")

fowleri_file <- find_inat_file("Anaxyrus_fowleri")
cat(sprintf("  fowleri -> %s\n", basename(fowleri_file)))
fowleri_raw  <- read.csv(fowleri_file) %>% parse_inat_year()

# Align columns before binding (guards against column-order/extra-column drift)
common_cols <- base::intersect(names(inat_spatial_raw$woodhousii), names(fowleri_raw))
if (length(common_cols) < ncol(inat_spatial_raw$woodhousii)) {
  warning(sprintf("woodhousii/fowleri share %d of %d columns; binding on shared set",
                  length(common_cols), ncol(inat_spatial_raw$woodhousii)))
}

inat_spatial_raw$woodhousii <- bind_rows(
  inat_spatial_raw$woodhousii[, common_cols, drop = FALSE],
  fowleri_raw[, common_cols, drop = FALSE]
)

cat(sprintf("  woodhousii now has %d records (added %d fowleri)\n",
            nrow(inat_spatial_raw$woodhousii), nrow(fowleri_raw)))

# --- Override public coords with private coords where available --------------
# Some iNat records carry both latitude/longitude and private_latitude/
# private_longitude (obscured public coords + true location). latitude is
# always populated; private_* only sometimes. Where private_* is present,
# use it in place of the public value. Coerced to numeric first so blank
# strings become NA and don't spuriously override.
apply_private_coords <- function(df) {
  if (all(c("private_latitude", "private_longitude") %in% names(df))) {
    priv_lat <- suppressWarnings(as.numeric(df$private_latitude))
    priv_lon <- suppressWarnings(as.numeric(df$private_longitude))
    df$latitude  <- ifelse(!is.na(priv_lat), priv_lat, df$latitude)
    df$longitude <- ifelse(!is.na(priv_lon), priv_lon, df$longitude)
  }
  df
}

inat_spatial_raw <- map(inat_spatial_raw, apply_private_coords)

cat("  Applied private-coordinate override where available\n")

cat("  Loaded", length(vertnet_spatial), "VertNet +",
    length(inat_spatial_raw), "iNat datasets\n\n")

# =============================================================================
# 2. Fetch Texas county boundaries and assign counties to iNat records
# =============================================================================
cat("--- Fetching county boundaries (tigris) ---\n")
options(tigris_use_cache = TRUE)

texas_counties <- tigris::counties(
  state        = "TX",
  cb           = FALSE,   # full-resolution TIGER/Line for precise point-in-polygon
  year         = 2020,
  progress_bar = FALSE
) %>%
  rename(county = NAME) %>%
  st_transform(3083)      # EPSG:3083 Texas Centric Albers Equal Area

cat("  ", nrow(texas_counties), "county polygons loaded (CRS:", st_crs(texas_counties)$epsg, ")\n\n")

cat("--- Assigning counties to iNat records ---\n")
inat_spatial <- mapply(
  assign_counties,
  inat_data    = inat_spatial_raw,
  species_name = unlist(SPECIES, use.names = FALSE),
  MoreArgs     = list(counties_sf = texas_counties),
  SIMPLIFY     = FALSE
) %>% set_names(names(SPECIES))
cat("\n")
# We expect a few dozen for each species to not be included - these seem to be records which actually did not occur within Texas
# But were included in Texas (possibly based on iNaturalists's point obscuring method)

# =============================================================================
# 3. Build urbanization_raw: all individual records in one long table
# =============================================================================
cat("--- Building urbanization_raw ---\n")

SHARED_COLS <- c("year", "county", "latitude", "longitude")

urbanization_raw <- map_dfr(names(SPECIES), function(tag) {
  sp <- SPECIES[[tag]]
  bind_rows(
    vertnet_spatial[[tag]] %>%
      select(any_of(SHARED_COLS)) %>%
      mutate(Dataset = "VertNet", Species = sp),
    inat_spatial[[tag]] %>%
      select(any_of(SHARED_COLS)) %>%
      mutate(Dataset = "iNaturalist", Species = sp)
  )
}) %>%
  filter(!is.na(year))

cat("  ", nrow(urbanization_raw), "total records\n")
cat("  ", n_distinct(urbanization_raw$Species), "species,",
    n_distinct(urbanization_raw$county), "counties\n\n")

# Stable per-record id, so point-extracted urban values can be joined back
urbanization_raw <- urbanization_raw %>%
  mutate(record_id = row_number())

cat("--- Building temporal_data ---\n")

VERTNET_COLLECTOR <- "recordedby"
INAT_COLLECTOR    <- "user_login"

temporal_data <- map_dfr(names(SPECIES), function(tag) {
  sp <- SPECIES[[tag]]
  
  vn <- vertnet_spatial[[tag]] %>%
    filter(!is.na(year)) %>%
    transmute(year,
              collector = .data[[VERTNET_COLLECTOR]],
              Dataset   = "VertNet",
              Species   = sp)
  
  inat <- inat_spatial[[tag]] %>%
    filter(!is.na(year)) %>%
    transmute(year,
              collector = .data[[INAT_COLLECTOR]],
              Dataset   = "iNaturalist",
              Species   = sp)
  
  bind_rows(vn, inat)
}) %>%
  group_by(Species, Dataset, year) %>%
  summarise(
    RecordCount      = n(),
    UniqueCollectors = n_distinct(collector[!is.na(collector) & collector != ""]),
    .groups          = "drop"
  )

cat("  ", nrow(temporal_data), "Species x Dataset x year rows\n")
cat("   Year range:", min(temporal_data$year), "-", max(temporal_data$year), "\n\n")

# =============================================================================
# 4. NLCD impervious surface: county means + point-level 1km neighborhood
# =============================================================================
cat("--- Extracting NLCD impervious surface ---\n")

NLCD_BUFFER_M <- 1000   # neighborhood radius for point extraction (metres)

nlcd_imperv         <- rast(NLCD_PATH)
texas_counties_proj <- st_transform(texas_counties, crs(nlcd_imperv))
nlcd_imperv_tx      <- crop(nlcd_imperv, ext(vect(texas_counties_proj)))

# --- 4a. County means (still used by county_counts / master aggregations) -----
county_imperv <- texas_counties_proj %>%
  mutate(pct_impervious = exact_extract(nlcd_imperv_tx, texas_counties_proj, "mean")) %>%
  st_transform(3083)

cat("  Top 10 most impervious counties:\n") #sanity check
county_imperv %>%
  st_drop_geometry() %>%
  select(county, pct_impervious) %>%
  arrange(desc(pct_impervious)) %>%
  head(10) %>%
  print()
cat("\n")

# --- 4b. Point-level 1km-neighborhood mean impervious -------------------------
# Every record is point-level (no-coordinate records were dropped upstream).
# For each record, mean impervious within NLCD_BUFFER_M. NLCD is a contemporary
# snapshot, so this is STATIC present-day urban cover at the record's location
# (the fixed-footprint definition: "is this spot urban TODAY").
cat("--- Point-level impervious (", NLCD_BUFFER_M, "m neighborhood) ---\n", sep = "")

pts_sf <- urbanization_raw %>%
  filter(!is.na(longitude), !is.na(latitude)) %>%
  select(record_id, longitude, latitude) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(crs(nlcd_imperv_tx)) %>%
  st_buffer(NLCD_BUFFER_M)

pts_sf$pct_impervious_local <-
  exact_extract(nlcd_imperv_tx, pts_sf, "mean", progress = FALSE)

urbanization_raw <- urbanization_raw %>%
  left_join(st_drop_geometry(pts_sf) %>% select(record_id, pct_impervious_local),
            by = "record_id")

cat("  Extracted", sum(!is.na(urbanization_raw$pct_impervious_local)),
    "point buffers. Range:",
    round(min(urbanization_raw$pct_impervious_local, na.rm = TRUE), 2), "-",
    round(max(urbanization_raw$pct_impervious_local, na.rm = TRUE), 2), "\n\n")

# Clean up large raster objects
rm(nlcd_imperv, nlcd_imperv_tx, texas_counties_proj, pts_sf)

# =============================================================================
# 5. County counts (wide format) with impervious cover
# =============================================================================
cat("--- Building county_counts ---\n")

county_counts <- urbanization_raw %>%
  filter(!is.na(county)) %>%
  group_by(Species, county, Dataset) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = Dataset, values_from = n, values_fill = 0L) %>%
  mutate(total = VertNet + iNaturalist) %>%
  filter(total > 0) %>%
  left_join(
    county_imperv %>% st_drop_geometry() %>% select(county, pct_impervious),
    by = "county"
  ) %>%
  mutate(prop_inat = iNaturalist / total)

# Standardise impervious cover, keeping scaling parameters for prediction grids
imperv_scaled                  <- scale(county_counts$pct_impervious)
county_counts$pct_impervious_z <- as.numeric(imperv_scaled)
imperv_center                  <- attr(imperv_scaled, "scaled:center")
imperv_scale                   <- attr(imperv_scaled, "scaled:scale")

cat("  ", nrow(county_counts), "species × county rows\n\n")

# =============================================================================
# 6. Master table: county × year × dataset × species (for GAMs)
# =============================================================================
cat("--- Building master ---\n")

master <- urbanization_raw %>%
  filter(!is.na(county), !is.na(year)) %>%
  group_by(county, year, Dataset, Species) %>%
  summarise(n_records = n(), .groups = "drop") %>%
  rename(county_name = county) %>%
  left_join(
    county_imperv %>% st_drop_geometry() %>%
      select(county, pct_impervious) %>%
      rename(county_name = county),
    by = "county_name"
  ) %>%
  mutate(
    year_c           = year - min(year),
    year_z           = as.numeric(scale(year)),
    pct_impervious_z = as.numeric(scale(pct_impervious)),
    Dataset          = as.factor(Dataset),
    county_f         = as.factor(county_name)
  )

cat("  ", nrow(master), "county × year × dataset × species rows\n\n")

# =============================================================================
# 7. County deltas and z-scores (for county delta analysis)
# =============================================================================
cat("--- Computing county deltas ---\n")

prepare_county_data <- function(vertnet_data, inat_data) {
  vertnet_counties <- vertnet_data %>%
    filter(!is.na(county)) %>%
    group_by(county) %>%
    summarise(vertnet_count = n(), .groups = "drop")
  
  inat_counties <- inat_data %>%
    filter(!is.na(county)) %>%
    group_by(county) %>%
    summarise(inat_count = n(), .groups = "drop")
  
  full_join(vertnet_counties, inat_counties, by = "county") %>%
    mutate(
      vertnet_count = replace_na(vertnet_count, 0),
      inat_count    = replace_na(inat_count, 0),
      total_count   = vertnet_count + inat_count
    ) %>%
    filter(total_count > 0)
}

county_data_all <- map2_dfr(
  vertnet_spatial, inat_spatial, prepare_county_data,
  .id = "tag"
) %>%
  mutate(Species = unlist(SPECIES[tag])) %>%
  select(-tag)

county_delta_combined <- county_data_all %>%
  mutate(
    delta_simple = inat_count - vertnet_count,
    prop_delta   = (inat_count - vertnet_count) / total_count
  ) %>%
  group_by(Species) %>%
  mutate(
    z_score_simple = (delta_simple - mean(delta_simple)) / sd(delta_simple),
    z_score_prop   = (prop_delta   - mean(prop_delta))   / sd(prop_delta)
  ) %>%
  ungroup()

cat("  ", nrow(county_delta_combined), "species × county delta rows\n")
cat("  ", n_distinct(county_delta_combined$county), "counties represented\n\n")

# =============================================================================
# 8. Time-filtered county data (VertNet pre-2008, iNat 2008+)
# =============================================================================
cat("--- Building time-filtered county data ---\n")

prepare_county_data_filtered <- function(vertnet_data, inat_data) {
  vertnet_counties <- vertnet_data %>%
    filter(!is.na(county), year < 2008) %>%
    group_by(county) %>%
    summarise(vertnet_count = n(), .groups = "drop")
  
  inat_counties <- inat_data %>%
    filter(!is.na(county), year >= 2008) %>%
    group_by(county) %>%
    summarise(inat_count = n(), .groups = "drop")
  
  full_join(vertnet_counties, inat_counties, by = "county") %>%
    mutate(
      vertnet_count = replace_na(vertnet_count, 0),
      inat_count    = replace_na(inat_count, 0),
      total_count   = vertnet_count + inat_count
    ) %>%
    filter(total_count > 0)
}

county_data_filtered <- map2_dfr(
  vertnet_spatial, inat_spatial, prepare_county_data_filtered,
  .id = "tag"
) %>%
  mutate(Species = unlist(SPECIES[tag])) %>%
  select(-tag)

cat("  ", nrow(county_data_filtered), "species × county rows (time-filtered)\n\n")


# =============================================================================
# 9. Standardize Species factor order across all objects
# =============================================================================
cat("--- Setting Species factor levels ---\n")

SPECIES_LEVELS <- unlist(SPECIES, use.names = FALSE)

factor_species <- function(df) {
  df %>% mutate(Species = factor(Species, levels = SPECIES_LEVELS))
}

urbanization_raw      <- factor_species(urbanization_raw)
county_counts         <- factor_species(county_counts)
master                <- factor_species(master)
county_delta_combined <- factor_species(county_delta_combined)
county_data_filtered  <- factor_species(county_data_filtered)
temporal_data         <- factor_species(temporal_data)

cat("  Factor levels:", paste(SPECIES_LEVELS, collapse = ", "), "\n\n")



# =============================================================================
# 10. Get total numbers of records and observers after filtering for county and year
# =============================================================================
county_year_summary <- map_dfr(names(SPECIES), function(tag) {
  sp <- SPECIES[[tag]]
  
  vn <- vertnet_spatial[[tag]] %>%
    filter(!is.na(county), !is.na(year))
  
  # Different species' files use different truncated column names
  # for the same field - check both, use whichever one actually exists.
  inst_col <- base::intersect(c("institutioncode", "institut_1"), names(vn))[1]
  
  vn_summary <- vn %>%
    summarise(
      Species            = sp,
      Dataset            = "VertNet",
      RecordCount        = n(),
      UniqueCollectors   = n_distinct(recordedby[!is.na(recordedby) & recordedby != ""]),
      InstitutionColumn  = ifelse(is.na(inst_col), NA_character_, inst_col),
      UniqueInstitutions = if (!is.na(inst_col)) {
        n_distinct(.data[[inst_col]][!is.na(.data[[inst_col]]) & .data[[inst_col]] != ""])
      } else {
        NA_integer_
      }
    )
  
  inat_data <- inat_spatial[[tag]] %>%
    filter(!is.na(county), !is.na(year))
  
  inat_summary <- inat_data %>%
    summarise(
      Species            = sp,
      Dataset            = "iNaturalist",
      RecordCount        = n(),
      UniqueCollectors   = n_distinct(user_login[!is.na(user_login) & user_login != ""]),
      InstitutionColumn  = NA_character_,
      UniqueInstitutions = NA_integer_
    )
  
  bind_rows(vn_summary, inat_summary)
})

# =============================================================================
# 11. Save all prepared objects
# =============================================================================
cat("--- Saving to", OUT_DIR, "---\n")

saveRDS(SPECIES,                file.path(OUT_DIR, "SPECIES.rds"))
saveRDS(SPECIES_LEVELS,         file.path(OUT_DIR, "SPECIES_LEVELS.rds"))
saveRDS(texas_counties,         file.path(OUT_DIR, "texas_counties.rds"))
saveRDS(county_imperv,          file.path(OUT_DIR, "county_imperv.rds"))
saveRDS(urbanization_raw,       file.path(OUT_DIR, "urbanization_raw.rds"))
saveRDS(county_counts,          file.path(OUT_DIR, "county_counts.rds"))
saveRDS(master,                 file.path(OUT_DIR, "master.rds"))
saveRDS(county_delta_combined,  file.path(OUT_DIR, "county_delta_combined.rds"))
saveRDS(county_data_filtered,   file.path(OUT_DIR, "county_data_filtered.rds"))
saveRDS(temporal_data,          file.path(OUT_DIR, "temporal_data.rds"))
saveRDS(vertnet_spatial,        file.path(OUT_DIR, "vertnet_spatial.rds"))
saveRDS(inat_spatial,           file.path(OUT_DIR, "inat_spatial.rds"))

# Save impervious scaling parameters (needed for prediction grids)
saveRDS(list(center = imperv_center, scale = imperv_scale),
        file.path(OUT_DIR, "imperv_scaling.rds"))

cat("  Done. Saved", length(list.files(OUT_DIR, pattern = "\\.rds$")), "files.\n\n")

# Print session info for reproducibility
cat("--- Session Info ---\n")
sessionInfo()