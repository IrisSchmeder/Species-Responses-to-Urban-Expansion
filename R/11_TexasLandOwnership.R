# ---------------------------------------------------------------------------
# 11_TexasLandOwnership.R
# Percent of Texas in public vs. private ownership, from PAD-US 4.1
# Data: USGS GAP (2024) PAD-US 4.1, state extract for TX
# Source data can be found at ScienceBase: https://www.sciencebase.gov/catalog/item/65294599d34e44db0e2ed7cf
#
# Paths resolve relative to the project root via here::here(); the PAD-US
# geodatabase lives in data/PADUS4_1_State_TX_GDB_KMZ/PADUS4_1_StateTX.gdb.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(here)                       # project-root-relative paths
  library(sf)
  library(dplyr)
})

sf_use_s2(FALSE)                      # planar ops after projecting
ALBERS   <- 3083                      # NAD83 / Texas Centric Albers Equal Area (m)
M2_ACRE  <- 4046.8564224

gdb <- here::here("data", "PADUS4_1_State_TX_GDB_KMZ", "PADUS4_1_StateTX.gdb")


LAYER <- "PADUS4_1Fee_State_TX"


# --- 1. Read, linearize curve geometries, validate -------------------------
fee_raw <- st_read(gdb, layer = LAYER, quiet = TRUE)

# ~300 features come in as MULTISURFACE (curves); cast to MULTIPOLYGON.
fee <- st_cast(fee_raw, "MULTIPOLYGON")

bad <- which(!st_is_valid(fee))
if (length(bad) > 0) {
  st_geometry(fee)[bad] <- st_make_valid(st_geometry(fee)[bad])
  still_bad <- which(!st_is_valid(fee))
  if (length(still_bad) > 0)
    st_geometry(fee)[still_bad] <- st_buffer(st_geometry(fee)[still_bad], 0)
}
stopifnot(all(st_is_valid(fee)))

# PAD-US TX native CRS is USGS Contiguous Albers Equal Area -> area is valid.


# --- 2. Naive acreage by owner type (overlaps NOT removed) -----------------
naive <- fee |>
  st_drop_geometry() |>
  group_by(Own_Type) |>
  summarise(acres_naive = sum(GIS_Acres, na.rm = TRUE), n = n()) |>
  arrange(desc(acres_naive))
print(as.data.frame(naive))


# --- 3. State boundary (denominator + clip mask) ---------------------------
# tigris needs network access for download
library(tigris)
options(tigris_use_cache = TRUE)

tx <- states(cb = FALSE, year = 2023) |>
  filter(STUSPS == "TX") |>
  st_transform(st_crs(fee))

tx_land     <- st_geometry(tx) |> st_make_valid()
aland_acres <- as.numeric(tx$ALAND) / M2_ACRE   # land only, excludes water


# --- 4. Public union, clipped to land (removes submerged SLB tracts) -------
public_types <- c("FED", "STAT", "LOC", "DIST", "JNT", "TRIB")

pub_union <- fee |>
  filter(Own_Type %in% public_types) |>
  st_geometry() |>
  st_union()

pub_clip  <- st_intersection(pub_union, tx_land)   # drops offshore/bay bottoms which appear in this dataset (These are included in the dataset becasue they are in public trust and some areas serve as oil/gas extraction sites which are owned by our state's School Land Board)
pub_acres <- as.numeric(st_area(pub_clip)) / M2_ACRE

# UNK, unioned and clipped separately, so it can be added or not.
unk_union <- fee |> filter(Own_Type == "UNK") |> st_geometry() |> st_union()
unk_clip  <- st_intersection(unk_union, tx_land)
# remove any area already inside public before counting, to avoid double count
unk_only  <- st_difference(unk_clip, pub_clip)
unk_acres <- as.numeric(st_area(unk_only)) / M2_ACRE


# --- 5. Public/private range -----------------------------------------------
# UNK is the designation for land which USGS cannot determine managing/owning entity - thus we do a calculation where all UNK is treated as private and all UNK is treated as public
# Lower bound on private: treat all UNK as public.
# Upper bound on private: treat all UNK as private.
pub_lo <- pub_acres                     # UNK counted as private
pub_hi <- pub_acres + unk_acres         # UNK counted as public

result <- data.frame(
  scenario     = c("UNK as private (min public)", "UNK as public (max public)"),
  public_acres = c(pub_lo, pub_hi),
  pct_public   = 100 * c(pub_lo, pub_hi) / aland_acres,
  pct_private  = 100 - 100 * c(pub_lo, pub_hi) / aland_acres
)

cat("\nState land area (ALAND):", round(aland_acres), "acres\n")
cat("UNK (clipped, non-overlapping):", round(unk_acres), "acres\n\n")
print(result, digits = 4, row.names = FALSE)

cat(sprintf(
  "\nTexas is approximately %.1f-%.1f%% privately owned land.\n",
  min(result$pct_private), max(result$pct_private)))