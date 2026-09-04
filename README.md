# Code for Schmeder et al.
This repository contains the code for the manuscript "Herpetofaunal Responses to Urban Expansion: Insights from Combining Two Texas Occurrence Datasets" submitted to Conservation Biology. It provides a pipeline to combine historical museum records with modern community science observations to model species occurrences in response to urbanization. The example demonstrated here uses eight herpetofauna species native to Texas.

## Abstract


## Data
The analysis uses downloaded data available from VertNet and iNaturalist. We also used 2024 Fractional Impervious Surface (CONUS) from NLCD (https://www.mrlc.gov/data) and an urbanization raster from Texas Department of Transportation (https://gis-txdot.opendata.arcgis.com/datasets/6aeee12f605d4b0b9fc74b31d2ea4ea5_0/about).

All paths are resolved relative to the project root using the [`here`](https://here.r-lib.org/) package, so the pipeline runs unchanged on any machine as long as the folder layout below is preserved. Open `TexasPairedSpecies.Rproj` in RStudio (or run from the repo root) and `here::here()` anchors to the project directory automatically.


### Repository layout

```
.
├── R/                          # all analysis scripts
├── data/
│   ├── raw/                    # occurrence CSVs (VertNet + iNaturalist) — see below
│   ├── NLCD/                   # NLCD fractional impervious raster (.tif)
│   ├── TexasBioticProvinces/   # Texas biotic provinces shapefile
│   ├── TxDOT_UrbanizedAreas/   # TxDOT urbanized-areas shapefile
│   └── prepared/               # generated .rds objects (created by 01)
├── figures/                    # generated figures (created by scripts)
├── results/                    # generated summary tables (created by scripts)
├── TexasPairedSpecies.Rproj
├── README.md
├── LICENSE
└── .gitignore
```

### Where to put the input data

- **`data/raw/`** — one VertNet CSV per species named `<Tag>_VertNet_Final.csv` (e.g. `Punctatus_VertNet_Final.csv`), and one iNaturalist CSV per species named with the `Genus_species` stem plus an observation count, e.g. `Incilius_nebulifer_30856obs_8March2025.csv` (matched by regex on the stem, so trailing counts/dates/spaces are tolerated). An `Anaxyrus_fowleri_*obs_*.csv` file is also required — its records are merged into *A. woodhousii*.
- **`data/NLCD/`** — `Annual_NLCD_FctImp_2024_CU_C1V1.tif`.
- **`data/TexasBioticProvinces/`** — `TexasBioticProvinces.shp` (+ sidecar files).
- **`data/TxDOT_UrbanizedAreas/`** — the TxDOT urbanized-areas shapefile (+ sidecars). Note: script 09 reads `Urbanized_Area.shp` and script 10 reads `TxDOT_Urbanized_Areas.shp` — place whichever your download provides and adjust the filename in the relevant script if it differs.

Texas county boundaries are downloaded automatically via `tigris`; the environmental layers above must be supplied.

### Running

Run scripts in numeric order. `01_DataPrep.R` must be run first — it writes the `.rds` objects in `data/prepared/` that every downstream script loads.

```r
source("R/01_DataPrep.R")   # then 03, 04, ... 10 in any order
```

`02_DatasetInfoTable.Rmd` is knit to PDF (Table 1).

## Scripts

### Data Prep & Info
- [R/01_DataPrep.R](R/01_DataPrep.R): Reads all csv files and builds global objects.
- [R/02_DatasetInfoTable.Rmd](R/02_DatasetInfoTable.Rmd): Makes table with dataset filtering information (table 1).

### Observational Bias
- [R/03_RecordAccumulation.R](R/03_RecordAccumulation.R): Makes record accumulation curves over time (figure 1).
- [R/04_TemporalBiasGAMs.R](R/04_TemporalBiasGAMs.R): Makes negative-binomial GAMs for record counts and unique observers over time (figure 2).
- [R/05_RatioDataset.R](R/05_RatioDataset.R): Calculates the proportion of records from each dataset within each county ranked by iNaturalist proportion (figure 3).

### Species Response
- [R/06_CountyDeltaAnalysis.R](R/06_CountyDeltaAnalysis.R): Makes County-level dataset delta analysis: paired t-tests and z-score maps (figure 4 and figure 5)
- [R/07_UrbanizationAnalysis.R](R/07_UrbanizationAnalysis.R): Calculates urbanization trends and impervious-surface GAMs (figure 6).

### Supplementary Figures
- [R/08_RecordsPerObserver.R](R/08_RecordsPerObserver.R): Calculates records per observer over time, by species and dataset (SI figure 1).
- [R/09_TexasProvincesMap.R](R/09_TexasProvincesMap.R): Produces map of level III ecoregions urbanized areas, major cities,  and the top four museums represented in VertNet herpetology records in Texas (SI figure 2).
- [R/10_TexasImperviousMap.R](R/10_TexasImperviousMap.R): Produces map of percent impervious cover by county, major cities, and the top four museums represented in VertNet herpetology records in Texas (SI figure 3).

## Notes

The required R packages include:
- here
- sf
- tidyverse
- conflicted
- tigris
- terra
- exactextractr
- scales
- knitr
- kableExtra
- mgcv
- patchwork
- cowplot
- rphylopic
- grid
- gridExtra
- lme4
- emmeans
- viridis
- broom
- ggspatial
- paletteer
- shadowtext
- units
- ggrepel

Install them in one call:

```r
install.packages(c(
  "here", "sf", "tidyverse", "conflicted", "tigris", "terra",
  "exactextractr", "scales", "knitr", "kableExtra", "mgcv", "patchwork",
  "cowplot", "rphylopic", "gridExtra", "lme4", "emmeans", "viridis",
  "broom", "ggspatial", "paletteer", "shadowtext", "units", "ggrepel"
))
```

(`grid` ships with base R.)
