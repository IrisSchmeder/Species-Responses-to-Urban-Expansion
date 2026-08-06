# Code for Schmeder et al.
This repository contains the code for the manuscript "Divergent Species Responses to 100 Years of Urban Expansion Across Texas" submitted to Conservation Biology. It provides a pipeline to combine historical museum records with modern community science observations to model species occurrences in response to urbanization. The example demonstrated here uses eight herpetofauna species native to Texas.

## Abstract


## Data
The analysis uses downloaded data available from VertNet and iNaturalist.

## Scripts
### Data Prep & Info
- [01_DataPrep.R](01_DataPrep.R): Reads all csv files and builds global objects.
- [02_DatasetInfoTable.Rmd](02_DatasetInfoTable.Rmd): Makes table with dataset filtering information (table 1).

### Observational Bias
- [03_RecordAccumulation.R](03_RecordAccumulation.R): Makes record accumulation curves over time (figure 1).
- [04_TemporalBiasGAMs.R](04_TemporalBiasGAMs.R): Makes negative-binomial GAMs for record counts and unique observers over time (figure 2).
- [05_RatioDataset.R](05_RatioDataset.R): Calculates the proportion of records from each dataset within each county ranked by iNaturalist proportion (figure 3).

### Species Response
- [06_CountyDeltaAnalysis.R](06_CountyDeltaAnalysis.R): Makes County-level dataset delta analysis: paired t-tests and z-score maps (figure 4 and figure 5)
- [07_UrbanizationAnalysis.R](07_UrbanizationAnalysis.R): Calculates urbanization trends and impervious-surface GAMs (figure 6).

### Supplementary Figures
- [08_RecordsPerObserver.R](08_RecordsPerObserver.R): Calculates records per observer over time, by species and dataset (SI figure 1).
- [09_TexasEcoregionsMap.R](09_TexasEcoregionsMap.R): Produces map of level III ecoregions urbanized areas, major cities,  and the top four museums represented in VertNet herpetology records in Texas (SI figure 2).
- [10_TexasImperviousMap.R](10_TexasImperviousMap.R): Produces map of percent impervious cover by county, major cities, and the top four museums represented in VertNet herpetology records in Texas (SI figure 3).

## Notes
The required R packages include:
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
- conflicted
- grid
- gridExtra
- lme4
- emmeans
- viridis
- scales
- broom
- ggspatial
- paletteer
- shadowtext
- units
- ggrepel
