# European-PM2.5-spatio-temporal-modelling

This repository contains R scripts for developing spatio-temporal PM2.5 exposure models across Europe (using the year 2019 as an example). The workflow is organized modularly according to the consecutive stages of the data and modeling pipeline.

---

## Workflow Overview

### 1. Collection and Processing of Daily PM2.5 Concentrations
* `1a. PM time-series.R`: Compiles a daily time series of particulate matter observations from European Environment Agency (EEA) monitoring stations (processed from raw, publicly available EEA downloads).

### 2. Collection and Processing of Spatio-Temporal Predictors

#### Aerosol Optical Depth (AOD) – MODIS Terra+Aqua (MCD19A2 V061)
* `2a. AOD MODIS step 1.R`: Preprocessing of downloaded MODIS GeoTIFF tiles.
* `2a. AOD MODIS step 2.R`: Imputation and spatial gap-filling of missing AOD observations.

#### Aerosol Optical Depth (AOD) – Copernicus Atmosphere Monitoring Service (CAMS)
* `2b. CAMS AOD Process.R`: Processing of CAMS NetCDF (`.nc`) files used as covariates to predict missing MODIS AOD values.
* `2b. CAMS AOD Process - monitors.R`: Extraction and alignment of CAMS AOD data at station coordinates.

#### Particulate Matter – CAMS European Air Quality Reanalysis
* `2c. CAMS AIR POLL Process - monitors.R`: Extraction and processing of CAMS regional air quality reanalysis fields at monitoring station locations.

#### Meteorological Variables – ECMWF ERA5 Reanalysis
* `2d. ERA5 Process - monitors.R`: Extraction of daily meteorological variables (temperature, boundary layer height, wind components, etc.) at monitoring station coordinates from ERA5 NetCDF files.

#### Vegetation Indices – MODIS Aqua (MYD13A2 V061)
* `2e. NDVI Process - monitors.R`: Extraction and temporal harmonization of monthly NDVI/EVI values at monitoring station locations.

#### Training Dataset
* `2f. Creation training dataset.R`: Integration of ground monitors, spatial (from EXPANSE Project) and spatio-temporal covariates into a single modeling dataset, including PM2.5 imputation via co-located PM10 stations where applicable.

### 3. Learner Training
* `3a. ModelPM25_RF.R`: Fitting and hyperparameter tuning of the Random Forest model.
* `3b. ModelPM25_XGB.R`: Fitting and hyperparameter tuning of the Extreme Gradient Boosting (XGBoost) model.

### 4. Ensemble Modeling
* `4. ModelPM25_ENSEMBLE.R`: Implementation of the ensemble model combining learner-base model predictions.

### 5. Cross-Validation
* `5a. Spatial-blocked cv scheme.R`: Partitioning of data (by monitors) into spatially independent clusters/folds using spatial blocked cross-validation.
* `5b. Cross-validation ensemble by group.R`: Performance evaluation (RMSE, R^2, inter, slope) of base learners and the ensemble model.

### 6. Spatial Prediction & Mapping
* `6a. Prediction RF model_Europe.R`: Continental-scale gridded prediction using Random Forest.
* `6b. Prediction XGB model_Europe.R`: Continental-scale gridded prediction using XGBoost.
* `6c. Prediction ENSEMBLE model_Europe.R`: Final continuous spatial prediction surface generated via the ensemble model.

---

## Requirements

The scripts are written in **R** and primarily rely on the following packages:
* **Spatial & Raster Data:** `terra`, `sf`, `ncdf4`, `exactextractr`
* **Data Wrangling:** `data.table`, `tidyverse` (`dplyr`, `lubridate`, `tidyr`)
* **Machine Learning & CV:** `ranger`, `xgboost`, `caret`
---

## Citation & License

If you use these scripts or adapt this workflow for research, please cite the corresponding publication:
> XXXXXX WILL BE ADDED

Distributed under the MIT License. See `LICENSE` for more information.
