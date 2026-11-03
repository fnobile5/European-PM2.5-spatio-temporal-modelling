#### European-PM2.5-spatio-temporal-modelling

#### Author: Federica Nobile
#### Affiliation: Institute of Environmental Medicine, Karolinska Institutet, Stockholm, Sweden
####              Department of Epidemiology, Lazio Region Health Service/ASL Roma 1, Rome, Italy

library(dplyr)
library(lubridate)
library(splitstackshape)
# library(caret)
library(ranger)
library(data.table)
library(lme4)
library(reshape)
library(reshape2)
library(foreign) 
library(ggplot2)
library(data.table)
library(mgcv)
library(gdata)
library(car)
# library(ggmap)
library(broom)
library(splines)
library(DataCombine)
library(gtools)
library(party)
library(Metrics)
library(writexl)
library(corrplot)
library(caret)
library(xgboost)
options(scipen = 999) 

data <- readRDS("./stage1.PM.final.RDS")
dim(data)
names(data)

#### Select variables of interest

data1 <- data %>%
          dplyr::select(day, dow, mm,
                        Air.Quality.Station.EoI.Code, CountryName, Country, Air.Quality.Network,
                        Air.Quality.Network.Name, Air.Quality.Station.Nat.Code, Air.Quality.Station.Name,
                        Altitude, Air.Quality.Station.Area, Air.Quality.Station.Type, year,
                        Longitude, Latitude,
                        pm25, final.pm25, pm10, 
                        allRoads_50, allRoads_100, allRoads_200, allRoads_300, allRoads_400, allRoads_500,
                        majorRoads_50, majorRoads_100, majorRoads_200, majorRoads_300, majorRoads_400, majorRoads_500,
                        imd_100, imd_200, imd_300, imd_400, imd_500,
                        ind_100, ind_200, ind_300, ind_400, ind_500,
                        nat_100, nat_200, nat_300, nat_400, nat_500,
                        por_100, por_200, por_300, por_400, por_500,
                        res_100, res_200, res_300, res_400, res_500,
                        tbu_100, tbu_200, tbu_300, tbu_400, tbu_500,
                        ugr_100, ugr_200, ugr_300, ugr_400, ugr_500,
                        pop, zoneID, alt10_enh,
                        u10, v10, dewp_2m, temp_2m, mslp, tp, blh, leaf_high, leaf_low, tcc, tco,
                        pm25_cams, pm10_cams,
                        aod1240_cams, aod469_cams, aod550_cams, aod670_cams, aod865_cams,
                        NDVI, EVI,
                        final.aod047)

hist(data1$final.pm25)
ks.test(data1$final.pm25, 'pnorm') 

data1 <- data1 %>% dplyr::mutate(ln.pm25 = log(final.pm25 + 0.5),
                          mm      = as.factor(mm),
                          dow     = as.factor(dow))
summary(data1$ln.pm25)
hist(data1$ln.pm25)
ks.test(data1$ln.pm25, 'pnorm') 

#### SETTING HYPERPARAMETERS [...] ####

### Hyperparameters from 2019 (as explained in the Supplementary Materials)

formula <- as.formula(ln.pm25 ~ day + dow + mm +
                        Longitude + Latitude +
                        allRoads_50 + allRoads_100 + allRoads_200 + allRoads_300 + 
                        majorRoads_100 + majorRoads_300 + 
                        imd_200 + imd_300 + 
                        nat_300 + 
                        tbu_300 + 
                        pop + alt10_enh +
                        u10 + v10 + dewp_2m + temp_2m + mslp + tp + blh + leaf_high + leaf_low + tcc + tco +
                        pm25_cams +
                        NDVI + EVI +
                        final.aod047)

fitxgb <- trainControl(method = "cv",
  number = 5)

xgbGriddef <- expand.grid(nrounds = 1000, 
                          max_depth = 10,
                          colsample_bytree = 0.8,
                          eta = 0.1,
                          gamma = 0,
                          min_child_weight = 5,
                          subsample = 1)

set.seed(1)
print(Sys.time()) 
mod.xgb <- train(formula,
                 data = data1,
                 method = "xgbTree",
                 trControl = fitxgb,
                 verbose = TRUE,
                 tuneGrid = xgbGriddef)
print(Sys.time()) 

mod.xgb

importBoost <- varImp(mod.xgb)
importBoost$importance
