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
library(gdata)
options(scipen = 999) 

data <- readRDS("./stage1.PM.final 18112024.RDS")
dim(data) 
names(data)

#### RUN THE DEFINITIVE MODEL FROM SCRIPT 2. ModelPM25_ExtremeGradientBoosting.R ####

data1 <- data %>%
  dplyr::select(day, dow, mm,
                Cntrycd, ArQltyN, ArQltyS, ArQSEIC, ArQlSNC, Altitud, ArQltST, ArQltSA, Longitude, Latitude,
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
                pm25_dm.2, pm10_dm.2,
                aod1240_cams, aod469_cams, aod550_cams, aod670_cams, aod865_cams,
                NDVI, EVI,
                final.aod047, final.aod055) %>%
  dplyr::rename(pm25_cams = pm25_dm.2,
                pm10_cams = pm10_dm.2)

rm(data)

data1 <- data1 %>% dplyr::mutate(ln.pm25 = log(final.pm25 + 0.5),
                                 mm      = as.factor(mm),
                                 dow     = as.factor(dow))

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

fitxgb <- trainControl(## No default cross validation
  method = "none")

xgbGriddef <- expand.grid(nrounds = 1000, 
                          max_depth = 10,
                          colsample_bytree = 0.8, 
                          eta = 0.1,
                          gamma = 0,
                          min_child_weight = 5,
                          subsample = 1)

set.seed(1)
mod2 <- train(formula,
              data = data1,
              method = "xgbTree",
              trControl = fitxgb,
              verbose = FALSE,
              tuneGrid = xgbGriddef)

varImp(mod2)$importance

keep(mod2, sure = T)
gc()


#### PREDICTION 5 KM ####

print("!!!!!!!!!!!!!!!!!!! PREDICTION 5 KM !!!!!!!!!!!!!!!!!!!")

##### Call stage3 dataset
stage3 <- readRDS("./stage3.2019.PM.all.grid5km.RDS")
names(stage3)
dim(stage3)
gc()

head(stage3)

# Define factor variables
stage3[,mm:=as.factor(mm)]
stage3[,dow:=as.factor(dow)]
stage3 <- as.data.frame(stage3)
gc()

##### Predict STAGE 1 model on Stage 3 dataset, by month

months <- c ("01", "02", "03", "04", "05", "06",
             "07", "08", "09", "10", "11", "12")

pred <- list()

for(i in 1:12){
  
  month  <- months[i]
  print(paste0(month, " ", Sys.time()))
  
  data             <- stage3 %>% dplyr::filter(mm == i)
  data$pred.rf.st3 <- predict(mod2,  data)
  data             <- data %>% dplyr::select(idcell, mm, day, Latitude, Longitude,
                                             pred.rf.st3)
  print(head(data))
  
  pred[[i]]        <- data
  
  rm(data)
  gc()
  
  print(paste0(month, " ", Sys.time(), " predicted!"))
  
}

final <- bind_rows(pred)
head(final)
dim(final)

setnames(final, "pred.rf.st3", "ln.pm25.pred.xgb")
final$pm25.pred.xgb <- exp(final$ln.pm25.pred.xgb)

head(final)
dim(final)

summary(final$pm25.pred.xgb)

tapply(final$pm25.pred.xgb, final$mm,
       function(x) format(summary(x)))

final %>% 
  dplyr::group_by(mm) %>% 
  dplyr::summarize(q01 = quantile(pm25.pred.xgb, 0.01),
                   q05 = quantile(pm25.pred.xgb, 0.05),
                   q25 = quantile(pm25.pred.xgb, 0.25),
                   q50 = quantile(pm25.pred.xgb, 0.50),
                   q75 = quantile(pm25.pred.xgb, 0.75),
                   q95 = quantile(pm25.pred.xgb, 0.95),
                   q99 = quantile(pm25.pred.xgb, 0.99))

saveRDS(final, "./stage3.2019...PM25.grid5km.XGB.RDS")

keep(mod2, sure = T)
gc()


#### PREDICTION 1 KM ####

print("!!!!!!!!!!!!!!!!!!! PREDICTION 1 KM !!!!!!!!!!!!!!!!!!!")

#### 1_1 ####

##### Call stage3 dataset
stage3 <- readRDS("./stage3.2019.PM.all.grid1km_1_1.RDS")
names(stage3)
dim(stage3)
gc()

head(stage3)

# Define factor variables
stage3[,mm:=as.factor(mm)]
stage3[,dow:=as.factor(dow)]
stage3 <- as.data.frame(stage3)
gc()

##### Predict STAGE 1 model on Stage 3 dataset, by month

months <- c ("01", "02", "03", "04", "05", "06",
             "07", "08", "09", "10", "11", "12")

pred <- list()

for(i in 1:12){
  
  month  <- months[i]
  print(paste0(month, " ", Sys.time()))
  
  data             <- stage3 %>% dplyr::filter(mm == i)
  data$pred.rf.st3 <- predict(mod2,  data)
  data             <- data %>% dplyr::select(idcell, mm, day, Latitude, Longitude,
                                             pred.rf.st3)
  print(head(data))
  
  pred[[i]]        <- data
  
  rm(data)
  gc()
  
  print(paste0(month, " ", Sys.time(), " predicted!"))
  
}


final <- bind_rows(pred)
head(final)
dim(final)

setnames(final, "pred.rf.st3", "ln.pm25.pred.xgb")
final$pm25.pred.xgb <- exp(final$ln.pm25.pred.xgb)

head(final)
dim(final)

summary(final$pm25.pred.xgb)

tapply(final$pm25.pred.xgb, final$mm,
       function(x) format(summary(x)))

final %>% 
  dplyr::group_by(mm) %>% 
  dplyr::summarize(q01 = quantile(pm25.pred.xgb, 0.01),
                   q05 = quantile(pm25.pred.xgb, 0.05),
                   q25 = quantile(pm25.pred.xgb, 0.25),
                   q50 = quantile(pm25.pred.xgb, 0.50),
                   q75 = quantile(pm25.pred.xgb, 0.75),
                   q95 = quantile(pm25.pred.xgb, 0.95),
                   q99 = quantile(pm25.pred.xgb, 0.99))

saveRDS(final, "./stage3.2019...PM25.grid1km.XGB_1_1.RDS")

keep(mod2, sure = T)
gc()

#### 1_2 ####

##### Call stage3 dataset
stage3 <- readRDS("./stage3.2019.PM.all.grid1km_1_2.RDS")
names(stage3)
dim(stage3)
gc()

head(stage3)

# Define factor variables
stage3[,mm:=as.factor(mm)]
stage3[,dow:=as.factor(dow)]
stage3 <- as.data.frame(stage3)
gc()

##### Predict STAGE 1 model on Stage 3 dataset, by month

months <- c ("01", "02", "03", "04", "05", "06",
             "07", "08", "09", "10", "11", "12")

pred <- list()

for(i in 1:12){
  
  month  <- months[i]
  print(paste0(month, " ", Sys.time()))
  
  data             <- stage3 %>% dplyr::filter(mm == i)
  data$pred.rf.st3 <- predict(mod2,  data)
  data             <- data %>% dplyr::select(idcell, mm, day, Latitude, Longitude,
                                             pred.rf.st3)
  print(head(data))
  
  pred[[i]]        <- data
  
  rm(data)
  gc()
  
  print(paste0(month, " ", Sys.time(), " predicted!"))
  
}


final <- bind_rows(pred)
head(final)
dim(final)

setnames(final, "pred.rf.st3", "ln.pm25.pred.xgb")
final$pm25.pred.xgb <- exp(final$ln.pm25.pred.xgb)

head(final)
dim(final)

summary(final$pm25.pred.xgb)

tapply(final$pm25.pred.xgb, final$mm,
       function(x) format(summary(x)))

final %>% 
  dplyr::group_by(mm) %>% 
  dplyr::summarize(q01 = quantile(pm25.pred.xgb, 0.01),
                   q05 = quantile(pm25.pred.xgb, 0.05),
                   q25 = quantile(pm25.pred.xgb, 0.25),
                   q50 = quantile(pm25.pred.xgb, 0.50),
                   q75 = quantile(pm25.pred.xgb, 0.75),
                   q95 = quantile(pm25.pred.xgb, 0.95),
                   q99 = quantile(pm25.pred.xgb, 0.99))

saveRDS(final, "./stage3.2019...PM25.grid1km.XGB_1_2.RDS")

keep(mod2, sure = T)
gc()

#### 2_1 ####

##### Call stage3 dataset
stage3 <- readRDS("./stage3.2019.PM.all.grid1km_2_1.RDS")
names(stage3)
dim(stage3)
gc()

head(stage3)

# Define factor variables
stage3[,mm:=as.factor(mm)]
stage3[,dow:=as.factor(dow)]
stage3 <- as.data.frame(stage3)
gc()

##### Predict STAGE 1 model on Stage 3 dataset, by month

months <- c ("01", "02", "03", "04", "05", "06",
             "07", "08", "09", "10", "11", "12")

pred <- list()

for(i in 1:12){
  
  month  <- months[i]
  print(paste0(month, " ", Sys.time()))
  
  data             <- stage3 %>% dplyr::filter(mm == i)
  data$pred.rf.st3 <- predict(mod2,  data)
  data             <- data %>% dplyr::select(idcell, mm, day, Latitude, Longitude,
                                             pred.rf.st3)
  print(head(data))
  
  pred[[i]]        <- data
  
  rm(data)
  gc()
  
  print(paste0(month, " ", Sys.time(), " predicted!"))
  
}


final <- bind_rows(pred)
head(final)
dim(final)

setnames(final, "pred.rf.st3", "ln.pm25.pred.xgb")
final$pm25.pred.xgb <- exp(final$ln.pm25.pred.xgb)

head(final)
dim(final)

summary(final$pm25.pred.xgb)

tapply(final$pm25.pred.xgb, final$mm,
       function(x) format(summary(x)))

final %>% 
  dplyr::group_by(mm) %>% 
  dplyr::summarize(q01 = quantile(pm25.pred.xgb, 0.01),
                   q05 = quantile(pm25.pred.xgb, 0.05),
                   q25 = quantile(pm25.pred.xgb, 0.25),
                   q50 = quantile(pm25.pred.xgb, 0.50),
                   q75 = quantile(pm25.pred.xgb, 0.75),
                   q95 = quantile(pm25.pred.xgb, 0.95),
                   q99 = quantile(pm25.pred.xgb, 0.99))

saveRDS(final, "./stage3.2019...PM25.grid1km.XGB_2_1.RDS")

keep(mod2, sure = T)
gc()

#### 2_2 ####

##### Call stage3 dataset
stage3 <- readRDS("./stage3.2019.PM.all.grid1km_2_2.RDS")
names(stage3)
dim(stage3)
gc()

head(stage3)

# Define factor variables
stage3[,mm:=as.factor(mm)]
stage3[,dow:=as.factor(dow)]
stage3 <- as.data.frame(stage3)
gc()

##### Predict STAGE 1 model on Stage 3 dataset, by month

months <- c ("01", "02", "03", "04", "05", "06",
             "07", "08", "09", "10", "11", "12")

pred <- list()

for(i in 1:12){
  
  month  <- months[i]
  print(paste0(month, " ", Sys.time()))
  
  data             <- stage3 %>% dplyr::filter(mm == i)
  data$pred.rf.st3 <- predict(mod2,  data)
  data             <- data %>% dplyr::select(idcell, mm, day, Latitude, Longitude,
                                             pred.rf.st3)
  print(head(data))
  
  pred[[i]]        <- data
  
  rm(data)
  gc()
  
  print(paste0(month, " ", Sys.time(), " predicted!"))
  
}


final <- bind_rows(pred)
head(final)
dim(final)

setnames(final, "pred.rf.st3", "ln.pm25.pred.xgb")
final$pm25.pred.xgb <- exp(final$ln.pm25.pred.xgb)

head(final)
dim(final)

summary(final$pm25.pred.xgb)

tapply(final$pm25.pred.xgb, final$mm,
       function(x) format(summary(x)))

final %>% 
  dplyr::group_by(mm) %>% 
  dplyr::summarize(q01 = quantile(pm25.pred.xgb, 0.01),
                   q05 = quantile(pm25.pred.xgb, 0.05),
                   q25 = quantile(pm25.pred.xgb, 0.25),
                   q50 = quantile(pm25.pred.xgb, 0.50),
                   q75 = quantile(pm25.pred.xgb, 0.75),
                   q95 = quantile(pm25.pred.xgb, 0.95),
                   q99 = quantile(pm25.pred.xgb, 0.99))

saveRDS(final, "./stage3.2019...PM25.grid1km.XGB_2_2.RDS")

keep(mod2, sure = T)
gc()

#### 3_1 ####

##### Call stage3 dataset
stage3 <- readRDS("./stage3.2019.PM.all.grid1km_3_1.RDS")
names(stage3)
dim(stage3)
gc()

head(stage3)

# Define factor variables
stage3[,mm:=as.factor(mm)]
stage3[,dow:=as.factor(dow)]
stage3 <- as.data.frame(stage3)
gc()

##### Predict STAGE 1 model on Stage 3 dataset, by month

months <- c ("01", "02", "03", "04", "05", "06",
             "07", "08", "09", "10", "11", "12")

pred <- list()

for(i in 1:12){
  
  month  <- months[i]
  print(paste0(month, " ", Sys.time()))
  
  data             <- stage3 %>% dplyr::filter(mm == i)
  data$pred.rf.st3 <- predict(mod2,  data)
  data             <- data %>% dplyr::select(idcell, mm, day, Latitude, Longitude,
                                             pred.rf.st3)
  print(head(data))
  
  pred[[i]]        <- data
  
  rm(data)
  gc()
  
  print(paste0(month, " ", Sys.time(), " predicted!"))
  
}


final <- bind_rows(pred)
head(final)
dim(final)

setnames(final, "pred.rf.st3", "ln.pm25.pred.xgb")
final$pm25.pred.xgb <- exp(final$ln.pm25.pred.xgb)

head(final)
dim(final)

summary(final$pm25.pred.xgb)

tapply(final$pm25.pred.xgb, final$mm,
       function(x) format(summary(x)))

final %>% 
  dplyr::group_by(mm) %>% 
  dplyr::summarize(q01 = quantile(pm25.pred.xgb, 0.01),
                   q05 = quantile(pm25.pred.xgb, 0.05),
                   q25 = quantile(pm25.pred.xgb, 0.25),
                   q50 = quantile(pm25.pred.xgb, 0.50),
                   q75 = quantile(pm25.pred.xgb, 0.75),
                   q95 = quantile(pm25.pred.xgb, 0.95),
                   q99 = quantile(pm25.pred.xgb, 0.99))

saveRDS(final, "./stage3.2019...PM25.grid1km.XGB_3_1.RDS")

keep(mod2, sure = T)
gc()

#### 3_2 ####

##### Call stage3 dataset
stage3 <- readRDS("./stage3.2019.PM.all.grid1km_3_2.RDS")
names(stage3)
dim(stage3)
gc()

head(stage3)

# Define factor variables
stage3[,mm:=as.factor(mm)]
stage3[,dow:=as.factor(dow)]
stage3 <- as.data.frame(stage3)
gc()

##### Predict STAGE 1 model on Stage 3 dataset, by month

months <- c ("01", "02", "03", "04", "05", "06",
             "07", "08", "09", "10", "11", "12")

pred <- list()

for(i in 1:12){
  
  month  <- months[i]
  print(paste0(month, " ", Sys.time()))
  
  data             <- stage3 %>% dplyr::filter(mm == i)
  data$pred.rf.st3 <- predict(mod2,  data)
  data             <- data %>% dplyr::select(idcell, mm, day, Latitude, Longitude,
                                             pred.rf.st3)
  print(head(data))
  
  pred[[i]]        <- data
  
  rm(data)
  gc()
  
  print(paste0(month, " ", Sys.time(), " predicted!"))
  
}


final <- bind_rows(pred)
head(final)
dim(final)

setnames(final, "pred.rf.st3", "ln.pm25.pred.xgb")
final$pm25.pred.xgb <- exp(final$ln.pm25.pred.xgb)

head(final)
dim(final)

summary(final$pm25.pred.xgb)

tapply(final$pm25.pred.xgb, final$mm,
       function(x) format(summary(x)))

final %>% 
  dplyr::group_by(mm) %>% 
  dplyr::summarize(q01 = quantile(pm25.pred.xgb, 0.01),
                   q05 = quantile(pm25.pred.xgb, 0.05),
                   q25 = quantile(pm25.pred.xgb, 0.25),
                   q50 = quantile(pm25.pred.xgb, 0.50),
                   q75 = quantile(pm25.pred.xgb, 0.75),
                   q95 = quantile(pm25.pred.xgb, 0.95),
                   q99 = quantile(pm25.pred.xgb, 0.99))

saveRDS(final, "./stage3.2019...PM25.grid1km.XGB_3_2.RDS")

keep(mod2, sure = T)
gc()


#### PREDICTION 100M ####

print("!!!!!!!!!!!!!!!!!!! PREDICTION 100 M !!!!!!!!!!!!!!!!!!!")

#### 1_1 ####

##### Call stage3 dataset
stage3 <- readRDS("./stage3.2019.PM.all.grid100m_1_1.RDS")
names(stage3)
dim(stage3)
gc()

head(stage3)

# Define factor variables
stage3[,mm:=as.factor(mm)]
stage3[,dow:=as.factor(dow)]
stage3 <- as.data.frame(stage3)
gc()

##### Predict STAGE 1 model on Stage 3 dataset, by month

months <- c ("01", "02", "03", "04", "05", "06",
             "07", "08", "09", "10", "11", "12")

pred <- list()

for(i in 1:12){
  
  month  <- months[i]
  print(paste0(month, " ", Sys.time()))
  
  data             <- stage3 %>% dplyr::filter(mm == i)
  data$pred.rf.st3 <- predict(mod2,  data)
  data             <- data %>% dplyr::select(idcell, mm, day, Latitude, Longitude,
                                             pred.rf.st3)
  print(head(data))
  
  pred[[i]]        <- data
  
  rm(data)
  gc()
  
  print(paste0(month, " ", Sys.time(), " predicted!"))
  
}


final <- bind_rows(pred)
head(final)
dim(final)

setnames(final, "pred.rf.st3", "ln.pm25.pred.xgb")
final$pm25.pred.xgb <- exp(final$ln.pm25.pred.xgb)

head(final)
dim(final)

summary(final$pm25.pred.xgb)

tapply(final$pm25.pred.xgb, final$mm,
       function(x) format(summary(x)))

final %>% 
  dplyr::group_by(mm) %>% 
  dplyr::summarize(q01 = quantile(pm25.pred.xgb, 0.01),
                   q05 = quantile(pm25.pred.xgb, 0.05),
                   q25 = quantile(pm25.pred.xgb, 0.25),
                   q50 = quantile(pm25.pred.xgb, 0.50),
                   q75 = quantile(pm25.pred.xgb, 0.75),
                   q95 = quantile(pm25.pred.xgb, 0.95),
                   q99 = quantile(pm25.pred.xgb, 0.99))

saveRDS(final, "./stage3.2019...PM25.grid100m.XGB_1_1.RDS")

keep(mod2, sure = T)
gc()

#### 1_2 ####

##### Call stage3 dataset
stage3 <- readRDS("./stage3.2019.PM.all.grid100m_1_2.RDS")
names(stage3)
dim(stage3)
gc()

head(stage3)

# Define factor variables
stage3[,mm:=as.factor(mm)]
stage3[,dow:=as.factor(dow)]
stage3 <- as.data.frame(stage3)
gc()

##### Predict STAGE 1 model on Stage 3 dataset, by month

months <- c ("01", "02", "03", "04", "05", "06",
             "07", "08", "09", "10", "11", "12")

pred <- list()

for(i in 1:12){
  
  month  <- months[i]
  print(paste0(month, " ", Sys.time()))
  
  data             <- stage3 %>% dplyr::filter(mm == i)
  data$pred.rf.st3 <- predict(mod2,  data)
  data             <- data %>% dplyr::select(idcell, mm, day, Latitude, Longitude,
                                             pred.rf.st3)
  print(head(data))
  
  pred[[i]]        <- data
  
  rm(data)
  gc()
  
  print(paste0(month, " ", Sys.time(), " predicted!"))
  
}


final <- bind_rows(pred)
head(final)
dim(final)

setnames(final, "pred.rf.st3", "ln.pm25.pred.xgb")
final$pm25.pred.xgb <- exp(final$ln.pm25.pred.xgb)

head(final)
dim(final)

summary(final$pm25.pred.xgb)

tapply(final$pm25.pred.xgb, final$mm,
       function(x) format(summary(x)))

final %>% 
  dplyr::group_by(mm) %>% 
  dplyr::summarize(q01 = quantile(pm25.pred.xgb, 0.01),
                   q05 = quantile(pm25.pred.xgb, 0.05),
                   q25 = quantile(pm25.pred.xgb, 0.25),
                   q50 = quantile(pm25.pred.xgb, 0.50),
                   q75 = quantile(pm25.pred.xgb, 0.75),
                   q95 = quantile(pm25.pred.xgb, 0.95),
                   q99 = quantile(pm25.pred.xgb, 0.99))

saveRDS(final, "./stage3.2019...PM25.grid100m.XGB_1_2.RDS")

keep(mod2, sure = T)
gc()

#### 2_1 ####

##### Call stage3 dataset
stage3 <- readRDS("./stage3.2019.PM.all.grid100m_2_1.RDS")
names(stage3)
dim(stage3)
gc()

head(stage3)

# Define factor variables
stage3[,mm:=as.factor(mm)]
stage3[,dow:=as.factor(dow)]
stage3 <- as.data.frame(stage3)
gc()

##### Predict STAGE 1 model on Stage 3 dataset, by month

months <- c ("01", "02", "03", "04", "05", "06",
             "07", "08", "09", "10", "11", "12")

pred <- list()

for(i in 1:12){
  
  month  <- months[i]
  print(paste0(month, " ", Sys.time()))
  
  data             <- stage3 %>% dplyr::filter(mm == i)
  data$pred.rf.st3 <- predict(mod2,  data)
  data             <- data %>% dplyr::select(idcell, mm, day, Latitude, Longitude,
                                             pred.rf.st3)
  print(head(data))
  
  pred[[i]]        <- data
  
  rm(data)
  gc()
  
  print(paste0(month, " ", Sys.time(), " predicted!"))
  
}


final <- bind_rows(pred)
head(final)
dim(final)

setnames(final, "pred.rf.st3", "ln.pm25.pred.xgb")
final$pm25.pred.xgb <- exp(final$ln.pm25.pred.xgb)

head(final)
dim(final)

summary(final$pm25.pred.xgb)

tapply(final$pm25.pred.xgb, final$mm,
       function(x) format(summary(x)))

final %>% 
  dplyr::group_by(mm) %>% 
  dplyr::summarize(q01 = quantile(pm25.pred.xgb, 0.01),
                   q05 = quantile(pm25.pred.xgb, 0.05),
                   q25 = quantile(pm25.pred.xgb, 0.25),
                   q50 = quantile(pm25.pred.xgb, 0.50),
                   q75 = quantile(pm25.pred.xgb, 0.75),
                   q95 = quantile(pm25.pred.xgb, 0.95),
                   q99 = quantile(pm25.pred.xgb, 0.99))

saveRDS(final, "./stage3.2019...PM25.grid100m.XGB_2_1.RDS")

keep(mod2, sure = T)
gc()

#### 2_2 ####

##### Call stage3 dataset
stage3 <- readRDS("./stage3.2019.PM.all.grid100m_2_2.RDS")
names(stage3)
dim(stage3)
gc()

head(stage3)

# Define factor variables
stage3[,mm:=as.factor(mm)]
stage3[,dow:=as.factor(dow)]
stage3 <- as.data.frame(stage3)
gc()

##### Predict STAGE 1 model on Stage 3 dataset, by month

months <- c ("01", "02", "03", "04", "05", "06",
             "07", "08", "09", "10", "11", "12")

pred <- list()

for(i in 1:12){
  
  month  <- months[i]
  print(paste0(month, " ", Sys.time()))
  
  data             <- stage3 %>% dplyr::filter(mm == i)
  data$pred.rf.st3 <- predict(mod2,  data)
  data             <- data %>% dplyr::select(idcell, mm, day, Latitude, Longitude,
                                             pred.rf.st3)
  print(head(data))
  
  pred[[i]]        <- data
  
  rm(data)
  gc()
  
  print(paste0(month, " ", Sys.time(), " predicted!"))
  
}


final <- bind_rows(pred)
head(final)
dim(final)

setnames(final, "pred.rf.st3", "ln.pm25.pred.xgb")
final$pm25.pred.xgb <- exp(final$ln.pm25.pred.xgb)

head(final)
dim(final)

summary(final$pm25.pred.xgb)

tapply(final$pm25.pred.xgb, final$mm,
       function(x) format(summary(x)))

final %>% 
  dplyr::group_by(mm) %>% 
  dplyr::summarize(q01 = quantile(pm25.pred.xgb, 0.01),
                   q05 = quantile(pm25.pred.xgb, 0.05),
                   q25 = quantile(pm25.pred.xgb, 0.25),
                   q50 = quantile(pm25.pred.xgb, 0.50),
                   q75 = quantile(pm25.pred.xgb, 0.75),
                   q95 = quantile(pm25.pred.xgb, 0.95),
                   q99 = quantile(pm25.pred.xgb, 0.99))

saveRDS(final, "./stage3.2019...PM25.grid100m.XGB_2_2.RDS")

keep(mod2, sure = T)
gc()

#### 3_1 ####

##### Call stage3 dataset
stage3 <- readRDS("./stage3.2019.PM.all.grid100m_3_1.RDS")
names(stage3)
dim(stage3)
gc()

head(stage3)

# Define factor variables
stage3[,mm:=as.factor(mm)]
stage3[,dow:=as.factor(dow)]
stage3 <- as.data.frame(stage3)
gc()

##### Predict STAGE 1 model on Stage 3 dataset, by month

months <- c ("01", "02", "03", "04", "05", "06",
             "07", "08", "09", "10", "11", "12")

pred <- list()

for(i in 1:12){
  
  month  <- months[i]
  print(paste0(month, " ", Sys.time()))
  
  data             <- stage3 %>% dplyr::filter(mm == i)
  data$pred.rf.st3 <- predict(mod2,  data)
  data             <- data %>% dplyr::select(idcell, mm, day, Latitude, Longitude,
                                             pred.rf.st3)
  print(head(data))
  
  pred[[i]]        <- data
  
  rm(data)
  gc()
  
  print(paste0(month, " ", Sys.time(), " predicted!"))
  
}


final <- bind_rows(pred)
head(final)
dim(final)

setnames(final, "pred.rf.st3", "ln.pm25.pred.xgb")
final$pm25.pred.xgb <- exp(final$ln.pm25.pred.xgb)

head(final)
dim(final)

summary(final$pm25.pred.xgb)

tapply(final$pm25.pred.xgb, final$mm,
       function(x) format(summary(x)))

final %>% 
  dplyr::group_by(mm) %>% 
  dplyr::summarize(q01 = quantile(pm25.pred.xgb, 0.01),
                   q05 = quantile(pm25.pred.xgb, 0.05),
                   q25 = quantile(pm25.pred.xgb, 0.25),
                   q50 = quantile(pm25.pred.xgb, 0.50),
                   q75 = quantile(pm25.pred.xgb, 0.75),
                   q95 = quantile(pm25.pred.xgb, 0.95),
                   q99 = quantile(pm25.pred.xgb, 0.99))

saveRDS(final, "./stage3.2019...PM25.grid100m.XGB_3_1.RDS")

keep(mod2, sure = T)
gc()

#### 3_2 ####

##### Call stage3 dataset
stage3 <- readRDS("./stage3.2019.PM.all.grid100m_3_2.RDS")
names(stage3)
dim(stage3)
gc()

head(stage3)

# Define factor variables
stage3[,mm:=as.factor(mm)]
stage3[,dow:=as.factor(dow)]
stage3 <- as.data.frame(stage3)
gc()

##### Predict STAGE 1 model on Stage 3 dataset, by month

months <- c ("01", "02", "03", "04", "05", "06",
             "07", "08", "09", "10", "11", "12")

pred <- list()

for(i in 1:12){
  
  month  <- months[i]
  print(paste0(month, " ", Sys.time()))
  
  data             <- stage3 %>% dplyr::filter(mm == i)
  data$pred.rf.st3 <- predict(mod2,  data)
  data             <- data %>% dplyr::select(idcell, mm, day, Latitude, Longitude,
                                             pred.rf.st3)
  print(head(data))
  
  pred[[i]]        <- data
  
  rm(data)
  gc()
  
  print(paste0(month, " ", Sys.time(), " predicted!"))
  
}


final <- bind_rows(pred)
head(final)
dim(final)

setnames(final, "pred.rf.st3", "ln.pm25.pred.xgb")
final$pm25.pred.xgb <- exp(final$ln.pm25.pred.xgb)

head(final)
dim(final)

summary(final$pm25.pred.xgb)

tapply(final$pm25.pred.xgb, final$mm,
       function(x) format(summary(x)))

final %>% 
  dplyr::group_by(mm) %>% 
  dplyr::summarize(q01 = quantile(pm25.pred.xgb, 0.01),
                   q05 = quantile(pm25.pred.xgb, 0.05),
                   q25 = quantile(pm25.pred.xgb, 0.25),
                   q50 = quantile(pm25.pred.xgb, 0.50),
                   q75 = quantile(pm25.pred.xgb, 0.75),
                   q95 = quantile(pm25.pred.xgb, 0.95),
                   q99 = quantile(pm25.pred.xgb, 0.99))

saveRDS(final, "./stage3.2019...PM25.grid100m.XGB_3_2.RDS")

keep(mod2, sure = T)
gc()

#### 4_1 ####

##### Call stage3 dataset
stage3 <- readRDS("./stage3.2019.PM.all.grid100m_4_1.RDS")
names(stage3)
dim(stage3)
gc()

head(stage3)

# Define factor variables
stage3[,mm:=as.factor(mm)]
stage3[,dow:=as.factor(dow)]
stage3 <- as.data.frame(stage3)
gc()

##### Predict STAGE 1 model on Stage 3 dataset, by month

months <- c ("01", "02", "03", "04", "05", "06",
             "07", "08", "09", "10", "11", "12")

pred <- list()

for(i in 1:12){
  
  month  <- months[i]
  print(paste0(month, " ", Sys.time()))
  
  data             <- stage3 %>% dplyr::filter(mm == i)
  data$pred.rf.st3 <- predict(mod2,  data)
  data             <- data %>% dplyr::select(idcell, mm, day, Latitude, Longitude,
                                             pred.rf.st3)
  print(head(data))
  
  pred[[i]]        <- data
  
  rm(data)
  gc()
  
  print(paste0(month, " ", Sys.time(), " predicted!"))
  
}


final <- bind_rows(pred)
head(final)
dim(final)

setnames(final, "pred.rf.st3", "ln.pm25.pred.xgb")
final$pm25.pred.xgb <- exp(final$ln.pm25.pred.xgb)

head(final)
dim(final)

summary(final$pm25.pred.xgb)

tapply(final$pm25.pred.xgb, final$mm,
       function(x) format(summary(x)))

final %>% 
  dplyr::group_by(mm) %>% 
  dplyr::summarize(q01 = quantile(pm25.pred.xgb, 0.01),
                   q05 = quantile(pm25.pred.xgb, 0.05),
                   q25 = quantile(pm25.pred.xgb, 0.25),
                   q50 = quantile(pm25.pred.xgb, 0.50),
                   q75 = quantile(pm25.pred.xgb, 0.75),
                   q95 = quantile(pm25.pred.xgb, 0.95),
                   q99 = quantile(pm25.pred.xgb, 0.99))

saveRDS(final, "./stage3.2019...PM25.grid100m.XGB_4_1.RDS")

keep(mod2, sure = T)
gc()

#### 4_2 ####

##### Call stage3 dataset
stage3 <- readRDS("./stage3.2019.PM.all.grid100m_4_2.RDS")
names(stage3)
dim(stage3)
gc()

head(stage3)

# Define factor variables
stage3[,mm:=as.factor(mm)]
stage3[,dow:=as.factor(dow)]
stage3 <- as.data.frame(stage3)
gc()

##### Predict STAGE 1 model on Stage 3 dataset, by month

months <- c ("01", "02", "03", "04", "05", "06",
             "07", "08", "09", "10", "11", "12")

pred <- list()

for(i in 1:12){
  
  month  <- months[i]
  print(paste0(month, " ", Sys.time()))
  
  data             <- stage3 %>% dplyr::filter(mm == i)
  data$pred.rf.st3 <- predict(mod2,  data)
  data             <- data %>% dplyr::select(idcell, mm, day, Latitude, Longitude,
                                             pred.rf.st3)
  print(head(data))
  
  pred[[i]]        <- data
  
  rm(data)
  gc()
  
  print(paste0(month, " ", Sys.time(), " predicted!"))
  
}


final <- bind_rows(pred)
head(final)
dim(final)

setnames(final, "pred.rf.st3", "ln.pm25.pred.xgb")
final$pm25.pred.xgb <- exp(final$ln.pm25.pred.xgb)

head(final)
dim(final)

summary(final$pm25.pred.xgb)

tapply(final$pm25.pred.xgb, final$mm,
       function(x) format(summary(x)))

final %>% 
  dplyr::group_by(mm) %>% 
  dplyr::summarize(q01 = quantile(pm25.pred.xgb, 0.01),
                   q05 = quantile(pm25.pred.xgb, 0.05),
                   q25 = quantile(pm25.pred.xgb, 0.25),
                   q50 = quantile(pm25.pred.xgb, 0.50),
                   q75 = quantile(pm25.pred.xgb, 0.75),
                   q95 = quantile(pm25.pred.xgb, 0.95),
                   q99 = quantile(pm25.pred.xgb, 0.99))

saveRDS(final, "./stage3.2019...PM25.grid100m.XGB_4_2.RDS")

keep(mod2, sure = T)
gc()

#### 5_1 ####

##### Call stage3 dataset
stage3 <- readRDS("./stage3.2019.PM.all.grid100m_5_1.RDS")
names(stage3)
dim(stage3)
gc()

head(stage3)

# Define factor variables
stage3[,mm:=as.factor(mm)]
stage3[,dow:=as.factor(dow)]
stage3 <- as.data.frame(stage3)
gc()

##### Predict STAGE 1 model on Stage 3 dataset, by month

months <- c ("01", "02", "03", "04", "05", "06",
             "07", "08", "09", "10", "11", "12")

pred <- list()

for(i in 1:12){
  
  month  <- months[i]
  print(paste0(month, " ", Sys.time()))
  
  data             <- stage3 %>% dplyr::filter(mm == i)
  data$pred.rf.st3 <- predict(mod2,  data)
  data             <- data %>% dplyr::select(idcell, mm, day, Latitude, Longitude,
                                             pred.rf.st3)
  print(head(data))
  
  pred[[i]]        <- data
  
  rm(data)
  gc()
  
  print(paste0(month, " ", Sys.time(), " predicted!"))
  
}


final <- bind_rows(pred)
head(final)
dim(final)

setnames(final, "pred.rf.st3", "ln.pm25.pred.xgb")
final$pm25.pred.xgb <- exp(final$ln.pm25.pred.xgb)

head(final)
dim(final)

summary(final$pm25.pred.xgb)

tapply(final$pm25.pred.xgb, final$mm,
       function(x) format(summary(x)))

final %>% 
  dplyr::group_by(mm) %>% 
  dplyr::summarize(q01 = quantile(pm25.pred.xgb, 0.01),
                   q05 = quantile(pm25.pred.xgb, 0.05),
                   q25 = quantile(pm25.pred.xgb, 0.25),
                   q50 = quantile(pm25.pred.xgb, 0.50),
                   q75 = quantile(pm25.pred.xgb, 0.75),
                   q95 = quantile(pm25.pred.xgb, 0.95),
                   q99 = quantile(pm25.pred.xgb, 0.99))

saveRDS(final, "./stage3.2019...PM25.grid100m.XGB_5_1.RDS")

keep(mod2, sure = T)
gc()

#### 5_2 ####

##### Call stage3 dataset
stage3 <- readRDS("./stage3.2019.PM.all.grid100m_5_2.RDS")
names(stage3)
dim(stage3)
gc()

head(stage3)

# Define factor variables
stage3[,mm:=as.factor(mm)]
stage3[,dow:=as.factor(dow)]
stage3 <- as.data.frame(stage3)
gc()

##### Predict STAGE 1 model on Stage 3 dataset, by month

months <- c ("01", "02", "03", "04", "05", "06",
             "07", "08", "09", "10", "11", "12")

pred <- list()

for(i in 1:12){
  
  month  <- months[i]
  print(paste0(month, " ", Sys.time()))
  
  data             <- stage3 %>% dplyr::filter(mm == i)
  data$pred.rf.st3 <- predict(mod2,  data)
  data             <- data %>% dplyr::select(idcell, mm, day, Latitude, Longitude,
                                             pred.rf.st3)
  print(head(data))
  
  pred[[i]]        <- data
  
  rm(data)
  gc()
  
  print(paste0(month, " ", Sys.time(), " predicted!"))
  
}


final <- bind_rows(pred)
head(final)
dim(final)

setnames(final, "pred.rf.st3", "ln.pm25.pred.xgb")
final$pm25.pred.xgb <- exp(final$ln.pm25.pred.xgb)

head(final)
dim(final)

summary(final$pm25.pred.xgb)

tapply(final$pm25.pred.xgb, final$mm,
       function(x) format(summary(x)))

final %>% 
  dplyr::group_by(mm) %>% 
  dplyr::summarize(q01 = quantile(pm25.pred.xgb, 0.01),
                   q05 = quantile(pm25.pred.xgb, 0.05),
                   q25 = quantile(pm25.pred.xgb, 0.25),
                   q50 = quantile(pm25.pred.xgb, 0.50),
                   q75 = quantile(pm25.pred.xgb, 0.75),
                   q95 = quantile(pm25.pred.xgb, 0.95),
                   q99 = quantile(pm25.pred.xgb, 0.99))

saveRDS(final, "./stage3.2019...PM25.grid100m.XGB_5_2.RDS")

keep(mod2, sure = T)
gc()

#### 6_1 ####

##### Call stage3 dataset
stage3 <- readRDS("./stage3.2019.PM.all.grid100m_6_1.RDS")
names(stage3)
dim(stage3)
gc()

head(stage3)

# Define factor variables
stage3[,mm:=as.factor(mm)]
stage3[,dow:=as.factor(dow)]
stage3 <- as.data.frame(stage3)
gc()

##### Predict STAGE 1 model on Stage 3 dataset, by month

months <- c ("01", "02", "03", "04", "05", "06",
             "07", "08", "09", "10", "11", "12")

pred <- list()

for(i in 1:12){
  
  month  <- months[i]
  print(paste0(month, " ", Sys.time()))
  
  data             <- stage3 %>% dplyr::filter(mm == i)
  data$pred.rf.st3 <- predict(mod2,  data)
  data             <- data %>% dplyr::select(idcell, mm, day, Latitude, Longitude,
                                             pred.rf.st3)
  print(head(data))
  
  pred[[i]]        <- data
  
  rm(data)
  gc()
  
  print(paste0(month, " ", Sys.time(), " predicted!"))
  
}


final <- bind_rows(pred)
head(final)
dim(final)

setnames(final, "pred.rf.st3", "ln.pm25.pred.xgb")
final$pm25.pred.xgb <- exp(final$ln.pm25.pred.xgb)

head(final)
dim(final)

summary(final$pm25.pred.xgb)

tapply(final$pm25.pred.xgb, final$mm,
       function(x) format(summary(x)))

final %>% 
  dplyr::group_by(mm) %>% 
  dplyr::summarize(q01 = quantile(pm25.pred.xgb, 0.01),
                   q05 = quantile(pm25.pred.xgb, 0.05),
                   q25 = quantile(pm25.pred.xgb, 0.25),
                   q50 = quantile(pm25.pred.xgb, 0.50),
                   q75 = quantile(pm25.pred.xgb, 0.75),
                   q95 = quantile(pm25.pred.xgb, 0.95),
                   q99 = quantile(pm25.pred.xgb, 0.99))

saveRDS(final, "./stage3.2019...PM25.grid100m.XGB_6_1.RDS")

keep(mod2, sure = T)
gc()

#### 6_2 ####

##### Call stage3 dataset
stage3 <- readRDS("./stage3.2019.PM.all.grid100m_6_2.RDS")
names(stage3)
dim(stage3)
gc()

head(stage3)
# stage3 <- stage3 %>% dplyr::mutate(dewp_2m = dewp_2m - 273.15,
#                                    temp_2m = temp_2m - 273.15,
#                                    mslp = mslp/100)

# Define factor variables
stage3[,mm:=as.factor(mm)]
stage3[,dow:=as.factor(dow)]
# stage3 <- stage3 %>% dplyr::select(-c(geometry))
stage3 <- as.data.frame(stage3)
gc()

##### Predict STAGE 1 model on Stage 3 dataset, by month

months <- c ("01", "02", "03", "04", "05", "06",
             "07", "08", "09", "10", "11", "12")

pred <- list()

for(i in 1:12){
  
  month  <- months[i]
  print(paste0(month, " ", Sys.time()))
  
  data             <- stage3 %>% dplyr::filter(mm == i)
  data$pred.rf.st3 <- predict(mod2,  data)
  data             <- data %>% dplyr::select(idcell, mm, day, Latitude, Longitude,
                                             pred.rf.st3)
  print(head(data))
  
  pred[[i]]        <- data
  
  rm(data)
  gc()
  
  print(paste0(month, " ", Sys.time(), " predicted!"))
  
}


final <- bind_rows(pred)
head(final)
dim(final)

setnames(final, "pred.rf.st3", "ln.pm25.pred.xgb")
final$pm25.pred.xgb <- exp(final$ln.pm25.pred.xgb)

head(final)
dim(final)

summary(final$pm25.pred.xgb)

tapply(final$pm25.pred.xgb, final$mm,
       function(x) format(summary(x)))

final %>% 
  dplyr::group_by(mm) %>% 
  dplyr::summarize(q01 = quantile(pm25.pred.xgb, 0.01),
                   q05 = quantile(pm25.pred.xgb, 0.05),
                   q25 = quantile(pm25.pred.xgb, 0.25),
                   q50 = quantile(pm25.pred.xgb, 0.50),
                   q75 = quantile(pm25.pred.xgb, 0.75),
                   q95 = quantile(pm25.pred.xgb, 0.95),
                   q99 = quantile(pm25.pred.xgb, 0.99))

saveRDS(final, "./stage3.2019...PM25.grid100m.XGB_6_2.RDS")

keep(mod2, sure = T)
gc()