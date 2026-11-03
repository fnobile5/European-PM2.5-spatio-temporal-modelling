#### European-PM2.5-spatio-temporal-modelling

#### Author: Federica Nobile
#### Affiliation: Institute of Environmental Medicine, Karolinska Institutet, Stockholm, Sweden
####              Department of Epidemiology, Lazio Region Health Service/ASL Roma 1, Rome, Italy

library(dplyr)
library(lubridate)
library(splitstackshape)
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

#### Call the dataset for stage 1

data <- readRDS("./stage1.PM.final.RDS")
dim(data) # 
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

data1 <- data1 %>% dplyr::mutate(ln.pm25 = log(final.pm25 + 0.5),
                                 mm      = as.factor(mm),
                                 dow     = as.factor(dow))
summary(data1)


##########################################################################
##########################################################################
#####                      RANDOM FOREST                             #####
##########################################################################
##########################################################################

# Define the formula

formula.rf <- as.formula(ln.pm25 ~ day + mm +
                           Longitude + Latitude +
                           allRoads_50 + allRoads_100 + allRoads_200 + allRoads_300 + 
                           imd_200 + imd_300 + 
                           tbu_300 + 
                           pop + alt10_enh +
                           u10 + v10 + dewp_2m + temp_2m + mslp + tp + blh +
                           leaf_high + leaf_low + tcc + tco +
                           pm25_cams +
                           NDVI + EVI +
                           final.aod047)

set.seed(1)
mod.rf <- ranger(formula.rf, data=data1, num.trees=1250, mtry=10, importance="impurity", keep.inbag=TRUE)
mod.rf$r.squared
mod.rf$prediction.error
mod_imp <- mod.rf
importance <- sort(importance(mod_imp), decreasing = TRUE)
importance <- importance/sum(importance)*100
cbind(importance)

saveRDS(mod.rf,"./mod_rf.rds")

data1$pred.rf <- predict(mod.rf, data1)$predictions

##########################################################################
##########################################################################
#####                 EXTREME GRADIENT BOOSTING                      #####
##########################################################################
##########################################################################

formula.xgb <- as.formula(ln.pm25 ~ day + dow + mm +
                            Longitude + Latitude +
                            allRoads_50 + allRoads_100 + allRoads_200 + allRoads_300 + 
                            majorRoads_100 + majorRoads_300 + 
                            imd_200 + imd_300 + 
                            nat_300 + 
                            tbu_300 + 
                            pop + alt10_enh +
                            u10 + v10 + dewp_2m + temp_2m + mslp + 
                            tp + blh + leaf_high + leaf_low + tcc + tco +
                            pm25_cams +
                            NDVI + EVI +
                            final.aod047)

set.seed(1)
mod.xgb <- train(formula.xgb,
                 data = data1,
                 method = "xgbTree",
                 trControl = fitxgb,
                 verbose = TRUE,
                 tuneGrid = xgbGrid)

mod.xgb
 
importBoost <- varImp(mod.xgb)
importBoost$importance

saveRDS(mod.xgb,"./mod_xgb.rds")

data1$pred.xgb <- predict(mod.xgb, data1)$predictions

##########################################################################
##########################################################################
#####                       ENSEMBLE MODEL                           #####
##########################################################################
##########################################################################

mod.en <- bam(ln.pm25 ~ s(Longitude, Latitude, by=pred.rf)+s(Longitude, Latitude, by=pred.xgb), data = data1)
saveRDS(mod.en,"./mod_rf_xgb.rds")

# Needed elements
formula.rf

formula.xgb

fitxgb <- trainControl(method = "none")

xgbGrid <- expand.grid(nrounds = 1000, 
                       max_depth = 10,
                       colsample_bytree = 0.8, 
                       eta = 0.1,
                       gamma = 0,
                       min_child_weight = 5,
                       subsample = 1)

##################################
##################################
#### 10-FOLD CROSS-VALIDATION ####
####      NESTED SCHEME       ####
##################################
##################################

sites.pm           <- read.csv("R code/spatial-blocked cv folds.csv")
names(sites.pm)

data1              <- as.data.table(data1)
sites.pm           <- as.data.table(sites.pm)
setnames(data1, "Air.Quality.Station.EoI.Code", "site")
setnames(sites.pm, "Air.Quality.Station.EoI.Code", "site")
setnames(sites.pm, "group", "split")
setkey(data1,    site, CountryName)
setkey(sites.pm, site, CountryName)
data2              <- merge(data1, sites.pm, all.x = T) 
head(data2)
summary(data2)
table(data2$split)

data2$site         <- as.character(data2$site) 
names(data2)
dim(data2) 
table(data2$split)
table(data2$split, data2$CountryName)


final.en <- data.frame()
train.en <- list()
valid.en   <- list()

for (i in 1:10){
  
  print(i)
  print(paste0(Sys.time(), " start..."))
  
  valid   <- subset(data2, split==i)
  portion <- subset(data2, split!=i)
  
  final.test  <- data.frame()
  final.valid <- data.frame()
  
  for(j in 1:9){
    
    n <- unique(portion$split)
    
    print(n)
    
    test                  <- subset(portion, split == n[j])
    train                 <- subset(portion, split != n[j])
    
    print(paste0("N (valid + test + train) = ", dim(valid)[1] + dim(test)[1] + dim(train)[1]))
    
    print(paste0(Sys.time(), " RF starts..."))
    mod.rf                <- ranger(formula.rf, data=train, num.trees=1250, mtry=10, importance="impurity", seed=1)
    test$pred.rf          <- predict(mod.rf, test)$predictions
    valid$pred.rf         <- predict(mod.rf, valid)$predictions
    print(paste0(Sys.time(), " RF done"))
    
    set.seed(1)
    print(paste0(Sys.time(), " XGB starts..."))
    mod.xgb <- train(formula.xgb,
                     data = train,
                     method = "xgbTree",
                     trControl = fitxgb,
                     verbose = TRUE,
                     tuneGrid = xgbGrid)
    test$pred.xgb        <- predict(mod.xgb, test)
    valid$pred.xgb       <- predict(mod.xgb, valid)
    print(paste0(Sys.time(), " XGB done"))
    
    test$itercv <- n[j]
    test        <- test %>% dplyr::select(day, CountryName, site, Longitude, Latitude,
                                          pm25, final.pm25, ln.pm25, pred.rf, pred.xgb, itercv)
    final.test <- rbind(final.test, test)
    
    valid$itercv <- i
    valid2       <- valid %>% dplyr::select(day, CountryName, site, Longitude, Latitude,
                                            pm25, final.pm25, ln.pm25, pred.rf, pred.xgb, itercv)
    final.valid  <- rbind(final.valid, valid2)
    
    valid        <- valid %>% dplyr::select(!c(pred.rf, pred.xgb, itercv))
    
    print(paste0(Sys.time(), " j = ", j, " -> ", n[j], " COMPLETED!"))
    
    rm(train, test, mod.rf, valid2, mod.xgb)
    gc()
    
  }
  
  final.valid2   <- final.valid %>% dplyr::group_by(CountryName, site, Longitude, Latitude, itercv, day, pm25, final.pm25, ln.pm25) %>%
    dplyr::summarise(pred.rf = mean(pred.rf, na.rm = T),
              pred.xgb = mean(pred.xgb, na.rm = T))
  final.valid2   <- as.data.frame(final.valid2)
  
  valid.en[[i]] <- final.valid
  saveRDS(final.valid, paste0("./Ensemble/valid_en_", i, ".RDS"))
  
  train.en[[i]] <- final.test
  saveRDS(final.test, paste0("./Ensemble/train_en_", i, ".RDS"))
  
  rm(final.valid, final.valid2, final.test)
  gc()
  
  mod.en <- bam(ln.pm25 ~ s(Longitude, Latitude, by=pred.rf)+s(Longitude, Latitude, by=pred.xgb), data = train.en[[i]])
  
  valid.en[[i]]$pred.en <- predict(object=mod.en, newdata = valid.en[[i]])
  valid.en[[i]]$itercv  <- i
  valid.en[[i]]         <- valid.en[[i]] %>% dplyr::select(day, CountryName, site, Longitude, Latitude, pm25, final.pm25, ln.pm25, pred.rf, pred.xgb, pred.en, itercv)
  
  print(paste0(Sys.time(), " i = ", i, " COMPLETED!"))
  
  saveRDS(valid.en[[i]], paste0("./Ensemble/valid_en_def_", i, ".RDS"))
  
  rm(mod.en)
  gc()
  
}

cv10.regcal.rf <- function(dataset, exp) 
{
  for (i in 1:10)
  {
    print(i)
    test                  <- subset(dataset, itercv==i)
    train                 <- subset(dataset, itercv!=i)
    train$appoggio        <- train[[exp]]
    out_train             <- lm(appoggio~pred.rf.cv,data=train)
    test$pred.rf.regcal   <- predict.lm(object=out_train, newdata=test, allow.new.levels=TRUE, re.form=NULL )
    assign(paste("test_s",i,sep=""),test)
    rm("test","train","out_train")
  }
  
  #BIND 1 dataset
  mod3.cv <<- data.table(rbind(test_s1,test_s2,test_s3,test_s4,test_s5,test_s6,test_s7,test_s8,test_s9,test_s10))
  mod3.cv <- as.data.table(mod3.cv)
  setkey(mod3.cv,site)
  rm(list = ls(pattern = "test_"))
  
  #table updates
  mod3.cv$buttare <- mod3.cv[[exp]]
  m2.fit.all.cv   <- lm(buttare~pred.rf.regcal, data=mod3.cv)
  
  cv.r2.all <- summary(m2.fit.all.cv)$r.squared
  # cv.rmspe  <- rmse(residuals(m2.fit.all.cv))
  cv.rmspe  <- rmse(m2.fit.all.cv$model$buttare, m2.fit.all.cv$fitted.values)
  cv.inter  <- summary(m2.fit.all.cv)$coef[1,1]
  cv.slope  <- summary(m2.fit.all.cv)$coef[2,1]
  
  #spatial
  spatialall.cv     <- data.frame(site=unique(mod3.cv$site), barpm=as.numeric(tapply(mod3.cv[[exp]],mod3.cv$site,mean,na.rm=T)),
                                  barpred=as.numeric(tapply(mod3.cv$pred.rf.regcal,mod3.cv$site,mean,na.rm=T)))
  spatialall.cv$site <- as.character(spatialall.cv$site)
  spatialall.cv <- as.data.table(spatialall.cv)
  setkey(spatialall.cv,site)
  
  m2.fit.all.cv.s   <- lm(barpm ~ barpred, data=spatialall.cv)
  
  cv.r2.spatial     <- summary(m2.fit.all.cv.s)$r.squared
  # cv.rmse.spatial   <- rmse(residuals(m2.fit.all.cv.s))
  cv.rmse.spatial   <- rmse(m2.fit.all.cv.s$model$barpm, m2.fit.all.cv.s$fitted.values)
  cv.spatial.inter  <- summary(m2.fit.all.cv.s)$coef[1,1]
  cv.spatial.slope  <- summary(m2.fit.all.cv.s)$coef[2,1]
  
  #temporal
  tempoall.cv         <- merge(mod3.cv,spatialall.cv)
  tempoall.cv$delpm   <- tempoall.cv[[exp]]-tempoall.cv$barpm
  tempoall.cv$delpred <- tempoall.cv$pred.rf.regcal-tempoall.cv$barpred
  mod_temporal.cv     <- lm(delpm ~ delpred, data=tempoall.cv)
  
  cv.r2.temporal     <- summary(mod_temporal.cv)$r.squared
  # cv.rmse.temporal   <- rmse(residuals(mod_temporal.cv))
  cv.rmse.temporal   <- rmse(mod_temporal.cv$model$delpm, mod_temporal.cv$fitted.values)
  cv.temporal.inter  <- summary(mod_temporal.cv)$coef[1,1]
  cv.temporal.slope  <- summary(mod_temporal.cv)$coef[2,1]
  
  #results
  options(scipen=999)
  results <<- rbind(cv.r2.all,      cv.rmspe,         cv.inter,          cv.slope,
                    cv.r2.spatial,  cv.rmse.spatial,  cv.spatial.inter,  cv.spatial.slope, 
                    cv.r2.temporal, cv.rmse.temporal, cv.temporal.inter, cv.temporal.slope)
  results <- round(results,digits=3)
  print(results)
}

data.learners <- list()
data.en       <- list()

for (i in 1:10){
  
  print(paste0(i, " starts...", Sys.time()))
  
  train <- readRDS(paste0("train_en_", i, ".RDS"))
  head(train)
  data.learners[[i]] <- train
  
  valid <- readRDS(paste0("valid_en_def_", i, ".RDS"))
  head(valid)
  data.en[[i]] <- valid
  
  rm(train, valid)
  print(paste0(i, " DONE ", Sys.time()))
  gc()
  
}

#### Performance - learners ####

perf.rf.log    <- list()
perf.rf.pm     <- list()
perf.rf.log.rc <- list()
perf.rf.pm.rc  <- list()

perf.xgb.log    <- list()
perf.xgb.pm     <- list()
perf.xgb.log.rc <- list()
perf.xgb.pm.rc  <- list()

for (i in 1:10){
  
  print(paste0(i, " starts...", Sys.time()))
  
  #### RANDOM FOREST
  
  mod2.cv            <- data.learners[[i]]
  mod2.cv$pred.rf.cv <- mod2.cv$pred.rf
  
  
  ### Check on the log.scale
  
  # Overall
  m1.fit.all  <- lm(ln.pm25~pred.rf.cv, data=mod2.cv)
  r2.all      <- summary(m1.fit.all)$r.squared
  rmspe       <- rmse(m1.fit.all$model$ln.pm25, m1.fit.all$fitted.values)
  inter       <- summary(m1.fit.all)$coef[1,1]
  slope       <- summary(m1.fit.all)$coef[2,1]
  
  #spatial
  spatialall      <- data.frame(site=unique(mod2.cv$site), barpm=as.numeric(tapply(mod2.cv$ln.pm25,mod2.cv$site,mean,na.rm=T)),
                                barpred=as.numeric(tapply(mod2.cv$pred.rf.cv,mod2.cv$site,mean,na.rm=T)))
  spatialall$site <- as.character(spatialall$site)
  m1.fit.all.s    <- lm(barpm ~ barpred, data=spatialall)
  r2.spatial      <- summary(m1.fit.all.s)$r.squared
  rmse.spatial       <- rmse(m1.fit.all.s$model$barpm, m1.fit.all.s$fitted.values)
  spatial.inter   <- summary(m1.fit.all.s)$coef[1,1]
  spatial.slope   <- summary(m1.fit.all.s)$coef[2,1]
  
  #temporal
  tempoall        <- left_join(mod2.cv,spatialall)
  delpm           <- tempoall$ln.pm25-tempoall$barpm
  delpred         <- tempoall$pred.rf.cv-tempoall$barpred
  mod_temporal    <- lm(delpm ~ delpred)
  r2.temporal     <- summary(mod_temporal)$r.squared
  rmse.temporal   <- rmse(mod_temporal$model$delpm, mod_temporal$fitted.values)
  temporal.inter  <- summary(mod_temporal)$coef[1,1]
  temporal.slope  <- summary(mod_temporal)$coef[2,1]
  
  #results
  options(scipen=999)
  cv.results <<- rbind(r2.all,      rmspe,         inter,          slope,
                       r2.spatial,  rmse.spatial,  spatial.inter,  spatial.slope, 
                       r2.temporal, rmse.temporal, temporal.inter, temporal.slope)
  cv.results <- round(cv.results,digits=3)
  print(cv.results)
  perf.rf.log[[i]]    <- cv.results
  rm(cv.results)
  
  
  ### Check on the PM scale
  
  mod2.cv$exp.pred.rf.cv <- exp(mod2.cv$pred.rf.cv)
  
  # Overall
  m1.fit.all  <- lm(final.pm25~exp.pred.rf.cv, data=mod2.cv)
  r2.all      <- summary(m1.fit.all)$r.squared
  rmspe       <- rmse(m1.fit.all$model$final.pm25, m1.fit.all$fitted.values)
  inter2      <- summary(m1.fit.all)$coef[1,1]
  slope2      <- summary(m1.fit.all)$coef[2,1]
  
  #spatial
  spatialall  <- data.frame(site=unique(mod2.cv$site), barpm=as.numeric(tapply(mod2.cv$final.pm25,mod2.cv$site,mean,na.rm=T)),
                            barpred=as.numeric(tapply(mod2.cv$exp.pred.rf.cv,mod2.cv$site,mean,na.rm=T)))
  spatialall$site <- as.character(spatialall$site)
  m1.fit.all.s    <- lm(barpm ~ barpred, data=spatialall)
  r2.spatial      <- summary(m1.fit.all.s)$r.squared
  rmse.spatial       <- rmse(m1.fit.all.s$model$barpm, m1.fit.all.s$fitted.values)
  spatial.inter   <- summary(m1.fit.all.s)$coef[1,1]
  spatial.slope   <- summary(m1.fit.all.s)$coef[2,1]
  
  #temporal
  tempoall        <- left_join(mod2.cv,spatialall)
  delpm           <- tempoall$final.pm25-tempoall$barpm
  delpred         <- tempoall$exp.pred.rf.cv-tempoall$barpred
  mod_temporal    <- lm(delpm ~ delpred)
  r2.temporal     <- summary(mod_temporal)$r.squared
  rmse.temporal   <- rmse(mod_temporal$model$delpm, mod_temporal$fitted.values)
  temporal.inter  <- summary(mod_temporal)$coef[1,1]
  temporal.slope  <- summary(mod_temporal)$coef[2,1]
  
  #results
  options(scipen=999)
  cv.results <<- rbind(r2.all,      rmspe,         inter2,         slope2,
                       r2.spatial,  rmse.spatial,  spatial.inter,  spatial.slope, 
                       r2.temporal, rmse.temporal, temporal.inter, temporal.slope)
  cv.results <- round(cv.results,digits=3)
  print(cv.results)
  perf.rf.pm[[i]]     <- cv.results
  rm(cv.results)
  
  
  ### Regression calibration (on the log scale)
  
  cv.results <- cv10.regcal.rf(mod2.cv, "ln.pm25")
  perf.rf.log.rc[[i]] <- cv.results
  rm(cv.results)
  
  
  ### Check regression calibration results on the PM scale
  
  mod3.cv$exp.pred.rf.regcal <- exp(mod3.cv$pred.rf.regcal)
  
  # Overall
  m1.fit.all  <- lm(final.pm25~exp.pred.rf.regcal, data=mod3.cv)
  r2.all      <- summary(m1.fit.all)$r.squared
  rmspe       <- rmse(m1.fit.all$model$final.pm25, m1.fit.all$fitted.values)
  inter2      <- summary(m1.fit.all)$coef[1,1]
  slope2      <- summary(m1.fit.all)$coef[2,1]
  
  #spatial
  spatialall  <- data.frame(site=unique(mod3.cv$site), barpm=as.numeric(tapply(mod3.cv$final.pm25,mod3.cv$site,mean,na.rm=T)),
                            barpred=as.numeric(tapply(mod3.cv$exp.pred.rf.regcal,mod3.cv$site,mean,na.rm=T)))
  spatialall$site <- as.character(spatialall$site)
  m1.fit.all.s    <- lm(barpm ~ barpred, data=spatialall)
  r2.spatial      <- summary(m1.fit.all.s)$r.squared
  rmse.spatial    <- rmse(m1.fit.all.s$model$barpm, m1.fit.all.s$fitted.values)
  spatial.inter   <- summary(m1.fit.all.s)$coef[1,1]
  spatial.slope   <- summary(m1.fit.all.s)$coef[2,1]
  
  #temporal
  tempoall        <- left_join(mod3.cv,spatialall)
  delpm           <- tempoall$final.pm25-tempoall$barpm
  delpred         <- tempoall$exp.pred.rf.regcal-tempoall$barpred
  mod_temporal    <- lm(delpm ~ delpred)
  r2.temporal     <- summary(mod_temporal)$r.squared
  rmse.temporal   <- rmse(mod_temporal$model$delpm, mod_temporal$fitted.values)
  temporal.inter  <- summary(mod_temporal)$coef[1,1]
  temporal.slope  <- summary(mod_temporal)$coef[2,1]
  
  #results
  options(scipen=999)
  cv.results <<- rbind(r2.all,      rmspe,         inter2,         slope2,
                       r2.spatial,  rmse.spatial,  spatial.inter,  spatial.slope, 
                       r2.temporal, rmse.temporal, temporal.inter, temporal.slope)
  cv.results <- round(cv.results,digits=3)
  print(cv.results)
  perf.rf.pm.rc[[i]]  <- cv.results
  rm(cv.results)
  
  print(paste0(i, " RF DONE ", Sys.time()))
  
  
  
  #### EXTREME GRADIENT BOOSTING
  
  mod2.cv            <- data.learners[[i]]
  mod2.cv$pred.rf.cv <- mod2.cv$pred.xgb
  
  
  ### Check on the log.scale
  
  # Overall
  m1.fit.all  <- lm(ln.pm25~pred.rf.cv, data=mod2.cv)
  r2.all      <- summary(m1.fit.all)$r.squared
  rmspe       <- rmse(m1.fit.all$model$ln.pm25, m1.fit.all$fitted.values)
  inter       <- summary(m1.fit.all)$coef[1,1]
  slope       <- summary(m1.fit.all)$coef[2,1]
  
  #spatial
  spatialall      <- data.frame(site=unique(mod2.cv$site), barpm=as.numeric(tapply(mod2.cv$ln.pm25,mod2.cv$site,mean,na.rm=T)),
                                barpred=as.numeric(tapply(mod2.cv$pred.rf.cv,mod2.cv$site,mean,na.rm=T)))
  spatialall$site <- as.character(spatialall$site)
  m1.fit.all.s    <- lm(barpm ~ barpred, data=spatialall)
  r2.spatial      <- summary(m1.fit.all.s)$r.squared
  rmse.spatial       <- rmse(m1.fit.all.s$model$barpm, m1.fit.all.s$fitted.values)
  spatial.inter   <- summary(m1.fit.all.s)$coef[1,1]
  spatial.slope   <- summary(m1.fit.all.s)$coef[2,1]
  
  #temporal
  tempoall        <- left_join(mod2.cv,spatialall)
  delpm           <- tempoall$ln.pm25-tempoall$barpm
  delpred         <- tempoall$pred.rf.cv-tempoall$barpred
  mod_temporal    <- lm(delpm ~ delpred)
  r2.temporal     <- summary(mod_temporal)$r.squared
  rmse.temporal   <- rmse(mod_temporal$model$delpm, mod_temporal$fitted.values)
  temporal.inter  <- summary(mod_temporal)$coef[1,1]
  temporal.slope  <- summary(mod_temporal)$coef[2,1]
  
  #results
  options(scipen=999)
  cv.results <<- rbind(r2.all,      rmspe,         inter,          slope,
                       r2.spatial,  rmse.spatial,  spatial.inter,  spatial.slope, 
                       r2.temporal, rmse.temporal, temporal.inter, temporal.slope)
  cv.results <- round(cv.results,digits=3)
  print(cv.results)
  perf.xgb.log[[i]]    <- cv.results
  rm(cv.results)
  
  
  ### Check on the PM scale
  
  mod2.cv$exp.pred.rf.cv <- exp(mod2.cv$pred.rf.cv)
  
  # Overall
  m1.fit.all  <- lm(final.pm25~exp.pred.rf.cv, data=mod2.cv)
  r2.all      <- summary(m1.fit.all)$r.squared
  rmspe       <- rmse(m1.fit.all$model$final.pm25, m1.fit.all$fitted.values)
  inter2      <- summary(m1.fit.all)$coef[1,1]
  slope2      <- summary(m1.fit.all)$coef[2,1]
  
  #spatial
  spatialall  <- data.frame(site=unique(mod2.cv$site), barpm=as.numeric(tapply(mod2.cv$final.pm25,mod2.cv$site,mean,na.rm=T)),
                            barpred=as.numeric(tapply(mod2.cv$exp.pred.rf.cv,mod2.cv$site,mean,na.rm=T)))
  spatialall$site <- as.character(spatialall$site)
  m1.fit.all.s    <- lm(barpm ~ barpred, data=spatialall)
  r2.spatial      <- summary(m1.fit.all.s)$r.squared
  rmse.spatial       <- rmse(m1.fit.all.s$model$barpm, m1.fit.all.s$fitted.values)
  spatial.inter   <- summary(m1.fit.all.s)$coef[1,1]
  spatial.slope   <- summary(m1.fit.all.s)$coef[2,1]
  
  #temporal
  tempoall        <- left_join(mod2.cv,spatialall)
  delpm           <- tempoall$final.pm25-tempoall$barpm
  delpred         <- tempoall$exp.pred.rf.cv-tempoall$barpred
  mod_temporal    <- lm(delpm ~ delpred)
  r2.temporal     <- summary(mod_temporal)$r.squared
  rmse.temporal   <- rmse(mod_temporal$model$delpm, mod_temporal$fitted.values)
  temporal.inter  <- summary(mod_temporal)$coef[1,1]
  temporal.slope  <- summary(mod_temporal)$coef[2,1]
  
  #results
  options(scipen=999)
  cv.results <<- rbind(r2.all,      rmspe,         inter2,         slope2,
                       r2.spatial,  rmse.spatial,  spatial.inter,  spatial.slope, 
                       r2.temporal, rmse.temporal, temporal.inter, temporal.slope)
  cv.results <- round(cv.results,digits=3)
  print(cv.results)
  perf.xgb.pm[[i]]     <- cv.results
  rm(cv.results)
  
  
  ### Regression calibration (on the log scale)
  
  cv.results <- cv10.regcal.rf(mod2.cv, "ln.pm25")
  perf.xgb.log.rc[[i]] <- cv.results
  rm(cv.results)
  
  
  ### Check regression calibration results on the PM scale
  
  mod3.cv$exp.pred.rf.regcal <- exp(mod3.cv$pred.rf.regcal)
  
  # Overall
  m1.fit.all  <- lm(final.pm25~exp.pred.rf.regcal, data=mod3.cv)
  r2.all      <- summary(m1.fit.all)$r.squared
  rmspe       <- rmse(m1.fit.all$model$final.pm25, m1.fit.all$fitted.values)
  inter2      <- summary(m1.fit.all)$coef[1,1]
  slope2      <- summary(m1.fit.all)$coef[2,1]
  
  #spatial
  spatialall  <- data.frame(site=unique(mod3.cv$site), barpm=as.numeric(tapply(mod3.cv$final.pm25,mod3.cv$site,mean,na.rm=T)),
                            barpred=as.numeric(tapply(mod3.cv$exp.pred.rf.regcal,mod3.cv$site,mean,na.rm=T)))
  spatialall$site <- as.character(spatialall$site)
  m1.fit.all.s    <- lm(barpm ~ barpred, data=spatialall)
  r2.spatial      <- summary(m1.fit.all.s)$r.squared
  rmse.spatial    <- rmse(m1.fit.all.s$model$barpm, m1.fit.all.s$fitted.values)
  spatial.inter   <- summary(m1.fit.all.s)$coef[1,1]
  spatial.slope   <- summary(m1.fit.all.s)$coef[2,1]
  
  #temporal
  tempoall        <- left_join(mod3.cv,spatialall)
  delpm           <- tempoall$final.pm25-tempoall$barpm
  delpred         <- tempoall$exp.pred.rf.regcal-tempoall$barpred
  mod_temporal    <- lm(delpm ~ delpred)
  r2.temporal     <- summary(mod_temporal)$r.squared
  rmse.temporal   <- rmse(mod_temporal$model$delpm, mod_temporal$fitted.values)
  temporal.inter  <- summary(mod_temporal)$coef[1,1]
  temporal.slope  <- summary(mod_temporal)$coef[2,1]
  
  #results
  options(scipen=999)
  cv.results <<- rbind(r2.all,      rmspe,         inter2,         slope2,
                       r2.spatial,  rmse.spatial,  spatial.inter,  spatial.slope, 
                       r2.temporal, rmse.temporal, temporal.inter, temporal.slope)
  cv.results <- round(cv.results,digits=3)
  print(cv.results)
  perf.xgb.pm.rc[[i]]  <- cv.results
  rm(cv.results)
  
  print(paste0(i, " XGB DONE ", Sys.time()))
  
}

table.perf.rf.log <- cbind(as.data.frame(perf.rf.log[[1]]),
                           as.data.frame(perf.rf.log[[2]]),
                           as.data.frame(perf.rf.log[[3]]),
                           as.data.frame(perf.rf.log[[4]]),
                           as.data.frame(perf.rf.log[[5]]),
                           as.data.frame(perf.rf.log[[6]]),
                           as.data.frame(perf.rf.log[[7]]),
                           as.data.frame(perf.rf.log[[8]]),
                           as.data.frame(perf.rf.log[[9]]))
colnames(table.perf.rf.log) <- c("iter1", "iter2", "iter3", "iter4", "iter5",
                                 "iter6", "iter7", "iter8", "iter9")
write_xlsx(table.perf.rf.log, "CV performance rf - log scale.xlsx")

table.perf.rf.pm <- cbind(as.data.frame(perf.rf.pm[[1]]),
                          as.data.frame(perf.rf.pm[[2]]),
                          as.data.frame(perf.rf.pm[[3]]),
                          as.data.frame(perf.rf.pm[[4]]),
                          as.data.frame(perf.rf.pm[[5]]),
                          as.data.frame(perf.rf.pm[[6]]),
                          as.data.frame(perf.rf.pm[[7]]),
                          as.data.frame(perf.rf.pm[[8]]),
                          as.data.frame(perf.rf.pm[[9]]))
colnames(table.perf.rf.pm) <- c("iter1", "iter2", "iter3", "iter4", "iter5",
                                "iter6", "iter7", "iter8", "iter9")
write_xlsx(table.perf.rf.pm, "CV performance rf - pm scale.xlsx")

table.perf.rf.log.rc <- cbind(as.data.frame(perf.rf.log.rc[[1]]),
                              as.data.frame(perf.rf.log.rc[[2]]),
                              as.data.frame(perf.rf.log.rc[[3]]),
                              as.data.frame(perf.rf.log.rc[[4]]),
                              as.data.frame(perf.rf.log.rc[[5]]),
                              as.data.frame(perf.rf.log.rc[[6]]),
                              as.data.frame(perf.rf.log.rc[[7]]),
                              as.data.frame(perf.rf.log.rc[[8]]),
                              as.data.frame(perf.rf.log.rc[[9]]))
colnames(table.perf.rf.log.rc) <- c("iter1", "iter2", "iter3", "iter4", "iter5",
                                    "iter6", "iter7", "iter8", "iter9")
write_xlsx(table.perf.rf.log.rc, "CV performance rf - log scale regcal.xlsx")

table.perf.rf.pm.rc <- cbind(as.data.frame(perf.rf.pm.rc[[1]]),
                             as.data.frame(perf.rf.pm.rc[[2]]),
                             as.data.frame(perf.rf.pm.rc[[3]]),
                             as.data.frame(perf.rf.pm.rc[[4]]),
                             as.data.frame(perf.rf.pm.rc[[5]]),
                             as.data.frame(perf.rf.pm.rc[[6]]),
                             as.data.frame(perf.rf.pm.rc[[7]]),
                             as.data.frame(perf.rf.pm.rc[[8]]),
                             as.data.frame(perf.rf.pm.rc[[9]]))
colnames(table.perf.rf.pm.rc) <- c("iter1", "iter2", "iter3", "iter4", "iter5",
                                   "iter6", "iter7", "iter8", "iter9")
write_xlsx(table.perf.rf.pm.rc, "CV performance rf - pm scale regcal.xlsx")


table.perf.xgb.log <- cbind(as.data.frame(perf.xgb.log[[1]]),
                            as.data.frame(perf.xgb.log[[2]]),
                            as.data.frame(perf.xgb.log[[3]]),
                            as.data.frame(perf.xgb.log[[4]]),
                            as.data.frame(perf.xgb.log[[5]]),
                            as.data.frame(perf.xgb.log[[6]]),
                            as.data.frame(perf.xgb.log[[7]]),
                            as.data.frame(perf.xgb.log[[8]]),
                            as.data.frame(perf.xgb.log[[9]]))
colnames(table.perf.xgb.log) <- c("iter1", "iter2", "iter3", "iter4", "iter5",
                                  "iter6", "iter7", "iter8", "iter9")
write_xlsx(table.perf.xgb.log, "CV performance xgb - log scale.xlsx")

table.perf.xgb.pm <- cbind(as.data.frame(perf.xgb.pm[[1]]),
                           as.data.frame(perf.xgb.pm[[2]]),
                           as.data.frame(perf.xgb.pm[[3]]),
                           as.data.frame(perf.xgb.pm[[4]]),
                           as.data.frame(perf.xgb.pm[[5]]),
                           as.data.frame(perf.xgb.pm[[6]]),
                           as.data.frame(perf.xgb.pm[[7]]),
                           as.data.frame(perf.xgb.pm[[8]]),
                           as.data.frame(perf.xgb.pm[[9]]))
colnames(table.perf.xgb.pm) <- c("iter1", "iter2", "iter3", "iter4", "iter5",
                                 "iter6", "iter7", "iter8", "iter9")
write_xlsx(table.perf.xgb.pm, "CV performance xgb - pm scale.xlsx")

table.perf.xgb.log.rc <- cbind(as.data.frame(perf.xgb.log.rc[[1]]),
                               as.data.frame(perf.xgb.log.rc[[2]]),
                               as.data.frame(perf.xgb.log.rc[[3]]),
                               as.data.frame(perf.xgb.log.rc[[4]]),
                               as.data.frame(perf.xgb.log.rc[[5]]),
                               as.data.frame(perf.xgb.log.rc[[6]]),
                               as.data.frame(perf.xgb.log.rc[[7]]),
                               as.data.frame(perf.xgb.log.rc[[8]]),
                               as.data.frame(perf.xgb.log.rc[[9]]))
colnames(table.perf.xgb.log.rc) <- c("iter1", "iter2", "iter3", "iter4", "iter5",
                                     "iter6", "iter7", "iter8", "iter9")
write_xlsx(table.perf.xgb.log.rc, "CV performance xgb - log scale regcal.xlsx")

table.perf.xgb.pm.rc <- cbind(as.data.frame(perf.xgb.pm.rc[[1]]),
                              as.data.frame(perf.xgb.pm.rc[[2]]),
                              as.data.frame(perf.xgb.pm.rc[[3]]),
                              as.data.frame(perf.xgb.pm.rc[[4]]),
                              as.data.frame(perf.xgb.pm.rc[[5]]),
                              as.data.frame(perf.xgb.pm.rc[[6]]),
                              as.data.frame(perf.xgb.pm.rc[[7]]),
                              as.data.frame(perf.xgb.pm.rc[[8]]),
                              as.data.frame(perf.xgb.pm.rc[[9]]))
colnames(table.perf.xgb.pm.rc) <- c("iter1", "iter2", "iter3", "iter4", "iter5",
                                    "iter6", "iter7", "iter8", "iter9")
write_xlsx(table.perf.xgb.pm.rc, "CV performance xgb - pm scale regcal.xlsx")


#### Performance - ensemble rf+xgb ####

mod2.cv            <- bind_rows(data.en)
dim(mod2.cv)
table(mod2.cv$itercv)
mod2.cv$pred.rf.cv <- mod2.cv$pred.en
head(mod2.cv)


### Check on the log.scale

# Overall
m1.fit.all  <- lm(ln.pm25~pred.rf.cv, data=mod2.cv)
r2.all      <- summary(m1.fit.all)$r.squared
rmspe       <- rmse(m1.fit.all$model$ln.pm25, m1.fit.all$fitted.values)
inter       <- summary(m1.fit.all)$coef[1,1]
slope       <- summary(m1.fit.all)$coef[2,1]

#spatial
spatialall      <- data.frame(site=unique(mod2.cv$site), barpm=as.numeric(tapply(mod2.cv$ln.pm25,mod2.cv$site,mean,na.rm=T)),
                              barpred=as.numeric(tapply(mod2.cv$pred.rf.cv,mod2.cv$site,mean,na.rm=T)))
spatialall$site <- as.character(spatialall$site)
m1.fit.all.s    <- lm(barpm ~ barpred, data=spatialall)
r2.spatial      <- summary(m1.fit.all.s)$r.squared
rmse.spatial       <- rmse(m1.fit.all.s$model$barpm, m1.fit.all.s$fitted.values)
spatial.inter   <- summary(m1.fit.all.s)$coef[1,1]
spatial.slope   <- summary(m1.fit.all.s)$coef[2,1]

#temporal
tempoall        <- left_join(mod2.cv,spatialall)
delpm           <- tempoall$ln.pm25-tempoall$barpm
delpred         <- tempoall$pred.rf.cv-tempoall$barpred
mod_temporal    <- lm(delpm ~ delpred)
r2.temporal     <- summary(mod_temporal)$r.squared
rmse.temporal   <- rmse(mod_temporal$model$delpm, mod_temporal$fitted.values)
temporal.inter  <- summary(mod_temporal)$coef[1,1]
temporal.slope  <- summary(mod_temporal)$coef[2,1]

#results
options(scipen=999)
cv.results <<- rbind(r2.all,      rmspe,         inter,          slope,
                     r2.spatial,  rmse.spatial,  spatial.inter,  spatial.slope, 
                     r2.temporal, rmse.temporal, temporal.inter, temporal.slope)
cv.results <- round(cv.results,digits=3)
print(cv.results)
# r2.all          0.780
# rmspe           0.295
# inter           0.002
# slope           0.999
# r2.spatial      0.783
# rmse.spatial    0.178
# spatial.inter   0.044
# spatial.slope   0.983
# r2.temporal     0.785
# rmse.temporal   0.345
# temporal.inter -0.002
# temporal.slope  0.996


### Check on the PM scale

mod2.cv$exp.pred.rf.cv <- exp(mod2.cv$pred.rf.cv)

# Overall
m1.fit.all  <- lm(final.pm25~exp.pred.rf.cv, data=mod2.cv)
r2.all      <- summary(m1.fit.all)$r.squared
rmspe       <- rmse(m1.fit.all$model$final.pm25, m1.fit.all$fitted.values)
inter2      <- summary(m1.fit.all)$coef[1,1]
slope2      <- summary(m1.fit.all)$coef[2,1]

#spatial
spatialall  <- data.frame(site=unique(mod2.cv$site), barpm=as.numeric(tapply(mod2.cv$final.pm25,mod2.cv$site,mean,na.rm=T)),
                          barpred=as.numeric(tapply(mod2.cv$exp.pred.rf.cv,mod2.cv$site,mean,na.rm=T)))
spatialall$site <- as.character(spatialall$site)
m1.fit.all.s    <- lm(barpm ~ barpred, data=spatialall)
r2.spatial      <- summary(m1.fit.all.s)$r.squared
rmse.spatial       <- rmse(m1.fit.all.s$model$barpm, m1.fit.all.s$fitted.values)
spatial.inter   <- summary(m1.fit.all.s)$coef[1,1]
spatial.slope   <- summary(m1.fit.all.s)$coef[2,1]

#temporal
tempoall        <- left_join(mod2.cv,spatialall)
delpm           <- tempoall$final.pm25-tempoall$barpm
delpred         <- tempoall$exp.pred.rf.cv-tempoall$barpred
mod_temporal    <- lm(delpm ~ delpred)
r2.temporal     <- summary(mod_temporal)$r.squared
rmse.temporal   <- rmse(mod_temporal$model$delpm, mod_temporal$fitted.values)
temporal.inter  <- summary(mod_temporal)$coef[1,1]
temporal.slope  <- summary(mod_temporal)$coef[2,1]

#results
options(scipen=999)
cv.results <<- rbind(r2.all,      rmspe,         inter2,         slope2,
                     r2.spatial,  rmse.spatial,  spatial.inter,  spatial.slope, 
                     r2.temporal, rmse.temporal, temporal.inter, temporal.slope)
cv.results <- round(cv.results,digits=3)
print(cv.results)

### Regression calibration (on the log scale)

cv.results <- cv10.regcal.rf(mod2.cv, "ln.pm25")

### Check regression calibration results on the PM scale

mod3.cv$exp.pred.rf.regcal <- exp(mod3.cv$pred.rf.regcal)

# Overall
m1.fit.all  <- lm(final.pm25~exp.pred.rf.regcal, data=mod3.cv)
r2.all      <- summary(m1.fit.all)$r.squared
rmspe       <- rmse(m1.fit.all$model$final.pm25, m1.fit.all$fitted.values)
inter2      <- summary(m1.fit.all)$coef[1,1]
slope2      <- summary(m1.fit.all)$coef[2,1]

#spatial
spatialall  <- data.frame(site=unique(mod3.cv$site), barpm=as.numeric(tapply(mod3.cv$final.pm25,mod3.cv$site,mean,na.rm=T)),
                          barpred=as.numeric(tapply(mod3.cv$exp.pred.rf.regcal,mod3.cv$site,mean,na.rm=T)))
spatialall$site <- as.character(spatialall$site)
m1.fit.all.s    <- lm(barpm ~ barpred, data=spatialall)
r2.spatial      <- summary(m1.fit.all.s)$r.squared
rmse.spatial    <- rmse(m1.fit.all.s$model$barpm, m1.fit.all.s$fitted.values)
spatial.inter   <- summary(m1.fit.all.s)$coef[1,1]
spatial.slope   <- summary(m1.fit.all.s)$coef[2,1]

#temporal
tempoall        <- left_join(mod3.cv,spatialall)
delpm           <- tempoall$final.pm25-tempoall$barpm
delpred         <- tempoall$exp.pred.rf.regcal-tempoall$barpred
mod_temporal    <- lm(delpm ~ delpred)
r2.temporal     <- summary(mod_temporal)$r.squared
rmse.temporal   <- rmse(mod_temporal$model$delpm, mod_temporal$fitted.values)
temporal.inter  <- summary(mod_temporal)$coef[1,1]
temporal.slope  <- summary(mod_temporal)$coef[2,1]

#results
options(scipen=999)
cv.results <<- rbind(r2.all,      rmspe,         inter2,         slope2,
                     r2.spatial,  rmse.spatial,  spatial.inter,  spatial.slope, 
                     r2.temporal, rmse.temporal, temporal.inter, temporal.slope)
cv.results <- round(cv.results,digits=3)
print(cv.results)
