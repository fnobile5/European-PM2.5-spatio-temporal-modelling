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

#### RUN

mod.rf.xgb <- readRDS("./Ensemble/mod_rf_xgb.rds")

#### PREDICTION 5 KM ####

print("!!!!!!!!!!!!!!!!!!! PREDICTION 5 KM !!!!!!!!!!!!!!!!!!!")

grid5km.rf  <- readRDS("./stage3.2019.FINALPRED.PM25.grid5km.RF.RDS")
head(grid5km.rf)
dim(grid5km.rf)
setDT(grid5km.rf)
setkey(grid5km.rf, "idcell", "day")

grid5km.xgb <- readRDS("./stage3.2019.FINALPRED.PM25.grid5km.XGB.RDS")
head(grid5km.xgb)
dim(grid5km.xgb)
setDT(grid5km.xgb)
setkey(grid5km.xgb, "idcell", "day")
grid5km <- merge(grid5km.rf, grid5km.xgb[,c("idcell", "day", "ln.pm25.pred.xgb", "pm25.pred.xgb")])
head(grid5km)

rm(grid5km.rf, grid5km.xgb)
gc()

months <- c ("01", "02", "03", "04", "05", "06",
             "07", "08", "09", "10", "11", "12")

pred <- list()

for(i in 1:12){
  
  month  <- months[i]
  print(paste0(month, " ", Sys.time()))
  
  data                        <- grid5km %>% dplyr::filter(mm == i)
  data$ln.pm25.pred.rf.xgb    <- predict(mod.rf.xgb, data)
  data                        <- data %>% dplyr::select(idcell, mm, day, Latitude, Longitude,
                                                        ln.pm25.pred.rf.xgb, ln.pm25.pred.rf, ln.pm25.pred.xgb)
  print(head(data))
  
  pred[[i]]        <- data
  
  rm(data)
  gc()
  
  print(paste0(month, " ", Sys.time(), " predicted!"))
  
}

final <- bind_rows(pred)
head(final)
dim(final)
final$pm25.pred.rf.xgb <- exp(final$ln.pm25.pred.rf.xgb)
final$pm25.pred.rf     <- exp(final$ln.pm25.pred.rf)
final$pm25.pred.xgb    <- exp(final$ln.pm25.pred.xgb)
summary(final)

tapply(final$pm25.pred.rf.xgb, final$mm,
       function(x) format(summary(x)))

saveRDS(final, "./stage3.2019.FINALPRED.PM25.grid5km.ENSEMBLE.RF.XGB.RDS")

keep(mod.rf.xgb, sure = T)
gc()


#### PREDICTION 1 KM ####

print("!!!!!!!!!!!!!!!!!!! PREDICTION 1 KM !!!!!!!!!!!!!!!!!!!")

#### 1_1 ####

print("!!!!!!!!!!!!!!!!!!! PREDICTION 1 KM -> 1_1 !!!!!!!!!!!!!!!!!!!")

grid1km.rf  <- readRDS("./stage3.2019.FINALPRED.PM25.grid1km.RF_1_1.RDS")
head(grid1km.rf)
dim(grid1km.rf)
setDT(grid1km.rf)
setkey(grid1km.rf, "idcell", "day")

grid1km.xgb <- readRDS("./stage3.2019.FINALPRED.PM25.grid1km.XGB_1_1.RDS")
head(grid1km.xgb)
dim(grid1km.xgb)
setDT(grid1km.xgb)
setkey(grid1km.xgb, "idcell", "day")
grid1km <- merge(grid1km.rf, grid1km.xgb[,c("idcell", "day", "ln.pm25.pred.xgb", "pm25.pred.xgb")])
head(grid1km)

rm(grid1km.rf, grid1km.xgb)
gc()

months <- c ("01", "02", "03", "04", "05", "06",
             "07", "08", "09", "10", "11", "12")

pred <- list()

for(i in 1:12){
  
  month  <- months[i]
  print(paste0(month, " ", Sys.time()))
  
  data                        <- grid1km %>% dplyr::filter(mm == i)
  data$ln.pm25.pred.rf.xgb    <- predict(mod.rf.xgb, data)
  data                        <- data %>% dplyr::select(idcell, mm, day, Latitude, Longitude,
                                                        ln.pm25.pred.rf.xgb, ln.pm25.pred.rf, ln.pm25.pred.xgb)
  print(head(data))
  
  pred[[i]]        <- data
  
  rm(data)
  gc()
  
  print(paste0(month, " ", Sys.time(), " predicted!"))
  
}

final <- bind_rows(pred)
head(final)
dim(final)
final$pm25.pred.rf.xgb <- exp(final$ln.pm25.pred.rf.xgb)
final$pm25.pred.rf     <- exp(final$ln.pm25.pred.rf)
final$pm25.pred.xgb    <- exp(final$ln.pm25.pred.xgb)
summary(final)

tapply(final$pm25.pred.rf.xgb, final$mm,
       function(x) format(summary(x)))

saveRDS(final, "./stage3.2019.FINALPRED.PM25.grid1km.ENSEMBLE.RF.XGB_1_1.RDS")

keep(mod.rf.xgb, sure = T)
gc()


#### 1_2 ####

print("!!!!!!!!!!!!!!!!!!! PREDICTION 1 KM -> 1_2 !!!!!!!!!!!!!!!!!!!")

grid1km.rf  <- readRDS("./stage3.2019.FINALPRED.PM25.grid1km.RF_1_2.RDS")
head(grid1km.rf)
dim(grid1km.rf)
setDT(grid1km.rf)
setkey(grid1km.rf, "idcell", "day")

grid1km.xgb <- readRDS("./stage3.2019.FINALPRED.PM25.grid1km.XGB_1_2.RDS")
head(grid1km.xgb)
dim(grid1km.xgb)
setDT(grid1km.xgb)
setkey(grid1km.xgb, "idcell", "day")
grid1km <- merge(grid1km.rf, grid1km.xgb[,c("idcell", "day", "ln.pm25.pred.xgb", "pm25.pred.xgb")])
head(grid1km)

rm(grid1km.rf, grid1km.xgb)
gc()

months <- c ("01", "02", "03", "04", "05", "06",
             "07", "08", "09", "10", "11", "12")

pred <- list()

for(i in 1:12){
  
  month  <- months[i]
  print(paste0(month, " ", Sys.time()))
  
  data                        <- grid1km %>% dplyr::filter(mm == i)
  data$ln.pm25.pred.rf.xgb    <- predict(mod.rf.xgb, data)
  data                        <- data %>% dplyr::select(idcell, mm, day, Latitude, Longitude,
                                                        ln.pm25.pred.rf.xgb, ln.pm25.pred.rf, ln.pm25.pred.xgb)
  print(head(data))
  
  pred[[i]]        <- data
  
  rm(data)
  gc()
  
  print(paste0(month, " ", Sys.time(), " predicted!"))
  
}

final <- bind_rows(pred)
head(final)
dim(final)
final$pm25.pred.rf.xgb <- exp(final$ln.pm25.pred.rf.xgb)
final$pm25.pred.rf     <- exp(final$ln.pm25.pred.rf)
final$pm25.pred.xgb    <- exp(final$ln.pm25.pred.xgb)
summary(final)

tapply(final$pm25.pred.rf.xgb, final$mm,
       function(x) format(summary(x)))

saveRDS(final, "./stage3.2019.FINALPRED.PM25.grid1km.ENSEMBLE.RF.XGB_1_2.RDS")

keep(mod.rf.xgb, sure = T)
gc()


#### 2_1 ####

print("!!!!!!!!!!!!!!!!!!! PREDICTION 1 KM -> 2_1 !!!!!!!!!!!!!!!!!!!")

grid1km.rf  <- readRDS("./stage3.2019.FINALPRED.PM25.grid1km.RF_2_1.RDS")
head(grid1km.rf)
dim(grid1km.rf)
setDT(grid1km.rf)
setkey(grid1km.rf, "idcell", "day")

grid1km.xgb <- readRDS("./stage3.2019.FINALPRED.PM25.grid1km.XGB_2_1.RDS")
head(grid1km.xgb)
dim(grid1km.xgb)
setDT(grid1km.xgb)
setkey(grid1km.xgb, "idcell", "day")
grid1km <- merge(grid1km.rf, grid1km.xgb[,c("idcell", "day", "ln.pm25.pred.xgb", "pm25.pred.xgb")])
head(grid1km)

rm(grid1km.rf, grid1km.xgb)
gc()

months <- c ("01", "02", "03", "04", "05", "06",
             "07", "08", "09", "10", "11", "12")

pred <- list()

for(i in 1:12){
  
  month  <- months[i]
  print(paste0(month, " ", Sys.time()))
  
  data                        <- grid1km %>% dplyr::filter(mm == i)
  data$ln.pm25.pred.rf.xgb    <- predict(mod.rf.xgb, data)
  data                        <- data %>% dplyr::select(idcell, mm, day, Latitude, Longitude,
                                                        ln.pm25.pred.rf.xgb, ln.pm25.pred.rf, ln.pm25.pred.xgb)
  print(head(data))
  
  pred[[i]]        <- data
  
  rm(data)
  gc()
  
  print(paste0(month, " ", Sys.time(), " predicted!"))
  
}

final <- bind_rows(pred)
head(final)
dim(final)
final$pm25.pred.rf.xgb <- exp(final$ln.pm25.pred.rf.xgb)
final$pm25.pred.rf     <- exp(final$ln.pm25.pred.rf)
final$pm25.pred.xgb    <- exp(final$ln.pm25.pred.xgb)
summary(final)

tapply(final$pm25.pred.rf.xgb, final$mm,
       function(x) format(summary(x)))

saveRDS(final, "./stage3.2019.FINALPRED.PM25.grid1km.ENSEMBLE.RF.XGB_2_1.RDS")

keep(mod.rf.xgb, sure = T)
gc()


#### 2_2 ####

print("!!!!!!!!!!!!!!!!!!! PREDICTION 1 KM -> 2_2 !!!!!!!!!!!!!!!!!!!")

grid1km.rf  <- readRDS("./stage3.2019.FINALPRED.PM25.grid1km.RF_2_2.RDS")
head(grid1km.rf)
dim(grid1km.rf)
setDT(grid1km.rf)
setkey(grid1km.rf, "idcell", "day")

grid1km.xgb <- readRDS("./stage3.2019.FINALPRED.PM25.grid1km.XGB_2_2.RDS")
head(grid1km.xgb)
dim(grid1km.xgb)
setDT(grid1km.xgb)
setkey(grid1km.xgb, "idcell", "day")
grid1km <- merge(grid1km.rf, grid1km.xgb[,c("idcell", "day", "ln.pm25.pred.xgb", "pm25.pred.xgb")])
head(grid1km)

rm(grid1km.rf, grid1km.xgb)
gc()

months <- c ("01", "02", "03", "04", "05", "06",
             "07", "08", "09", "10", "11", "12")

pred <- list()

for(i in 1:12){
  
  month  <- months[i]
  print(paste0(month, " ", Sys.time()))
  
  data                        <- grid1km %>% dplyr::filter(mm == i)
  data$ln.pm25.pred.rf.xgb    <- predict(mod.rf.xgb, data)
  data                        <- data %>% dplyr::select(idcell, mm, day, Latitude, Longitude,
                                                        ln.pm25.pred.rf.xgb, ln.pm25.pred.rf, ln.pm25.pred.xgb)
  print(head(data))
  
  pred[[i]]        <- data
  
  rm(data)
  gc()
  
  print(paste0(month, " ", Sys.time(), " predicted!"))
  
}

final <- bind_rows(pred)
head(final)
dim(final)
final$pm25.pred.rf.xgb <- exp(final$ln.pm25.pred.rf.xgb)
final$pm25.pred.rf     <- exp(final$ln.pm25.pred.rf)
final$pm25.pred.xgb    <- exp(final$ln.pm25.pred.xgb)
summary(final)

tapply(final$pm25.pred.rf.xgb, final$mm,
       function(x) format(summary(x)))

saveRDS(final, "./stage3.2019.FINALPRED.PM25.grid1km.ENSEMBLE.RF.XGB_2_2.RDS")

keep(mod.rf.xgb, sure = T)
gc()


#### 3_1 ####

print("!!!!!!!!!!!!!!!!!!! PREDICTION 1 KM -> 3_1 !!!!!!!!!!!!!!!!!!!")

grid1km.rf  <- readRDS("./stage3.2019.FINALPRED.PM25.grid1km.RF_3_1.RDS")
head(grid1km.rf)
dim(grid1km.rf)
setDT(grid1km.rf)
setkey(grid1km.rf, "idcell", "day")

grid1km.xgb <- readRDS("./stage3.2019.FINALPRED.PM25.grid1km.XGB_3_1.RDS")
head(grid1km.xgb)
dim(grid1km.xgb)
setDT(grid1km.xgb)
setkey(grid1km.xgb, "idcell", "day")
grid1km <- merge(grid1km.rf, grid1km.xgb[,c("idcell", "day", "ln.pm25.pred.xgb", "pm25.pred.xgb")])
head(grid1km)

rm(grid1km.rf, grid1km.xgb)
gc()

months <- c ("01", "02", "03", "04", "05", "06",
             "07", "08", "09", "10", "11", "12")

pred <- list()

for(i in 1:12){
  
  month  <- months[i]
  print(paste0(month, " ", Sys.time()))
  
  data                        <- grid1km %>% dplyr::filter(mm == i)
  data$ln.pm25.pred.rf.xgb    <- predict(mod.rf.xgb, data)
  data                        <- data %>% dplyr::select(idcell, mm, day, Latitude, Longitude,
                                                        ln.pm25.pred.rf.xgb, ln.pm25.pred.rf, ln.pm25.pred.xgb)
  print(head(data))
  
  pred[[i]]        <- data
  
  rm(data)
  gc()
  
  print(paste0(month, " ", Sys.time(), " predicted!"))
  
}

final <- bind_rows(pred)
head(final)
dim(final)
final$pm25.pred.rf.xgb <- exp(final$ln.pm25.pred.rf.xgb)
final$pm25.pred.rf     <- exp(final$ln.pm25.pred.rf)
final$pm25.pred.xgb    <- exp(final$ln.pm25.pred.xgb)
summary(final)

tapply(final$pm25.pred.rf.xgb, final$mm,
       function(x) format(summary(x)))

saveRDS(final, "./stage3.2019.FINALPRED.PM25.grid1km.ENSEMBLE.RF.XGB_3_1.RDS")

keep(mod.rf.xgb, sure = T)
gc()


#### 3_2 ####

print("!!!!!!!!!!!!!!!!!!! PREDICTION 1 KM -> 3_2 !!!!!!!!!!!!!!!!!!!")

grid1km.rf  <- readRDS("./stage3.2019.FINALPRED.PM25.grid1km.RF_3_2.RDS")
head(grid1km.rf)
dim(grid1km.rf)
setDT(grid1km.rf)
setkey(grid1km.rf, "idcell", "day")

grid1km.xgb <- readRDS("./stage3.2019.FINALPRED.PM25.grid1km.XGB_3_2.RDS")
head(grid1km.xgb)
dim(grid1km.xgb)
setDT(grid1km.xgb)
setkey(grid1km.xgb, "idcell", "day")
grid1km <- merge(grid1km.rf, grid1km.xgb[,c("idcell", "day", "ln.pm25.pred.xgb", "pm25.pred.xgb")])
head(grid1km)

rm(grid1km.rf, grid1km.xgb)
gc()

months <- c ("01", "02", "03", "04", "05", "06",
             "07", "08", "09", "10", "11", "12")

pred <- list()

for(i in 1:12){
  
  month  <- months[i]
  print(paste0(month, " ", Sys.time()))
  
  data                        <- grid1km %>% dplyr::filter(mm == i)
  data$ln.pm25.pred.rf.xgb    <- predict(mod.rf.xgb, data)
  data                        <- data %>% dplyr::select(idcell, mm, day, Latitude, Longitude,
                                                        ln.pm25.pred.rf.xgb, ln.pm25.pred.rf, ln.pm25.pred.xgb)
  print(head(data))
  
  pred[[i]]        <- data
  
  rm(data)
  gc()
  
  print(paste0(month, " ", Sys.time(), " predicted!"))
  
}

final <- bind_rows(pred)
head(final)
dim(final)
final$pm25.pred.rf.xgb <- exp(final$ln.pm25.pred.rf.xgb)
final$pm25.pred.rf     <- exp(final$ln.pm25.pred.rf)
final$pm25.pred.xgb    <- exp(final$ln.pm25.pred.xgb)
summary(final)

tapply(final$pm25.pred.rf.xgb, final$mm,
       function(x) format(summary(x)))

saveRDS(final, "./stage3.2019.FINALPRED.PM25.grid1km.ENSEMBLE.RF.XGB_3_2.RDS")

keep(mod.rf.xgb, sure = T)
gc()


#### PREDICTION 100 M ####

print("!!!!!!!!!!!!!!!!!!! PREDICTION 100 M !!!!!!!!!!!!!!!!!!!")

#### 1_1 ####

print("!!!!!!!!!!!!!!!!!!! PREDICTION 100 M -> 1_1 !!!!!!!!!!!!!!!!!!!")

grid100m.rf  <- readRDS("./stage3.2019.FINALPRED.PM25.grid100m.RF_1_1.RDS")
head(grid100m.rf)
dim(grid100m.rf)
setDT(grid100m.rf)
setkey(grid100m.rf, "idcell", "day")

grid100m.xgb <- readRDS("./stage3.2019.FINALPRED.PM25.grid100m.XGB_1_1.RDS")
head(grid100m.xgb)
dim(grid100m.xgb)
setDT(grid100m.xgb)
setkey(grid100m.xgb, "idcell", "day")
grid100m <- merge(grid100m.rf, grid100m.xgb[,c("idcell", "day", "ln.pm25.pred.xgb", "pm25.pred.xgb")])
head(grid100m)

rm(grid100m.rf, grid100m.xgb)
gc()

months <- c ("01", "02", "03", "04", "05", "06",
             "07", "08", "09", "10", "11", "12")

pred <- list()

for(i in 1:12){
  
  month  <- months[i]
  print(paste0(month, " ", Sys.time()))
  
  data                        <- grid100m %>% dplyr::filter(mm == i)
  data$ln.pm25.pred.rf.xgb    <- predict(mod.rf.xgb, data)
  data                        <- data %>% dplyr::select(idcell, mm, day, Latitude, Longitude,
                                                        ln.pm25.pred.rf.xgb, ln.pm25.pred.rf, ln.pm25.pred.xgb)
  print(head(data))
  
  pred[[i]]        <- data
  
  rm(data)
  gc()
  
  print(paste0(month, " ", Sys.time(), " predicted!"))
  
}

final <- bind_rows(pred)
head(final)
dim(final)
final$pm25.pred.rf.xgb <- exp(final$ln.pm25.pred.rf.xgb)
final$pm25.pred.rf     <- exp(final$ln.pm25.pred.rf)
final$pm25.pred.xgb    <- exp(final$ln.pm25.pred.xgb)
summary(final)

tapply(final$pm25.pred.rf.xgb, final$mm,
       function(x) format(summary(x)))

saveRDS(final, "./stage3.2019.FINALPRED.PM25.grid100m.ENSEMBLE.RF.XGB_1_1.RDS")

keep(mod.rf.xgb, sure = T)
gc()


#### 1_2 ####

print("!!!!!!!!!!!!!!!!!!! PREDICTION 100 M -> 1_2 !!!!!!!!!!!!!!!!!!!")

grid100m.rf  <- readRDS("./stage3.2019.FINALPRED.PM25.grid100m.RF_1_2.RDS")
head(grid100m.rf)
dim(grid100m.rf)
setDT(grid100m.rf)
setkey(grid100m.rf, "idcell", "day")

grid100m.xgb <- readRDS("./stage3.2019.FINALPRED.PM25.grid100m.XGB_1_2.RDS")
head(grid100m.xgb)
dim(grid100m.xgb)
setDT(grid100m.xgb)
setkey(grid100m.xgb, "idcell", "day")
grid100m <- merge(grid100m.rf, grid100m.xgb[,c("idcell", "day", "ln.pm25.pred.xgb", "pm25.pred.xgb")])
head(grid100m)

rm(grid100m.rf, grid100m.xgb)
gc()

months <- c ("01", "02", "03", "04", "05", "06",
             "07", "08", "09", "10", "11", "12")

pred <- list()

for(i in 1:12){
  
  month  <- months[i]
  print(paste0(month, " ", Sys.time()))
  
  data                        <- grid100m %>% dplyr::filter(mm == i)
  data$ln.pm25.pred.rf.xgb    <- predict(mod.rf.xgb, data)
  data                        <- data %>% dplyr::select(idcell, mm, day, Latitude, Longitude,
                                                        ln.pm25.pred.rf.xgb, ln.pm25.pred.rf, ln.pm25.pred.xgb)
  print(head(data))
  
  pred[[i]]        <- data
  
  rm(data)
  gc()
  
  print(paste0(month, " ", Sys.time(), " predicted!"))
  
}

final <- bind_rows(pred)
head(final)
dim(final)
final$pm25.pred.rf.xgb <- exp(final$ln.pm25.pred.rf.xgb)
final$pm25.pred.rf     <- exp(final$ln.pm25.pred.rf)
final$pm25.pred.xgb    <- exp(final$ln.pm25.pred.xgb)
summary(final)

tapply(final$pm25.pred.rf.xgb, final$mm,
       function(x) format(summary(x)))

saveRDS(final, "./stage3.2019.FINALPRED.PM25.grid100m.ENSEMBLE.RF.XGB_1_2.RDS")

keep(mod.rf.xgb, sure = T)
gc()


#### 2_1 ####

print("!!!!!!!!!!!!!!!!!!! PREDICTION 100 M -> 2_1 !!!!!!!!!!!!!!!!!!!")

grid100m.rf  <- readRDS("./stage3.2019.FINALPRED.PM25.grid100m.RF_2_1.RDS")
head(grid100m.rf)
dim(grid100m.rf)
setDT(grid100m.rf)
setkey(grid100m.rf, "idcell", "day")

grid100m.xgb <- readRDS("./stage3.2019.FINALPRED.PM25.grid100m.XGB_2_1.RDS")
head(grid100m.xgb)
dim(grid100m.xgb)
setDT(grid100m.xgb)
setkey(grid100m.xgb, "idcell", "day")
grid100m <- merge(grid100m.rf, grid100m.xgb[,c("idcell", "day", "ln.pm25.pred.xgb", "pm25.pred.xgb")])
head(grid100m)

rm(grid100m.rf, grid100m.xgb)
gc()

months <- c ("01", "02", "03", "04", "05", "06",
             "07", "08", "09", "10", "11", "12")

pred <- list()

for(i in 1:12){
  
  month  <- months[i]
  print(paste0(month, " ", Sys.time()))
  
  data                        <- grid100m %>% dplyr::filter(mm == i)
  data$ln.pm25.pred.rf.xgb    <- predict(mod.rf.xgb, data)
  data                        <- data %>% dplyr::select(idcell, mm, day, Latitude, Longitude,
                                                        ln.pm25.pred.rf.xgb, ln.pm25.pred.rf, ln.pm25.pred.xgb)
  print(head(data))
  
  pred[[i]]        <- data
  
  rm(data)
  gc()
  
  print(paste0(month, " ", Sys.time(), " predicted!"))
  
}

final <- bind_rows(pred)
head(final)
dim(final)
final$pm25.pred.rf.xgb <- exp(final$ln.pm25.pred.rf.xgb)
final$pm25.pred.rf     <- exp(final$ln.pm25.pred.rf)
final$pm25.pred.xgb    <- exp(final$ln.pm25.pred.xgb)
summary(final)

tapply(final$pm25.pred.rf.xgb, final$mm,
       function(x) format(summary(x)))

saveRDS(final, "./stage3.2019.FINALPRED.PM25.grid100m.ENSEMBLE.RF.XGB_2_1.RDS")

keep(mod.rf.xgb, sure = T)
gc()


#### 2_2 ####

print("!!!!!!!!!!!!!!!!!!! PREDICTION 100 M -> 2_2 !!!!!!!!!!!!!!!!!!!")

grid100m.rf  <- readRDS("./stage3.2019.FINALPRED.PM25.grid100m.RF_2_2.RDS")
head(grid100m.rf)
dim(grid100m.rf)
setDT(grid100m.rf)
setkey(grid100m.rf, "idcell", "day")

grid100m.xgb <- readRDS("./stage3.2019.FINALPRED.PM25.grid100m.XGB_2_2.RDS")
head(grid100m.xgb)
dim(grid100m.xgb)
setDT(grid100m.xgb)
setkey(grid100m.xgb, "idcell", "day")
grid100m <- merge(grid100m.rf, grid100m.xgb[,c("idcell", "day", "ln.pm25.pred.xgb", "pm25.pred.xgb")])
head(grid100m)

rm(grid100m.rf, grid100m.xgb)
gc()

months <- c ("01", "02", "03", "04", "05", "06",
             "07", "08", "09", "10", "11", "12")

pred <- list()

for(i in 1:12){
  
  month  <- months[i]
  print(paste0(month, " ", Sys.time()))
  
  data                        <- grid100m %>% dplyr::filter(mm == i)
  data$ln.pm25.pred.rf.xgb    <- predict(mod.rf.xgb, data)
  data                        <- data %>% dplyr::select(idcell, mm, day, Latitude, Longitude,
                                                        ln.pm25.pred.rf.xgb, ln.pm25.pred.rf, ln.pm25.pred.xgb)
  print(head(data))
  
  pred[[i]]        <- data
  
  rm(data)
  gc()
  
  print(paste0(month, " ", Sys.time(), " predicted!"))
  
}

final <- bind_rows(pred)
head(final)
dim(final)
final$pm25.pred.rf.xgb <- exp(final$ln.pm25.pred.rf.xgb)
final$pm25.pred.rf     <- exp(final$ln.pm25.pred.rf)
final$pm25.pred.xgb    <- exp(final$ln.pm25.pred.xgb)
summary(final)

tapply(final$pm25.pred.rf.xgb, final$mm,
       function(x) format(summary(x)))

saveRDS(final, "./stage3.2019.FINALPRED.PM25.grid100m.ENSEMBLE.RF.XGB_2_2.RDS")

keep(mod.rf.xgb, sure = T)
gc()


#### 3_1 ####

print("!!!!!!!!!!!!!!!!!!! PREDICTION 100 M -> 3_1 !!!!!!!!!!!!!!!!!!!")

grid100m.rf  <- readRDS("./stage3.2019.FINALPRED.PM25.grid100m.RF_3_1.RDS")
head(grid100m.rf)
dim(grid100m.rf)
setDT(grid100m.rf)
setkey(grid100m.rf, "idcell", "day")

grid100m.xgb <- readRDS("./stage3.2019.FINALPRED.PM25.grid100m.XGB_3_1.RDS")
head(grid100m.xgb)
dim(grid100m.xgb)
setDT(grid100m.xgb)
setkey(grid100m.xgb, "idcell", "day")
grid100m <- merge(grid100m.rf, grid100m.xgb[,c("idcell", "day", "ln.pm25.pred.xgb", "pm25.pred.xgb")])
head(grid100m)

rm(grid100m.rf, grid100m.xgb)
gc()

months <- c ("01", "02", "03", "04", "05", "06",
             "07", "08", "09", "10", "11", "12")

pred <- list()

for(i in 1:12){
  
  month  <- months[i]
  print(paste0(month, " ", Sys.time()))
  
  data                        <- grid100m %>% dplyr::filter(mm == i)
  data$ln.pm25.pred.rf.xgb    <- predict(mod.rf.xgb, data)
  data                        <- data %>% dplyr::select(idcell, mm, day, Latitude, Longitude,
                                                        ln.pm25.pred.rf.xgb, ln.pm25.pred.rf, ln.pm25.pred.xgb)
  print(head(data))
  
  pred[[i]]        <- data
  
  rm(data)
  gc()
  
  print(paste0(month, " ", Sys.time(), " predicted!"))
  
}

final <- bind_rows(pred)
head(final)
dim(final)
final$pm25.pred.rf.xgb <- exp(final$ln.pm25.pred.rf.xgb)
final$pm25.pred.rf     <- exp(final$ln.pm25.pred.rf)
final$pm25.pred.xgb    <- exp(final$ln.pm25.pred.xgb)
summary(final)

tapply(final$pm25.pred.rf.xgb, final$mm,
       function(x) format(summary(x)))

saveRDS(final, "./stage3.2019.FINALPRED.PM25.grid100m.ENSEMBLE.RF.XGB_3_1.RDS")

keep(mod.rf.xgb, sure = T)
gc()


#### 3_2 ####

print("!!!!!!!!!!!!!!!!!!! PREDICTION 100 M -> 3_2 !!!!!!!!!!!!!!!!!!!")

grid100m.rf  <- readRDS("./stage3.2019.FINALPRED.PM25.grid100m.RF_3_2.RDS")
head(grid100m.rf)
dim(grid100m.rf)
setDT(grid100m.rf)
setkey(grid100m.rf, "idcell", "day")

grid100m.xgb <- readRDS("./stage3.2019.FINALPRED.PM25.grid100m.XGB_3_2.RDS")
head(grid100m.xgb)
dim(grid100m.xgb)
setDT(grid100m.xgb)
setkey(grid100m.xgb, "idcell", "day")
grid100m <- merge(grid100m.rf, grid100m.xgb[,c("idcell", "day", "ln.pm25.pred.xgb", "pm25.pred.xgb")])
head(grid100m)

rm(grid100m.rf, grid100m.xgb)
gc()

months <- c ("01", "02", "03", "04", "05", "06",
             "07", "08", "09", "10", "11", "12")

pred <- list()

for(i in 1:12){
  
  month  <- months[i]
  print(paste0(month, " ", Sys.time()))
  
  data                        <- grid100m %>% dplyr::filter(mm == i)
  data$ln.pm25.pred.rf.xgb    <- predict(mod.rf.xgb, data)
  data                        <- data %>% dplyr::select(idcell, mm, day, Latitude, Longitude,
                                                        ln.pm25.pred.rf.xgb, ln.pm25.pred.rf, ln.pm25.pred.xgb)
  print(head(data))
  
  pred[[i]]        <- data
  
  rm(data)
  gc()
  
  print(paste0(month, " ", Sys.time(), " predicted!"))
  
}

final <- bind_rows(pred)
head(final)
dim(final)
final$pm25.pred.rf.xgb <- exp(final$ln.pm25.pred.rf.xgb)
final$pm25.pred.rf     <- exp(final$ln.pm25.pred.rf)
final$pm25.pred.xgb    <- exp(final$ln.pm25.pred.xgb)
summary(final)

tapply(final$pm25.pred.rf.xgb, final$mm,
       function(x) format(summary(x)))

saveRDS(final, "./stage3.2019.FINALPRED.PM25.grid100m.ENSEMBLE.RF.XGB_3_2.RDS")

keep(mod.rf.xgb, sure = T)
gc()


#### 4_1 ####

print("!!!!!!!!!!!!!!!!!!! PREDICTION 100 M -> 4_1 !!!!!!!!!!!!!!!!!!!")

grid100m.rf  <- readRDS("./stage3.2019.FINALPRED.PM25.grid100m.RF_4_1.RDS")
head(grid100m.rf)
dim(grid100m.rf)
setDT(grid100m.rf)
setkey(grid100m.rf, "idcell", "day")

grid100m.xgb <- readRDS("./stage3.2019.FINALPRED.PM25.grid100m.XGB_4_1.RDS")
head(grid100m.xgb)
dim(grid100m.xgb)
setDT(grid100m.xgb)
setkey(grid100m.xgb, "idcell", "day")
grid100m <- merge(grid100m.rf, grid100m.xgb[,c("idcell", "day", "ln.pm25.pred.xgb", "pm25.pred.xgb")])
head(grid100m)

rm(grid100m.rf, grid100m.xgb)
gc()

months <- c ("01", "02", "03", "04", "05", "06",
             "07", "08", "09", "10", "11", "12")

pred <- list()

for(i in 1:12){
  
  month  <- months[i]
  print(paste0(month, " ", Sys.time()))
  
  data                        <- grid100m %>% dplyr::filter(mm == i)
  data$ln.pm25.pred.rf.xgb    <- predict(mod.rf.xgb, data)
  data                        <- data %>% dplyr::select(idcell, mm, day, Latitude, Longitude,
                                                        ln.pm25.pred.rf.xgb, ln.pm25.pred.rf, ln.pm25.pred.xgb)
  print(head(data))
  
  pred[[i]]        <- data
  
  rm(data)
  gc()
  
  print(paste0(month, " ", Sys.time(), " predicted!"))
  
}

final <- bind_rows(pred)
head(final)
dim(final)
final$pm25.pred.rf.xgb <- exp(final$ln.pm25.pred.rf.xgb)
final$pm25.pred.rf     <- exp(final$ln.pm25.pred.rf)
final$pm25.pred.xgb    <- exp(final$ln.pm25.pred.xgb)
summary(final)

tapply(final$pm25.pred.rf.xgb, final$mm,
       function(x) format(summary(x)))

saveRDS(final, "./stage3.2019.FINALPRED.PM25.grid100m.ENSEMBLE.RF.XGB_4_1.RDS")

keep(mod.rf.xgb, sure = T)
gc()


#### 4_2 ####

print("!!!!!!!!!!!!!!!!!!! PREDICTION 100 M -> 4_2 !!!!!!!!!!!!!!!!!!!")

grid100m.rf  <- readRDS("./stage3.2019.FINALPRED.PM25.grid100m.RF_4_2.RDS")
head(grid100m.rf)
dim(grid100m.rf)
setDT(grid100m.rf)
setkey(grid100m.rf, "idcell", "day")

grid100m.xgb <- readRDS("./stage3.2019.FINALPRED.PM25.grid100m.XGB_4_2.RDS")
head(grid100m.xgb)
dim(grid100m.xgb)
setDT(grid100m.xgb)
setkey(grid100m.xgb, "idcell", "day")
grid100m <- merge(grid100m.rf, grid100m.xgb[,c("idcell", "day", "ln.pm25.pred.xgb", "pm25.pred.xgb")])
head(grid100m)

rm(grid100m.rf, grid100m.xgb)
gc()

months <- c ("01", "02", "03", "04", "05", "06",
             "07", "08", "09", "10", "11", "12")

pred <- list()

for(i in 1:12){
  
  month  <- months[i]
  print(paste0(month, " ", Sys.time()))
  
  data                        <- grid100m %>% dplyr::filter(mm == i)
  data$ln.pm25.pred.rf.xgb    <- predict(mod.rf.xgb, data)
  data                        <- data %>% dplyr::select(idcell, mm, day, Latitude, Longitude,
                                                        ln.pm25.pred.rf.xgb, ln.pm25.pred.rf, ln.pm25.pred.xgb)
  print(head(data))
  
  pred[[i]]        <- data
  
  rm(data)
  gc()
  
  print(paste0(month, " ", Sys.time(), " predicted!"))
  
}

final <- bind_rows(pred)
head(final)
dim(final)
final$pm25.pred.rf.xgb <- exp(final$ln.pm25.pred.rf.xgb)
final$pm25.pred.rf     <- exp(final$ln.pm25.pred.rf)
final$pm25.pred.xgb    <- exp(final$ln.pm25.pred.xgb)
summary(final)

tapply(final$pm25.pred.rf.xgb, final$mm,
       function(x) format(summary(x)))

saveRDS(final, "./stage3.2019.FINALPRED.PM25.grid100m.ENSEMBLE.RF.XGB_4_2.RDS")

keep(mod.rf.xgb, sure = T)
gc()


#### 5_1 ####

print("!!!!!!!!!!!!!!!!!!! PREDICTION 100 M -> 5_1 !!!!!!!!!!!!!!!!!!!")

grid100m.rf  <- readRDS("./stage3.2019.FINALPRED.PM25.grid100m.RF_5_1.RDS")
head(grid100m.rf)
dim(grid100m.rf)
setDT(grid100m.rf)
setkey(grid100m.rf, "idcell", "day")

grid100m.xgb <- readRDS("./stage3.2019.FINALPRED.PM25.grid100m.XGB_5_1.RDS")
head(grid100m.xgb)
dim(grid100m.xgb)
setDT(grid100m.xgb)
setkey(grid100m.xgb, "idcell", "day")
grid100m <- merge(grid100m.rf, grid100m.xgb[,c("idcell", "day", "ln.pm25.pred.xgb", "pm25.pred.xgb")])
head(grid100m)

rm(grid100m.rf, grid100m.xgb)
gc()

months <- c ("01", "02", "03", "04", "05", "06",
             "07", "08", "09", "10", "11", "12")

pred <- list()

for(i in 1:12){
  
  month  <- months[i]
  print(paste0(month, " ", Sys.time()))
  
  data                        <- grid100m %>% dplyr::filter(mm == i)
  data$ln.pm25.pred.rf.xgb    <- predict(mod.rf.xgb, data)
  data                        <- data %>% dplyr::select(idcell, mm, day, Latitude, Longitude,
                                                        ln.pm25.pred.rf.xgb, ln.pm25.pred.rf, ln.pm25.pred.xgb)
  print(head(data))
  
  pred[[i]]        <- data
  
  rm(data)
  gc()
  
  print(paste0(month, " ", Sys.time(), " predicted!"))
  
}

final <- bind_rows(pred)
head(final)
dim(final)
final$pm25.pred.rf.xgb <- exp(final$ln.pm25.pred.rf.xgb)
final$pm25.pred.rf     <- exp(final$ln.pm25.pred.rf)
final$pm25.pred.xgb    <- exp(final$ln.pm25.pred.xgb)
summary(final)

tapply(final$pm25.pred.rf.xgb, final$mm,
       function(x) format(summary(x)))

saveRDS(final, "./stage3.2019.FINALPRED.PM25.grid100m.ENSEMBLE.RF.XGB_5_1.RDS")

keep(mod.rf.xgb, sure = T)
gc()


#### 5_2 ####

print("!!!!!!!!!!!!!!!!!!! PREDICTION 100 M -> 5_2 !!!!!!!!!!!!!!!!!!!")

grid100m.rf  <- readRDS("./stage3.2019.FINALPRED.PM25.grid100m.RF_5_2.RDS")
head(grid100m.rf)
dim(grid100m.rf)
setDT(grid100m.rf)
setkey(grid100m.rf, "idcell", "day")

grid100m.xgb <- readRDS("./stage3.2019.FINALPRED.PM25.grid100m.XGB_5_2.RDS")
head(grid100m.xgb)
dim(grid100m.xgb)
setDT(grid100m.xgb)
setkey(grid100m.xgb, "idcell", "day")
grid100m <- merge(grid100m.rf, grid100m.xgb[,c("idcell", "day", "ln.pm25.pred.xgb", "pm25.pred.xgb")])
head(grid100m)

rm(grid100m.rf, grid100m.xgb)
gc()

months <- c ("01", "02", "03", "04", "05", "06",
             "07", "08", "09", "10", "11", "12")

pred <- list()

for(i in 1:12){
  
  month  <- months[i]
  print(paste0(month, " ", Sys.time()))
  
  data                        <- grid100m %>% dplyr::filter(mm == i)
  data$ln.pm25.pred.rf.xgb    <- predict(mod.rf.xgb, data)
  data                        <- data %>% dplyr::select(idcell, mm, day, Latitude, Longitude,
                                                        ln.pm25.pred.rf.xgb, ln.pm25.pred.rf, ln.pm25.pred.xgb)
  print(head(data))
  
  pred[[i]]        <- data
  
  rm(data)
  gc()
  
  print(paste0(month, " ", Sys.time(), " predicted!"))
  
}

final <- bind_rows(pred)
head(final)
dim(final)
final$pm25.pred.rf.xgb <- exp(final$ln.pm25.pred.rf.xgb)
final$pm25.pred.rf     <- exp(final$ln.pm25.pred.rf)
final$pm25.pred.xgb    <- exp(final$ln.pm25.pred.xgb)
summary(final)

tapply(final$pm25.pred.rf.xgb, final$mm,
       function(x) format(summary(x)))

saveRDS(final, "./stage3.2019.FINALPRED.PM25.grid100m.ENSEMBLE.RF.XGB_5_2.RDS")

keep(mod.rf.xgb, sure = T)
gc()


#### 6_1 ####

print("!!!!!!!!!!!!!!!!!!! PREDICTION 100 M -> 6_1 !!!!!!!!!!!!!!!!!!!")

grid100m.rf  <- readRDS("./stage3.2019.FINALPRED.PM25.grid100m.RF_6_1.RDS")
head(grid100m.rf)
dim(grid100m.rf)
setDT(grid100m.rf)
setkey(grid100m.rf, "idcell", "day")

grid100m.xgb <- readRDS("./stage3.2019.FINALPRED.PM25.grid100m.XGB_6_1.RDS")
head(grid100m.xgb)
dim(grid100m.xgb)
setDT(grid100m.xgb)
setkey(grid100m.xgb, "idcell", "day")
grid100m <- merge(grid100m.rf, grid100m.xgb[,c("idcell", "day", "ln.pm25.pred.xgb", "pm25.pred.xgb")])
head(grid100m)

rm(grid100m.rf, grid100m.xgb)
gc()

months <- c ("01", "02", "03", "04", "05", "06",
             "07", "08", "09", "10", "11", "12")

pred <- list()

for(i in 1:12){
  
  month  <- months[i]
  print(paste0(month, " ", Sys.time()))
  
  data                        <- grid100m %>% dplyr::filter(mm == i)
  data$ln.pm25.pred.rf.xgb    <- predict(mod.rf.xgb, data)
  data                        <- data %>% dplyr::select(idcell, mm, day, Latitude, Longitude,
                                                        ln.pm25.pred.rf.xgb, ln.pm25.pred.rf, ln.pm25.pred.xgb)
  print(head(data))
  
  pred[[i]]        <- data
  
  rm(data)
  gc()
  
  print(paste0(month, " ", Sys.time(), " predicted!"))
  
}

final <- bind_rows(pred)
head(final)
dim(final)
final$pm25.pred.rf.xgb <- exp(final$ln.pm25.pred.rf.xgb)
final$pm25.pred.rf     <- exp(final$ln.pm25.pred.rf)
final$pm25.pred.xgb    <- exp(final$ln.pm25.pred.xgb)
summary(final)

tapply(final$pm25.pred.rf.xgb, final$mm,
       function(x) format(summary(x)))

saveRDS(final, "./stage3.2019.FINALPRED.PM25.grid100m.ENSEMBLE.RF.XGB_6_1.RDS")

keep(mod.rf.xgb, sure = T)
gc()


#### 6_2 ####

print("!!!!!!!!!!!!!!!!!!! PREDICTION 100 M -> 6_2 !!!!!!!!!!!!!!!!!!!")

grid100m.rf  <- readRDS("./stage3.2019.FINALPRED.PM25.grid100m.RF_6_2.RDS")
head(grid100m.rf)
dim(grid100m.rf)
setDT(grid100m.rf)
setkey(grid100m.rf, "idcell", "day")

grid100m.xgb <- readRDS("./stage3.2019.FINALPRED.PM25.grid100m.XGB_6_2.RDS")
head(grid100m.xgb)
dim(grid100m.xgb)
setDT(grid100m.xgb)
setkey(grid100m.xgb, "idcell", "day")
grid100m <- merge(grid100m.rf, grid100m.xgb[,c("idcell", "day", "ln.pm25.pred.xgb", "pm25.pred.xgb")])
head(grid100m)

rm(grid100m.rf, grid100m.xgb)
gc()

months <- c ("01", "02", "03", "04", "05", "06",
             "07", "08", "09", "10", "11", "12")

pred <- list()

for(i in 1:12){
  
  month  <- months[i]
  print(paste0(month, " ", Sys.time()))
  
  data                        <- grid100m %>% dplyr::filter(mm == i)
  data$ln.pm25.pred.rf.xgb    <- predict(mod.rf.xgb, data)
  data                        <- data %>% dplyr::select(idcell, mm, day, Latitude, Longitude,
                                                        ln.pm25.pred.rf.xgb, ln.pm25.pred.rf, ln.pm25.pred.xgb)
  print(head(data))
  
  pred[[i]]        <- data
  
  rm(data)
  gc()
  
  print(paste0(month, " ", Sys.time(), " predicted!"))
  
}

final <- bind_rows(pred)
head(final)
dim(final)
final$pm25.pred.rf.xgb <- exp(final$ln.pm25.pred.rf.xgb)
final$pm25.pred.rf     <- exp(final$ln.pm25.pred.rf)
final$pm25.pred.xgb    <- exp(final$ln.pm25.pred.xgb)
summary(final)

tapply(final$pm25.pred.rf.xgb, final$mm,
       function(x) format(summary(x)))

saveRDS(final, "./stage3.2019.FINALPRED.PM25.grid100m.ENSEMBLE.RF.XGB_6_2.RDS")

keep(mod.rf.xgb, sure = T)
gc()