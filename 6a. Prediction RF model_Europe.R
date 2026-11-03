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


mod2 <- readRDS("./Ensemble/mod_rf.rds")
mod2$r.squared
mod2$prediction.error
mod_imp <- mod2
importance <- sort(importance(mod_imp), decreasing = TRUE)
importance <- importance/sum(importance)*100
cbind(importance)

rm(mod_imp)
gc()


#### PREDICTION 5 KM ####

print("!!!!!!!!!!!!!!!!!!! PREDICTION 5 KM !!!!!!!!!!!!!!!!!!!")

##### Call stage3 dataset (complete on grid cells)
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
  data$pred.rf.st3 <- predict(mod2,  data)$predictions
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

setnames(final, "pred.rf.st3", "ln.pm25.pred.rf")
final$pm25.pred.rf <- exp(final$ln.pm25.pred.rf)

head(final)
dim(final)

summary(final$pm25.pred.rf)

tapply(final$pm25.pred.rf, final$mm,
       function(x) format(summary(x)))
final %>% 
  dplyr::group_by(mm) %>% 
  dplyr::summarize(q01 = quantile(pm25.pred.rf, 0.01),
                   q05 = quantile(pm25.pred.rf, 0.05),
                   q25 = quantile(pm25.pred.rf, 0.25),
                   q50 = quantile(pm25.pred.rf, 0.50),
                   q75 = quantile(pm25.pred.rf, 0.75),
                   q95 = quantile(pm25.pred.rf, 0.95),
                   q99 = quantile(pm25.pred.rf, 0.99))

saveRDS(final, ".l/stage3.2019...PM25.grid5km.RF.RDS")

keep(mod2, sure = T)
gc()


#### PREDICTION 1 KM ####

print("!!!!!!!!!!!!!!!!!!! PREDICTION 1 KM !!!!!!!!!!!!!!!!!!!")

#### 1_1 ####

##### Call stage3 dataset (complete on grid cells)
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
  data$pred.rf.st3 <- predict(mod2,  data)$predictions
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

setnames(final, "pred.rf.st3", "ln.pm25.pred.rf")
final$pm25.pred.rf <- exp(final$ln.pm25.pred.rf)

head(final)
dim(final) 

summary(final$pm25.pred.rf)

tapply(final$pm25.pred.rf, final$mm,
       function(x) format(summary(x)))

final %>% 
  dplyr::group_by(mm) %>% 
  dplyr::summarize(q01 = quantile(pm25.pred.rf, 0.01),
                   q05 = quantile(pm25.pred.rf, 0.05),
                   q25 = quantile(pm25.pred.rf, 0.25),
                   q50 = quantile(pm25.pred.rf, 0.50),
                   q75 = quantile(pm25.pred.rf, 0.75),
                   q95 = quantile(pm25.pred.rf, 0.95),
                   q99 = quantile(pm25.pred.rf, 0.99))

saveRDS(final, "./stage3.2019...PM25.grid1km.RF_1_1.RDS")

keep(mod2, sure = T)
gc()

#### 1_2 ####

##### Call stage3 dataset (complete on grid cells)
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
  data$pred.rf.st3 <- predict(mod2,  data)$predictions
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

setnames(final, "pred.rf.st3", "ln.pm25.pred.rf")
final$pm25.pred.rf <- exp(final$ln.pm25.pred.rf)

head(final)

summary(final$pm25.pred.rf)


tapply(final$pm25.pred.rf, final$mm,
       function(x) format(summary(x)))

final %>% 
  dplyr::group_by(mm) %>% 
  dplyr::summarize(q01 = quantile(pm25.pred.rf, 0.01),
                   q05 = quantile(pm25.pred.rf, 0.05),
                   q25 = quantile(pm25.pred.rf, 0.25),
                   q50 = quantile(pm25.pred.rf, 0.50),
                   q75 = quantile(pm25.pred.rf, 0.75),
                   q95 = quantile(pm25.pred.rf, 0.95),
                   q99 = quantile(pm25.pred.rf, 0.99))

saveRDS(final, "./stage3.2019...PM25.grid1km.RF_1_2.RDS")

keep(mod2, sure = T)
gc()

#### 2_1 ####

##### Call stage3 dataset (complete on grid cells)
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
  data$pred.rf.st3 <- predict(mod2,  data)$predictions
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

setnames(final, "pred.rf.st3", "ln.pm25.pred.rf")
final$pm25.pred.rf <- exp(final$ln.pm25.pred.rf)

head(final)
dim(final)

summary(final$pm25.pred.rf)


tapply(final$pm25.pred.rf, final$mm,
       function(x) format(summary(x)))

final %>% 
  dplyr::group_by(mm) %>% 
  dplyr::summarize(q01 = quantile(pm25.pred.rf, 0.01),
                   q05 = quantile(pm25.pred.rf, 0.05),
                   q25 = quantile(pm25.pred.rf, 0.25),
                   q50 = quantile(pm25.pred.rf, 0.50),
                   q75 = quantile(pm25.pred.rf, 0.75),
                   q95 = quantile(pm25.pred.rf, 0.95),
                   q99 = quantile(pm25.pred.rf, 0.99))

saveRDS(final, "./stage3.2019...PM25.grid1km.RF_2_1.RDS")

keep(mod2, sure = T)
gc()

#### 2_2 ####

##### Call stage3 dataset (complete on grid cells)
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
  data$pred.rf.st3 <- predict(mod2,  data)$predictions
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

setnames(final, "pred.rf.st3", "ln.pm25.pred.rf")
final$pm25.pred.rf <- exp(final$ln.pm25.pred.rf)

head(final)
dim(final)

summary(final$pm25.pred.rf)

tapply(final$pm25.pred.rf, final$mm,
       function(x) format(summary(x)))

final %>% 
  dplyr::group_by(mm) %>% 
  dplyr::summarize(q01 = quantile(pm25.pred.rf, 0.01),
                   q05 = quantile(pm25.pred.rf, 0.05),
                   q25 = quantile(pm25.pred.rf, 0.25),
                   q50 = quantile(pm25.pred.rf, 0.50),
                   q75 = quantile(pm25.pred.rf, 0.75),
                   q95 = quantile(pm25.pred.rf, 0.95),
                   q99 = quantile(pm25.pred.rf, 0.99))

saveRDS(final, "./stage3.2019...PM25.grid1km.RF_2_2.RDS")

keep(mod2, sure = T)
gc()

#### 3_1 ####

##### Call stage3 dataset (complete on grid cells)
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
  data$pred.rf.st3 <- predict(mod2,  data)$predictions
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

setnames(final, "pred.rf.st3", "ln.pm25.pred.rf")
final$pm25.pred.rf <- exp(final$ln.pm25.pred.rf)

head(final)
dim(final)

summary(final$pm25.pred.rf)

tapply(final$pm25.pred.rf, final$mm,
       function(x) format(summary(x)))

final %>% 
  dplyr::group_by(mm) %>% 
  dplyr::summarize(q01 = quantile(pm25.pred.rf, 0.01),
                   q05 = quantile(pm25.pred.rf, 0.05),
                   q25 = quantile(pm25.pred.rf, 0.25),
                   q50 = quantile(pm25.pred.rf, 0.50),
                   q75 = quantile(pm25.pred.rf, 0.75),
                   q95 = quantile(pm25.pred.rf, 0.95),
                   q99 = quantile(pm25.pred.rf, 0.99))

saveRDS(final, "./stage3.2019...PM25.grid1km.RF_3_1.RDS")

keep(mod2, sure = T)
gc()

#### 3_2 ####

##### Call stage3 dataset (complete on grid cells)
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
  data$pred.rf.st3 <- predict(mod2,  data)$predictions
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

setnames(final, "pred.rf.st3", "ln.pm25.pred.rf")
final$pm25.pred.rf <- exp(final$ln.pm25.pred.rf)

head(final)
dim(final)

summary(final$pm25.pred.rf)

tapply(final$pm25.pred.rf, final$mm,
       function(x) format(summary(x)))

final %>% 
  dplyr::group_by(mm) %>% 
  dplyr::summarize(q01 = quantile(pm25.pred.rf, 0.01),
                   q05 = quantile(pm25.pred.rf, 0.05),
                   q25 = quantile(pm25.pred.rf, 0.25),
                   q50 = quantile(pm25.pred.rf, 0.50),
                   q75 = quantile(pm25.pred.rf, 0.75),
                   q95 = quantile(pm25.pred.rf, 0.95),
                   q99 = quantile(pm25.pred.rf, 0.99))

saveRDS(final, "./stage3.2019...PM25.grid1km.RF_3_2.RDS")

keep(mod2, sure = T)
gc()


#### PREDICTION 100M ####

print("!!!!!!!!!!!!!!!!!!! PREDICTION 100 M !!!!!!!!!!!!!!!!!!!")

#### 1_1 ####

##### Call stage3 dataset (complete on grid cells)
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
  data$pred.rf.st3 <- predict(mod2,  data)$predictions
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

setnames(final, "pred.rf.st3", "ln.pm25.pred.rf")
final$pm25.pred.rf <- exp(final$ln.pm25.pred.rf)

head(final)
dim(final)

summary(final$pm25.pred.rf)

tapply(final$pm25.pred.rf, final$mm,
       function(x) format(summary(x)))

final %>% 
  dplyr::group_by(mm) %>% 
  dplyr::summarize(q01 = quantile(pm25.pred.rf, 0.01),
                   q05 = quantile(pm25.pred.rf, 0.05),
                   q25 = quantile(pm25.pred.rf, 0.25),
                   q50 = quantile(pm25.pred.rf, 0.50),
                   q75 = quantile(pm25.pred.rf, 0.75),
                   q95 = quantile(pm25.pred.rf, 0.95),
                   q99 = quantile(pm25.pred.rf, 0.99))

saveRDS(final, "./stage3.2019...PM25.grid100m.RF_1_1.RDS")

keep(mod2, sure = T)
gc()

#### 1_2 ####

##### Call stage3 dataset (complete on grid cells)
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
  data$pred.rf.st3 <- predict(mod2,  data)$predictions
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

setnames(final, "pred.rf.st3", "ln.pm25.pred.rf")
final$pm25.pred.rf <- exp(final$ln.pm25.pred.rf)

head(final)
dim(final)

summary(final$pm25.pred.rf)

tapply(final$pm25.pred.rf, final$mm,
       function(x) format(summary(x)))

final %>% 
  dplyr::group_by(mm) %>% 
  dplyr::summarize(q01 = quantile(pm25.pred.rf, 0.01),
                   q05 = quantile(pm25.pred.rf, 0.05),
                   q25 = quantile(pm25.pred.rf, 0.25),
                   q50 = quantile(pm25.pred.rf, 0.50),
                   q75 = quantile(pm25.pred.rf, 0.75),
                   q95 = quantile(pm25.pred.rf, 0.95),
                   q99 = quantile(pm25.pred.rf, 0.99))

saveRDS(final, "./stage3.2019...PM25.grid100m.RF_1_2.RDS")

keep(mod2, sure = T)
gc()

#### 2_1 ####

##### Call stage3 dataset (complete on grid cells)
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
  data$pred.rf.st3 <- predict(mod2,  data)$predictions
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

setnames(final, "pred.rf.st3", "ln.pm25.pred.rf")
final$pm25.pred.rf <- exp(final$ln.pm25.pred.rf)

head(final)
dim(final)

summary(final$pm25.pred.rf)

tapply(final$pm25.pred.rf, final$mm,
       function(x) format(summary(x)))

final %>% 
  dplyr::group_by(mm) %>% 
  dplyr::summarize(q01 = quantile(pm25.pred.rf, 0.01),
                   q05 = quantile(pm25.pred.rf, 0.05),
                   q25 = quantile(pm25.pred.rf, 0.25),
                   q50 = quantile(pm25.pred.rf, 0.50),
                   q75 = quantile(pm25.pred.rf, 0.75),
                   q95 = quantile(pm25.pred.rf, 0.95),
                   q99 = quantile(pm25.pred.rf, 0.99))

saveRDS(final, "./stage3.2019...PM25.grid100m.RF_2_1.RDS")

keep(mod2, sure = T)
gc()

#### 2_2 ####

##### Call stage3 dataset (complete on grid cells)
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
  data$pred.rf.st3 <- predict(mod2,  data)$predictions
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

setnames(final, "pred.rf.st3", "ln.pm25.pred.rf")
final$pm25.pred.rf <- exp(final$ln.pm25.pred.rf)

head(final)
dim(final)

summary(final$pm25.pred.rf)

tapply(final$pm25.pred.rf, final$mm,
       function(x) format(summary(x)))

final %>% 
  dplyr::group_by(mm) %>% 
  dplyr::summarize(q01 = quantile(pm25.pred.rf, 0.01),
                   q05 = quantile(pm25.pred.rf, 0.05),
                   q25 = quantile(pm25.pred.rf, 0.25),
                   q50 = quantile(pm25.pred.rf, 0.50),
                   q75 = quantile(pm25.pred.rf, 0.75),
                   q95 = quantile(pm25.pred.rf, 0.95),
                   q99 = quantile(pm25.pred.rf, 0.99))

saveRDS(final, "./stage3.2019...PM25.grid100m.RF_2_2.RDS")

keep(mod2, sure = T)
gc()

#### 3_1 ####

##### Call stage3 dataset (complete on grid cells)
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
  data$pred.rf.st3 <- predict(mod2,  data)$predictions
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

setnames(final, "pred.rf.st3", "ln.pm25.pred.rf")
final$pm25.pred.rf <- exp(final$ln.pm25.pred.rf)

head(final)
dim(final)

summary(final$pm25.pred.rf)

tapply(final$pm25.pred.rf, final$mm,
       function(x) format(summary(x)))

final %>% 
  dplyr::group_by(mm) %>% 
  dplyr::summarize(q01 = quantile(pm25.pred.rf, 0.01),
                   q05 = quantile(pm25.pred.rf, 0.05),
                   q25 = quantile(pm25.pred.rf, 0.25),
                   q50 = quantile(pm25.pred.rf, 0.50),
                   q75 = quantile(pm25.pred.rf, 0.75),
                   q95 = quantile(pm25.pred.rf, 0.95),
                   q99 = quantile(pm25.pred.rf, 0.99))

saveRDS(final, "./stage3.2019...PM25.grid100m.RF_3_1.RDS")

keep(mod2, sure = T)
gc()

#### 3_2 ####

##### Call stage3 dataset (complete on grid cells)
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
  data$pred.rf.st3 <- predict(mod2,  data)$predictions
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

setnames(final, "pred.rf.st3", "ln.pm25.pred.rf")
final$pm25.pred.rf <- exp(final$ln.pm25.pred.rf)

head(final)
dim(final)

summary(final$pm25.pred.rf)

tapply(final$pm25.pred.rf, final$mm,
       function(x) format(summary(x)))

final %>% 
  dplyr::group_by(mm) %>% 
  dplyr::summarize(q01 = quantile(pm25.pred.rf, 0.01),
                   q05 = quantile(pm25.pred.rf, 0.05),
                   q25 = quantile(pm25.pred.rf, 0.25),
                   q50 = quantile(pm25.pred.rf, 0.50),
                   q75 = quantile(pm25.pred.rf, 0.75),
                   q95 = quantile(pm25.pred.rf, 0.95),
                   q99 = quantile(pm25.pred.rf, 0.99))

saveRDS(final, "./stage3.2019...PM25.grid100m.RF_3_2.RDS")

keep(mod2, sure = T)
gc()

#### 4_1 ####

##### Call stage3 dataset (complete on grid cells)
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
  data$pred.rf.st3 <- predict(mod2,  data)$predictions
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

setnames(final, "pred.rf.st3", "ln.pm25.pred.rf")
final$pm25.pred.rf <- exp(final$ln.pm25.pred.rf)

head(final)
dim(final)

summary(final$pm25.pred.rf)

tapply(final$pm25.pred.rf, final$mm,
       function(x) format(summary(x)))

final %>% 
  dplyr::group_by(mm) %>% 
  dplyr::summarize(q01 = quantile(pm25.pred.rf, 0.01),
                   q05 = quantile(pm25.pred.rf, 0.05),
                   q25 = quantile(pm25.pred.rf, 0.25),
                   q50 = quantile(pm25.pred.rf, 0.50),
                   q75 = quantile(pm25.pred.rf, 0.75),
                   q95 = quantile(pm25.pred.rf, 0.95),
                   q99 = quantile(pm25.pred.rf, 0.99))

saveRDS(final, "./stage3.2019...PM25.grid100m.RF_4_1.RDS")

keep(mod2, sure = T)
gc()

#### 4_2 ####

##### Call stage3 dataset (complete on grid cells)
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
  data$pred.rf.st3 <- predict(mod2,  data)$predictions
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

setnames(final, "pred.rf.st3", "ln.pm25.pred.rf")
final$pm25.pred.rf <- exp(final$ln.pm25.pred.rf)

head(final)
dim(final)

summary(final$pm25.pred.rf)

tapply(final$pm25.pred.rf, final$mm,
       function(x) format(summary(x)))

final %>% 
  dplyr::group_by(mm) %>% 
  dplyr::summarize(q01 = quantile(pm25.pred.rf, 0.01),
                   q05 = quantile(pm25.pred.rf, 0.05),
                   q25 = quantile(pm25.pred.rf, 0.25),
                   q50 = quantile(pm25.pred.rf, 0.50),
                   q75 = quantile(pm25.pred.rf, 0.75),
                   q95 = quantile(pm25.pred.rf, 0.95),
                   q99 = quantile(pm25.pred.rf, 0.99))

saveRDS(final, "./stage3.2019...PM25.grid100m.RF_4_2.RDS")

keep(mod2, sure = T)
gc()

#### 5_1 ####

##### Call stage3 dataset (complete on grid cells)
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
  data$pred.rf.st3 <- predict(mod2,  data)$predictions
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

setnames(final, "pred.rf.st3", "ln.pm25.pred.rf")
final$pm25.pred.rf <- exp(final$ln.pm25.pred.rf)

head(final)
dim(final)

summary(final$pm25.pred.rf)

tapply(final$pm25.pred.rf, final$mm,
       function(x) format(summary(x)))

final %>% 
  dplyr::group_by(mm) %>% 
  dplyr::summarize(q01 = quantile(pm25.pred.rf, 0.01),
                   q05 = quantile(pm25.pred.rf, 0.05),
                   q25 = quantile(pm25.pred.rf, 0.25),
                   q50 = quantile(pm25.pred.rf, 0.50),
                   q75 = quantile(pm25.pred.rf, 0.75),
                   q95 = quantile(pm25.pred.rf, 0.95),
                   q99 = quantile(pm25.pred.rf, 0.99))

saveRDS(final, "./stage3.2019...PM25.grid100m.RF_5_1.RDS")

keep(mod2, sure = T)
gc()

#### 5_2 ####

##### Call stage3 dataset (complete on grid cells)
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
  data$pred.rf.st3 <- predict(mod2,  data)$predictions
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

setnames(final, "pred.rf.st3", "ln.pm25.pred.rf")
final$pm25.pred.rf <- exp(final$ln.pm25.pred.rf)

head(final)
dim(final)

summary(final$pm25.pred.rf)

tapply(final$pm25.pred.rf, final$mm,
       function(x) format(summary(x)))

final %>% 
  dplyr::group_by(mm) %>% 
  dplyr::summarize(q01 = quantile(pm25.pred.rf, 0.01),
                   q05 = quantile(pm25.pred.rf, 0.05),
                   q25 = quantile(pm25.pred.rf, 0.25),
                   q50 = quantile(pm25.pred.rf, 0.50),
                   q75 = quantile(pm25.pred.rf, 0.75),
                   q95 = quantile(pm25.pred.rf, 0.95),
                   q99 = quantile(pm25.pred.rf, 0.99))

saveRDS(final, "./stage3.2019...PM25.grid100m.RF_5_2.RDS")

keep(mod2, sure = T)
gc()

#### 6_1 ####

##### Call stage3 dataset (complete on grid cells)
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
  data$pred.rf.st3 <- predict(mod2,  data)$predictions
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

setnames(final, "pred.rf.st3", "ln.pm25.pred.rf")
final$pm25.pred.rf <- exp(final$ln.pm25.pred.rf)

head(final)
dim(final)

summary(final$pm25.pred.rf)

tapply(final$pm25.pred.rf, final$mm,
       function(x) format(summary(x)))

final %>% 
  dplyr::group_by(mm) %>% 
  dplyr::summarize(q01 = quantile(pm25.pred.rf, 0.01),
                   q05 = quantile(pm25.pred.rf, 0.05),
                   q25 = quantile(pm25.pred.rf, 0.25),
                   q50 = quantile(pm25.pred.rf, 0.50),
                   q75 = quantile(pm25.pred.rf, 0.75),
                   q95 = quantile(pm25.pred.rf, 0.95),
                   q99 = quantile(pm25.pred.rf, 0.99))

saveRDS(final, "./stage3.2019...PM25.grid100m.RF_6_1.RDS")

keep(mod2, sure = T)
gc()

#### 6_2 ####

##### Call stage3 dataset (complete on grid cells)
stage3 <- readRDS("./stage3.2019.PM.all.grid100m_6_2.RDS")
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
  data$pred.rf.st3 <- predict(mod2,  data)$predictions
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

setnames(final, "pred.rf.st3", "ln.pm25.pred.rf")
final$pm25.pred.rf <- exp(final$ln.pm25.pred.rf)

head(final)
dim(final)

summary(final$pm25.pred.rf)

tapply(final$pm25.pred.rf, final$mm,
       function(x) format(summary(x)))

final %>% 
  dplyr::group_by(mm) %>% 
  dplyr::summarize(q01 = quantile(pm25.pred.rf, 0.01),
                   q05 = quantile(pm25.pred.rf, 0.05),
                   q25 = quantile(pm25.pred.rf, 0.25),
                   q50 = quantile(pm25.pred.rf, 0.50),
                   q75 = quantile(pm25.pred.rf, 0.75),
                   q95 = quantile(pm25.pred.rf, 0.95),
                   q99 = quantile(pm25.pred.rf, 0.99))

saveRDS(final, "./stage3.2019...PM25.grid100m.RF_6_2.RDS")

keep(mod2, sure = T)
gc()