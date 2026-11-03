#### European-PM2.5-spatio-temporal-modelling

#### Author: Federica Nobile
#### Affiliation: Institute of Environmental Medicine, Karolinska Institutet, Stockholm, Sweden
####              Department of Epidemiology, Lazio Region Health Service/ASL Roma 1, Rome, Italy


# Load necessary libraries
# library(raster)
# library(sf)
library(tidyr)
library(lubridate)
library(dplyr)
library(data.table)
library(caret)
library(ranger)
library(splitstackshape)
library(writexl)

## 2.5 I create a loop by month (I save a monthly dataset with predicted 047 and 055 values after each iteration)

months <- c ("01", "02", "03", "04", "05", "06",
             "07", "08", "09", "10", "11", "12")

# countries <- c("Norway", "Sweden", "Finland", "United Kingdom",
#                "Estonia", "Latvia", "Denmark", "Lithuania",
#                "Ireland", "Poland", "Germany", "Netherlands",
#                "Belgium", "France", "Czechia", "Luxembourg",
#                "Slovakia", "Romania", "Hungary", "Austria",
#                "Switzerland", "Liechtenstein", "Slovenia", "Italy",
#                "Croatia", "Serbia", "Bulgaria", "Bosnia and Herz.",
#                "Spain", "Montenegro", "Portugal", "Monaco",
#                "Kosovo", "Greece", "Macedonia", "Andorra",
#                "Albania", "Malta")

r2.047         <- list()
prederr.047    <- list()
importance.047 <- list()


formula047  <- as.formula(aod047 ~ aod550_cams + aod469_cams + aod670_cams + aod865_cams + aod1240_cams + 
                         jd + x_3035 + y_3035)

points_eu_coord <- readRDS("./Europe_grid_1000_3035_centroids_permaiac_country_coord.RDS")
setDT(points_eu_coord)
setkey(points_eu_coord, layer)

for (i in 1:12){

  month <- months[i]
  print(paste0(month, " ", Sys.time()))
  
  data    <- readRDS(paste0("daily_aod047_cams_maiac_1km_2019", month, ".RDS"))
  data    <- data %>%  
    dplyr::mutate(aod047 = aod047*0.001)
  
  setDT(data)
  setkey(data, layer)
  data    <- merge(data, points_eu_coord, all.x = T)
  data$jd <- yday(data$day)
  data$mm <- month(data$day)
  gc()
  
  
  #### Training dataset for 047 model
  
  data_complete047 <- data %>% dplyr::filter(!is.na(aod047)) %>%
    dplyr::select(layer, day, aod1240_cams, aod469_cams, aod550_cams, aod670_cams, aod865_cams, aod047, NAME, x_3035, y_3035, jd, mm)  
  setDT(data_complete047)
  setkey(data_complete047, layer)
  print(paste0("047 ",dim(data_complete047)[1]))
  print(head(data_complete047))

  gc()
  
  
  #### Training - 047 model
  
  set.seed(1)
  subdata047          <- stratified(data_complete047, c("mm", "NAME"), 0.01)
  print(dim(subdata047))
  
  print("Running 047 model....")
  set.seed(1)
  mod047              <- ranger(formula047, data=subdata047, num.trees=1000, mtry=6, importance="impurity")
  print(cbind(mod047$variable.importance))
  importance.047[[i]] <- cbind(mod047$variable.importance)
  print(mod047$prediction.error)
  prederr.047[[i]]    <- mod047$prediction.error
  print(mod047$r.squared)
  r2.047[[i]]         <- mod047$r.squared
  
  gc()
  
  #### Saving to avoid to lost everything

  writexl::write_xlsx(as.data.frame(importance.047[[i]]), paste0("Final/importance_2019", month, ".xlsx"))
  writexl::write_xlsx(as.data.frame(prederr.047[[i]]), paste0("Final/prederr_2019", month, ".xlsx"))
  writexl::write_xlsx(as.data.frame(r2.047[[i]]), paste0("Final/r2_2019", month, ".xlsx"))
  
  gc()
  
  
  #### Prediction - 047
  
  pred047 <- data.frame()
  
  countries <- unique(data$NAME)
    
  for(j in countries){

    cc             <- subset(data, NAME == j)
    cc$dd          <- day(cc$day)
    print(paste0(j, " ", dim(cc)[1]))
    
    cc1            <- cc %>% dplyr::filter(dd < 15)
    cc2            <- cc %>% dplyr::filter(dd >= 15)
    rm(cc)
    gc()
    
    print("Predicting 047 - 1....")
    cc1$pred.aod047 <- predict(mod047, cc1)$predictions
    gc()

    print("Predicting 047 - 2....")
    cc2$pred.aod047 <- predict(mod047, cc2)$predictions
    gc()

    cc <- rbind(cc1, cc2)
    rm(cc1, cc2)
    gc()
    
    cc$final.aod047 <- cc$pred.aod047
    cc <- cc %>%
            dplyr::mutate(final.aod047 = ifelse(!is.na(aod047), aod047, final.aod047))
                    
    print(summary(cc))
    cc             <- cc %>% dplyr::select(!(dd))
    pred047    <- rbind(pred047, cc)
    gc()
    
  }
  
  saveRDS(pred047, paste0("./daily_pred_aod_maiac_2019", month, ".RDS"))
  
  print(paste0(month, " ", Sys.time(), " completed!"))
  rm(cc, data, data_complete047,
     subdata047,
     mod047,
     pred047)
  gc()
  
}
