#### European-PM2.5-spatio-temporal-modelling

#### Author: Federica Nobile
#### Affiliation: Institute of Environmental Medicine, Karolinska Institutet, Stockholm, Sweden
####              Department of Epidemiology, Lazio Region Health Service/ASL Roma 1, Rome, Italy

library(sf)
library(dplyr)
library(tidyverse)
library(raster)
library(lubridate)
library(geos)
library(splitstackshape)
library(caret)
library(ranger)
library(data.table)
library(corrplot)

#### MONITORS ####

monitors      <- read_sf("./monitors_2019.shp")
monitors      <- unique(monitors)

setnames(monitors, c("Smp_P_I", "Smplngp", "CntryNm", "Country", "Ar_Ql_N", "A_Q_N_N", "A_Q_S_E", "A_Q_S_N_", "Ar_Q_S_N", "Altitud", "A_Q_S_A", "A_Q_S_T"), c("Sampling.Point.Id", "Samplingpoint", "CountryName", "Country", "Air.Quality.Network", "Air.Quality.Network.Name", "Air.Quality.Station.EoI.Code", "Air.Quality.Station.Nat.Code", "Air.Quality.Station.Name", "Altitude", "Air.Quality.Station.Area", "Air.Quality.Station.Type"))

length(unique(monitors$Air.Quality.Station.EoI.Code))
length(unique(monitors$geometry))
split <- split(monitors, monitors$Air.Quality.Station.EoI.Code)

monitors <- monitors %>% dplyr::select(!c(Sampling.Point.Id, Samplingpoint))
monitors <- unique(monitors)

# For the others' year experience
monitors <- monitors %>%
  arrange(year, CountryName, Country, Air.Quality.Station.EoI.Code,
          Air.Quality.Station.Area, desc(Air.Quality.Station.Type)) %>%
  group_by(year, CountryName, Country, Air.Quality.Station.EoI.Code) %>%
  slice(1) %>%
  ungroup()

table(monitors$CountryName)

coordinates <- st_coordinates(monitors)

# If you want to add these as new columns in your sf object
monitors$Longitude <- coordinates[, "X"]
monitors$Latitude  <- coordinates[, "Y"]

### Expand daily dataset
dates                  <- seq(as.Date("2019-01-01"), as.Date("2019-12-31"), by = "day")
expanded_monitors      <- monitors[rep(1:nrow(monitors), times = 365), ]

expanded_monitors$day  <- rep(dates, each = nrow(monitors))
head(expanded_monitors)

#### PM OBSERVATIONS ####

obs_pm25 <- readRDS("AIR POLL OBS/timeseries_2019_pm25.RDS")
head(obs_pm25)
obs_pm25 <- obs_pm25 %>% mutate(day = as.Date(as.numeric(date), "1970-01-01"))

head(expanded_monitors)
head(obs_pm25)

obs_pm25 <- obs_pm25 %>% dplyr::select(Air.Quality.Station.EoI.Code, day, pm25) 
obs_pm25 <- unique(obs_pm25) 
obs_pm25 <- obs_pm25 %>% group_by(Air.Quality.Station.EoI.Code, day) %>% summarise(pm25 = mean(pm25, na.rm = T)) 
summary(obs_pm25)

data <- merge(expanded_monitors, obs_pm25, by = c("Air.Quality.Station.EoI.Code", "day"), all.x = TRUE)


obs_pm10 <- readRDS("AIR POLL OBS/timeseries_2019_pm10.RDS")
obs_pm10 <- obs_pm10 %>% mutate(day = as.Date(as.numeric(date), "1970-01-01"))

obs_pm10 <- obs_pm10 %>% dplyr::select(Air.Quality.Station.EoI.Code, day, pm10)
obs_pm10 <- unique(obs_pm10) 
obs_pm10 <- obs_pm10 %>% group_by(Air.Quality.Station.EoI.Code, day) %>% summarise(pm10 = mean(pm10, na.rm = T))
summary(obs_pm10)

head(data)
head(obs_pm10)

data <- merge(data, obs_pm10, by = c("Air.Quality.Station.EoI.Code", "day"), all.x = TRUE)

head(data)
dim(data)

# I consider only the monitors with at least one of the 1 air pollutants
data2 <- data %>% filter(!is.na(pm25) | !is.na(pm10))
dim(data2)
length(unique(data2$Air.Quality.Station.EoI.Code))

#### SPATIAL PREDICTORS ####

sp_predictors <- read.csv("X:/NOBILE/PhD KI - Study II/2019/SPATIAL PREDICTORS/pred_pm25All.csv")
head(sp_predictors)
sp_predictors <- sp_predictors[,-c(1)]
sp_predictors <- unique(sp_predictors)

setnames(sp_predictors, c("Smp_P_I", "Smplngp", "CntryNm", "Country", "Ar_Ql_N", "A_Q_N_N", "A_Q_S_E", "A_Q_S_N_", "Ar_Q_S_N", "Altitud", "A_Q_S_A", "A_Q_S_T"), c("Sampling.Point.Id", "Samplingpoint", "CountryName", "Country", "Air.Quality.Network", "Air.Quality.Network.Name", "Air.Quality.Station.EoI.Code", "Air.Quality.Station.Nat.Code", "Air.Quality.Station.Name", "Altitude", "Air.Quality.Station.Area", "Air.Quality.Station.Type"))

sp_predictors <- sp_predictors %>% dplyr::select(Air.Quality.Station.EoI.Code, allRoads_100:zoneID)
sp_predictors <- unique(sp_predictors)
sp_predictors <- sp_predictors %>% dplyr::group_by(Air.Quality.Station.EoI.Code) %>% slice(1) %>% ungroup()

head(data2)
head(sp_predictors)

data3 <- merge(data2, sp_predictors, by = c("Air.Quality.Station.EoI.Code"), all.x = TRUE)
dim(data3)   
head(data3)


#### SPATIO-TEMPORAL PREDICTORS ####

#### Era5 Cams ####

era5_cams <- readRDS("ERA5/era5_2019_monitors.RDS")
head(data3)
head(era5_cams)

setnames(era5_cams, c("Smp_P_I", "Smplngp", "CntryNm", "Country", "Ar_Ql_N", "A_Q_N_N", "A_Q_S_E", "A_Q_S_N_", "Ar_Q_S_N", "Altitud", "A_Q_S_A", "A_Q_S_T"), c("Sampling.Point.Id", "Samplingpoint", "CountryName", "Country", "Air.Quality.Network", "Air.Quality.Network.Name", "Air.Quality.Station.EoI.Code", "Air.Quality.Station.Nat.Code", "Air.Quality.Station.Name", "Altitude", "Air.Quality.Station.Area", "Air.Quality.Station.Type"))

era5_cams <- era5_cams %>% dplyr::select(Air.Quality.Station.EoI.Code, day, u10:tco)
era5_cams <- unique(era5_cams)

data4 <- merge(data3, era5_cams, by = c("day", "Air.Quality.Station.EoI.Code"), all.x = TRUE)
dim(data4)  
head(data4)

#### Dispersion models Cams ####

pm25_cams <- readRDS("CAMS AIR POLL/pm25_cams_2019_monitors.RDS")
head(data4)
data4 <- data4 %>% dplyr::select(!c(geometry))
head(pm25_cams)
pm25_cams <- pm25_cams %>% ungroup() %>% st_drop_geometry() %>% dplyr::select(!c(geometry))
pm25_cams <- pm25_cams %>% dplyr::select(Air.Quality.Station.EoI.Code, day, pm25_cams)
pm25_cams <- unique(pm25_cams)

data5 <- merge(data4, pm25_cams, by = c("Air.Quality.Station.EoI.Code", "day"), all.x = TRUE)
dim(data5)   
head(data5)


pm10_cams <- readRDS("CAMS AIR POLL/pm10_cams_2019_monitors.RDS")
head(data5)
head(pm10_cams)
pm10_cams <- pm10_cams %>% ungroup() %>% st_drop_geometry() %>% dplyr::select(!c(geometry))
pm10_cams <- pm10_cams %>% dplyr::select(Air.Quality.Station.EoI.Code, day, pm10_cams)
pm10_cams <- unique(pm10_cams)

data6 <- merge(data5, pm10_cams, by = c("Air.Quality.Station.EoI.Code", "day"), all.x = TRUE)
dim(data6)
head(data6)

# cor(data6$pm25, data6$pm25_cams, use = "pairwise.complete.obs")
# cor(data6$pm10, data6$pm10_cams, use = "pairwise.complete.obs")

#### AOD Cams ####

aod_cams <- readRDS("CAMS AOD/aod_cams_2019_monitors.RDS")
head(data6)
head(aod_cams)

setnames(aod_cams, c("Smp_P_I", "Smplngp", "CntryNm", "Country", "Ar_Ql_N", "A_Q_N_N", "A_Q_S_E", "A_Q_S_N_", "Ar_Q_S_N", "Altitud", "A_Q_S_A", "A_Q_S_T"), c("Sampling.Point.Id", "Samplingpoint", "CountryName", "Country", "Air.Quality.Network", "Air.Quality.Network.Name", "Air.Quality.Station.EoI.Code", "Air.Quality.Station.Nat.Code", "Air.Quality.Station.Name", "Altitude", "Air.Quality.Station.Area", "Air.Quality.Station.Type"))

aod_cams <- aod_cams %>% dplyr::select(Air.Quality.Station.EoI.Code, day, aod1240_cams:aod865_cams)
aod_cams <- unique(aod_cams)

data7 <- merge(data6, aod_cams, by = c("day", "Air.Quality.Station.EoI.Code"), all.x = TRUE)
dim(data7)
head(data7)

#### AOD Modis ####

grid_aod <- read_sf("./AOD MODIS/AOD Raw/Europe_grid_1000_3035_centroids_permaiac.shp") # the 2019 dataset is my reference!
st_crs(monitors) == st_crs(grid_aod)

monitors                <- st_transform(monitors, st_crs(grid_aod))
nearest                 <- st_nearest_feature(monitors, grid_aod)
monitors$id_aod_1km     <- grid_aod$layer[nearest]
dist                    <- geos_distance(monitors, grid_aod[nearest,])
monitors$dist_to_points <- dist

names(monitors)
monitors <- monitors %>% dplyr::select(!c(geometry)) %>% st_drop_geometry()
monitors <- monitors %>% dplyr::select(Air.Quality.Station.EoI.Code, dist_to_points, id_aod_1km)
monitors <- unique(monitors)

data8 <- merge(data7, monitors, by = c("Air.Quality.Station.EoI.Code"), all.x = TRUE)
names(data8)

summary(monitors$dist_to_points)

months <- c ("01", "02", "03", "04", "05", "06",
             "07", "08", "09", "10", "11", "12")

data9 <- data.frame()

data8$mm <- month(data8$day)

for(i in 1:12){
  
  month  <- months[i]
  print(paste0(month, " ", Sys.time()))
  
  data1  <- data8 %>% dplyr::filter(mm == i)
  
  data2  <- readRDS(paste0("AOD MODIS/Final/daily_pred_aod_maiac_2019", month,".RDS"))
  data2  <- data2 %>% dplyr::filter(layer %in% data1$id_aod_1km) %>% dplyr::select(layer, day, final.aod047)
  
  data1  <- merge(data1, data2, by.x = c("day", "id_aod_1km"), by.y = c("day", "layer"), all.x = T)
  
  data9  <- rbind(data9, data1)
  
  rm(data1, data2)
  gc()
  
  print(paste0(month, " ", Sys.time(), " completed!"))
  
}


#### Ndvi Modis ####

ndvi <- readRDS("NDVI/daily_ndvi_evi_monitors_2019.RDS")
ndvi <- ndvi %>% ungroup() %>% st_drop_geometry() %>% dplyr::select(!c(geometry))
head(ndvi)

setnames(ndvi, c("Smp_P_I", "Smplngp", "CntryNm", "Country", "Ar_Ql_N", "A_Q_N_N", "A_Q_S_E", "A_Q_S_N_", "Ar_Q_S_N", "Altitud", "A_Q_S_A", "A_Q_S_T"), c("Sampling.Point.Id", "Samplingpoint", "CountryName", "Country", "Air.Quality.Network", "Air.Quality.Network.Name", "Air.Quality.Station.EoI.Code", "Air.Quality.Station.Nat.Code", "Air.Quality.Station.Name", "Altitude", "Air.Quality.Station.Area", "Air.Quality.Station.Type"))

ndvi <- ndvi %>% dplyr::select(Air.Quality.Station.EoI.Code, mm, EVI, NDVI)
ndvi <- unique(ndvi)

head(data9)
data9 <- data9 %>% mutate(mm = month(day))

data10 <- merge(data9, ndvi, by = c("mm", "Air.Quality.Station.EoI.Code"), all.x = TRUE)
dim(data10) 
head(data10)

# data10 <- data10 %>% mutate(EVI  = EVI * 0.0001,
#                           NDVI = NDVI * 0.0001) # no needed now

summary(data10$EVI)
summary(data10$NDVI)



summary(data10[,1:120])
summary(data10[,121:224])
saveRDS(data10, "./stage1.2019.PM.all.RDS")


#### Check that I'm considering the monitors of interest and whether there is a need of filling some missing values ####

names(data10)
table(data10$CountryName)

data12 <- data10 %>%
          dplyr::filter(CountryName != "IS") %>% # only Continental Europe
          dplyr::filter(CountryName != "CY") %>% # we don't have the spatial predictors
          rename(year = year.x) %>%
          dplyr::select(!c("year.y"))
dim(data12)
names(data12)
head(data12)

europe <- read_sf("./Europe.shp")
plot(st_geometry(europe))

monitors_def <- data12 %>% 
                dplyr::select(CountryName, Country, Air.Quality.Network, Air.Quality.Network.Name, Air.Quality.Station.EoI.Code, Air.Quality.Station.Nat.Code, Air.Quality.Station.Name, Altitude, Air.Quality.Station.Area, Air.Quality.Station.Type, Latitude, Longitude, geometry)
monitors_def <- unique(monitors_def)
plot(st_geometry(monitors_def), add = TRUE)

# I exclude the monitors related to far islands (this is specific for each year)

xmin <- -10
xmax <- 40
ymin <- 35
ymax <- 72

monitors_def2 <- monitors_def[
  st_coordinates(monitors_def)[,1] >= xmin &
    st_coordinates(monitors_def)[,1] <= xmax &
    st_coordinates(monitors_def)[,2] >= ymin &
    st_coordinates(monitors_def)[,2] <= ymax,
]

plot(st_geometry(monitors_def2))
plot(st_geometry(europe), add = T)
plot(st_geometry(monitors_def), col = "red", add = T)


data13 <- data12 %>%
  dplyr::filter(Air.Quality.Station.EoI.Code %in% monitors_def2$Air.Quality.Station.EoI.Code)
dim(data13)

plot(st_geometry(europe))

monitors_def <- data13 %>% 
                  dplyr::select(CountryName, Country, Air.Quality.Network, Air.Quality.Network.Name, Air.Quality.Station.EoI.Code, Air.Quality.Station.Nat.Code, Air.Quality.Station.Name, Altitude, Air.Quality.Station.Area, Air.Quality.Station.Type, Latitude, Longitude, geometry)
monitors_def <- unique(monitors_def)
plot(st_geometry(monitors_def))
plot(st_geometry(europe), add = T)

summary(data13[,1:110])
summary(data13[,111:223])
missing_values <- sapply(data13, function(x) sum(is.na(x)))
print(missing_values)

vedi1 <- data13 %>% filter(is.na(ind_200)) # borders -> negligible
plot(st_geometry(vedi1), col = "blue")

vedi2 <- data13 %>% filter(is.na(final.aod047)) # none

vedi3 <- data13 %>% filter(is.na(NDVI)) # 0 because solved in the specific script
plot(st_geometry(vedi3), col = "green", add = TRUE)

vedi4 <- data13 %>% filter(is.na(temp_2m)) # borders but the percentage of missing is not negligible!
plot(st_geometry(vedi4), col = "red", add = TRUE) 

### FILLING MISSING VALUES OF 9KM TEMP USING THE NEAREST POINT

data.a <- data13 %>% dplyr::filter(!is.na(temp_2m))

data.b <- data13 %>% dplyr::filter(is.na(temp_2m))
names(data.b)

cams <- read_sf("ERA5/centroids_era5_9km.shp") 
st_crs(monitors_def) == st_crs(cams)

nearest                 <- st_nearest_feature(monitors_def, cams)
monitors_def$id_era5_9  <- cams$id_era5_9[nearest]
dist                    <- geos_distance(monitors_def, cams[nearest,])
monitors_def$dist       <- dist

names(monitors_def)
head(monitors_def)
monitors_def <- monitors_def %>% dplyr::select(c(Air.Quality.Station.EoI.Code, id_era5_9, geometry)) %>% st_drop_geometry()

data.b2 <- merge(data.b, monitors_def, by = c("Air.Quality.Station.EoI.Code"), all.x = TRUE)
names(data.b2)
data.b2 <- data.b2 %>% dplyr::select(!c(u10, v10, dewp_2m, temp_2m, leaf_high, leaf_low))

tm_cams <- readRDS("ERA5/daily_era5_9km_2019.RDS")
head(tm_cams)
tm_cams <- tm_cams %>% dplyr:: filter(id_era5_9 %in% data.b2$id_era5_9)
data.b3 <- merge(data.b2, tm_cams, by.x = c("id_era5_9", "day", "year"), by.y = c("id_era5_9", "day", "year"), all.x = TRUE)
names(data.b3)
head(data.b3)
data.b3 <- data.b3 %>%
              dplyr::select(!c(id_era5_9, x, y))

summary(data.b3[,111:223])

names(data.a)
names(data.b3)
data14 <- rbind(data.a, data.b3)


#### Predicting PM2.5 where we have PM10 observed values ####

names(data14)
table(data14$Air.Quality.Station.Type)
table(data14$Air.Quality.Station.Area)
data14$dow          <- wday(data14$day)
data14$mm           <- as.factor(data14$mm)
data14$dow          <- as.factor(data14$dow)
data14$tipo_zona    <- as.factor(data14$Air.Quality.Station.Type)

summary(data14$pm25)
summary(data14$pm10)
vedi <- data14 %>% filter(is.na(pm25) & is.na(pm10)) 

data14 <- data14 %>% mutate(pm25 = ifelse(pm25 < 0, 0, pm25),
                            pm10 = ifelse(pm10 < 0, 0, pm10))

data <- data14 %>% dplyr::filter(!is.na(pm25) & !is.na(pm10)) 
summary(data$pm25)

summary(data$pm10)


# Set hyperparameters [.........]

### I use hyperparameters from 2019 (as explained in the Supplementary Materials) ###

formula  <- as.formula(pm25 ~ pm10 + tipo_zona + mm + dow + Longitude + Latitude)

data <- as.data.frame(data)
set.seed(1)
mod <- ranger(formula, data=data, num.trees=1000, mtry=5, importance="impurity")
cbind(mod$variable.importance)

mod$prediction.error 
mod$r.squared 

data14           <- as.data.frame(data14)
data_nopm10      <- subset(data14, is.na(pm10))
data_sipm10      <- subset(data14, !is.na(pm10))
data_1           <- subset(data_sipm10, mm == 1 | mm == 2 | mm == 3 | mm == 4)
data_1$pred.pm25 <- predict(mod, data_1)$predictions

data_2           <- subset(data_sipm10, mm == 5 | mm == 6 | mm == 7 | mm == 8)
data_2$pred.pm25 <- predict(mod, data_2)$predictions

data_3           <- subset(data_sipm10, mm == 9 | mm == 10 | mm == 11 | mm == 12)
data_3$pred.pm25 <- predict(mod, data_3)$predictions

data_nopm10$pred.pm25 <- NA

data15 <- rbind(data_1, data_2, data_3, data_nopm10)
dim(data15)

data15 <- data15 %>%
  dplyr::mutate(final.pm25   = pred.pm25,
                final.pm25   = ifelse(!is.na(pm25), pm25, final.pm25))
summary(data15$pm25)
summary(data15$pm10)
summary(data15$pred.pm25)
summary(data15$final.pm25)


#### Final checks ####

missing_values <- sapply(data15, function(x) sum(is.na(x)))
print(missing_values)

data15 <- data15 %>% dplyr::filter(!is.na(imd_8000))
missing_values <- sapply(data15, function(x) sum(is.na(x)))
print(missing_values)

data15 <- data15 %>% dplyr::filter(!is.na(EVI))
missing_values <- sapply(data15, function(x) sum(is.na(x)))
print(missing_values)

summary(data15)
data15 <- data15 %>% dplyr::mutate(dewp_2m = dewp_2m - 273.15,
                                   temp_2m = temp_2m - 273.15,
                                   mslp = mslp/100)

names(data15)
data15 <- data15 %>% dplyr::select(!c(id_aod_1km, dist_to_points, geometry)) %>% st_drop_geometry
summary(data15[,1:130])
summary(data15[,131:224])

length(unique(data15$Air.Quality.Station.EoI.Code)) 

saveRDS(data15, "./stage1.PM.final.RDS")
