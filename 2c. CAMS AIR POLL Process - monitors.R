#### European-PM2.5-spatio-temporal-modelling

#### Author: Federica Nobile
#### Affiliation: Institute of Environmental Medicine, Karolinska Institutet, Stockholm, Sweden
####              Department of Epidemiology, Lazio Region Health Service/ASL Roma 1, Rome, Italy

# Load necessary libraries
library(raster)
library(sf)
library(sp)
library(exactextractr)
# library(ggplot2)
# library(rnaturalearth)
# library(rnaturalearthdata)
library(tidyverse)
library(tidyr)
library(ncdf4)
# library(dplyr)
library(data.table)

# 1. Read monitors shapefile
monitors      <- read_sf("./monitors_2019.shp")
monitors      <- unique(monitors) 
length(unique(monitors$A_Q_S_E)) 

#### PM2.5 ####

# 2. Read the AOD maps (nc and then raster)
### Open file and get its variables
cams_raster <- raster::brick("RAW DATA PM25/particulate_matter_25um_2019_01/cams.eaq.vra.ENSa.pm2p5.l0.2019-01.nc")
plot(cams_raster$X2019.01.01.00.00.00)

# 3. Check and transform CRS if needed
monitors_crs  <- st_crs(monitors)
monitors_repr <- st_transform(monitors, crs = st_crs(cams_raster$X2019.01.01.00.00.00))
month         <- as.data.frame(monitors_repr)

plot(st_geometry(monitors_repr), add = TRUE)

rm(cams_raster)


final <- list()
months <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")

# Start with a loop by month 

for(i in 1:11){
  
  print(paste0(months[i], " ", Sys.time()))
  
  # Read the AOD maps (nc and then raster) - aod nc already read
  cams_raster <- raster::brick(paste0("RAW DATA PM25/particulate_matter_25um_2019_", months[i], "/cams.eaq.vra.ENSa.pm2p5.l0.2019-", months[i], ".nc"))
  
  # Perform spatial join and assign temperature values to monitor points (extract to point)
  cams_data <- raster::extract(cams_raster, monitors_repr)
  month     <- cbind(month, cams_data)
  
  # From wide to long format
  month1 <- month %>% tidyr::pivot_longer(!c(Smp_P_I, Smplngp, CntryNm, Country, Ar_Ql_N, A_Q_N_N, A_Q_S_E, A_Q_S_N_, Ar_Q_S_N, Altitud, A_Q_S_A, A_Q_S_T, year, geometry), names_to = "date", values_to = "pm25_cams")
  
  day          <- rep(seq(as.Date(paste0("2019-", months[i], "-01")), as.Date(paste0("2019-", months[i+1], "-01")), by = "day"), each = 24)
  day          <- day[c(1:(length(day)-24))]
  layers       <- month1 %>% dplyr::select(date)
  layers       <- unique(layers)
  layers$layer <- substr(layers$date, 2, 11)
  layers       <- layers[order(layers$layer),]
  dates        <- cbind(as.data.frame(day), as.data.frame(layers))
  
  print(paste0(months[i], " ", Sys.time(), " completed!"))
  
  month2       <- merge(month1, dates[,c(1:2)], by = c("date"), all.x = TRUE)

  print(paste0(months[i], " ", Sys.time(), " join completed!"))
  
  final[[i]] <- month2
  rm(month, month1, month2, day, layers, dates)
  month      <- as.data.frame(monitors_repr)
  
  gc()
  
}

for(i in 12){
  
  print(paste0(months[i], " ", Sys.time()))
  
  # Read the AOD maps (nc and then raster) - aod nc already read
  cams_raster <- raster::brick(paste0("RAW DATA PM25/particulate_matter_25um_2019_", months[i], "/cams.eaq.vra.ENSa.pm2p5.l0.2019-", months[i], ".nc"))
  
  # Perform spatial join and assign temperature values to monitor points (extract to point)
  cams_data <- raster::extract(cams_raster, monitors_repr)
  month     <- cbind(month, cams_data)
  
  # From wide to long format
  month1 <- month %>% tidyr::pivot_longer(!c(Smp_P_I, Smplngp, CntryNm, Country, Ar_Ql_N, A_Q_N_N, A_Q_S_E, A_Q_S_N_, Ar_Q_S_N, Altitud, A_Q_S_A, A_Q_S_T, year, geometry), names_to = "date", values_to = "pm25_cams")
  
  day          <- rep(seq(as.Date(paste0("2019-12-01")), as.Date(paste0("2019-01-01")), by = "day"), each = 24)
  day          <- day[c(1:(length(day)-24))]
  layers       <- month1 %>% dplyr::select(date)
  layers       <- unique(layers)
  layers$layer <- substr(layers$date, 2, 11)
  layers       <- layers[order(layers$layer),]
  dates        <- cbind(as.data.frame(day), as.data.frame(layers))
  
  print(paste0(months[i], " ", Sys.time(), " completed!"))
  
  month2       <- merge(month1, dates[,c(1:2)], by = c("date"), all.x = TRUE)
  
  print(paste0(months[i], " ", Sys.time(), " join completed!"))
  
  final[[i]] <- month2
  rm(month, month1, month2, day, layers, dates)
  month      <- as.data.frame(monitors_repr)
  
  gc()
  
}

final2 <-  bind_rows(final)
head(final2)
setnames(final2, c("Smp_P_I", "Smplngp", "CntryNm", "Country", "Ar_Ql_N", "A_Q_N_N", "A_Q_S_E", "A_Q_S_N_", "Ar_Q_S_N", "Altitud", "A_Q_S_A", "A_Q_S_T"), c("Sampling.Point.Id", "Samplingpoint", "CountryName", "Country", "Air.Quality.Network", "Air.Quality.Network.Name", "Air.Quality.Station.EoI.Code", "Air.Quality.Station.Nat.Code", "Air.Quality.Station.Name", "Altitude", "Air.Quality.Station.Area", "Air.Quality.Station.Type"))
head(final2)

final3 <- final2 %>% 
  group_by(Sampling.Point.Id, Samplingpoint, CountryName, Country, Air.Quality.Network, Air.Quality.Network.Name, Air.Quality.Station.EoI.Code, Air.Quality.Station.Nat.Code, Air.Quality.Station.Name, Altitude, Air.Quality.Station.Area, Air.Quality.Station.Type, geometry, day) %>%
                      summarise(pm25_cams = mean(pm25_cams, na.rm = TRUE))
dim(final3) # 2011515      15 (5511 monitors * 365 days)

# Save
saveRDS(final3, "pm25_cams_2019_monitors.RDS")


#### PM10 ####

# 2. Read the AOD maps (nc and then raster)
### Open file and get its variables
cams_raster <- raster::brick("RAW DATA PM10/particulate_matter_10um_2019_01/cams.eaq.vra.ENSa.pm10.l0.2019-01.nc")

plot(cams_raster$X2019.01.01.00.00.00)
plot(cams_raster$X2019.01.31.23.00.00, add= TRUE)

# 3. Check and transform CRS if needed
monitors_crs  <- st_crs(monitors)
monitors_repr <- st_transform(monitors, crs = st_crs(cams_raster$X2019.01.01.00.00.00))
month         <- as.data.frame(monitors_repr)

plot(st_geometry(monitors_repr), add = TRUE)

rm(cams_raster)


final <- list()
months <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")

# Start with a loop by month 

for(i in 1:11){
  
  print(paste0(months[i], " ", Sys.time()))
  
  # Read the AOD maps (nc and then raster) - aod nc already read
  cams_raster <- raster::brick(paste0("RAW DATA PM10/particulate_matter_10um_2019_", months[i], "/cams.eaq.vra.ENSa.pm10.l0.2019-", months[i], ".nc"))
  
  # Perform spatial join and assign temperature values to monitor points (extract to point)
  cams_data <- raster::extract(cams_raster, monitors_repr)
  month     <- cbind(month, cams_data)
  
  # From wide to long format
  month1 <- month %>% tidyr::pivot_longer(!c(Smp_P_I, Smplngp, CntryNm, Country, Ar_Ql_N, A_Q_N_N, A_Q_S_E, A_Q_S_N_, Ar_Q_S_N, Altitud, A_Q_S_A, A_Q_S_T, year, geometry), names_to = "date", values_to = "pm10_cams")
  
  day          <- rep(seq(as.Date(paste0("2019-", months[i], "-01")), as.Date(paste0("2019-", months[i+1], "-01")), by = "day"), each = 24)
  day          <- day[c(1:(length(day)-24))]
  layers       <- month1 %>% dplyr::select(date)
  layers       <- unique(layers)
  layers$layer <- substr(layers$date, 2, 11)
  layers       <- layers[order(layers$layer),]
  dates        <- cbind(as.data.frame(day), as.data.frame(layers))
  
  print(paste0(months[i], " ", Sys.time(), " completed!"))
  
  month2       <- merge(month1, dates[,c(1:2)], by = c("date"), all.x = TRUE)
  
  print(paste0(months[i], " ", Sys.time(), " join completed!"))
  
  final[[i]] <- month2
  rm(month, month1, month2, day, layers, dates)
  month      <- as.data.frame(monitors_repr)
  
  gc()
  
}

for(i in 12){
  
  print(paste0(months[i], " ", Sys.time()))
  
  # Read the AOD maps (nc and then raster) - aod nc already read
  cams_raster <- raster::brick(paste0("RAW DATA PM10/particulate_matter_10um_2019_", months[i], "/cams.eaq.vra.ENSa.pm10.l0.2019-", months[i], ".nc"))
  
  # Perform spatial join and assign temperature values to monitor points (extract to point)
  cams_data <- raster::extract(cams_raster, monitors_repr)
  month     <- cbind(month, cams_data)
  
  # From wide to long format
  month1 <- month %>% tidyr::pivot_longer(!c(Smp_P_I, Smplngp, CntryNm, Country, Ar_Ql_N, A_Q_N_N, A_Q_S_E, A_Q_S_N_, Ar_Q_S_N, Altitud, A_Q_S_A, A_Q_S_T, year, geometry), names_to = "date", values_to = "pm10_cams")
  
  day          <- rep(seq(as.Date(paste0("2019-12-01")), as.Date(paste0("2019-01-01")), by = "day"), each = 24)
  day          <- day[c(1:(length(day)-24))]
  layers       <- month1 %>% dplyr::select(date)
  layers       <- unique(layers)
  layers$layer <- substr(layers$date, 2, 11)
  layers       <- layers[order(layers$layer),]
  dates        <- cbind(as.data.frame(day), as.data.frame(layers))
  
  print(paste0(months[i], " ", Sys.time(), " completed!"))
  
  month2       <- merge(month1, dates[,c(1:2)], by = c("date"), all.x = TRUE)
  
  print(paste0(months[i], " ", Sys.time(), " join completed!"))
  
  final[[i]] <- month2
  rm(month, month1, month2, day, layers, dates)
  month      <- as.data.frame(monitors_repr)
  
  gc()
  
}

final2 <-  bind_rows(final)
head(final2)
setnames(final2, c("Smp_P_I", "Smplngp", "CntryNm", "Country", "Ar_Ql_N", "A_Q_N_N", "A_Q_S_E", "A_Q_S_N_", "Ar_Q_S_N", "Altitud", "A_Q_S_A", "A_Q_S_T"), c("Sampling.Point.Id", "Samplingpoint", "CountryName", "Country", "Air.Quality.Network", "Air.Quality.Network.Name", "Air.Quality.Station.EoI.Code", "Air.Quality.Station.Nat.Code", "Air.Quality.Station.Name", "Altitude", "Air.Quality.Station.Area", "Air.Quality.Station.Type"))
head(final2)

final3 <- final2 %>% 
  group_by(Sampling.Point.Id, Samplingpoint, CountryName, Country, Air.Quality.Network, Air.Quality.Network.Name, Air.Quality.Station.EoI.Code, Air.Quality.Station.Nat.Code, Air.Quality.Station.Name, Altitude, Air.Quality.Station.Area, Air.Quality.Station.Type, geometry, day) %>%
  summarise(pm10_cams = mean(pm10_cams, na.rm = TRUE))
dim(final3) 

# Save
saveRDS(final3, "pm10_cams_2019_monitors.RDS")
