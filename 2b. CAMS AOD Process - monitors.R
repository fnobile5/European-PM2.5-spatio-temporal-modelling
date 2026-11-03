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
library(tidyr)
library(lubridate)
library(dplyr)
library(ncdf4)
library(data.table)

variables <- c("1240", 
               "469",
               "550",
               "670",
               "865")                

# 1. Read monitors shapefile
monitors      <- read_sf("./monitors_2019.shp")
monitors      <- unique(monitors) 

# 2. Read the AOD maps (nc and then raster)
aod_raster   <- raster::brick("total_aerosol_optical_depth_1240nm_2019.nc")

plot(aod_raster$X1514764800)
plot(aod_raster$X1536354000, add= TRUE)

# 3. Check and transform CRS if needed
monitors_crs  <- st_crs(monitors)
aod_crs       <- crs(aod_raster)
monitors_repr <- st_transform(monitors, crs = st_crs(aod_raster$X1514764800))
month         <- as.data.frame(monitors_repr)

plot(st_geometry(monitors_repr), add = TRUE)

rm(i, var_name, aod_raster, lon, lat, var_data)


final <- list()


for(i in 1:5){
  
  print(paste0(variables[i], " ", Sys.time()))
  
  # Read the AOD maps (nc and then raster) 
  aod_raster   <- raster::brick(paste0("total_aerosol_optical_depth_", variables[i], "nm_2019.nc"))
  
  # Perform spatial join and assign temperature values to monitor points (extract to point)
  aod_data <- raster::extract(aod_raster, monitors_repr)
  month    <- cbind(month, aod_data)
      
  print(paste0(variables[i], "...", Sys.time(), " completed!"))
  
  final[[i]] <- month
  rm(month)
  month      <- as.data.frame(monitors_repr)
  
  gc()
  
}

# From wide to long format

variables
names(final[[1]])
colnames(final[[1]])[15:2934] <- as.character(rep(seq(as.Date("2019-01-01"), as.Date("2019-12-31"), by = "days"), each = 8))
names(final[[1]]) <- make.unique(names(final[[1]]))
data1 <- final[[1]] %>% pivot_longer(!c(Smp_P_I, Smplngp, CntryNm, Country, Ar_Ql_N, A_Q_N_N, A_Q_S_E, A_Q_S_N_, Ar_Q_S_N, Altitud, A_Q_S_A, A_Q_S_T, year, geometry), names_to = "date", values_to = "aod1240_cams") %>%
  dplyr::mutate(day = as.Date(date)) %>%
  st_drop_geometry() %>%
  dplyr::group_by(Smp_P_I, Smplngp, CntryNm, Country, Ar_Ql_N, A_Q_N_N, A_Q_S_E, A_Q_S_N_, Ar_Q_S_N, Altitud, A_Q_S_A, A_Q_S_T, day) %>%
  dplyr::summarise(aod1240_cams = mean(aod1240_cams, na.rn = T)) 

names(final[[2]])
colnames(final[[2]])[15:2934] <- as.character(rep(seq(as.Date("2019-01-01"), as.Date("2019-12-31"), by = "days"), each = 8))
names(final[[2]]) <- make.unique(names(final[[2]]))
data2 <- final[[2]] %>% pivot_longer(!c(Smp_P_I, Smplngp, CntryNm, Country, Ar_Ql_N, A_Q_N_N, A_Q_S_E, A_Q_S_N_, Ar_Q_S_N, Altitud, A_Q_S_A, A_Q_S_T, year, geometry), names_to = "date", values_to = "aod469_cams")  %>%
  dplyr::mutate(day = as.Date(date)) %>%
  st_drop_geometry() %>%
  dplyr::group_by(Smp_P_I, Smplngp, CntryNm, Country, Ar_Ql_N, A_Q_N_N, A_Q_S_E, A_Q_S_N_, Ar_Q_S_N, Altitud, A_Q_S_A, A_Q_S_T, day) %>%
  dplyr::summarise(aod469_cams = mean(aod469_cams, na.rn = T)) 

names(final[[3]])
colnames(final[[3]])[15:2934] <- as.character(rep(seq(as.Date("2019-01-01"), as.Date("2019-12-31"), by = "days"), each = 8))
names(final[[3]]) <- make.unique(names(final[[3]]))
data3 <- final[[3]] %>% pivot_longer(!c(Smp_P_I, Smplngp, CntryNm, Country, Ar_Ql_N, A_Q_N_N, A_Q_S_E, A_Q_S_N_, Ar_Q_S_N, Altitud, A_Q_S_A, A_Q_S_T, year, geometry), names_to = "date", values_to = "aod550_cams")  %>%
  dplyr::mutate(day = as.Date(date)) %>%
  st_drop_geometry() %>%
  dplyr::group_by(Smp_P_I, Smplngp, CntryNm, Country, Ar_Ql_N, A_Q_N_N, A_Q_S_E, A_Q_S_N_, Ar_Q_S_N, Altitud, A_Q_S_A, A_Q_S_T, day) %>%
  dplyr::summarise(aod550_cams = mean(aod550_cams, na.rn = T)) 

names(final[[4]])
colnames(final[[4]])[15:2934] <- as.character(rep(seq(as.Date("2019-01-01"), as.Date("2019-12-31"), by = "days"), each = 8))
names(final[[4]]) <- make.unique(names(final[[4]]))
data4 <- final[[4]] %>% pivot_longer(!c(Smp_P_I, Smplngp, CntryNm, Country, Ar_Ql_N, A_Q_N_N, A_Q_S_E, A_Q_S_N_, Ar_Q_S_N, Altitud, A_Q_S_A, A_Q_S_T, year, geometry), names_to = "date", values_to = "aod670_cams")  %>%
  dplyr::mutate(day = as.Date(date)) %>%
  st_drop_geometry() %>%
  dplyr::group_by(Smp_P_I, Smplngp, CntryNm, Country, Ar_Ql_N, A_Q_N_N, A_Q_S_E, A_Q_S_N_, Ar_Q_S_N, Altitud, A_Q_S_A, A_Q_S_T, day) %>%
  dplyr::summarise(aod670_cams = mean(aod670_cams, na.rn = T))   

names(final[[5]])
colnames(final[[5]])[15:2934] <- as.character(rep(seq(as.Date("2019-01-01"), as.Date("2019-12-31"), by = "days"), each = 8))
names(final[[5]]) <- make.unique(names(final[[5]]))
data5 <- final[[5]] %>% pivot_longer(!c(Smp_P_I, Smplngp, CntryNm, Country, Ar_Ql_N, A_Q_N_N, A_Q_S_E, A_Q_S_N_, Ar_Q_S_N, Altitud, A_Q_S_A, A_Q_S_T, year, geometry), names_to = "date", values_to = "aod865_cams")  %>%
  dplyr::mutate(day = as.Date(date)) %>%
  st_drop_geometry() %>%
  dplyr::group_by(Smp_P_I, Smplngp, CntryNm, Country, Ar_Ql_N, A_Q_N_N, A_Q_S_E, A_Q_S_N_, Ar_Q_S_N, Altitud, A_Q_S_A, A_Q_S_T, day) %>%
  dplyr::summarise(aod865_cams = mean(aod865_cams, na.rn = T))  

names(data1)
names(data2)
names(data3)
names(data4)
names(data5)

# Final dataset
setDT(data1)
setkey(data1, Smp_P_I, Smplngp, CntryNm, Country, Ar_Ql_N, A_Q_N_N, A_Q_S_E, A_Q_S_N_, Ar_Q_S_N, Altitud, A_Q_S_A, A_Q_S_T, day)
setDT(data2)
setkey(data2, Smp_P_I, Smplngp, CntryNm, Country, Ar_Ql_N, A_Q_N_N, A_Q_S_E, A_Q_S_N_, Ar_Q_S_N, Altitud, A_Q_S_A, A_Q_S_T, day)
final2 <- merge(data1, data2, all.x = TRUE)

setDT(final2)
setkey(final2, Smp_P_I, Smplngp, CntryNm, Country, Ar_Ql_N, A_Q_N_N, A_Q_S_E, A_Q_S_N_, Ar_Q_S_N, Altitud, A_Q_S_A, A_Q_S_T, day)
setDT(data3)
setkey(data3, Smp_P_I, Smplngp, CntryNm, Country, Ar_Ql_N, A_Q_N_N, A_Q_S_E, A_Q_S_N_, Ar_Q_S_N, Altitud, A_Q_S_A, A_Q_S_T, day)
final2 <- merge(final2, data3, all.x = TRUE)

setDT(data4)
setkey(data4, Smp_P_I, Smplngp, CntryNm, Country, Ar_Ql_N, A_Q_N_N, A_Q_S_E, A_Q_S_N_, Ar_Q_S_N, Altitud, A_Q_S_A, A_Q_S_T, day)
final2 <- merge(final2, data4, all.x = TRUE)

setDT(data5)
setkey(data5, Smp_P_I, Smplngp, CntryNm, Country, Ar_Ql_N, A_Q_N_N, A_Q_S_E, A_Q_S_N_, Ar_Q_S_N, Altitud, A_Q_S_A, A_Q_S_T, day)
final2 <- merge(final2, data5, all.x = TRUE)

names(final2)
head(final2)
final3 <- final2 %>% 
  dplyr::mutate(year = 2019) %>%
  dplyr::select(Smp_P_I, Smplngp, CntryNm, Country, Ar_Ql_N, A_Q_N_N, A_Q_S_E, A_Q_S_N_, Ar_Q_S_N, Altitud, A_Q_S_A, A_Q_S_T, day, year,
                            aod1240_cams, aod469_cams, aod550_cams, aod670_cams, aod865_cams)

summary(final3)

# Save
saveRDS(final3, "aod_cams_2019_monitors.RDS")
