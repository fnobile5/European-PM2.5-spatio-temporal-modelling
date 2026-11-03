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
library(data.table)

variables <- c("10m_u_component_of_wind",
               "10m_v_component_of_wind",
               "2m_dewpoint_temperature",
               "2m_temperature",
               "mean_sea_level_pressure",            # not available for 9km
               "total_precipitation",
               "boundary_layer_height",              # not available for 9km
               "leaf_area_index_high_vegetation",
               "leaf_area_index_low_vegetation",
               "total_cloud_cover",                  # not available for 9km
               "total_column_ozone")                 # not available for 9km

m <- "01"
i <- 4

# Read temperature data
start_date <- paste0("2019-", m, "-01")
end_date   <- paste0("2019-", m, "-31")

# 1. Read monitors shapefile
monitors      <- read_sf("./monitors_2019.shp")
monitors      <- unique(monitors) 
length(unique(monitors$A_Q_S_E)) 

# 2. Read the temperature map as a raster
temp_raster <- raster::brick(paste("ERA5_data/", variables[i], "_2019_", m, ".nc", sep = ""))

# 3. Check and transform CRS if needed
monitors_crs  <- st_crs(monitors)
temp_crs      <- crs(temp_raster)
monitors_crs == temp_crs 
monitors_repr <- st_transform(monitors, crs = st_crs(temp_raster$X0))
month         <- as.data.frame(monitors_repr)

plot(temp_raster$X0)
plot(st_geometry(monitors_repr), add = TRUE)

rm(i, m)
                      

final <- list()


for(i in 1:11){
  
  print(paste0(variables[i], " ", Sys.time()))
  
  for(m in c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")){
    
    if(m == "02"){

      # Read the temperature map as a raster
      temp_raster <- raster::brick(paste("ERA5_data/", variables[i], "_2019_", m, ".nc", sep = ""))
      
      # Perform spatial join and assign temperature values to monitor points (extract to point)
      temperature_data <- raster::extract(temp_raster, monitors_repr)
      month            <- cbind(month, temperature_data)
      
      gc()
    }
    
    if(m == "04" | m == "06" | m == "09" | m == "11"){

      # Read the temperature map as a raster
      temp_raster <- raster::brick(paste("ERA5_data/", variables[i], "_2019_", m, ".nc", sep = ""))
      
      # Perform spatial join and assign temperature values to monitor points (extract to point)
      temperature_data <- raster::extract(temp_raster, monitors_repr)
      month            <- cbind(month, temperature_data)
      
      gc()
    }
    
    if(m == "01" | m == "03" | m == "05" | m == "07" | m == "08" | m == "10" | m == "12"){

      # Read the temperature map as a raster
      temp_raster <- raster::brick(paste("ERA5_data/", variables[i], "_2019_", m, ".nc", sep = ""))
      
      # Perform spatial join and assign temperature values to monitor points (extract to point)
      temperature_data <- raster::extract(temp_raster, monitors_repr)
      month            <- cbind(month, temperature_data)
      
      gc()
    }
    
    print(paste0(variables[i], m, "...", Sys.time(), " done! "))
    
  }
  
  print(paste0(variables[i], "...", Sys.time(), " completed!"))
  
  final[[i]] <- month
  rm(month)
  month      <- as.data.frame(monitors_repr)
  
}

# From wide to long format

variables
names(final[[1]])
colnames(final[[1]])[15:379] <- as.character(seq(as.Date("2019-01-01"), as.Date("2019-12-31"), by = "days"))
data1 <- final[[1]] %>% 
                dplyr::mutate(year = 2019) %>%
                pivot_longer(!c(Smp_P_I, Smplngp, CntryNm, Country, Ar_Ql_N, A_Q_N_N, A_Q_S_E, A_Q_S_N_, Ar_Q_S_N, Altitud, A_Q_S_A, A_Q_S_T, year, geometry), names_to = "date", values_to = "u10") %>%
                dplyr::mutate(day = as.Date(date))

names(final[[2]])
colnames(final[[2]])[15:379] <- as.character(seq(as.Date("2019-01-01"), as.Date("2019-12-31"), by = "days"))
data2 <- final[[2]]  %>% 
  dplyr::mutate(year = 2019) %>%
  pivot_longer(!c(Smp_P_I, Smplngp, CntryNm, Country, Ar_Ql_N, A_Q_N_N, A_Q_S_E, A_Q_S_N_, Ar_Q_S_N, Altitud, A_Q_S_A, A_Q_S_T, year, geometry), names_to = "date", values_to = "v10") %>%
  dplyr::mutate(day = as.Date(date))

names(final[[3]])
colnames(final[[3]])[15:379] <- as.character(seq(as.Date("2019-01-01"), as.Date("2019-12-31"), by = "days"))
data3 <- final[[3]] %>% 
  dplyr::mutate(year = 2019) %>%
  pivot_longer(!c(Smp_P_I, Smplngp, CntryNm, Country, Ar_Ql_N, A_Q_N_N, A_Q_S_E, A_Q_S_N_, Ar_Q_S_N, Altitud, A_Q_S_A, A_Q_S_T, year, geometry), names_to = "date", values_to = "dewp_2m") %>%
  dplyr::mutate(day = as.Date(date))

names(final[[4]])
colnames(final[[4]])[15:379] <- as.character(seq(as.Date("2019-01-01"), as.Date("2019-12-31"), by = "days"))
data4 <- final[[4]] %>% 
  dplyr::mutate(year = 2019) %>%
  pivot_longer(!c(Smp_P_I, Smplngp, CntryNm, Country, Ar_Ql_N, A_Q_N_N, A_Q_S_E, A_Q_S_N_, Ar_Q_S_N, Altitud, A_Q_S_A, A_Q_S_T, year, geometry), names_to = "date", values_to = "temp_2m") %>%
  dplyr::mutate(day = as.Date(date))

names(final[[5]])
colnames(final[[5]])[15:8774] <- as.character(rep(seq(as.Date("2019-01-01"), as.Date("2019-12-31"), by = "days"), each = 24))
names(final[[5]]) <- make.unique(names(final[[5]]))
data5 <- final[[5]] %>% 
  dplyr::mutate(year = 2019) %>%
  pivot_longer(!c(Smp_P_I, Smplngp, CntryNm, Country, Ar_Ql_N, A_Q_N_N, A_Q_S_E, A_Q_S_N_, Ar_Q_S_N, Altitud, A_Q_S_A, A_Q_S_T, year, geometry), names_to = "date", values_to = "mslp") %>%
  dplyr::mutate(day = as.Date(date)) %>%
  st_drop_geometry() %>%
  dplyr::group_by(Smp_P_I, Smplngp, CntryNm, Country, Ar_Ql_N, A_Q_N_N, A_Q_S_E, A_Q_S_N_, Ar_Q_S_N, Altitud, A_Q_S_A, A_Q_S_T, year, day) %>%
  dplyr::summarise(mslp = mean(mslp, na.rm = T))

names(final[[6]])
colnames(final[[6]])[15:8774] <- as.character(rep(seq(as.Date("2019-01-01"), as.Date("2019-12-31"), by = "days"), each = 24))
names(final[[6]]) <- make.unique(names(final[[6]]))
data6 <- final[[6]] %>% 
  dplyr::mutate(year = 2019) %>%
  pivot_longer(!c(Smp_P_I, Smplngp, CntryNm, Country, Ar_Ql_N, A_Q_N_N, A_Q_S_E, A_Q_S_N_, Ar_Q_S_N, Altitud, A_Q_S_A, A_Q_S_T, year, geometry), names_to = "date", values_to = "tp") %>%
  dplyr::mutate(day = as.Date(date)) %>%
  st_drop_geometry() %>%
  dplyr::group_by(Smp_P_I, Smplngp, CntryNm, Country, Ar_Ql_N, A_Q_N_N, A_Q_S_E, A_Q_S_N_, Ar_Q_S_N, Altitud, A_Q_S_A, A_Q_S_T, year, day) %>%
  dplyr::summarise(tp = mean(tp, na.rm = T)) 

names(final[[7]])
colnames(final[[7]])[15:8774] <- as.character(rep(seq(as.Date("2019-01-01"), as.Date("2019-12-31"), by = "days"), each = 24))
names(final[[7]]) <- make.unique(names(final[[7]]))
data7 <- final[[7]] %>% 
  dplyr::mutate(year = 2019) %>%
  pivot_longer(!c(Smp_P_I, Smplngp, CntryNm, Country, Ar_Ql_N, A_Q_N_N, A_Q_S_E, A_Q_S_N_, Ar_Q_S_N, Altitud, A_Q_S_A, A_Q_S_T, year, geometry), names_to = "date", values_to = "blh") %>%
  dplyr::mutate(day = as.Date(date)) %>%
  st_drop_geometry() %>%
  dplyr::group_by(Smp_P_I, Smplngp, CntryNm, Country, Ar_Ql_N, A_Q_N_N, A_Q_S_E, A_Q_S_N_, Ar_Q_S_N, Altitud, A_Q_S_A, A_Q_S_T, year, day) %>%
  dplyr::summarise(blh = mean(blh, na.rm = T))  

names(final[[8]])
colnames(final[[8]])[15:379] <- as.character(seq(as.Date("2019-01-01"), as.Date("2019-12-31"), by = "days"))
data8 <- final[[8]] %>% 
  dplyr::mutate(year = 2019) %>%
  pivot_longer(!c(Smp_P_I, Smplngp, CntryNm, Country, Ar_Ql_N, A_Q_N_N, A_Q_S_E, A_Q_S_N_, Ar_Q_S_N, Altitud, A_Q_S_A, A_Q_S_T, year, geometry), names_to = "date", values_to = "leaf_high") %>%
  dplyr::mutate(day = as.Date(date))

names(final[[9]])
colnames(final[[9]])[15:379] <- as.character(seq(as.Date("2019-01-01"), as.Date("2019-12-31"), by = "days"))
data9 <- final[[9]] %>% 
  dplyr::mutate(year = 2019) %>%
  pivot_longer(!c(Smp_P_I, Smplngp, CntryNm, Country, Ar_Ql_N, A_Q_N_N, A_Q_S_E, A_Q_S_N_, Ar_Q_S_N, Altitud, A_Q_S_A, A_Q_S_T, year, geometry), names_to = "date", values_to = "leaf_low") %>%
  dplyr::mutate(day = as.Date(date))

names(final[[10]])
colnames(final[[10]])[15:8774] <- as.character(rep(seq(as.Date("2019-01-01"), as.Date("2019-12-31"), by = "days"), each = 24))
names(final[[10]]) <- make.unique(names(final[[10]]))
data10 <- final[[10]] %>% 
  dplyr::mutate(year = 2019) %>%
  pivot_longer(!c(Smp_P_I, Smplngp, CntryNm, Country, Ar_Ql_N, A_Q_N_N, A_Q_S_E, A_Q_S_N_, Ar_Q_S_N, Altitud, A_Q_S_A, A_Q_S_T, year, geometry), names_to = "date", values_to = "tcc") %>%
  dplyr::mutate(day = as.Date(date)) %>%
  st_drop_geometry() %>%
  dplyr::group_by(Smp_P_I, Smplngp, CntryNm, Country, Ar_Ql_N, A_Q_N_N, A_Q_S_E, A_Q_S_N_, Ar_Q_S_N, Altitud, A_Q_S_A, A_Q_S_T, year, day) %>%
  dplyr::summarise(tcc = mean(tcc, na.rm = T)) 

names(final[[11]])
colnames(final[[11]])[15:8774] <- as.character(rep(seq(as.Date("2019-01-01"), as.Date("2019-12-31"), by = "days"), each = 24))
names(final[[11]]) <- make.unique(names(final[[11]]))
data11 <- final[[11]] %>% 
  dplyr::mutate(year = 2019) %>%
  pivot_longer(!c(Smp_P_I, Smplngp, CntryNm, Country, Ar_Ql_N, A_Q_N_N, A_Q_S_E, A_Q_S_N_, Ar_Q_S_N, Altitud, A_Q_S_A, A_Q_S_T, year, geometry), names_to = "date", values_to = "tco") %>%
  dplyr::mutate(day = as.Date(date)) %>%
  st_drop_geometry() %>%
  dplyr::group_by(Smp_P_I, Smplngp, CntryNm, Country, Ar_Ql_N, A_Q_N_N, A_Q_S_E, A_Q_S_N_, Ar_Q_S_N, Altitud, A_Q_S_A, A_Q_S_T, year, day) %>%
  dplyr::summarise(tco = mean(tco, na.rm = T)) 

# Final dataset
names(data1)
data1 <- data1 %>% dplyr::select(-c(date)) %>% st_drop_geometry()
data1 <- data1 %>% dplyr::select(-c(geometry)) 
names(data2)
data2 <- data2 %>% dplyr::select(-c(date)) %>% st_drop_geometry()
data2 <- data2 %>% dplyr::select(-c(geometry))
names(data3)
data3 <- data3 %>% dplyr::select(-c(date)) %>% st_drop_geometry()
data3 <- data3 %>% dplyr::select(-c(geometry))
names(data4)
data4 <- data4 %>% dplyr::select(-c(date)) %>% st_drop_geometry()
data4 <- data4 %>% dplyr::select(-c(geometry))
names(data5)
names(data6)
names(data7)
names(data8)
data8 <- data8 %>% dplyr::select(-c(date)) %>% st_drop_geometry()
data8 <- data8 %>% dplyr::select(-c(geometry))
names(data9)
data9 <- data9 %>% dplyr::select(-c(date)) %>% st_drop_geometry()
data9 <- data9 %>% dplyr::select(-c(geometry)) 
names(data10)
names(data11)

rm(final)
gc()

setDT(data1)
setkey(data1, Smp_P_I, Smplngp, CntryNm, Country, Ar_Ql_N, A_Q_N_N, A_Q_S_E, A_Q_S_N_, Ar_Q_S_N, Altitud, A_Q_S_A, A_Q_S_T, year, day)
setDT(data2)
setkey(data2, Smp_P_I, Smplngp, CntryNm, Country, Ar_Ql_N, A_Q_N_N, A_Q_S_E, A_Q_S_N_, Ar_Q_S_N, Altitud, A_Q_S_A, A_Q_S_T, year, day)
final2 <- merge(data1, data2, all.x = TRUE)

setDT(final2)
setkey(final2, Smp_P_I, Smplngp, CntryNm, Country, Ar_Ql_N, A_Q_N_N, A_Q_S_E, A_Q_S_N_, Ar_Q_S_N, Altitud, A_Q_S_A, A_Q_S_T, year, day)
setDT(data3)
setkey(data3, Smp_P_I, Smplngp, CntryNm, Country, Ar_Ql_N, A_Q_N_N, A_Q_S_E, A_Q_S_N_, Ar_Q_S_N, Altitud, A_Q_S_A, A_Q_S_T, year, day)
final2 <- merge(final2, data3, all.x = TRUE)

setDT(data4)
setkey(data4, Smp_P_I, Smplngp, CntryNm, Country, Ar_Ql_N, A_Q_N_N, A_Q_S_E, A_Q_S_N_, Ar_Q_S_N, Altitud, A_Q_S_A, A_Q_S_T, year, day)
final2 <- merge(final2, data4, all.x = TRUE)

setDT(data5)
setkey(data5, Smp_P_I, Smplngp, CntryNm, Country, Ar_Ql_N, A_Q_N_N, A_Q_S_E, A_Q_S_N_, Ar_Q_S_N, Altitud, A_Q_S_A, A_Q_S_T, year, day)
final2 <- merge(final2, data5, all.x = TRUE)

setDT(data6)
setkey(data6, Smp_P_I, Smplngp, CntryNm, Country, Ar_Ql_N, A_Q_N_N, A_Q_S_E, A_Q_S_N_, Ar_Q_S_N, Altitud, A_Q_S_A, A_Q_S_T, year, day)
final2 <- merge(final2, data6, all.x = TRUE)

setDT(data7)
setkey(data7, Smp_P_I, Smplngp, CntryNm, Country, Ar_Ql_N, A_Q_N_N, A_Q_S_E, A_Q_S_N_, Ar_Q_S_N, Altitud, A_Q_S_A, A_Q_S_T, year, day)
final2 <- merge(final2, data7, all.x = TRUE)

setDT(data8)
setkey(data8, Smp_P_I, Smplngp, CntryNm, Country, Ar_Ql_N, A_Q_N_N, A_Q_S_E, A_Q_S_N_, Ar_Q_S_N, Altitud, A_Q_S_A, A_Q_S_T, year, day)
final2 <- merge(final2, data8, all.x = TRUE)

setDT(data9)
setkey(data9, Smp_P_I, Smplngp, CntryNm, Country, Ar_Ql_N, A_Q_N_N, A_Q_S_E, A_Q_S_N_, Ar_Q_S_N, Altitud, A_Q_S_A, A_Q_S_T, year, day)
final2 <- merge(final2, data9, all.x = TRUE)

setDT(data10)
setkey(data10, Smp_P_I, Smplngp, CntryNm, Country, Ar_Ql_N, A_Q_N_N, A_Q_S_E, A_Q_S_N_, Ar_Q_S_N, Altitud, A_Q_S_A, A_Q_S_T, year, day)
final2 <- merge(final2, data10, all.x = TRUE)

setDT(data11)
setkey(data11, Smp_P_I, Smplngp, CntryNm, Country, Ar_Ql_N, A_Q_N_N, A_Q_S_E, A_Q_S_N_, Ar_Q_S_N, Altitud, A_Q_S_A, A_Q_S_T, year, day)
final2 <- merge(final2, data11, all.x = TRUE)

names(final2)
final2 <- final2 %>% dplyr::select(Smp_P_I, Smplngp, CntryNm, Country, Ar_Ql_N, A_Q_N_N, A_Q_S_E, A_Q_S_N_, Ar_Q_S_N, Altitud, A_Q_S_A, A_Q_S_T, year, day,
                            u10, v10, dewp_2m, temp_2m, mslp,
                            tp, blh, leaf_high, leaf_low, tcc, tco)

dim(final2)

# Save
saveRDS(final2, "era5_2019_monitors.RDS")
