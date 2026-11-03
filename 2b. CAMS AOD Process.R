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
library(imputeTS)
library(ncdf4)

#### VARIABLES 75 KM ####

variables <- c("1240", 
               "469",
               "550",
               "670",
               "865") 

# 1. Read the temperature map as a raster

raster   <- raster::brick("total_aerosol_optical_depth_1240nm_2019.nc")

plot(raster$X1514764800)
plot(raster$X1536354000, add= TRUE)

# 2. Extract point values -> aim 1: create the dataset structure

centroids             <- as.data.frame(rasterToPoints(raster, spatial = FALSE))
centroids$id_cams_aod <- 1:length(centroids$x)
head(centroids)
centroids_points  <- st_as_sf(centroids, coords = c("x","y"), crs = st_crs(raster))

plot(raster$X1514764800)
plot(st_geometry(centroids_points), add = TRUE, cex = 0.2)

data <- centroids %>%
            select(x, y, id_cams_aod) %>%
            st_drop_geometry()

write_sf(centroids_points[,c("id_cams_aod", "geometry")], "centroids_cams_aod_75km.shp")
                
rm(i)


# 3. Create a loop to extract all point values - 2019 

final <- list()
data0 <- data

for(i in 1:5){
  
  print(paste0(variables[i], " ", Sys.time()))
  
  # Read the AOD maps (nc and then raster) 
  raster   <- raster::brick(paste0("total_aerosol_optical_depth_", variables[i], "nm_2019.nc"))
  
  # Extract point values
  centroids           <- as.data.frame(rasterToPoints(raster, spatial = FALSE))
  centroids$id_cams_aod <- 1:length(centroids$x)
  
  # Merge in a dataset
  data0 <- merge(data0, centroids, by = c("x", "y", "id_cams_aod"), all.x = TRUE)
  
  final[[i]] <- data0
  
  print(paste0(variables[i], "...", Sys.time(), " completed!"))
  
  rm(data0)
  data0 <- data
  
}


# 4. From wide to long format

variables
names(final[[1]])
colnames(final[[1]])[4:2923] <- as.character(rep(seq(as.Date("2019-01-01"), as.Date("2019-12-31"), by = "days"), each = 8))
names(final[[1]]) <- make.unique(names(final[[1]]))
data1 <- final[[1]] %>% pivot_longer(!c(x, y, id_cams_aod), names_to = "date", values_to = "aod1240_cams") %>%
  mutate(day = as.Date(date)) %>%
  st_drop_geometry() %>%
  group_by(x, y, id_cams_aod, day) %>%
  summarise(aod1240_cams = mean(aod1240_cams, na.rn = T)) 

names(final[[2]])
colnames(final[[2]])[4:2923] <- as.character(rep(seq(as.Date("2019-01-01"), as.Date("2019-12-31"), by = "days"), each = 8))
names(final[[2]]) <- make.unique(names(final[[2]]))
data2 <- final[[2]] %>% pivot_longer(!c(x, y, id_cams_aod), names_to = "date", values_to = "aod469_cams")  %>%
  mutate(day = as.Date(date)) %>%
  st_drop_geometry() %>%
  group_by(x, y, id_cams_aod, day) %>%
  summarise(aod469_cams = mean(aod469_cams, na.rn = T)) 

names(final[[3]])
colnames(final[[3]])[4:2923] <- as.character(rep(seq(as.Date("2019-01-01"), as.Date("2019-12-31"), by = "days"), each = 8))
names(final[[3]]) <- make.unique(names(final[[3]]))
data3 <- final[[3]] %>% pivot_longer(!c(x, y, id_cams_aod), names_to = "date", values_to = "aod550_cams")  %>%
  mutate(day = as.Date(date)) %>%
  st_drop_geometry() %>%
  group_by(x, y, id_cams_aod, day) %>%
  summarise(aod550_cams = mean(aod550_cams, na.rn = T)) 

names(final[[4]])
colnames(final[[4]])[4:2923] <- as.character(rep(seq(as.Date("2019-01-01"), as.Date("2019-12-31"), by = "days"), each = 8))
names(final[[4]]) <- make.unique(names(final[[4]]))
data4 <- final[[4]] %>% pivot_longer(!c(x, y, id_cams_aod), names_to = "date", values_to = "aod670_cams")  %>%
  mutate(day = as.Date(date)) %>%
  st_drop_geometry() %>%
  group_by(x, y, id_cams_aod, day) %>%
  summarise(aod670_cams = mean(aod670_cams, na.rn = T))   

names(final[[5]])
colnames(final[[5]])[4:2923] <- as.character(rep(seq(as.Date("2019-01-01"), as.Date("2019-12-31"), by = "days"), each = 8))
names(final[[5]]) <- make.unique(names(final[[5]]))
data5 <- final[[5]] %>% pivot_longer(!c(x, y, id_cams_aod), names_to = "date", values_to = "aod865_cams")  %>%
  mutate(day = as.Date(date)) %>%
  st_drop_geometry() %>%
  group_by(x, y, id_cams_aod, day) %>%
  summarise(aod865_cams = mean(aod865_cams, na.rn = T))  

names(data1)
names(data2)
names(data3)
names(data4)
names(data5)

dim(data1)
dim(data2)
dim(data3)
dim(data4)
dim(data5)

rm(final)
gc()


# 5. Put all togheter

final2 <- merge(data1, data2, by = c("id_cams_aod", "x", "y", "day"))
gc()
final2 <- merge(final2, as.data.frame(data3), by = c("id_cams_aod", "x", "y", "day"))
gc()
final2 <- merge(final2, as.data.frame(data4), by = c("id_cams_aod", "x", "y", "day"))
gc()
final2 <- merge(final2, as.data.frame(data5), by = c("id_cams_aod", "x", "y", "day"))
gc()

names(final2)
final3 <- final2 %>% select(id_cams_aod, x, y, day, aod1240_cams, aod469_cams, aod550_cams, aod670_cams, aod865_cams)
summary(final3) 


saveRDS(final3, "daily_cams_aod_75km_2019.RDS")
gc()
