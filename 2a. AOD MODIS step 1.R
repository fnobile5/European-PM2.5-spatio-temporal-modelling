#### European-PM2.5-spatio-temporal-modelling

#### Author: Federica Nobile
#### Affiliation: Institute of Environmental Medicine, Karolinska Institutet, Stockholm, Sweden
####              Department of Epidemiology, Lazio Region Health Service/ASL Roma 1, Rome, Italy

# Load necessary libraries
library(terra)
library(sf)
library(tidyr)
library(lubridate)
library(data.table)
# library(dplyr)

#### 1a. I ASSIGN AOD MODIS VALUES TO 1 KM GRID - AOD047 ####
# Raster to point is to be preferred to the nearest (due to missing values)

# 1. Read points shapefile (centroids 1 km grid)

points <- read_sf("./Europe_grid_1000_3035_centroids_permaiac.shp")
head(points)
gc()

# 2. Read the rasters
# Since now I'm using GEE for the download, I need to re-construct the complete layers (downloaded in pieces) - different from 2019 approach

files <- list.files(
  path = "AOD Raw/",
  pattern = "^Daily_OpticalDepth047_2019-0000000000.*\\.tif$", 
  full.names = TRUE
)

rasters             <- lapply(files, rast)
europa_mosaic_047_1 <- do.call(mosaic, rasters)

# plot(europa_mosaic_047_1$'0_Optical_Depth_047')
# plot(europa_mosaic_047_1$'364_Optical_Depth_047')
gc()

files <- list.files(
  path = "AOD Raw/",
  pattern = "^Daily_OpticalDepth047_2019-0000001280.*\\.tif$", 
  full.names = TRUE
)

rasters             <- lapply(files, rast)
europa_mosaic_047_2 <- do.call(mosaic, rasters)

# plot(europa_mosaic_047_2$'0_Optical_Depth_047')
# plot(europa_mosaic_047_2$'364_Optical_Depth_047')
gc()

files <- list.files(
  path = "AOD Raw/",
  pattern = "^Daily_OpticalDepth047_2019-0000002560.*\\.tif$", 
  full.names = TRUE
)

rasters             <- lapply(files, rast)
europa_mosaic_047_3 <- do.call(mosaic, rasters)

# plot(europa_mosaic_047_3$'0_Optical_Depth_047')
# plot(europa_mosaic_047_3$'364_Optical_Depth_047')
gc()

files <- list.files(
  path = "AOD Raw/",
  pattern = "^Daily_OpticalDepth047_2019-0000003840.*\\.tif$",
  full.names = TRUE
)

rasters             <- lapply(files, rast)
europa_mosaic_047_4 <- do.call(mosaic, rasters)

# plot(europa_mosaic_047_4$'0_Optical_Depth_047')
# plot(europa_mosaic_047_4$'364_Optical_Depth_047')
gc()

files <- list.files(
  path = "AOD Raw/",
  pattern = "^Daily_OpticalDepth047_2019-0000005120.*\\.tif$",  
  full.names = TRUE
)

rasters             <- lapply(files, rast)
europa_mosaic_047_5 <- do.call(mosaic, rasters)

# plot(europa_mosaic_047_5$'0_Optical_Depth_047')
# plot(europa_mosaic_047_5$'364_Optical_Depth_047')
gc()


# 3. Check and transform CRS if needed

points_crs  <- st_crs(points)
raster_crs  <- crs(europa_mosaic_047_2$'0_Optical_Depth_047')
points_crs == raster_crs
points_repr <- st_transform(points, crs = st_crs(europa_mosaic_047_2$'0_Optical_Depth_047'))

day <- points %>% dplyr::select(layer) %>% st_drop_geometry()

plot(europa_mosaic_047_2$'0_Optical_Depth_047')
# plot(st_geometry(points_repr), add = TRUE) # It can take a long time

rm(i)

columns <- names(europa_mosaic_047_1)

# 4. Loop for computing 30 days together

for(i in c(1, 31, 61, 91, 121, 151, 181, 211, 241, 271, 301)){
  
  print(paste0(i, " ", Sys.time()))
  
  j = i + 29
  print(paste0(j, " ", Sys.time()))
  
  print(columns[i:j])
  
  # Read the map as a raster
  raster <- mosaic(europa_mosaic_047_1[[columns[i:j]]], europa_mosaic_047_2[[columns[i:j]]], europa_mosaic_047_3[[columns[i:j]]], europa_mosaic_047_4[[columns[i:j]]], europa_mosaic_047_5[[columns[i:j]]])
      
  # Perform spatial join and assign raster values to points of the grid (extract to point)
  data              <- raster::extract(raster, points_repr)
  data              <- as.data.frame(data)
  # colnames(data)[2] <- paste0("day_", columns[i])
  
  day <- cbind(day, data)
  
  print(paste0(i, " ", Sys.time(), " completed!"))
  gc()
  
} 

head(day)
names(day)

# 12th group

i = 331
print(paste0(i, " ", Sys.time()))
j = 365
print(paste0(j, " ", Sys.time()))

print(columns[i:j])

# Read the map as a raster
raster <- mosaic(europa_mosaic_047_1[[columns[i:j]]], europa_mosaic_047_2[[columns[i:j]]], europa_mosaic_047_3[[columns[i:j]]], europa_mosaic_047_4[[columns[i:j]]], europa_mosaic_047_5[[columns[i:j]]])

# Perform spatial join and assign raster values to points of the grid (extract to point)
data              <- raster::extract(raster, points_repr)
data              <- as.data.frame(data)
# colnames(data)[2] <- paste0("day_", columns[i])

day <- cbind(day, data)

print(paste0(i, " ", Sys.time(), " completed!"))
gc()


names(day)
day <- day[,-c(2, 33, 64, 95, 126, 157, 188, 219, 250, 281, 312, 343)]
names(day)
head(day)
saveRDS(day, "./daily_aod47_maiac_1km_2019.RDS")
rm(list = c(
  "europa_mosaic_047",
  "europa_mosaic_047_1",
  "europa_mosaic_047_2",
  "europa_mosaic_047_3",
  "europa_mosaic_047_4",
  "europa_mosaic_047_5"
))
gc()


# 5. To long format - by month

dates <- as.data.frame(c(1:365))
dates$`c(1:365)` <- paste0("day", dates$`c(1:365)`)
dates$day <- c(seq(as.Date("2019-01-01"), as.Date("2019-12-31"), 1))
colnames(dates)[1] <- "jd"
dates <- cbind(dates, columns)
colnames(dates)[3] <- "label"
head(dates)



day1_1 <- day %>% dplyr::select(layer, '0_Optical_Depth_047':'30_Optical_Depth_047')
colnames(day1_1)[2:32] <- c(dates$day[1:31]) 
head(day1_1)
day1_1 <- day1_1 %>% pivot_longer(!c(layer), names_to = "date", values_to = "aod047")
day1_1 <- day1_1 %>% dplyr::mutate(day = as.Date(as.numeric(date), "1970-01-01"))
head(day1_1)
dim(day1_1) 
saveRDS(day1_1, "./daily_aod_maiac_1km_201901.RDS")
rm(day1_1)
gc()

day1_2 <- day %>% dplyr::select(layer, '31_Optical_Depth_047':'58_Optical_Depth_047')
colnames(day1_2)[2:29] <- c(dates$day[32:59]) 
head(day1_2)
day1_2 <- day1_2 %>% pivot_longer(!c(layer), names_to = "date", values_to = "aod047")
day1_2 <- day1_2 %>% dplyr::mutate(day = as.Date(as.numeric(date), "1970-01-01"))
head(day1_2)
dim(day1_2) 
saveRDS(day1_2, "./daily_aod_maiac_1km_201902.RDS")
rm(day1_2)
gc()

day1_3 <- day %>% dplyr::select(layer, '59_Optical_Depth_047':'89_Optical_Depth_047')
colnames(day1_3)[2:32] <- c(dates$day[60:90])
head(day1_3)
day1_3 <- day1_3 %>% pivot_longer(!c(layer), names_to = "date", values_to = "aod047")
day1_3 <- day1_3 %>% dplyr::mutate(day = as.Date(as.numeric(date), "1970-01-01"))
head(day1_3)
dim(day1_3) 
saveRDS(day1_3, "./daily_aod_maiac_1km_201903.RDS")
rm(day1_3)
gc()

day1_4 <- day %>% dplyr::select(layer, '90_Optical_Depth_047':'119_Optical_Depth_047')
colnames(day1_4)[2:31] <- c(dates$day[91:120]) 
head(day1_4)
day1_4 <- day1_4 %>% pivot_longer(!c(layer), names_to = "date", values_to = "aod047")
day1_4 <- day1_4 %>% dplyr::mutate(day = as.Date(as.numeric(date), "1970-01-01"))
head(day1_4)
dim(day1_4) 
saveRDS(day1_4, "./daily_aod_maiac_1km_201904.RDS")
rm(day1_4)
gc()

day1_5 <- day %>% dplyr::select(layer, '120_Optical_Depth_047':'150_Optical_Depth_047')
colnames(day1_5)[2:32] <- c(dates$day[121:151]) #
head(day1_5)
day1_5 <- day1_5 %>% pivot_longer(!c(layer), names_to = "date", values_to = "aod047")
day1_5 <- day1_5 %>% dplyr::mutate(day = as.Date(as.numeric(date), "1970-01-01"))
head(day1_5)
dim(day1_5) 
saveRDS(day1_5, "./daily_aod_maiac_1km_201905.RDS")
rm(day1_5)
gc()

day1_6 <- day %>% dplyr::select(layer, '151_Optical_Depth_047':'180_Optical_Depth_047')
colnames(day1_6)[2:31] <- c(dates$day[152:181]) 
head(day1_6)
day1_6 <- day1_6 %>% pivot_longer(!c(layer), names_to = "date", values_to = "aod047")
day1_6 <- day1_6 %>% dplyr::mutate(day = as.Date(as.numeric(date), "1970-01-01"))
head(day1_6)
dim(day1_6) 
saveRDS(day1_6, "./daily_aod_maiac_1km_201906.RDS")
rm(day1_6)
gc()

day1_7 <- day %>% dplyr::select(layer, '181_Optical_Depth_047':'211_Optical_Depth_047')
colnames(day1_7)[2:32] <- c(dates$day[182:212]) 
head(day1_7)
day1_7 <- day1_7 %>% pivot_longer(!c(layer), names_to = "date", values_to = "aod047")
day1_7 <- day1_7 %>% dplyr::mutate(day = as.Date(as.numeric(date), "1970-01-01"))
head(day1_7)
dim(day1_7)
saveRDS(day1_7, "./daily_aod_maiac_1km_201907.RDS")
rm(day1_7)
gc()

day1_8 <- day %>% dplyr::select(layer, '212_Optical_Depth_047':'242_Optical_Depth_047')
colnames(day1_8)[2:32] <- c(dates$day[213:243]) 
head(day1_8)
day1_8 <- day1_8 %>% pivot_longer(!c(layer), names_to = "date", values_to = "aod047")
day1_8 <- day1_8 %>% dplyr::mutate(day = as.Date(as.numeric(date), "1970-01-01"))
head(day1_8)
dim(day1_8)  
saveRDS(day1_8, "./daily_aod_maiac_1km_201908.RDS")
rm(day1_8)
gc()

day2_9 <- day %>% dplyr::select(layer, '243_Optical_Depth_047':'272_Optical_Depth_047')
colnames(day2_9)[2:31] <- c(dates$day[244:273]) 
head(day2_9)
day2_9 <- day2_9 %>% pivot_longer(!c(layer), names_to = "date", values_to = "aod047")
day2_9 <- day2_9 %>% dplyr::mutate(day = as.Date(as.numeric(date), "1970-01-01"))
head(day2_9)
dim(day2_9) 
saveRDS(day2_9, "./daily_aod_maiac_1km_201909.RDS")
rm(day2_9)
gc()

day2_10 <- day %>% dplyr::select(layer, '273_Optical_Depth_047':'303_Optical_Depth_047')
colnames(day2_10)[2:32] <- c(dates$day[274:304]) 
head(day2_10)
day2_10 <- day2_10 %>% pivot_longer(!c(layer), names_to = "date", values_to = "aod047")
day2_10 <- day2_10 %>% dplyr::mutate(day = as.Date(as.numeric(date), "1970-01-01"))
head(day2_10)
dim(day2_10) 
saveRDS(day2_10, "./daily_aod_maiac_1km_201910.RDS")
rm(day2_10)
gc()

day2_11 <- day %>% dplyr::select(layer, '304_Optical_Depth_047':'333_Optical_Depth_047')
colnames(day2_11)[2:31] <- c(dates$day[305:334]) 
head(day2_11)
day2_11 <- day2_11 %>% pivot_longer(!c(layer), names_to = "date", values_to = "aod047")
day2_11 <- day2_11 %>% dplyr::mutate(day = as.Date(as.numeric(date), "1970-01-01"))
head(day2_11)
dim(day2_11) 
saveRDS(day2_11, "./daily_aod_maiac_1km_201911.RDS")
rm(day2_11)
gc()

day2_12 <- day %>% dplyr::select(layer, '334_Optical_Depth_047':'364_Optical_Depth_047')
colnames(day2_12)[2:32] <- c(dates$day[335:365]) 
head(day2_12)
day2_12 <- day2_12 %>% pivot_longer(!c(layer), names_to = "date", values_to = "aod047")
day2_12 <- day2_12 %>% dplyr::mutate(day = as.Date(as.numeric(date), "1970-01-01"))
head(day2_12)
dim(day2_12) 
saveRDS(day2_12, "./daily_aod_maiac_1km_201912.RDS")
rm(day2_12)
gc()


#### 2. I ASSIGN AOD CAMS VALUES TO 1 KM GRID #### (I start again from here after updating of AOD from CAMS)
# The nearest is fine

# 1. Read points shapefile (centroids 1 km grid)

points <- read_sf("./Europe_grid_1000_3035_centroids_permaiac.shp") 
head(points)

# 2. Read centroids shapefile (75 km grid - cams)

points_75 <- read_sf("./centroids_cams_aod_75km.shp") 
head(points_75)
st_crs(points) == st_crs(points_75) 
points_75_repr <- st_transform(points_75, crs = st_crs(points))

# 3. Assign the nearest 75 km point to the 1 km point

nearest             <- st_nearest_feature(points, points_75_repr)
points$id_cams_75km <- points_75_repr$id_cms_[nearest]

# To check the distances
library(geos)
dist                  <- geos_distance(points, points_75_repr[nearest,])
points$dist_to_points <- dist
summary(points$dist_to_points)

points_final <- points %>% dplyr::filter(!(layer %in% vedi$layer)) %>%# 4,946,280
  dplyr::select(!(geometry)) %>%
  st_drop_geometry()
points_final <- as.data.frame(points_final)
head(points_final)

# 4. Merge with the AOD CAMS dataset - divided by month

aod_cams_data <- readRDS("./daily_cams_aod_75km_2019.RDS")
head(aod_cams_data)
aod_cams_data <- aod_cams_data %>% dplyr::mutate(mm = substr(day, 6, 7)) %>%
                                    dplyr::rename(id_cams_75km = id_cams_aod)

cams1 <- subset(aod_cams_data, mm == "01")
setDT(cams1)
setkey(cams1, id_cams_75km, day)
points_final1 <- points_final %>%
  tidyr::crossing(day = seq.Date(from = as.Date("2019-01-01"),
                                 to   = as.Date("2019-01-31"),
                                 by   = "day"))
setDT(points_final1)
setkey(points_final1, id_cams_75km, day)
data1 <- merge(points_final1, cams1, all.x = TRUE)
dim(data1) 
head(data1)
rm(cams1)
saveRDS(data1, "./daily_aod_cams_1km_201901.RDS")
rm(data1)
gc()

cams2 <- subset(aod_cams_data, mm == "02")
setDT(cams2)
setkey(cams2, id_cams_75km, day)
points_final2 <- points_final %>%
  tidyr::crossing(day = seq.Date(from = as.Date("2019-02-01"),
                                 to   = as.Date("2019-02-28"),
                                 by   = "day"))
setDT(points_final2)
setkey(points_final2, id_cams_75km, day)
data2 <- merge(points_final2, cams2, all.x = TRUE)
dim(data2) 
head(data2)
rm(cams2)
saveRDS(data2, "./daily_aod_cams_1km_201902.RDS")
rm(data2)
gc()

cams3 <- subset(aod_cams_data, mm == "03")
setDT(cams3)
setkey(cams3, id_cams_75km, day)
points_final3 <- points_final %>%
  tidyr::crossing(day = seq.Date(from = as.Date("2019-03-01"),
                                 to   = as.Date("2019-03-31"),
                                 by   = "day"))
setDT(points_final3)
setkey(points_final3, id_cams_75km, day)
data3 <- merge(points_final3, cams3, all.x = TRUE)
dim(data3) 
head(data3)
rm(cams3)
saveRDS(data3, "./daily_aod_cams_1km_201903.RDS")
rm(data3)
gc()

cams4 <- subset(aod_cams_data, mm == "04")
setDT(cams4)
setkey(cams4, id_cams_75km, day)
points_final4 <- points_final %>%
  tidyr::crossing(day = seq.Date(from = as.Date("2019-04-01"),
                                 to   = as.Date("2019-04-30"),
                                 by   = "day"))
setDT(points_final4)
setkey(points_final4, id_cams_75km, day)
data4 <- merge(points_final4, cams4, all.x = TRUE)
dim(data4) 
head(data4)
rm(cams4)
saveRDS(data4, "./daily_aod_cams_1km_201904.RDS")
rm(data4)
gc()

cams5 <- subset(aod_cams_data, mm == "05")
setDT(cams5)
setkey(cams5, id_cams_75km, day)
points_final5 <- points_final %>%
  tidyr::crossing(day = seq.Date(from = as.Date("2019-05-01"),
                                 to   = as.Date("2019-05-31"),
                                 by   = "day"))
setDT(points_final5)
setkey(points_final5, id_cams_75km, day)
data5 <- merge(points_final5, cams5, all.x = TRUE)
dim(data5) 
head(data5)
rm(cams5)
saveRDS(data5, "./daily_aod_cams_1km_201905.RDS")
rm(data5)
gc()

cams6 <- subset(aod_cams_data, mm == "06")
setDT(cams6)
setkey(cams6, id_cams_75km, day)
points_final6 <- points_final %>%
  tidyr::crossing(day = seq.Date(from = as.Date("2019-06-01"),
                                 to   = as.Date("2019-06-30"),
                                 by   = "day"))
setDT(points_final6)
setkey(points_final6, id_cams_75km, day)
data6 <- merge(points_final6, cams6, all.x = TRUE)
dim(data6) 
head(data6)
rm(cams6)
saveRDS(data6, "./daily_aod_cams_1km_201906.RDS")
rm(data6)
gc()

cams7 <- subset(aod_cams_data, mm == "07")
setDT(cams7)
setkey(cams7, id_cams_75km, day)
points_final7 <- points_final %>%
  tidyr::crossing(day = seq.Date(from = as.Date("2019-07-01"),
                                 to   = as.Date("2019-07-31"),
                                 by   = "day"))
setDT(points_final7)
setkey(points_final7, id_cams_75km, day)
data7 <- merge(points_final7, cams7, all.x = TRUE)
dim(data7) 
head(data7)
rm(cams7)
saveRDS(data7, "./daily_aod_cams_1km_201907.RDS")
rm(data7)
gc()

cams8 <- subset(aod_cams_data, mm == "08")
setDT(cams8)
setkey(cams8, id_cams_75km, day)
points_final8 <- points_final %>%
  tidyr::crossing(day = seq.Date(from = as.Date("2019-08-01"),
                                 to   = as.Date("2019-08-31"),
                                 by   = "day"))
setDT(points_final8)
setkey(points_final8, id_cams_75km, day)
data8 <- merge(points_final8, cams8, all.x = TRUE)
dim(data8) 
head(data8)
rm(cams8)
saveRDS(data8, "./daily_aod_cams_1km_201908.RDS")
rm(data8)
gc()

cams9 <- subset(aod_cams_data, mm == "09")
setDT(cams9)
setkey(cams9, id_cams_75km, day)
points_final9 <- points_final %>%
  tidyr::crossing(day = seq.Date(from = as.Date("2019-09-01"),
                                 to   = as.Date("2019-09-30"),
                                 by   = "day"))
setDT(points_final9)
setkey(points_final9, id_cams_75km, day)
data9 <- merge(points_final9, cams9, all.x = TRUE)
dim(data9) 
head(data9)
rm(cams9)
saveRDS(data9, "./daily_aod_cams_1km_201909.RDS")
rm(data9)
gc()

cams10 <- subset(aod_cams_data, mm == "10")
setDT(cams10)
setkey(cams10, id_cams_75km, day)
points_final10 <- points_final %>%
  tidyr::crossing(day = seq.Date(from = as.Date("2019-10-01"),
                                 to   = as.Date("2019-10-31"),
                                 by   = "day"))
setDT(points_final10)
setkey(points_final10, id_cams_75km, day)
data10 <- merge(points_final10, cams10, all.x = TRUE)
dim(data10) 
head(data10)
rm(cams10)
saveRDS(data10, "./daily_aod_cams_1km_201910.RDS")
rm(data10)
gc()

cams11 <- subset(aod_cams_data, mm == "11")
setDT(cams11)
setkey(cams11, id_cams_75km, day)
points_final11 <- points_final %>%
  tidyr::crossing(day = seq.Date(from = as.Date("2019-11-01"),
                                 to   = as.Date("2019-11-30"),
                                 by   = "day"))
setDT(points_final11)
setkey(points_final11, id_cams_75km, day)
data11 <- merge(points_final11, cams11, all.x = TRUE)
dim(data11) 
head(data11)
rm(cams11)
saveRDS(data11, "./daily_aod_cams_1km_201911.RDS")
rm(data11)
gc()

cams12 <- subset(aod_cams_data, mm == "12")
setDT(cams12)
setkey(cams12, id_cams_75km, day)
points_final12 <- points_final %>%
  tidyr::crossing(day = seq.Date(from = as.Date("2019-12-01"),
                                 to   = as.Date("2019-12-31"),
                                 by   = "day"))
setDT(points_final12)
setkey(points_final12, id_cams_75km, day)
data12 <- merge(points_final12, cams12, all.x = TRUE)
dim(data12) 
head(data12)
rm(cams12)
saveRDS(data12, "./daily_aod_cams_1km_201912.RDS")
rm(data12)
gc()


rm(list = ls())
gc()

#### 3. I CREATE MONTHLY DATASETS WITH AOD CAMS VALUES AND AOD MAIAC VALUES (EVEN IF WITH MISSING VALUES) - 1 KM GRID ####

months <- c("01", "02", "03", "04", "05", "06",
             "07", "08", "09", "10", "11", "12")

for(i in 1:12){
  
  month <- months[i]
  print(paste0(month, " ", Sys.time()))
  
  cams     <- readRDS(paste0("./daily_aod_cams_1km_2019", month, ".RDS"))
  print(head(cams))
  
  maiac047 <- readRDS(paste0("./daily_aod_maiac_1km_2019", month, ".RDS"))
  print(head(maiac047))
  
  setDT(cams)
  setkey(cams, layer, day)
  setDT(maiac047)
  setkey(maiac047, layer, day)
  data <- merge(cams, maiac047, all.x = TRUE)
  print(names(data))
  data <- data %>% dplyr::select(layer, id_cams_75km, dist_to_points, day, aod1240_cams, aod469_cams, aod550_cams, aod670_cams, aod865_cams, aod047)
  print(head(data))
  
  saveRDS(data, paste0("./daily_aod047_cams_maiac_1km_2019", month, ".RDS"))
  
  rm(cams, maiac047)
  gc()
  
  print(paste0(month, " ", Sys.time(), " completed!"))
  
  gc()
  
} 