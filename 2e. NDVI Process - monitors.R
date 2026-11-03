#### European-PM2.5-spatio-temporal-modelling

#### Author: Federica Nobile
#### Affiliation: Institute of Environmental Medicine, Karolinska Institutet, Stockholm, Sweden
####              Department of Epidemiology, Lazio Region Health Service/ASL Roma 1, Rome, Italy

# Load necessary libraries
library(terra)
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
library(data.table)
library(gdata)
library(nngeo)

variables <- c("NDVI", 
               "EVI" )                 

# 1. Read monitors shapefile; but actually it doesn't be needed because I have to work on the complete grid for spatially imputing missings.
monitors      <- read_sf("./monitors_2019.shp")
monitors      <- unique(monitors)

# 2. Read the rasters

files <- list.files(
  path = "GEE Layers/",
  pattern = "^NDVI_Monthly_.*\\.tif$",  
  full.names = TRUE
)

rasters             <- lapply(files, rast)
files # check on the order

for (i in seq_along(rasters)) { 
  names(rasters[[i]]) <- paste0("layer_", i) # unique names
}

raster_ndvi         <- terra::rast(rasters)

raster_ndvi_raster  <- raster(raster_ndvi)
coordinates         <- as.data.frame(rasterToPoints(raster_ndvi_raster, spatial = FALSE))
coordinates         <- coordinates %>% dplyr::arrange(x, y) %>% dplyr::mutate(id = dplyr::row_number()) %>% dplyr::select(-c(layer_1))
coordinates_ndvi    <- st_as_sf(coordinates, coords = c("x", "y"), crs = 4326) 

plot(raster_ndvi$layer_1)
plot(raster_ndvi$layer_2)
plot(raster_ndvi$layer_3)
plot(raster_ndvi$layer_4)
plot(raster_ndvi$layer_5)
plot(raster_ndvi$layer_6)
plot(raster_ndvi$layer_7)
plot(raster_ndvi$layer_8)
plot(raster_ndvi$layer_9)
plot(raster_ndvi$layer_10)
plot(raster_ndvi$layer_11)
plot(raster_ndvi$layer_12)
# plot(st_geometry(coordinates_ndvi), add = T)
gc()



files <- list.files(
  path = "GEE Layers/",
  pattern = "^EVI_Monthly_.*\\.tif$",  
  full.names = TRUE
)

rasters             <- lapply(files, rast)
files # check on the order

for (i in seq_along(rasters)) { 
  names(rasters[[i]]) <- paste0("layer_", i) # unique names
}

raster_evi         <- terra::rast(rasters)
raster_evi_raster  <- raster(raster_evi)
coordinates        <- as.data.frame(rasterToPoints(raster_evi_raster, spatial = FALSE))
coordinates        <- coordinates %>% dplyr::arrange(x, y) %>% dplyr::mutate(id = dplyr::row_number()) %>% dplyr::select(-c(layer_1))
coordinates_evi    <- st_as_sf(coordinates, coords = c("x", "y"), crs = 4326) 

plot(raster_evi$layer_1)
plot(raster_evi$layer_2)
plot(raster_evi$layer_3)
plot(raster_evi$layer_4)
plot(raster_evi$layer_5)
plot(raster_evi$layer_6)
plot(raster_evi$layer_7)
plot(raster_evi$layer_8)
plot(raster_evi$layer_9)
plot(raster_evi$layer_10)
plot(raster_evi$layer_11)
plot(raster_evi$layer_12)
gc()

# 3. Check and transform CRS if needed
monitors_crs  <- st_crs(monitors)
monitors_repr <- st_transform(monitors, crs = st_crs(raster_evi$layer_1))
plot(st_geometry(monitors_repr), add = T)
# final <- as.data.frame(monitors_repr)

st_crs(coordinates_ndvi) == st_crs(raster_ndvi$layer_1) # TRUE
st_crs(coordinates_evi)  == st_crs(raster_evi$layer_1) # TRUE

final_ndvi <- as.data.frame(coordinates_ndvi)
final_evi  <- as.data.frame(coordinates_evi)


# 4. Create a loop to extract all point values - 2019 - for the 2 variables

for(i in 1:2){
  
  print(paste0(variables[i], " ", Sys.time()))
  
  if(variables[i] == "NDVI"){
    # Extract point values
    ndvi_data <- raster::extract(raster_ndvi, coordinates_ndvi)
    colnames(ndvi_data)[2:13] <- paste0("ndvi_", colnames(ndvi_data)[2:13])
    final_ndvi <- cbind(final_ndvi, ndvi_data[,-1])
    
    print(paste0(variables[i], "...", Sys.time(), " completed!"))
    
    
  } 
  
  
  if(variables[i] == "EVI"){
    
    # Extract point values
    evi_data <- raster::extract(raster_evi, coordinates_evi)
    colnames(evi_data)[2:13] <- paste0("evi_", colnames(evi_data)[2:13])
    final_evi <- cbind(final_evi, evi_data[,-1])
    
    print(paste0(variables[i], "...", Sys.time(), " completed!"))
    gc()
    
    
  }
  
  
} 

head(final_ndvi)
head(final_evi)

coordinates_ndvi <- unique(coordinates_ndvi)
head(coordinates_ndvi)
length(unique(coordinates_ndvi$geometry)) 
head(coordinates_ndvi) 

coordinates_evi <- unique(coordinates_evi)
head(coordinates_evi) 

write_sf(coordinates_ndvi, "centroids_1km.shp")

# 5. Long format

final_ndvi2 <- final_ndvi %>% 
  pivot_longer(!c(id, geometry), names_to = "label", values_to = "ndvi")
final_ndvi2 <- final_ndvi2 %>% dplyr::mutate(mm = substr(label, 12, nchar(label)))
final_ndvi2 <- final_ndvi2 %>% dplyr::select(id, mm, ndvi)

final_evi2 <- final_evi %>% 
  pivot_longer(!c(id, geometry), names_to = "label", values_to = "evi")
final_evi2 <- final_evi2 %>% dplyr::mutate(mm = substr(label, 11, nchar(label)))
final_evi2 <- final_evi2 %>% dplyr::select(id, mm, evi)

head(final_ndvi2)
setDT(final_ndvi2)
setkey(final_ndvi2, id, mm)
head(final_evi2)
setDT(final_evi2)
setkey(final_evi2, id, mm)

final2 <- merge(final_ndvi2, final_evi2)
head(final2)
gc()

summary(final2) 

final2 <- final2 %>% dplyr::mutate(EVI  = evi * 0.0001,
                                   NDVI = ndvi * 0.0001)
summary(final2) 

missing_ndvi <- final2 %>% dplyr::filter(is.na(NDVI))
missing_evi  <- final2 %>% dplyr::filter(is.na(EVI))

means <- final2 %>% dplyr::group_by(id) %>% dplyr::summarize(ndvi.final = mean(NDVI, na.rm=TRUE), evi.final = mean(EVI, na.rm=TRUE))
head(means)

setDT(final2)
setkey(final2, id)

setDT(means)
setkey(means, id)

final.2 <- merge(final2, means, all.x = T)
head(final.2)
gc()

final.2 <- final.2 %>% dplyr::mutate(ndvi_assigned = ifelse(is.na(NDVI), 1, 0),
                                     evi_assigned  = ifelse(is.na(EVI), 1, 0),
                                     NDVI = ifelse(is.na(NDVI), ndvi.final, NDVI),
                                     EVI  = ifelse(is.na(EVI),  evi.final,  EVI))
head(final.2)
final.2 <- final.2 %>% dplyr::select(!c(ndvi.final, evi.final, ndvi, evi))
summary(final.2)

final.2 <- final.2 %>% dplyr::mutate(mm = as.numeric(mm))
head(final.2)
gc()

saveRDS(final.2, "daily_ndvi_evi_1km_2019.RDS")

keep(monitors, monitors_repr, monitors_crs,
     final.2, coordinates, coordinates_ndvi, sure = T)
gc()


# Now I create the dataset specific for the monitors (nearest)

monitors_buffer <- st_buffer(monitors_repr, dist = 10000)
# plot(raster_evi$layer_1)
plot(st_geometry(monitors_buffer), fill = "red", add = T)

coordinates_ndvi <- st_filter(coordinates_ndvi, monitors_buffer) # for making the following command faster!

st_crs(monitors_repr) == st_crs(coordinates_ndvi)

nearest                    <- st_nearest_feature(monitors_repr, coordinates_ndvi)
monitors_repr$id           <- coordinates_ndvi$id[nearest]
head(monitors_repr)
head(monitors_repr)


# Complete dataset

head(final.2)
head(monitors_repr)

total <- monitors_repr[rep(seq_len(nrow(monitors_repr)), times = 12), ]
total$mm <- rep(1:12, each = nrow(monitors_repr))

setDT(final.2)
setkey(final.2, id, mm)
setDT(total)
setkey(total, id, mm)

total <- merge(total, final.2, all.x = T)
dim(total) 
head(total)

saveRDS(total, "daily_ndvi_evi_monitors_2019.RDS")
