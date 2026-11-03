#### European-PM2.5-spatio-temporal-modelling

#### Author: Federica Nobile
#### Affiliation: Institute of Environmental Medicine, Karolinska Institutet, Stockholm, Sweden
####              Department of Epidemiology, Lazio Region Health Service/ASL Roma 1, Rome, Italy

library(dplyr)
library(lubridate)
library(sf)
library(spatstat.geom)
library(ggplot2)
library(sp)

#### Call the dataset for stage 1

data <- readRDS("./stage1.PM.final.RDS")
dim(data) 
names(data)

#### Select variables of interest

data1 <- data %>%
  dplyr::select(day, dow, mm,
                Air.Quality.Station.EoI.Code, CountryName, Country, Air.Quality.Network,
                Air.Quality.Network.Name, Air.Quality.Station.Nat.Code, Air.Quality.Station.Name,
                Altitude, Air.Quality.Station.Area, Air.Quality.Station.Type, year,
                Longitude, Latitude,
                pm25, final.pm25, pm10, 
                allRoads_50, allRoads_100, allRoads_200, allRoads_300, allRoads_400, allRoads_500,
                majorRoads_50, majorRoads_100, majorRoads_200, majorRoads_300, majorRoads_400, majorRoads_500,
                imd_100, imd_200, imd_300, imd_400, imd_500,
                ind_100, ind_200, ind_300, ind_400, ind_500,
                nat_100, nat_200, nat_300, nat_400, nat_500,
                por_100, por_200, por_300, por_400, por_500,
                res_100, res_200, res_300, res_400, res_500,
                tbu_100, tbu_200, tbu_300, tbu_400, tbu_500,
                ugr_100, ugr_200, ugr_300, ugr_400, ugr_500,
                pop, zoneID, alt10_enh,
                u10, v10, dewp_2m, temp_2m, mslp, tp, blh, leaf_high, leaf_low, tcc, tco,
                pm25_cams, pm10_cams,
                aod1240_cams, aod469_cams, aod550_cams, aod670_cams, aod865_cams,
                NDVI, EVI,
                final.aod047)

data1 <- data1 %>% dplyr::mutate(ln.pm25 = log(final.pm25 + 0.5),
                                 mm      = as.factor(mm),
                                 dow     = as.factor(dow)) 
summary(data1)

names(data1)
sites.pm           <- data1 %>% dplyr::select(Air.Quality.Station.EoI.Code, CountryName, Longitude, Latitude)
sites.pm           <- unique(sites.pm)
table(sites.pm$CountryName)

sites.pm.shape <- st_as_sf(sites.pm, coords = c("Longitude", "Latitude"), crs = "EPSG:4326")
plot(sites.pm.shape$geometry)

# 1. Transform to a projected CRS (Europe Lambert Equal Area - EPSG: 3035) for distance calculations
sf_points_proj <- st_transform(sites.pm.shape, 3035)
plot(sf_points_proj$geometry)

# 2. Define function to create 10 balanced (in terms of size) groups ensuring spatial diversity
create_final_balanced_groups <- function(sf_points_proj, num_groups = 10, min_dist_km = 10, max_attempts_per_group = 500) {
  
  set.seed(42) # For reproducibility
  
  df <- sf_points_proj %>%
    mutate(x = st_coordinates(.)[,1], y = st_coordinates(.)[,2], 
           lat = st_coordinates(st_transform(sf_points_proj, 4326))[,2], 
           lon = st_coordinates(st_transform(sf_points_proj, 4326))[,1]) %>%
    st_drop_geometry()
  
  df$group <- NA  # Initialize group column
  
  # Divide Europe into latitude-longitude bins (stratified sampling)
  lat_bins <- cut(df$lat, breaks = 10)  # Divide latitude into bins
  lon_bins <- cut(df$lon, breaks = 10)  # Divide longitude into bins
  df$lat_bin <- as.numeric(lat_bins)
  df$lon_bin <- as.numeric(lon_bins)
  
  # Number of points per group
  points_per_group <- ceiling(nrow(df) / num_groups)
  
  for (g in 1:num_groups) {
    
    available_points <- df %>% filter(is.na(group))
    if (nrow(available_points) == 0) break
    
    # Select a first point ensuring spread across Europe
    first_point <- available_points %>%
      group_by(lat_bin, lon_bin) %>%
      sample_n(1, replace = TRUE) %>%
      ungroup() %>%
      sample_n(1)
    
    print(paste0(first_point$Cntrycd, " ", first_point$ArQSEIC, " ...first point -> ", Sys.time()))
    
    group_points <- first_point
    df$group[df$x == first_point$x & df$y == first_point$y] <- g
    
    attempt_count <- 0
    
    while (nrow(group_points) < points_per_group && nrow(available_points) > 0) {
      
      candidate <- available_points %>% sample_n(1)  # Pick randomly
      
      # Compute distances to all points in the group
      distances <- sqrt((group_points$x - candidate$x)^2 + (group_points$y - candidate$y)^2)
      
      if (all(distances > min_dist_km * 1000)) {  # Convert km to meters
        group_points <- rbind(group_points, candidate)
        df$group[df$x == candidate$x & df$y == candidate$y] <- g
      }
      
      available_points <- df %>% filter(is.na(group))  # Update available points
      
      attempt_count <- attempt_count + 1
      
      # EARLY STOPPING CONDITION: If we exceed max_attempts_per_group, break and reassign points later
      if (attempt_count >= max_attempts_per_group) {
        message(sprintf("Group %d reached max attempts (%d). Stopping search and reassigning points later.", g, max_attempts_per_group))
        break
      }
    }
    
    print(paste0("Group", g, " ...completed -> ", Sys.time()))
    
  }
  
  # If some points remain unassigned, redistribute them across existing groups
  unassigned_points <- df %>% filter(is.na(group))
  
  if (nrow(unassigned_points) > 0) {
    message(sprintf("%d points could not be assigned due to distance constraints. Redistributing...", nrow(unassigned_points)))
    
    for (i in 1:nrow(unassigned_points)) {
      random_group <- sample(1:num_groups, 1)  # Assign to a random group
      df$group[df$x == unassigned_points$x[i] & df$y == unassigned_points$y[i]] <- random_group
    }
  }
  
  # **Balance groups: Ensure each group has similar size**
  group_counts <- df %>% count(group)
  max_group_size <- max(group_counts$n)
  min_group_size <- min(group_counts$n)
  
  while (min_group_size < max_group_size - 2) {  # Allow only small variation
    
    # Find the smallest and largest groups
    small_group <- group_counts$group[which.min(group_counts$n)]
    large_group <- group_counts$group[which.max(group_counts$n)]
    
    # Find a point from the largest group that can be moved
    candidate_move <- df %>% filter(group == large_group) %>% sample_n(1)
    
    # Move the point
    df$group[df$x == candidate_move$x & df$y == candidate_move$y] <- small_group
    
    # Update group sizes
    group_counts <- df %>% count(group)
    max_group_size <- max(group_counts$n)
    min_group_size <- min(group_counts$n)
  }
  
  return(df)
}

# Apply function
grouped_points <- create_final_balanced_groups(sf_points_proj)

# 3. Merge groups back to sf object and convert back to WGS84 (lat/lon)
sf_points_proj$group <- grouped_points$group
sf_points_final <- st_transform(sf_points_proj, 4326)

# 4. Save results
write.csv(st_drop_geometry(sf_points_final), "./spatial-blocked cv folds.csv", row.names = FALSE)