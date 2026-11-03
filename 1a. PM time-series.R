#### CLEANING DATA ####

# YEAR: 2018
# LAST UPDATE: 23/03/2026

.libPaths("C:/Users/stafoggia/Documents/R/win-library/4.1")

# install.packages("arrow", repos = c("https://apache.r-universe.dev", "https://cloud.r-project.org"))
library(arrow)
library(dplyr)
library(purrr)
library(data.table)
library(sf)
library(tmap)
library(gdata)

#### PM2.5 - 2018 --------------------------------------------------------------

# # I need to run the following commands in local
# 
# files <- list.files("Z:/meteo_condivisa/NOBILE/PhD KI - Study II/ParquetFiles/ParquetFiles_PM25_Ea1/", full.names = TRUE) # 3147
# 
# files <- files[-1] # 3146
# 
# pm25_2018 <- map_dfr(files, function(f) {
# 
#   open_dataset(f) %>%
#     dplyr::mutate(year = year(Start),
#                   file_name = basename(f)) %>%
#     dplyr::filter(Validity == 1,
#                   year == 2018) %>%
#     collect()
# 
# })
# pm25_2018 <- pm25_2018 %>% filter(Validity == 1)
# pm25_2018 <- as.data.frame(pm25_2018)
# head(pm25_2018)
# saveRDS(pm25_2018, "Z:/meteo_condivisa/NOBILE/PhD KI - Study II/2018/AIR POLL OBS/pm25_2018_appoggio.RDS")

pm25_2018 <- readRDS("Z:/meteo_condivisa/NOBILE/PhD KI - Study II/2018/AIR POLL OBS/pm25_2018_appoggio.RDS")
pm25_2018 <- pm25_2018 %>% mutate(CountryName = substr(Samplingpoint, 1, 2))
length(unique(pm25_2018$file_name)) # 1615
length(unique(pm25_2018$Samplingpoint)) # 1615

#### Monitors ####

monitors <- read.csv("Z:/meteo_condivisa/NOBILE/PhD KI - Study II/ParquetFiles/DataExtract.csv") # 79115
head(monitors)
monitors <- monitors %>%
  dplyr::filter(Air.Pollutant == "PM2.5") %>%
  dplyr::mutate(date_start = as.Date(substr(Operational.Activity.Begin, 1, 10), format = "%d/%m/%Y"),
                date_end   = as.Date(substr(Operational.Activity.End, 1, 10), format = "%d/%m/%Y")) %>%
  # dplyr::filter(date_start <= as.Date("2018-12-31") & (date_end >= as.Date("2018-01-01") | is.na(date_end) | date_end == "")) %>%
  dplyr::select(Country, Air.Quality.Network, Air.Quality.Network.Name, Air.Quality.Station.EoI.Code, Air.Quality.Station.Nat.Code, Air.Quality.Station.Name, Sampling.Point.Id, Longitude, Latitude, Altitude, Air.Quality.Station.Area, Air.Quality.Station.Type, date_start, date_end)
monitors <- unique(monitors) # 3910
length(unique(monitors$Sampling.Point.Id)) # 3910
length(unique(monitors$Air.Quality.Station.EoI.Code)) # 3417

#### Linkage ####

head(pm25_2018)
pm25_2018 <- pm25_2018 %>% dplyr::mutate('Sampling.Point.Id' = substr(file_name, 1, nchar(file_name)-8))
setDT(pm25_2018)
setkey(pm25_2018, Sampling.Point.Id)

head(monitors)
setDT(monitors)
setkey(monitors, Sampling.Point.Id)

final_pm25 <- merge(pm25_2018, monitors, all.x = T)
head(final_pm25)
summary(final_pm25)
length(unique(final_pm25$Sampling.Point.Id)) # 1615
length(unique(final_pm25$Air.Quality.Station.EoI.Code)) # 1344

vedi <- final_pm25 %>% filter(is.na(Air.Quality.Station.EoI.Code)) # 238,240
length(unique(vedi$Sampling.Point.Id)) # 215
table(vedi$CountryName)
# IT 
# 238240

final_pm25_ok <- final_pm25 %>% filter(!is.na(Air.Quality.Station.EoI.Code)) # 7,472,096

# There are some issues specific to Italy

vedi_IT <- pm25_2018 %>% filter(Sampling.Point.Id %in% vedi$Sampling.Point.Id & CountryName == "IT")
vedi_IT <- vedi_IT %>% dplyr::mutate(Air.Quality.Station.EoI.Code = substr(Sampling.Point.Id, 5, 11))
setkey(vedi_IT,  Air.Quality.Station.EoI.Code)

monitors_IT <- monitors %>% dplyr::filter(Country == "Italy") # 499
split_IT    <- split(monitors_IT, monitors_IT$Air.Quality.Station.EoI.Code)
monitors_IT <- monitors_IT %>%
  group_by(across(Country:Air.Quality.Station.Name)) %>%
  summarise(Longitude                = first(Longitude),
            Latitude                 = first(Latitude),
            Altitude                 = first(Altitude),
            Air.Quality.Station.Area = first(Air.Quality.Station.Area),
            Air.Quality.Station.Type = first(Air.Quality.Station.Type),
            Sampling.Point.Id        = first(Sampling.Point.Id),
            date_start               = first(date_start),
            date_end                 = first(date_end)) # 405
setDT(monitors_IT)
setkey(monitors_IT, Air.Quality.Station.EoI.Code)

vedi_IT <- merge(vedi_IT, monitors_IT, all.x = T)
summary(vedi_IT) # OK!!!!! :)

final_pm25_it <- vedi_IT %>% rename(Sampling.Point.Id = Sampling.Point.Id.x) %>% select(!(Sampling.Point.Id.y)) # 238,240
final_pm25 <- rbind(final_pm25_ok, final_pm25_it)
summary(final_pm25)
final_pm25 <- final_pm25 %>% mutate(Value = ifelse(Value < 0, 0, Value))

length(unique(final_pm25$Sampling.Point.Id)) # 1615
length(unique(final_pm25$Air.Quality.Station.EoI.Code)) # 1558

#### Final aggregation ####

time_series_new_day      <- final_pm25 %>% filter(AggType == "day") %>% rename(pm25 = Value) %>% mutate(n = 1)
head(time_series_new_day)
time_series_new_day$date <- as.Date(time_series_new_day$Start)

time_series_new_hh       <- final_pm25 %>% filter(!(AggType == "day"))
time_series_new_hh$date  <- as.Date(time_series_new_hh$Start)
time_series_new_hh2      <- time_series_new_hh %>%
  group_by(Air.Quality.Station.Nat.Code, Air.Quality.Station.EoI.Code, date) %>%
  summarise(
    n = n(),
    pm25 = mean(Value, na.rm = TRUE),
    across(
      setdiff(names(.), c("Air.Quality.Station.EoI.Code", "Air.Quality.Station.Nat.Code", "date", "Value")),
      ~ unique(.)[1]
    ),
    .groups = "drop"
  )
time_series_new_hh2 <- as.data.frame(time_series_new_hh2)

head(time_series_new_day)
head(time_series_new_hh2)
time_series_pm25 <- rbind(time_series_new_day, time_series_new_hh2) # 354,051

length(unique(time_series_pm25$Sampling.Point.Id)) # 1615
length(unique(time_series_pm25$Air.Quality.Station.EoI.Code)) # 1558
summary(time_series_pm25$pm25)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.000   6.921  10.750  13.723  16.681 324.298

saveRDS(time_series_pm25, "Z:/meteo_condivisa/NOBILE/PhD KI - Study II/2018/AIR POLL OBS/timeseries_2018_pm25.RDS")

rm(list=ls())
gc()

#### PM10 - 2018 --------------------------------------------------------------
# 
# # I need to run the following commands in local
# 
# files <- list.files("Z:/meteo_condivisa/NOBILE/PhD KI - Study II/ParquetFiles/ParquetFiles_PM10_Ea1/", full.names = TRUE) # 5234
# 
# files <- files[-1] # 5233
# 
# pm10_2018 <- map_dfr(files, function(f) {
# 
#   open_dataset(f) %>%
#     dplyr::mutate(year = year(Start),
#                   file_name = basename(f)) %>%
#     dplyr::filter(Validity == 1,
#                   year == 2018) %>%
#     collect()
# 
# })
# pm10_2018 <- pm10_2018 %>% filter(Validity == 1)
# pm10_2018 <- as.data.frame(pm10_2018)
# head(pm10_2018)
# saveRDS(pm10_2018, "Z:/meteo_condivisa/NOBILE/PhD KI - Study II/2018/AIR POLL OBS/pm10_2018_appoggio.RDS")

pm10_2018 <- readRDS("X:/NOBILE/PhD KI - Study II/2018/AIR POLL OBS/pm10_2018_appoggio.RDS")
pm10_2018 <- pm10_2018 %>% mutate(CountryName = substr(Samplingpoint, 1, 2))
length(unique(pm10_2018$file_name)) # 3260
length(unique(pm10_2018$Samplingpoint)) # 3260

#### Monitors ####

monitors <- read.csv("X:/NOBILE/PhD KI - Study II/ParquetFiles/DataExtract.csv") # 79115
head(monitors)
monitors <- monitors %>%
  dplyr::filter(Air.Pollutant == "PM10") %>%
  dplyr::mutate(date_start = as.Date(substr(Operational.Activity.Begin, 1, 10), format = "%d/%m/%Y"),
                date_end   = as.Date(substr(Operational.Activity.End, 1, 10), format = "%d/%m/%Y")) %>%
  # dplyr::filter(date_start <= as.Date("2018-12-31") & (date_end >= as.Date("2018-01-01") | is.na(date_end) | date_end == "")) %>%
  dplyr::select(Country, Air.Quality.Network, Air.Quality.Network.Name, Air.Quality.Station.EoI.Code, Air.Quality.Station.Nat.Code, Air.Quality.Station.Name, Sampling.Point.Id, Longitude, Latitude, Altitude, Air.Quality.Station.Area, Air.Quality.Station.Type, date_start, date_end)
monitors <- unique(monitors) # 6416
length(unique(monitors$Sampling.Point.Id)) # 6416
length(unique(monitors$Air.Quality.Station.EoI.Code)) # 5237

#### Linkage ####

head(pm10_2018)
pm10_2018 <- pm10_2018 %>% dplyr::mutate('Sampling.Point.Id' = substr(file_name, 1, nchar(file_name)-8))
setDT(pm10_2018)
setkey(pm10_2018, Sampling.Point.Id)

head(monitors)
setDT(monitors)
setkey(monitors, Sampling.Point.Id)

final_pm10 <- merge(pm10_2018, monitors, all.x = T)
head(final_pm10)
summary(final_pm10)
length(unique(final_pm10$Sampling.Point.Id)) # 3260
length(unique(final_pm10$Air.Quality.Station.EoI.Code)) # 2601

vedi <- final_pm10 %>% filter(is.na(Air.Quality.Station.EoI.Code)) # 635322
length(unique(vedi$Sampling.Point.Id)) # 450
table(vedi$CountryName)
# IT 
# 635322

final_pm10_ok <- final_pm10 %>% filter(!is.na(Air.Quality.Station.EoI.Code)) # 15,370,048

# There are some issues specific to Italy

vedi_IT <- pm10_2018 %>% filter(Sampling.Point.Id %in% vedi$Sampling.Point.Id & CountryName == "IT")
vedi_IT <- vedi_IT %>% dplyr::mutate(Air.Quality.Station.EoI.Code = substr(Sampling.Point.Id, 5, 11))
setkey(vedi_IT,  Air.Quality.Station.EoI.Code)

monitors_IT <- monitors %>% dplyr::filter(Country == "Italy") # 936
split_IT    <- split(monitors_IT, monitors_IT$Air.Quality.Station.EoI.Code)
monitors_IT <- monitors_IT %>%
  group_by(across(Country:Air.Quality.Station.Name)) %>%
  summarise(Longitude                = first(Longitude),
            Latitude                 = first(Latitude),
            Altitude                 = first(Altitude),
            Air.Quality.Station.Area = first(Air.Quality.Station.Area),
            Air.Quality.Station.Type = first(Air.Quality.Station.Type),
            Sampling.Point.Id        = first(Sampling.Point.Id),
            date_start               = first(date_start),
            date_end                 = first(date_end)) # 733
setDT(monitors_IT)
setkey(monitors_IT, Air.Quality.Station.EoI.Code)

vedi_IT <- merge(vedi_IT, monitors_IT, all.x = T)
summary(vedi_IT) # OK!!!!! :)

final_pm10_it <- vedi_IT %>% rename(Sampling.Point.Id = Sampling.Point.Id.x) %>% select(!(Sampling.Point.Id.y)) # 635,322
final_pm10 <- rbind(final_pm10_ok, final_pm10_it)
summary(final_pm10)
final_pm10 <- final_pm10 %>% mutate(Value = ifelse(Value < 0, 0, Value))

length(unique(final_pm10$Sampling.Point.Id)) # 3260
length(unique(final_pm10$Air.Quality.Station.EoI.Code)) # 3049

#### Final aggregation ####

time_series_new_day      <- final_pm10 %>% filter(AggType == "day") %>% rename(pm10 = Value) %>% mutate(n = 1)
head(time_series_new_day)
time_series_new_day$date <- as.Date(time_series_new_day$Start)

time_series_new_hh       <- final_pm10 %>% filter(!(AggType == "day"))
time_series_new_hh$date  <- as.Date(time_series_new_hh$Start)
time_series_new_hh2      <- time_series_new_hh %>%
  group_by(Air.Quality.Station.Nat.Code, Air.Quality.Station.EoI.Code, date) %>%
  summarise(
    n = n(),
    pm10 = mean(Value, na.rm = TRUE),
    across(
      setdiff(names(.), c("Air.Quality.Station.EoI.Code", "Air.Quality.Station.Nat.Code", "date", "Value")),
      ~ unique(.)[1]
    ),
    .groups = "drop"
  )
time_series_new_hh2 <- as.data.frame(time_series_new_hh2)

head(time_series_new_day)
head(time_series_new_hh2)
time_series_pm10 <- rbind(time_series_new_day, time_series_new_hh2) # 16,005,370

length(unique(time_series_pm10$Sampling.Point.Id)) # 3260
length(unique(time_series_pm10$Air.Quality.Station.EoI.Code)) # 3049
summary(time_series_pm10$pm10)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00   12.66   18.67   22.00   26.85  852.30 

saveRDS(time_series_pm10, "X:/NOBILE/PhD KI - Study II/2018/AIR POLL OBS/timeseries_2018_pm10.RDS")

keep(time_series_pm10, sure = T)
gc()

time_series_pm25 <- readRDS("X:/NOBILE/PhD KI - Study II/2018/AIR POLL OBS/timeseries_2018_pm25.RDS")

#### MONITORS SHAPE ------------------------------------------------------------

names(time_series_pm25)
monitors_pm25 <- time_series_pm25 %>% select(!c(Pollutant, Start, End, pm25, Unit, AggType, Validity, Verification, ResultTime, DataCapture, FkObservationLog, file_name, date_start, date_end, n, date))

monitors_pm10 <- time_series_pm10 %>% select(!c(Pollutant, Start, End, pm10, Unit, AggType, Validity, Verification, ResultTime, DataCapture, FkObservationLog, file_name, date_start, date_end, n, date))

monitors <- rbind(monitors_pm25, monitors_pm10)
monitors <- unique(monitors) # 4875

monitors_sf <- st_as_sf(
  monitors,
  coords = c("Longitude", "Latitude"), 
  crs = 4326   
)

tmap_mode("view")
tm_shape(monitors_sf) +
  tm_dots(size = 0.05, col = "red") +
  tm_basemap("OpenStreetMap")

sf::write_sf(monitors_sf, "X:/NOBILE/PhD KI - Study II/2018/AIR POLL OBS/monitors_2018.shp")
