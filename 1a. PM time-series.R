#### European-PM2.5-spatio-temporal-modelling

#### Author: Federica Nobile
#### Affiliation: Institute of Environmental Medicine, Karolinska Institutet, Stockholm, Sweden
####              Department of Epidemiology, Lazio Region Health Service/ASL Roma 1, Rome, Italy

library(arrow)
library(dplyr)
library(purrr)
library(data.table)
library(sf)
library(tmap)
library(gdata)

#### PM2.5 - 2019 --------------------------------------------------------------

files <- list.files("./ParquetFiles_PM25_Ea1/", full.names = TRUE)

files <- files[-1]

pm25_2019 <- map_dfr(files, function(f) {

  open_dataset(f) %>%
    dplyr::mutate(year = year(Start),
                  file_name = basename(f)) %>%
    dplyr::filter(Validity == 1,
                  year == 2019) %>%
    collect()

})
pm25_2019 <- pm25_2019 %>% filter(Validity == 1)
pm25_2019 <- as.data.frame(pm25_2019)
head(pm25_2019)
saveRDS(pm25_2019, "./pm25_2019_appoggio.RDS")

pm25_2019 <- pm25_2019 %>% mutate(CountryName = substr(Samplingpoint, 1, 2))
length(unique(pm25_2019$file_name)) 
length(unique(pm25_2019$Samplingpoint))

#### Monitors ####

monitors <- read.csv("./DataExtract.csv")
head(monitors)
monitors <- monitors %>%
  dplyr::filter(Air.Pollutant == "PM2.5") %>%
  dplyr::mutate(date_start = as.Date(substr(Operational.Activity.Begin, 1, 10), format = "%d/%m/%Y"),
                date_end   = as.Date(substr(Operational.Activity.End, 1, 10), format = "%d/%m/%Y"))
  dplyr::select(Country, Air.Quality.Network, Air.Quality.Network.Name, Air.Quality.Station.EoI.Code, Air.Quality.Station.Nat.Code, Air.Quality.Station.Name, Sampling.Point.Id, Longitude, Latitude, Altitude, Air.Quality.Station.Area, Air.Quality.Station.Type, date_start, date_end)
monitors <- unique(monitors)
length(unique(monitors$Sampling.Point.Id)) 
length(unique(monitors$Air.Quality.Station.EoI.Code)) 

#### Linkage ####

head(pm25_2019)
pm25_2019 <- pm25_2019 %>% dplyr::mutate('Sampling.Point.Id' = substr(file_name, 1, nchar(file_name)-8))
setDT(pm25_2019)
setkey(pm25_2019, Sampling.Point.Id)

head(monitors)
setDT(monitors)
setkey(monitors, Sampling.Point.Id)

final_pm25 <- merge(pm25_2019, monitors, all.x = T)
head(final_pm25)
summary(final_pm25)
length(unique(final_pm25$Sampling.Point.Id)) 
length(unique(final_pm25$Air.Quality.Station.EoI.Code)) 

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

saveRDS(time_series_pm25, "./timeseries_2019_pm25.RDS")

rm(list=ls())
gc()

#### PM10 - 2019 --------------------------------------------------------------

files <- list.files("./ParquetFiles_PM10_Ea1/", full.names = TRUE)

files <- files[-1]

pm10_2019 <- map_dfr(files, function(f) {

  open_dataset(f) %>%
    dplyr::mutate(year = year(Start),
                  file_name = basename(f)) %>%
    dplyr::filter(Validity == 1,
                  year == 2019) %>%
    collect()

})
pm10_2019 <- pm10_2019 %>% filter(Validity == 1)
pm10_2019 <- as.data.frame(pm10_2019)
head(pm10_2019)
saveRDS(pm10_2019, "./pm10_2019_appoggio.RDS")

#### Monitors ####

monitors <- read.csv("./DataExtract.csv")
head(monitors)
monitors <- monitors %>%
  dplyr::filter(Air.Pollutant == "PM10") %>%
  dplyr::mutate(date_start = as.Date(substr(Operational.Activity.Begin, 1, 10), format = "%d/%m/%Y"),
                date_end   = as.Date(substr(Operational.Activity.End, 1, 10), format = "%d/%m/%Y"))
  dplyr::select(Country, Air.Quality.Network, Air.Quality.Network.Name, Air.Quality.Station.EoI.Code, Air.Quality.Station.Nat.Code, Air.Quality.Station.Name, Sampling.Point.Id, Longitude, Latitude, Altitude, Air.Quality.Station.Area, Air.Quality.Station.Type, date_start, date_end)
monitors <- unique(monitors)
length(unique(monitors$Sampling.Point.Id))
length(unique(monitors$Air.Quality.Station.EoI.Code))

#### Linkage ####

head(pm10_2019)
pm10_2019 <- pm10_2019 %>% dplyr::mutate('Sampling.Point.Id' = substr(file_name, 1, nchar(file_name)-8))
setDT(pm10_2019)
setkey(pm10_2019, Sampling.Point.Id)

head(monitors)
setDT(monitors)
setkey(monitors, Sampling.Point.Id)

final_pm10 <- merge(pm10_2019, monitors, all.x = T)
head(final_pm10)
summary(final_pm10)
length(unique(final_pm10$Sampling.Point.Id)) # 3260
length(unique(final_pm10$Air.Quality.Station.EoI.Code)) # 2601

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

saveRDS(time_series_pm10, "./timeseries_2019_pm10.RDS")

keep(time_series_pm10, sure = T)
gc()

time_series_pm25 <- readRDS("./timeseries_2019_pm25.RDS")

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

sf::write_sf(monitors_sf, "./monitors_2019.shp")
