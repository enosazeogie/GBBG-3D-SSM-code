## process drone data to estimate altitude error
##
require(tidyverse)
require(MASS)

fs <- list.files("data/drone", pattern = "UAV")
## remove apparently spurious data sets
fs <- fs[-c(3,5,6,9,10)]
dat <- lapply(1:length(fs), function(i) {
  x <- read_delim(paste0("data/drone/", fs[i]), delim = "\t", skip = 2) 
  if(nrow(x) > 0) {
    x |>
      mutate(date = lubridate::ymd_hms(paste(`#gmt_date`, `gmt_time`), tz="UTC")) |>
      select(date, lon=longitude, lat=latitude, galt=gps_altitude_m, lalt = laser_altitude_cm, num_sats) |>
      mutate(lon = as.numeric(lon), 
             lat = as.numeric(lat), 
             galt = as.numeric(galt), 
             lalt = as.numeric(lalt)) |>
      ## convert laser alt from cm to m
      mutate(lalt = lalt/100) |>
      ## correct GPS alt for height above geoid (local MSL)
      ##  value from Lato et al. 2022
      mutate(galtc = galt - 51.45) |>
      mutate(id = i) |>
      select(id, everything())
  }
  }) |> bind_rows() |>
  filter(date > ISOdate(2013,01,01))
  
## calculate raw GPS error, removing apparently spurious laser alt measurements 
raw.err <- with(subset(dat, lalt < 100), MASS::fitdistr(galt - lalt, "normal"))

## calculate corrected GPS error - where correction is for height above geoid
c.err <- with(subset(dat, lalt < 100), MASS::fitdistr(galtc - lalt, "normal"))

           
