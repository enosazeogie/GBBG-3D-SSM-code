library(TMB)
require(tidyverse)
require(sf)


setwd(file.path("C:/Users/User/Documents/2024_Gulls/GullSSM/gps_alt_Ian Jonsen_rw/gps_alt"))
#setwd((file.path("C:\\Users\\IanMaywar\\Documents\\gps_alt\\gps_alt\\")))

source("fn/ssm_filter_crw.R")
source("fn/ssm_filter.R")

d <- read_csv("data/pre_proc/cleaned_alltags_Mar11_24.csv")

names(d) <- tolower(names(d))
d <- d |>
  mutate(lc = "G") |>
  select(id=device_id,
         date=utc_datetime,
         lc,
         lon=longitude,
         lat=latitude,
         z=altitude_m,
         hdop) |>
  filter(lon <0, lat > 0) |>
  mutate(id = as.character(id))

prj <- "+proj=merc +units=km +datum=WGS84"

d <- d |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326) |>
  st_transform(crs = prj)

dxy <- d |> 
  st_coordinates() |>
  as.data.frame() |>
  rename(x = X, y = Y)

d <- bind_cols(d, dxy) |>
  select(id, date, x, y, z, hdop)
d$geometry <- NULL

## Fit to a single track with 20 min time step & t-distribution df = 2
fit.rw <- ssm_filter(subset(d, id == 233345 & z > -20), 
                  time.step = 0.5, 
                  df = 4)

fit.crw <- ssm_filter_crw(subset(d, id == 233345 & z > -20),
                     time.step = 0.5, 
                     df = 4)

ggplot(fit.rw$predicted, aes(x=date, y=z, ymin = z - z.se * 1.96, ymax = z + z.se * 1.96)) +
  geom_point(data = fit.rw$data, aes(date, z), inherit.aes = FALSE) +
  geom_ribbon(fill = "red", alpha = 0.25) +
  geom_line(linewidth = 1, col = 'red') +
  geom_ribbon(data = fit.crw$predicted, aes(x = date, y = z, min = z - z.se * 1.96, ymax = z + z.se * 1.96), fill = "blue", alpha = 0.25) +
  geom_line(data = fit.crw$predicted, aes(x = date, y = z), linewidth = 1, col = 'blue') +
  ylim(-20,50)


##########troubleshooting lack of time matches with raw and ssm

sum(d$id == "233345")  # Count rows where 'variable' equals 'value'
subset_data <- subset(d, id == "233345")

ssm_subset_rw <- data.frame(
  z = fit.rw$fitted$z, 
  date = as.POSIXct(fit.rw$fitted$date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
)

colnames(ssm_subset_rw) <- c("z_rw", "date")  # Rename columns

# Perform a full join to retain all dates from both datasets
merged_data <- full_join(subset_data, ssm_subset_rw, by = "date")

# Identify matches and gaps
matches <- merged_data[complete.cases(merged_data), ]  # Rows where both datasets have values
gaps <- merged_data[!complete.cases(merged_data), ]  # Rows where one dataset is missing values

# View results
head(matches)  # Preview matching dates
head(gaps)  # Preview gaps

subset_data <- subset_data[order(subset_data$date), ]  # Ensure data is sorted by date
subset_data$interval <- c(NA, diff(subset_data$date))  # Compute time differences

ssm_subset_rw <- ssm_subset_rw[order(ssm_subset_rw$date), ]  # Ensure data is sorted by date
ssm_subset_rw$interval <- c(NA, diff(ssm_subset_rw$date))  # Compute time differences


merged_data <- merged_data[order(merged_data$date), ]  # Ensure data is sorted by date
merged_data$interval <- c(NA, diff(merged_data$date))  # Compute time differences

####interpolation not working
library(dplyr)
library(zoo)

m_data <- ssm_subset_rw %>% 
  left_join(subset_data, by = "date") %>%  # Merge datasets by date
  arrange(date) %>%
  mutate(z_interp = na.approx(z_rw, na.rm = FALSE))  # Interpolate missing values
