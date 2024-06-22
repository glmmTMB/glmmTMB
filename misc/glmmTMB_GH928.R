fn <- "wind_data.csv"
if (!file.exists(fn)) {
    download.file("https://github.com/glmmTMB/glmmTMB/files/14537932/wind_data.csv", dest = fn)
}

library(tidyverse)

wind_data_full <- read_csv(fn)
wind_data <- wind_data_full %>%
  filter(datetime_UTC >= ymd_hms('2007-09-10 00:00:00') & datetime_UTC <= ymd_hms('2007-10-15 23:59:59'))
  
wind_data <- wind_data %>%
  arrange(datetime_UTC) %>% # Make sure the data is sorted by datetime_UTC
  mutate(time = as.numeric(difftime(datetime_UTC, min(datetime_UTC), units = "hours")) + 1)

wind_data$datetime_UTC <- NULL # date type variables sometimes cause problems with predict()

set.seed(123) 
wind_data <- wind_data[order(wind_data$time), ]
total_rows <- nrow(wind_data)
train_size <- floor(0.6 * total_rows)
valid_size <- floor(0.2 * total_rows)
test_size <- total_rows - train_size - valid_size 

train_data <- wind_data[1:train_size, ]
valid_data <- wind_data[(train_size + 1):(train_size + valid_size), ]
test_data <- wind_data[(train_size + valid_size + 1):total_rows, ]

test_mod <- glmmTMB(spd ~ s(time) + s(U,V), data = train_data) 
predictions <- predict(test_mod, valid_data, type="response")

