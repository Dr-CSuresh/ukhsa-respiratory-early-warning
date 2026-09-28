# 1. LOAD PACKAGES

library(dplyr)
library(readr)
library(lubridate)
library(forecast)


# 2. LOAD CLEAN RESPIRATORY DATA

respiratory <- read_csv(
  "data/processed/respiratory_surveillance_long.csv",
  show_col_types = FALSE
) |>
  mutate(
    date = as.Date(date)
  )


# 3. FUNCTION TO GENERATE CURRENT FORECAST

forecast_current <- function(data) {
  
  data <- data |>
    arrange(date)
  
  disease_name <- unique(data$disease)
  
  y <- ts(
    data$positivity,
    frequency = 52
  )
  
  model <- auto.arima(
    y,
    seasonal = TRUE,
    stepwise = TRUE,
    approximation = TRUE
  )
  
  fc <- forecast(
    model,
    h = 4,
    level = c(80, 95)
  )
  
  last_date <- max(data$date)
  
  forecast_dates <- seq(
    last_date + weeks(1),
    by = "1 week",
    length.out = 4
  )
  
  tibble(
    disease = disease_name,
    date = forecast_dates,
    horizon = 1:4,
    
    forecast = as.numeric(
      fc$mean
    ),
    
    lower_80 = as.numeric(
      fc$lower[, "80%"]
    ),
    
    upper_80 = as.numeric(
      fc$upper[, "80%"]
    ),
    
    lower_95 = as.numeric(
      fc$lower[, "95%"]
    ),
    
    upper_95 = as.numeric(
      fc$upper[, "95%"]
    )
  )
}


# 4. GENERATE FORECASTS FOR ALL PATHOGENS

current_forecasts <- respiratory |>
  group_by(disease) |>
  group_split() |>
  lapply(forecast_current) |>
  bind_rows()


# 5. CHECK RESULTS

print(current_forecasts)


# 6. SAVE CURRENT FORECASTS

write_csv(
  current_forecasts,
  "data/processed/current_forecasts.csv"
)


# 7. CONFIRM

cat(
  "\nCurrent forecasts saved successfully.\n",
  "Forecast origin:",
  as.character(max(respiratory$date)),
  "\nPathogens:",
  n_distinct(current_forecasts$disease),
  "\nForecast rows:",
  nrow(current_forecasts),
  "\n"
)