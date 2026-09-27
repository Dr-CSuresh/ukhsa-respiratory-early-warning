# UKHSA Respiratory Virus Surveillance
# Respiratory virus forecasting


# 1. PACKAGES

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(lubridate)
library(forecast)
library(purrr)


# 2. IMPORT CLEAN DATA

respiratory <- read_csv(
  "data/processed/respiratory_surveillance_long.csv",
  show_col_types = FALSE
)

glimpse(respiratory)


# 3. CREATE RESPIRATORY SEASON

respiratory <- respiratory |>
  mutate(
    month = month(date),
    
    calendar_year = year(date),
    
    season_start_year = if_else(
      month >= 7,
      calendar_year,
      calendar_year - 1
    ),
    
    season = paste0(
      season_start_year,
      "/",
      substr(
        season_start_year + 1,
        3,
        4
      )
    )
  )


# 4. DEFINE TRAINING AND TEST PERIODS

# We will train the models using data before 2025/26.
#
# The complete 2025/26 season will be hidden from the models
# and used to test how well they forecast unseen data.

training_data <- respiratory |>
  filter(
    date < as.Date("2025-07-01")
  )

test_data <- respiratory |>
  filter(
    date >= as.Date("2025-07-01"),
    date < as.Date("2026-07-01")
  )


# 5. CHECK TRAINING AND TEST DATA

training_summary <- training_data |>
  group_by(disease) |>
  summarise(
    observations = n(),
    first_date = min(date),
    last_date = max(date),
    .groups = "drop"
  )

print(training_summary)


test_summary <- test_data |>
  group_by(disease) |>
  summarise(
    observations = n(),
    first_date = min(date),
    last_date = max(date),
    .groups = "drop"
  )

print(test_summary)


# 6. CREATE LIST OF VIRUSES

viruses <- sort(
  unique(
    respiratory$disease
  )
)

print(viruses)


# 7. CREATE FORECASTING FUNCTION

forecast_virus <- function(
    virus_name
) {
  
  message(
    "\nForecasting: ",
    virus_name
  )
  
  
  # ----------------------------------------
  # Extract training data
  # ----------------------------------------
  
  train <- training_data |>
    filter(
      disease == virus_name
    ) |>
    arrange(date)
  
  
  # ----------------------------------------
  # Extract test data
  # ----------------------------------------
  
  test <- test_data |>
    filter(
      disease == virus_name
    ) |>
    arrange(date)
  
  
  # ----------------------------------------
  # Convert training data to time series
  # ----------------------------------------
  
  train_ts <- ts(
    train$positivity,
    frequency = 52
  )
  
  
  # ----------------------------------------
  # Forecast horizon
  # ----------------------------------------
  
  h <- nrow(test)
  
  
  # ----------------------------------------
  # MODEL 1:
  # Seasonal naïve
  # ----------------------------------------
  
  snaive_model <- snaive(
    train_ts,
    h = h
  )
  
  
  # ----------------------------------------
  # MODEL 2:
  # ETS
  # ----------------------------------------
  
  ets_fit <- ets(
    train_ts
  )
  
  ets_model <- forecast(
    ets_fit,
    h = h
  )
  
  
  # ----------------------------------------
  # MODEL 3:
  # ARIMA
  # ----------------------------------------
  
  arima_fit <- auto.arima(
    train_ts,
    seasonal = TRUE
  )
  
  arima_model <- forecast(
    arima_fit,
    h = h
  )
  
  
  # ----------------------------------------
  # Create forecast table
  # ----------------------------------------
  
  forecast_table <- tibble(
    disease = virus_name,
    
    date = test$date,
    
    observed = test$positivity,
    
    seasonal_naive =
      as.numeric(
        snaive_model$mean
      ),
    
    ETS =
      as.numeric(
        ets_model$mean
      ),
    
    ARIMA =
      as.numeric(
        arima_model$mean
      )
  )
  
  
  return(
    forecast_table
  )
}


# 8. RUN FORECASTS FOR ALL VIRUSES

forecast_results <- map_dfr(
  viruses,
  forecast_virus
)


# 9. INSPECT FORECAST RESULTS

glimpse(
  forecast_results
)

print(
  forecast_results,
  n = 30
)


# 10. CONVERT FORECASTS TO LONG FORMAT

forecast_long <- forecast_results |>
  pivot_longer(
    cols = c(
      seasonal_naive,
      ETS,
      ARIMA
    ),
    
    names_to = "model",
    
    values_to = "forecast"
  )


# 11. CALCULATE FORECAST ERRORS

forecast_long <- forecast_long |>
  mutate(
    error =
      observed - forecast,
    
    absolute_error =
      abs(error),
    
    squared_error =
      error^2
  )


# 12. CALCULATE MODEL ACCURACY

accuracy_results <- forecast_long |>
  group_by(
    disease,
    model
  ) |>
  summarise(
    MAE = mean(
      absolute_error,
      na.rm = TRUE
    ),
    
    RMSE = sqrt(
      mean(
        squared_error,
        na.rm = TRUE
      )
    ),
    
    .groups = "drop"
  ) |>
  arrange(
    disease,
    MAE
  )

print(
  accuracy_results,
  n = Inf
)


# 13. IDENTIFY LOWEST-MAE MODEL

# This is descriptive:
# it identifies which model had the lowest MAE
# for each virus in this particular held-out season.

lowest_mae_model <- accuracy_results |>
  group_by(disease) |>
  slice_min(
    MAE,
    n = 1,
    with_ties = FALSE
  ) |>
  ungroup()

print(
  lowest_mae_model,
  n = Inf
)


# 14. COMPARE MODEL ERRORS

model_comparison <- accuracy_results |>
  pivot_wider(
    names_from = model,
    
    values_from = c(
      MAE,
      RMSE
    )
  )

print(
  model_comparison,
  n = Inf
)


# 15. PLOT OBSERVED VS FORECASTS

forecast_plot <- ggplot(
  forecast_long,
  aes(
    x = date
  )
) +
  
  geom_line(
    aes(
      y = observed
    ),
    linewidth = 0.9
  ) +
  
  geom_line(
    aes(
      y = forecast,
      linetype = model
    ),
    linewidth = 0.7
  ) +
  
  facet_wrap(
    ~ disease,
    scales = "free_y",
    ncol = 2
  ) +
  
  labs(
    title =
      "Respiratory Virus Forecasting in England",
    
    subtitle =
      "Forecasts generated from pre-July 2025 data and evaluated against the 2025/26 season",
    
    x = NULL,
    
    y =
      "PCR positivity (%)",
    
    linetype =
      "Forecast model",
    
    caption =
      "Models: seasonal naïve, ETS and ARIMA\nSource: UK Health Security Agency"
  ) +
  
  theme_minimal()

forecast_plot


# 16. PLOT FORECAST ERRORS

accuracy_plot <- ggplot(
  accuracy_results,
  aes(
    x = model,
    y = MAE
  )
) +
  
  geom_col() +
  
  facet_wrap(
    ~ disease,
    scales = "free_y",
    ncol = 2
  ) +
  
  labs(
    title =
      "Forecast Accuracy by Respiratory Virus",
    
    subtitle =
      "Mean absolute error during the held-out 2025/26 respiratory season",
    
    x =
      "Forecast model",
    
    y =
      "Mean absolute error",
    
    caption =
      "Lower MAE indicates smaller forecast error\nSource: UK Health Security Agency"
  ) +
  
  theme_minimal()

accuracy_plot


# 17. CALCULATE ERROR BY FORECAST HORIZON

# This allows us to see whether forecasts become
# progressively worse further into the future.

horizon_results <- forecast_long |>
  group_by(
    disease,
    model
  ) |>
  arrange(date) |>
  mutate(
    forecast_week =
      row_number()
  ) |>
  ungroup()


# 18. PLOT ERROR BY FORECAST HORIZON

horizon_plot <- ggplot(
  horizon_results,
  aes(
    x = forecast_week,
    y = absolute_error,
    linetype = model
  )
) +
  
  geom_line() +
  
  facet_wrap(
    ~ disease,
    scales = "free_y",
    ncol = 2
  ) +
  
  labs(
    title =
      "Forecast Error by Forecast Horizon",
    
    subtitle =
      "Absolute error across the held-out 2025/26 respiratory season",
    
    x =
      "Weeks ahead",
    
    y =
      "Absolute forecast error",
    
    linetype =
      "Forecast model",
    
    caption =
      "Source: UK Health Security Agency"
  ) +
  
  theme_minimal()

horizon_plot


# 19. CREATE OUTPUT FOLDERS

dir.create(
  "outputs/figures",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "outputs/tables",
  recursive = TRUE,
  showWarnings = FALSE
)


# 20. SAVE FIGURES

ggsave(
  "outputs/figures/forecast_comparison.png",
  forecast_plot,
  width = 12,
  height = 10,
  dpi = 300
)

ggsave(
  "outputs/figures/forecast_accuracy.png",
  accuracy_plot,
  width = 12,
  height = 9,
  dpi = 300
)

ggsave(
  "outputs/figures/forecast_error_by_horizon.png",
  horizon_plot,
  width = 12,
  height = 9,
  dpi = 300
)


# 21. SAVE TABLES

write_csv(
  forecast_results,
  "outputs/tables/forecast_results.csv"
)

write_csv(
  accuracy_results,
  "outputs/tables/forecast_accuracy.csv"
)

write_csv(
  lowest_mae_model,
  "outputs/tables/lowest_mae_model.csv"
)

write_csv(
  horizon_results,
  "outputs/tables/forecast_horizon_errors.csv"
)


# 22. PRINT FORECASTING SUMMARY

cat(
  "\n========================================\n"
)

cat(
  "FORECAST EVALUATION COMPLETE\n"
)

cat(
  "========================================\n"
)

cat(
  "Viruses evaluated:",
  n_distinct(
    forecast_results$disease
  ),
  "\n"
)

cat(
  "Test period:",
  as.character(
    min(
      forecast_results$date
    )
  ),
  "to",
  as.character(
    max(
      forecast_results$date
    )
  ),
  "\n"
)

cat(
  "Forecast models:",
  n_distinct(
    forecast_long$model
  ),
  "\n"
)

cat(
  "Test observations:",
  nrow(
    forecast_results
  ),
  "\n"
)

cat(
  "========================================\n"
)


# END