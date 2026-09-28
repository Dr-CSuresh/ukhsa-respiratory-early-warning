# UKHSA Respiratory Virus Surveillance
# Rolling-origin forecasting


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


# 3. DEFINE FORECAST SETTINGS

forecast_start <- as.Date("2025-07-07")

forecast_end <- as.Date("2026-06-29")

maximum_horizon <- 4


# 4. CREATE ROLLING FORECAST FUNCTION

rolling_forecast_virus <- function(
    virus_name
) {
  
  message(
    "\n========================================"
  )
  
  message(
    "\nRolling forecasts: ",
    virus_name
  )
  
  message(
    "\n========================================"
  )
  
  
  # ----------------------------------------
  # Extract virus data
  # ----------------------------------------
  
  virus_data <- respiratory |>
    filter(
      disease == virus_name
    ) |>
    arrange(date)
  
  
  # ----------------------------------------
  # Identify forecast origins
  # ----------------------------------------
  
  forecast_origins <- virus_data |>
    filter(
      date >= forecast_start,
      date <= forecast_end
    ) |>
    pull(date)
  
  
  # ----------------------------------------
  # Run forecast from each origin
  # ----------------------------------------
  
  virus_results <- map_dfr(
    seq_along(forecast_origins),
    
    function(i) {
      
      origin_date <- forecast_origins[i]
      
      
      message(
        "  ",
        virus_name,
        " | origin ",
        i,
        "/",
        length(forecast_origins),
        " | ",
        origin_date
      )
      
      
      # ------------------------------------
      # Training data available at origin
      # ------------------------------------
      
      train <- virus_data |>
        filter(
          date <= origin_date
        )
      
      
      train_ts <- ts(
        train$positivity,
        frequency = 52
      )
      
      
      # ------------------------------------
      # Actual future observations
      # ------------------------------------
      
      future_data <- virus_data |>
        filter(
          date > origin_date
        ) |>
        slice_head(
          n = maximum_horizon
        )
      
      
      # Skip if fewer than 4 future weeks
      # are available
      
      if (
        nrow(future_data) <
        maximum_horizon
      ) {
        
        return(NULL)
        
      }
      
      
      # ------------------------------------
      # MODEL 1:
      # Seasonal naïve
      # ------------------------------------
      
      snaive_forecast <- snaive(
        train_ts,
        h = maximum_horizon
      )
      
      
      # ------------------------------------
      # MODEL 2:
      # ARIMA
      # ------------------------------------
      
      arima_fit <- auto.arima(
        train_ts,
        seasonal = TRUE,
        stepwise = TRUE,
        approximation = TRUE
      )
      
      
      arima_forecast <- forecast(
        arima_fit,
        h = maximum_horizon
      )
      
      
      # ------------------------------------
      # MODEL 3:
      # STL + ETS
      # ------------------------------------
      
      stl_forecast <- stlf(
        train_ts,
        h = maximum_horizon,
        method = "ets"
      )
      
      
      # ------------------------------------
      # Store results
      # ------------------------------------
      
      tibble(
        disease =
          virus_name,
        
        origin_date =
          origin_date,
        
        target_date =
          future_data$date,
        
        horizon =
          1:maximum_horizon,
        
        observed =
          future_data$positivity,
        
        seasonal_naive =
          as.numeric(
            snaive_forecast$mean
          ),
        
        ARIMA =
          as.numeric(
            arima_forecast$mean
          ),
        
        STL_ETS =
          as.numeric(
            stl_forecast$mean
          )
      )
      
    }
  )
  
  
  message(
    "Finished ",
    virus_name,
    " ✓"
  )
  
  
  return(
    virus_results
  )
}


# 5. CREATE LIST OF VIRUSES

viruses <- sort(
  unique(
    respiratory$disease
  )
)

print(viruses)


# 6. RUN ROLLING FORECASTS

rolling_results <- map_dfr(
  viruses,
  rolling_forecast_virus
)


# 7. INSPECT RESULTS

glimpse(
  rolling_results
)

print(
  rolling_results,
  n = 30
)


# 8. CONVERT TO LONG FORMAT

rolling_long <- rolling_results |>
  pivot_longer(
    cols = c(
      seasonal_naive,
      ARIMA,
      STL_ETS
    ),
    
    names_to =
      "model",
    
    values_to =
      "forecast"
  )


# 9. CALCULATE FORECAST ERRORS

rolling_long <- rolling_long |>
  mutate(
    error =
      observed - forecast,
    
    absolute_error =
      abs(error),
    
    squared_error =
      error^2
  )


# 10. CALCULATE ACCURACY BY HORIZON

accuracy_by_horizon <- rolling_long |>
  group_by(
    disease,
    model,
    horizon
  ) |>
  summarise(
    forecasts =
      n(),
    
    MAE =
      mean(
        absolute_error,
        na.rm = TRUE
      ),
    
    RMSE =
      sqrt(
        mean(
          squared_error,
          na.rm = TRUE
        )
      ),
    
    .groups =
      "drop"
  ) |>
  arrange(
    disease,
    horizon,
    MAE
  )

print(
  accuracy_by_horizon,
  n = Inf
)


# 11. CALCULATE OVERALL ACCURACY

overall_accuracy <- rolling_long |>
  group_by(
    disease,
    model
  ) |>
  summarise(
    forecasts =
      n(),
    
    MAE =
      mean(
        absolute_error,
        na.rm = TRUE
      ),
    
    RMSE =
      sqrt(
        mean(
          squared_error,
          na.rm = TRUE
        )
      ),
    
    .groups =
      "drop"
  ) |>
  arrange(
    disease,
    MAE
  )

print(
  overall_accuracy,
  n = Inf
)


# 12. IDENTIFY LOWEST MAE BY HORIZON

lowest_mae_by_horizon <- accuracy_by_horizon |>
  group_by(
    disease,
    horizon
  ) |>
  slice_min(
    MAE,
    n = 1,
    with_ties = FALSE
  ) |>
  ungroup()

print(
  lowest_mae_by_horizon,
  n = Inf
)


# 13. COUNT LOWEST-MAE RESULTS

lowest_mae_counts <- lowest_mae_by_horizon |>
  count(
    model,
    name = "lowest_mae_count"
  ) |>
  arrange(
    desc(
      lowest_mae_count
    )
  )

print(
  lowest_mae_counts
)

# 14. COMPARE MODELS WITH SEASONAL NAIVE BENCHMARK

benchmark_comparison <- accuracy_by_horizon |>
  select(
    disease,
    model,
    horizon,
    MAE
  ) |>
  group_by(
    disease,
    horizon
  ) |>
  mutate(
    seasonal_naive_MAE =
      MAE[
        model == "seasonal_naive"
      ],
    
    MAE_improvement =
      seasonal_naive_MAE - MAE,
    
    percent_MAE_improvement =
      100 *
      (
        seasonal_naive_MAE - MAE
      ) /
      seasonal_naive_MAE
  ) |>
  ungroup() |>
  filter(
    model != "seasonal_naive"
  ) |>
  arrange(
    disease,
    horizon,
    desc(
      percent_MAE_improvement
    )
  )

print(
  benchmark_comparison,
  n = Inf
)

write_csv(
  benchmark_comparison,
  "outputs/tables/forecast_benchmark_comparison.csv"
)

# 15. PLOT MAE BY FORECAST HORIZON

horizon_accuracy_plot <- ggplot(
  accuracy_by_horizon,
  aes(
    x = horizon,
    y = MAE,
    linetype = model
  )
) +
  
  geom_line(
    linewidth = 0.8
  ) +
  
  geom_point(
    size = 2
  ) +
  
  facet_wrap(
    ~ disease,
    scales = "free_y",
    ncol = 2
  ) +
  
  scale_x_continuous(
    breaks = 1:4
  ) +
  
  labs(
    title =
      "Rolling Forecast Accuracy by Forecast Horizon",
    
    subtitle =
      "Out-of-sample 1–4 week forecasts during the 2025/26 respiratory season",
    
    x =
      "Forecast horizon (weeks ahead)",
    
    y =
      "Mean absolute error",
    
    linetype =
      "Forecast model",
    
    caption =
      "Lower MAE indicates smaller forecast error\nSource: UK Health Security Agency"
  ) +
  
  theme_minimal()

horizon_accuracy_plot


# 16. PLOT ONE-WEEK-AHEAD FORECASTS

one_week_forecasts <- rolling_long |>
  filter(
    horizon == 1
  )


one_week_plot <- ggplot(
  one_week_forecasts,
  aes(
    x = target_date
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
      "One-Week-Ahead Respiratory Virus Forecasts",
    
    subtitle =
      "Rolling forecasts during the 2025/26 respiratory season",
    
    x =
      NULL,
    
    y =
      "PCR positivity (%)",
    
    linetype =
      "Forecast model",
    
    caption =
      "Source: UK Health Security Agency"
  ) +
  
  theme_minimal()

one_week_plot


# 17. PLOT FOUR-WEEK-AHEAD FORECASTS

four_week_forecasts <- rolling_long |>
  filter(
    horizon == 4
  )


four_week_plot <- ggplot(
  four_week_forecasts,
  aes(
    x = target_date
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
      "Four-Week-Ahead Respiratory Virus Forecasts",
    
    subtitle =
      "Rolling forecasts during the 2025/26 respiratory season",
    
    x =
      NULL,
    
    y =
      "PCR positivity (%)",
    
    linetype =
      "Forecast model",
    
    caption =
      "Source: UK Health Security Agency"
  ) +
  
  theme_minimal()

four_week_plot


# 18. CREATE OUTPUT FOLDERS

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


# 19. SAVE FIGURES

ggsave(
  "outputs/figures/rolling_forecast_accuracy.png",
  horizon_accuracy_plot,
  width = 12,
  height = 9,
  dpi = 300
)

ggsave(
  "outputs/figures/rolling_forecast_1_week.png",
  one_week_plot,
  width = 12,
  height = 10,
  dpi = 300
)

ggsave(
  "outputs/figures/rolling_forecast_4_week.png",
  four_week_plot,
  width = 12,
  height = 10,
  dpi = 300
)


# 20. SAVE TABLES

write_csv(
  rolling_results,
  "outputs/tables/rolling_forecast_results.csv"
)

write_csv(
  accuracy_by_horizon,
  "outputs/tables/rolling_accuracy_by_horizon.csv"
)

write_csv(
  overall_accuracy,
  "outputs/tables/rolling_overall_accuracy.csv"
)

write_csv(
  lowest_mae_by_horizon,
  "outputs/tables/rolling_lowest_mae_by_horizon.csv"
)


# 21. PRINT SUMMARY

cat(
  "\n========================================\n"
)

cat(
  "ROLLING FORECAST EVALUATION COMPLETE\n"
)

cat(
  "========================================\n"
)

cat(
  "Viruses evaluated:",
  n_distinct(
    rolling_results$disease
  ),
  "\n"
)

cat(
  "Forecast horizons:",
  maximum_horizon,
  "weeks\n"
)

cat(
  "Forecast origins:",
  n_distinct(
    rolling_results$origin_date
  ),
  "\n"
)

cat(
  "Forecast observations:",
  nrow(
    rolling_results
  ),
  "\n"
)

cat(
  "========================================\n"
)


# END