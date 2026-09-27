# UKHSA Respiratory Virus Surveillance
# Statistical signal detection


# 1. PACKAGES

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(lubridate)


# 2. IMPORT CLEAN DATA

respiratory <- read_csv(
  "data/processed/respiratory_surveillance_long.csv",
  show_col_types = FALSE
)

glimpse(respiratory)


# 3. CREATE RESPIRATORY-SEASON VARIABLES

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
    ),
    
    season_start_date = as.Date(
      paste0(
        season_start_year,
        "-07-01"
      )
    ),
    
    week_of_season = floor(
      as.numeric(
        date - season_start_date
      ) / 7
    ) + 1
  )


# 4. CHECK WEEK OF RESPIRATORY SEASON

respiratory |>
  select(
    date,
    disease,
    epiweek,
    season,
    week_of_season
  ) |>
  head(20) |>
  print()


# 5. REMOVE INCOMPLETE CURRENT SEASON FROM BASELINE

# 2026/27 is still incomplete.
#
# We do not want an incomplete season contributing to
# the historical expected values.

historical_data <- respiratory |>
  filter(
    season != "2026/27"
  )


# 6. CALCULATE SEASONAL BASELINE

# For every virus and week of the respiratory season,
# calculate the historical distribution of positivity.
#
# Median = expected activity
# 95th percentile = provisional upper surveillance threshold

seasonal_baseline <- historical_data |>
  group_by(
    disease,
    week_of_season
  ) |>
  summarise(
    expected_positivity = median(
      positivity,
      na.rm = TRUE
    ),
    
    upper_threshold = quantile(
      positivity,
      probs = 0.95,
      na.rm = TRUE
    ),
    
    historical_min = min(
      positivity,
      na.rm = TRUE
    ),
    
    historical_max = max(
      positivity,
      na.rm = TRUE
    ),
    
    observations = n(),
    
    .groups = "drop"
  )


# 7. INSPECT BASELINE

print(
  seasonal_baseline,
  n = 30
)


# 8. ADD BASELINE TO OBSERVED DATA

surveillance <- respiratory |>
  left_join(
    seasonal_baseline,
    by = c(
      "disease",
      "week_of_season"
    )
  )


# 9. CALCULATE DEVIATION FROM EXPECTED ACTIVITY

surveillance <- surveillance |>
  mutate(
    absolute_deviation =
      positivity - expected_positivity,
    
    relative_to_expected =
      if_else(
        expected_positivity > 0,
        positivity / expected_positivity,
        NA_real_
      )
  )


# 10. CREATE SIGNAL FLAG

# A provisional signal occurs when observed positivity
# exceeds the historical 95th percentile for that virus
# and week of the respiratory season.

surveillance <- surveillance |>
  mutate(
    signal = positivity > upper_threshold
  )


# 11. COUNT SIGNALS BY VIRUS

signal_summary <- surveillance |>
  group_by(disease) |>
  summarise(
    observations = n(),
    
    signals = sum(
      signal,
      na.rm = TRUE
    ),
    
    signal_percent = round(
      100 * mean(
        signal,
        na.rm = TRUE
      ),
      2
    ),
    
    .groups = "drop"
  )

print(signal_summary)


# 12. SHOW STRONGEST HISTORICAL SIGNALS

strongest_signals <- surveillance |>
  filter(
    signal == TRUE
  ) |>
  arrange(
    desc(relative_to_expected)
  ) |>
  select(
    date,
    season,
    disease,
    positivity,
    expected_positivity,
    upper_threshold,
    relative_to_expected
  )

print(
  strongest_signals,
  n = 30
)


# 13. PLOT OBSERVED ACTIVITY AGAINST EXPECTED BASELINE

signal_plot <- ggplot(
  surveillance,
  aes(
    x = date
  )
) +
  geom_line(
    aes(
      y = positivity
    )
  ) +
  geom_line(
    aes(
      y = upper_threshold
    ),
    linetype = "dashed"
  ) +
  geom_point(
    data = surveillance |>
      filter(
        signal == TRUE
      ),
    
    aes(
      y = positivity
    ),
    
    size = 1.5
  ) +
  facet_wrap(
    ~ disease,
    scales = "free_y",
    ncol = 2
  ) +
  labs(
    title = "Respiratory Virus Statistical Surveillance Signals",
    subtitle = "Observed PCR positivity compared with historical seasonal thresholds",
    x = NULL,
    y = "PCR positivity (%)",
    caption = "Dashed line = historical 95th percentile for corresponding week of respiratory season\nSource: UK Health Security Agency"
  ) +
  theme_minimal()

signal_plot


# 14. LOOK AT CURRENT SEASON

current_season <- surveillance |>
  filter(
    season == "2026/27"
  )


current_signals <- current_season |>
  filter(
    signal == TRUE
  ) |>
  select(
    date,
    disease,
    positivity,
    expected_positivity,
    upper_threshold,
    relative_to_expected
  )

print(
  current_signals,
  n = Inf
)


# 15. GET LATEST OBSERVATION FOR EACH VIRUS

latest_status <- surveillance |>
  group_by(disease) |>
  filter(
    date == max(date)
  ) |>
  ungroup() |>
  select(
    disease,
    date,
    positivity,
    expected_positivity,
    upper_threshold,
    signal
  )

print(latest_status)


# 16. CREATE SIMPLE STATUS LABEL

latest_status <- latest_status |>
  mutate(
    status = case_when(
      signal == TRUE ~
        "Above historical seasonal threshold",
      
      signal == FALSE ~
        "Within historical seasonal range",
      
      TRUE ~
        "Status unavailable"
    )
  )

print(latest_status)


# 17. CREATE OUTPUT FOLDERS

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

dir.create(
  "outputs/alerts",
  recursive = TRUE,
  showWarnings = FALSE
)


# 18. SAVE SIGNAL FIGURE

ggsave(
  "outputs/figures/surveillance_signals.png",
  signal_plot,
  width = 12,
  height = 10,
  dpi = 300
)


# 19. SAVE SURVEILLANCE DATA

write_csv(
  surveillance,
  "data/processed/respiratory_surveillance_signals.csv"
)


# 20. SAVE SIGNAL TABLES

write_csv(
  signal_summary,
  "outputs/tables/signal_summary.csv"
)

write_csv(
  strongest_signals,
  "outputs/tables/strongest_surveillance_signals.csv"
)

write_csv(
  latest_status,
  "outputs/alerts/latest_surveillance_status.csv"
)


# 21. PRINT CURRENT SURVEILLANCE STATUS

cat(
  "\n========================================\n"
)

cat(
  "CURRENT RESPIRATORY SURVEILLANCE STATUS\n"
)

cat(
  "========================================\n"
)

latest_status |>
  select(
    disease,
    positivity,
    expected_positivity,
    upper_threshold,
    status
  ) |>
  print(
    n = Inf
  )

cat(
  "========================================\n"
)


# END