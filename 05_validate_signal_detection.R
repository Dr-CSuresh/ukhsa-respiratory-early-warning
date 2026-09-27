# UKHSA Respiratory Virus Surveillance
# Validate statistical signal detection


# 1. PACKAGES

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(lubridate)
library(purrr)


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


# 4. CHECK SEASON COMPLETENESS

season_completeness <- respiratory |>
  group_by(
    disease,
    season,
    season_start_year
  ) |>
  summarise(
    weeks_observed = n(),
    .groups = "drop"
  )

print(
  season_completeness,
  n = Inf
)


# 5. IDENTIFY COMPLETE SEASONS

complete_seasons <- season_completeness |>
  filter(
    weeks_observed >= 50
  ) |>
  distinct(season) |>
  arrange(season)

print(
  complete_seasons,
  n = Inf
)


# 6. CREATE VALIDATION FUNCTION

# For each held-out season:
#
# 1. Remove that season from the training data
# 2. Calculate the historical seasonal baseline
# 3. Apply the threshold to the held-out season
# 4. Record which weeks would have generated signals


validate_season <- function(
    held_out_season
) {
  
  message(
    "Validating season: ",
    held_out_season
  )
  
  
  # Training data
  
  training_data <- respiratory |>
    filter(
      season != held_out_season,
      season != "2026/27"
    )
  
  
  # Held-out test data
  
  test_data <- respiratory |>
    filter(
      season == held_out_season
    )
  
  
  # Build baseline without held-out season
  
  baseline <- training_data |>
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
      
      baseline_observations = n(),
      
      .groups = "drop"
    )
  
  
  # Apply baseline to held-out season
  
  validated <- test_data |>
    left_join(
      baseline,
      by = c(
        "disease",
        "week_of_season"
      )
    ) |>
    mutate(
      absolute_deviation =
        positivity - expected_positivity,
      
      relative_to_expected =
        if_else(
          expected_positivity > 0,
          positivity / expected_positivity,
          NA_real_
        ),
      
      signal =
        positivity > upper_threshold,
      
      held_out_season =
        held_out_season
    )
  
  
  return(validated)
}


# 7. RUN LEAVE-ONE-SEASON-OUT VALIDATION

validation_results <- map_dfr(
  complete_seasons$season,
  validate_season
)


# 8. INSPECT VALIDATION RESULTS

glimpse(
  validation_results
)


validation_results |>
  select(
    date,
    disease,
    season,
    positivity,
    expected_positivity,
    upper_threshold,
    signal
  ) |>
  head(30) |>
  print()


# 9. SUMMARISE SIGNALS BY VIRUS

validation_summary <- validation_results |>
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

print(validation_summary)


# 10. SUMMARISE SIGNALS BY VIRUS AND SEASON

signals_by_season <- validation_results |>
  group_by(
    disease,
    season
  ) |>
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

print(
  signals_by_season,
  n = Inf
)


# 11. FIND STRONGEST VALIDATED SIGNALS

strongest_validated_signals <- validation_results |>
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
  strongest_validated_signals,
  n = 30
)


# 12. LOOK SPECIFICALLY AT RSV 2021/22

rsv_2021 <- validation_results |>
  filter(
    disease == "RSV",
    season == "2021/22",
    signal == TRUE
  ) |>
  select(
    date,
    positivity,
    expected_positivity,
    upper_threshold,
    relative_to_expected
  )

print(
  rsv_2021,
  n = Inf
)


# 13. COUNT CONSECUTIVE SIGNALS

# This helps distinguish isolated threshold crossings
# from sustained periods of unusual activity.

validation_results <- validation_results |>
  arrange(
    disease,
    season,
    date
  ) |>
  group_by(
    disease,
    season
  ) |>
  mutate(
    signal_start =
      signal == TRUE &
      lag(
        signal,
        default = FALSE
      ) == FALSE,
    
    signal_run =
      cumsum(signal_start)
  ) |>
  group_by(
    disease,
    season,
    signal_run
  ) |>
  mutate(
    consecutive_signal_weeks =
      if_else(
        signal == TRUE,
        sum(
          signal,
          na.rm = TRUE
        ),
        0L
      )
  ) |>
  ungroup()


# 14. IDENTIFY SUSTAINED SIGNALS

# For this exploratory analysis,
# define a sustained signal as at least
# 3 consecutive weeks above threshold.

sustained_signals <- validation_results |>
  filter(
    signal == TRUE,
    consecutive_signal_weeks >= 3
  )


# 15. SUMMARISE SUSTAINED SIGNAL EPISODES

sustained_episodes <- sustained_signals |>
  group_by(
    disease,
    season,
    signal_run
  ) |>
  summarise(
    signal_start_date =
      min(date),
    
    signal_end_date =
      max(date),
    
    duration_weeks =
      n(),
    
    maximum_positivity =
      max(
        positivity,
        na.rm = TRUE
      ),
    
    maximum_relative_to_expected =
      max(
        relative_to_expected,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  ) |>
  arrange(
    disease,
    signal_start_date
  )

print(
  sustained_episodes,
  n = Inf
)


# 16. PLOT VALIDATED SIGNALS

validation_plot <- ggplot(
  validation_results,
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
    data = validation_results |>
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
    title = "Retrospective Validation of Respiratory Surveillance Signals",
    subtitle = "Each season evaluated against a baseline constructed without that season",
    x = NULL,
    y = "PCR positivity (%)",
    caption = paste0(
      "Points indicate weeks above the held-out seasonal threshold\n",
      "Source: UK Health Security Agency"
    )
  ) +
  theme_minimal()

validation_plot


# 17. PLOT NUMBER OF SIGNALS BY SEASON

signal_count_plot <- ggplot(
  signals_by_season,
  aes(
    x = season,
    y = signals,
    group = disease
  )
) +
  geom_line() +
  geom_point(
    size = 2
  ) +
  facet_wrap(
    ~ disease,
    scales = "free_y",
    ncol = 2
  ) +
  labs(
    title = "Validated Surveillance Signals by Respiratory Season",
    subtitle = "Leave-one-season-out evaluation",
    x = "Respiratory season",
    y = "Number of signal weeks",
    caption = "Source: UK Health Security Agency"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

signal_count_plot


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
  "outputs/figures/validated_surveillance_signals.png",
  validation_plot,
  width = 12,
  height = 10,
  dpi = 300
)

ggsave(
  "outputs/figures/validated_signals_by_season.png",
  signal_count_plot,
  width = 12,
  height = 9,
  dpi = 300
)


# 20. SAVE TABLES

write_csv(
  validation_results,
  "outputs/tables/signal_validation_results.csv"
)

write_csv(
  validation_summary,
  "outputs/tables/signal_validation_summary.csv"
)

write_csv(
  signals_by_season,
  "outputs/tables/signals_by_season.csv"
)

write_csv(
  strongest_validated_signals,
  "outputs/tables/strongest_validated_signals.csv"
)

write_csv(
  sustained_episodes,
  "outputs/tables/sustained_signal_episodes.csv"
)


# 21. VALIDATION SUMMARY

cat(
  "\n========================================\n"
)

cat(
  "SIGNAL DETECTION VALIDATION COMPLETE\n"
)

cat(
  "========================================\n"
)

cat(
  "Seasons evaluated:",
  n_distinct(
    validation_results$season
  ),
  "\n"
)

cat(
  "Viruses evaluated:",
  n_distinct(
    validation_results$disease
  ),
  "\n"
)

cat(
  "Held-out observations:",
  nrow(
    validation_results
  ),
  "\n"
)

cat(
  "Total signals:",
  sum(
    validation_results$signal,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Sustained signal episodes:",
  nrow(
    sustained_episodes
  ),
  "\n"
)

cat(
  "========================================\n"
)


# END