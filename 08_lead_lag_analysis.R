# UKHSA Respiratory Virus Surveillance
# Lead-lag analysis


# 1. PACKAGES

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(lubridate)


# 2. IMPORT ALIGNED DATA

surveillance <- read_csv(
  "data/processed/syndromic_influenza_aligned.csv",
  show_col_types = FALSE
) |>
  mutate(
    week_start = as.Date(week_start)
  ) |>
  arrange(
    week_start
  )


glimpse(
  surveillance
)


# 3. CHECK DATA

cat(
  "\nWeeks available:",
  nrow(surveillance),
  "\n"
)

cat(
  "Date range:",
  format(min(surveillance$week_start)),
  "to",
  format(max(surveillance$week_start)),
  "\n"
)

cat(
  "Missing NHS 111 values:",
  sum(is.na(surveillance$weekly_calls)),
  "\n"
)

cat(
  "Missing influenza values:",
  sum(is.na(surveillance$influenza_positivity)),
  "\n"
)


# 4. STANDARDISE BOTH SERIES

surveillance <- surveillance |>
  mutate(
    nhs111_z =
      as.numeric(
        scale(weekly_calls)
      ),
    
    influenza_z =
      as.numeric(
        scale(influenza_positivity)
      )
  )


# 5. PLOT STANDARDISED SERIES

standardised_long <- surveillance |>
  select(
    week_start,
    nhs111_z,
    influenza_z
  ) |>
  pivot_longer(
    cols = c(
      nhs111_z,
      influenza_z
    ),
    names_to = "indicator",
    values_to = "z_score"
  ) |>
  mutate(
    indicator = recode(
      indicator,
      nhs111_z = "NHS 111 ARI calls",
      influenza_z = "Influenza PCR positivity"
    )
  )


standardised_plot <- ggplot(
  standardised_long,
  aes(
    x = week_start,
    y = z_score,
    linetype = indicator
  )
) +
  geom_line(
    linewidth = 0.8
  ) +
  labs(
    title =
      "NHS 111 Respiratory Calls and Influenza Activity",
    
    subtitle =
      "Standardised weekly surveillance indicators",
    
    x = NULL,
    
    y =
      "Standardised value (z-score)",
    
    linetype =
      "Indicator",
    
    caption =
      "Source: UK Health Security Agency"
  ) +
  theme_minimal()

standardised_plot


# 6. CREATE LAGGED CORRELATION FUNCTION

calculate_lag_correlation <- function(
    data,
    max_lag = 6
) {
  
  lags <- -max_lag:max_lag
  
  
  results <- lapply(
    lags,
    
    function(k) {
      
      if (k > 0) {
        
        # Positive lag:
        # NHS 111 occurs earlier than influenza
        
        nhs111_values <-
          data$nhs111_z[
            1:(nrow(data) - k)
          ]
        
        influenza_values <-
          data$influenza_z[
            (1 + k):nrow(data)
          ]
        
      } else if (k < 0) {
        
        # Negative lag:
        # Influenza occurs earlier than NHS 111
        
        lag_abs <- abs(k)
        
        nhs111_values <-
          data$nhs111_z[
            (1 + lag_abs):nrow(data)
          ]
        
        influenza_values <-
          data$influenza_z[
            1:(nrow(data) - lag_abs)
          ]
        
      } else {
        
        nhs111_values <-
          data$nhs111_z
        
        influenza_values <-
          data$influenza_z
        
      }
      
      
      tibble(
        lag_weeks = k,
        
        correlation = cor(
          nhs111_values,
          influenza_values,
          use = "complete.obs"
        ),
        
        paired_weeks =
          sum(
            complete.cases(
              nhs111_values,
              influenza_values
            )
          )
      )
      
    }
  )
  
  
  bind_rows(
    results
  )
}


# 7. RAW CROSS-CORRELATION

raw_lag_correlations <- calculate_lag_correlation(
  surveillance,
  max_lag = 6
)


print(
  raw_lag_correlations,
  n = Inf
)


# 8. IDENTIFY STRONGEST RAW CORRELATION

strongest_raw <- raw_lag_correlations |>
  slice_max(
    correlation,
    n = 1,
    with_ties = FALSE
  )


print(
  strongest_raw
)


# 9. PLOT RAW LAG CORRELATIONS

raw_lag_plot <- ggplot(
  raw_lag_correlations,
  aes(
    x = lag_weeks,
    y = correlation
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_line(
    linewidth = 0.8
  ) +
  geom_point(
    size = 2.5
  ) +
  scale_x_continuous(
    breaks = -6:6
  ) +
  labs(
    title =
      "Lagged Association Between NHS 111 Calls and Influenza Positivity",
    
    subtitle =
      "Positive lags indicate NHS 111 activity occurs earlier than influenza activity",
    
    x =
      "Lag (weeks)",
    
    y =
      "Pearson correlation",
    
    caption =
      "Raw surveillance series; shared seasonality may contribute to correlation"
  ) +
  theme_minimal()

raw_lag_plot


# 10. CALCULATE WEEK-TO-WEEK CHANGES

changes <- surveillance |>
  mutate(
    nhs111_change =
      nhs111_z - lag(nhs111_z),
    
    influenza_change =
      influenza_z - lag(influenza_z)
  ) |>
  filter(
    !is.na(nhs111_change),
    !is.na(influenza_change)
  )


# 11. STANDARDISE WEEKLY CHANGES

changes <- changes |>
  mutate(
    nhs111_change_z =
      as.numeric(
        scale(nhs111_change)
      ),
    
    influenza_change_z =
      as.numeric(
        scale(influenza_change)
      )
  )


# 12. PLOT WEEKLY CHANGES

changes_long <- changes |>
  select(
    week_start,
    nhs111_change_z,
    influenza_change_z
  ) |>
  pivot_longer(
    cols = c(
      nhs111_change_z,
      influenza_change_z
    ),
    names_to = "indicator",
    values_to = "standardised_change"
  ) |>
  mutate(
    indicator = recode(
      indicator,
      nhs111_change_z =
        "Change in NHS 111 ARI calls",
      
      influenza_change_z =
        "Change in influenza positivity"
    )
  )


changes_plot <- ggplot(
  changes_long,
  aes(
    x = week_start,
    y = standardised_change,
    linetype = indicator
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_line(
    linewidth = 0.7
  ) +
  labs(
    title =
      "Week-to-Week Changes in Respiratory Surveillance",
    
    subtitle =
      "Differencing reduces the influence of shared long-term seasonal patterns",
    
    x = NULL,
    
    y =
      "Standardised weekly change",
    
    linetype =
      "Indicator",
    
    caption =
      "Source: UK Health Security Agency"
  ) +
  theme_minimal()

changes_plot


# 13. PREPARE DIFFERENCED DATA

change_data <- changes |>
  transmute(
    week_start,
    
    nhs111_z =
      nhs111_change_z,
    
    influenza_z =
      influenza_change_z
  )


# 14. CROSS-CORRELATION OF WEEKLY CHANGES

change_lag_correlations <- calculate_lag_correlation(
  change_data,
  max_lag = 6
)


print(
  change_lag_correlations,
  n = Inf
)


# 15. IDENTIFY STRONGEST CHANGE CORRELATION

strongest_change <- change_lag_correlations |>
  slice_max(
    correlation,
    n = 1,
    with_ties = FALSE
  )


print(
  strongest_change
)


# 16. PLOT CHANGE CORRELATIONS

change_lag_plot <- ggplot(
  change_lag_correlations,
  aes(
    x = lag_weeks,
    y = correlation
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_line(
    linewidth = 0.8
  ) +
  geom_point(
    size = 2.5
  ) +
  scale_x_continuous(
    breaks = -6:6
  ) +
  labs(
    title =
      "Lagged Association of Week-to-Week Changes",
    
    subtitle =
      "Positive lags indicate changes in NHS 111 activity occur earlier than changes in influenza positivity",
    
    x =
      "Lag (weeks)",
    
    y =
      "Pearson correlation",
    
    caption =
      "Analysis of first-differenced standardised surveillance series"
  ) +
  theme_minimal()

change_lag_plot


# 17. CREATE RESPIRATORY SEASON

surveillance <- surveillance |>
  mutate(
    season_start_year =
      if_else(
        month(week_start) >= 7,
        year(week_start),
        year(week_start) - 1
      ),
    
    season =
      paste0(
        season_start_year,
        "/",
        substr(
          season_start_year + 1,
          3,
          4
        )
      )
  )


# 18. CHECK NUMBER OF WEEKS PER SEASON

season_coverage <- surveillance |>
  group_by(
    season
  ) |>
  summarise(
    first_week =
      min(week_start),
    
    last_week =
      max(week_start),
    
    weeks =
      n(),
    
    .groups =
      "drop"
  )


print(
  season_coverage,
  n = Inf
)


# 19. KEEP COMPLETE RESPIRATORY SEASONS

complete_seasons <- season_coverage |>
  filter(
    weeks >= 50
  ) |>
  pull(
    season
  )


cat(
  "\nComplete seasons available:\n"
)

print(
  complete_seasons
)


# 20. SEASON-SPECIFIC LAG ANALYSIS

season_lag_correlations <- lapply(
  complete_seasons,
  
  function(season_name) {
    
    season_data <- surveillance |>
      filter(
        season == season_name
      )
    
    
    calculate_lag_correlation(
      season_data,
      max_lag = 6
    ) |>
      mutate(
        season =
          season_name
      )
    
  }
) |>
  bind_rows()


print(
  season_lag_correlations,
  n = Inf
)


# 21. STRONGEST CORRELATION BY SEASON

strongest_by_season <- season_lag_correlations |>
  group_by(
    season
  ) |>
  slice_max(
    correlation,
    n = 1,
    with_ties = FALSE
  ) |>
  ungroup()


print(
  strongest_by_season,
  n = Inf
)


# 22. PLOT SEASON-SPECIFIC LAG CORRELATIONS

season_lag_plot <- ggplot(
  season_lag_correlations,
  aes(
    x = lag_weeks,
    y = correlation,
    linetype = season
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_line(
    linewidth = 0.8
  ) +
  geom_point(
    size = 2
  ) +
  scale_x_continuous(
    breaks = -6:6
  ) +
  labs(
    title =
      "Lagged Association by Respiratory Season",
    
    subtitle =
      "Positive lags indicate NHS 111 activity occurs earlier than influenza activity",
    
    x =
      "Lag (weeks)",
    
    y =
      "Pearson correlation",
    
    linetype =
      "Respiratory season",
    
    caption =
      "Only seasons with at least 50 observed weeks included"
  ) +
  theme_minimal()

season_lag_plot


# 23. CREATE OUTPUT FOLDERS

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


# 24. SAVE TABLES

write_csv(
  raw_lag_correlations,
  "outputs/tables/raw_lag_correlations.csv"
)

write_csv(
  change_lag_correlations,
  "outputs/tables/change_lag_correlations.csv"
)

write_csv(
  season_coverage,
  "outputs/tables/lead_lag_season_coverage.csv"
)

write_csv(
  season_lag_correlations,
  "outputs/tables/season_lag_correlations.csv"
)

write_csv(
  strongest_by_season,
  "outputs/tables/strongest_lag_by_season.csv"
)


# 25. SAVE FIGURES

ggsave(
  "outputs/figures/lead_lag_standardised_series.png",
  standardised_plot,
  width = 12,
  height = 7,
  dpi = 300
)

ggsave(
  "outputs/figures/raw_lag_correlation.png",
  raw_lag_plot,
  width = 10,
  height = 6,
  dpi = 300
)

ggsave(
  "outputs/figures/weekly_change_series.png",
  changes_plot,
  width = 12,
  height = 7,
  dpi = 300
)

ggsave(
  "outputs/figures/change_lag_correlation.png",
  change_lag_plot,
  width = 10,
  height = 6,
  dpi = 300
)

ggsave(
  "outputs/figures/season_lag_correlation.png",
  season_lag_plot,
  width = 10,
  height = 6,
  dpi = 300
)


# 26. FINAL SUMMARY

cat(
  "\n========================================\n"
)

cat(
  "LEAD-LAG ANALYSIS COMPLETE\n"
)

cat(
  "========================================\n"
)

cat(
  "Weeks analysed:",
  nrow(surveillance),
  "\n"
)

cat(
  "Raw strongest correlation:\n"
)

print(
  strongest_raw
)

cat(
  "\nWeekly-change strongest correlation:\n"
)

print(
  strongest_change
)

cat(
  "\nComplete seasons analysed:\n"
)

print(
  complete_seasons
)

cat(
  "========================================\n"
)


# END
