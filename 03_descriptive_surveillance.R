# UKHSA Respiratory Virus Surveillance
# Descriptive surveillance and seasonality


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


# 3. CREATE TIME VARIABLES

# Respiratory seasons run from July to June.
# For example:
# July 2024 to June 2025 = 2024/25 season.

respiratory <- respiratory |>
  mutate(
    month = month(date),
    
    month_name = month(
      date,
      label = TRUE,
      abbr = TRUE
    ),
    
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


# 4. CHECK RESPIRATORY SEASONS

respiratory |>
  distinct(season) |>
  arrange(season) |>
  print(n = Inf)


# 5. OVERALL SUMMARY BY VIRUS

virus_summary <- respiratory |>
  group_by(disease) |>
  summarise(
    observations = n(),
    
    mean_positivity = mean(
      positivity,
      na.rm = TRUE
    ),
    
    median_positivity = median(
      positivity,
      na.rm = TRUE
    ),
    
    sd_positivity = sd(
      positivity,
      na.rm = TRUE
    ),
    
    minimum = min(
      positivity,
      na.rm = TRUE
    ),
    
    maximum = max(
      positivity,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )

print(virus_summary)


# 6. CALCULATE TYPICAL WEEKLY SEASONALITY

weekly_seasonality <- respiratory |>
  group_by(
    disease,
    epiweek
  ) |>
  summarise(
    mean_positivity = mean(
      positivity,
      na.rm = TRUE
    ),
    
    median_positivity = median(
      positivity,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )


# 7. PLOT TYPICAL SEASONAL PATTERN

seasonality_plot <- ggplot(
  weekly_seasonality,
  aes(
    x = epiweek,
    y = median_positivity
  )
) +
  geom_line(
    linewidth = 0.8
  ) +
  facet_wrap(
    ~ disease,
    scales = "free_y",
    ncol = 2
  ) +
  labs(
    title = "Typical Seasonal Pattern of Respiratory Viruses in England",
    subtitle = "Median weekly PCR positivity, 2017–2026",
    x = "Epidemiological week",
    y = "Median PCR positivity (%)",
    caption = "Source: UK Health Security Agency"
  ) +
  theme_minimal()

seasonality_plot


# 8. CALCULATE MONTHLY SEASONALITY

monthly_seasonality <- respiratory |>
  group_by(
    disease,
    month,
    month_name
  ) |>
  summarise(
    median_positivity = median(
      positivity,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )


# 9. CREATE SEASONAL HEATMAP

heatmap_plot <- ggplot(
  monthly_seasonality,
  aes(
    x = month_name,
    y = disease,
    fill = median_positivity
  )
) +
  geom_tile() +
  labs(
    title = "Seasonal Activity of Respiratory Viruses in England",
    subtitle = "Median PCR positivity by calendar month, 2017–2026",
    x = NULL,
    y = NULL,
    fill = "Median\npositivity (%)",
    caption = "Source: UK Health Security Agency"
  ) +
  theme_minimal()

heatmap_plot


# 10. IDENTIFY PEAK WEEK FOR EACH RESPIRATORY SEASON

seasonal_peaks <- respiratory |>
  group_by(
    disease,
    season,
    season_start_year
  ) |>
  slice_max(
    order_by = positivity,
    n = 1,
    with_ties = FALSE
  ) |>
  ungroup() |>
  select(
    disease,
    season,
    season_start_year,
    peak_date = date,
    peak_epiweek = epiweek,
    peak_positivity = positivity
  ) |>
  arrange(
    disease,
    season_start_year
  )

print(
  seasonal_peaks,
  n = Inf
)


# 11. CHECK HOW MANY WEEKS ARE AVAILABLE PER SEASON

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


# 12. KEEP COMPLETE RESPIRATORY SEASONS

# A complete respiratory season should contain approximately
# 52 weeks.
#
# We use >= 50 weeks to allow for calendar variation.

complete_seasons <- season_completeness |>
  filter(
    weeks_observed >= 50
  ) |>
  select(
    disease,
    season
  )

seasonal_peaks_complete <- seasonal_peaks |>
  inner_join(
    complete_seasons,
    by = c(
      "disease",
      "season"
    )
  )


# 13. SUMMARISE SEASONAL PEAKS

peak_summary <- seasonal_peaks_complete |>
  group_by(disease) |>
  summarise(
    seasons_analysed = n(),
    
    median_peak_positivity = median(
      peak_positivity,
      na.rm = TRUE
    ),
    
    maximum_peak_positivity = max(
      peak_positivity,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )

print(peak_summary)


# 14. PLOT PEAK MAGNITUDE BY SEASON

peak_magnitude_plot <- ggplot(
  seasonal_peaks_complete,
  aes(
    x = season_start_year,
    y = peak_positivity
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
    title = "Seasonal Peak Respiratory Virus Activity",
    subtitle = "Maximum weekly PCR positivity within each complete respiratory season",
    x = "Respiratory season starting year",
    y = "Peak PCR positivity (%)",
    caption = "Source: UK Health Security Agency"
  ) +
  theme_minimal()

peak_magnitude_plot


# 15. PLOT TIMING OF SEASONAL PEAKS

peak_timing_plot <- ggplot(
  seasonal_peaks_complete,
  aes(
    x = season_start_year,
    y = peak_epiweek
  )
) +
  geom_line() +
  geom_point(
    size = 2
  ) +
  facet_wrap(
    ~ disease,
    ncol = 2
  ) +
  scale_y_continuous(
    breaks = seq(
      1,
      53,
      by = 4
    )
  ) +
  labs(
    title = "Timing of Seasonal Respiratory Virus Peaks",
    subtitle = "Epidemiological week containing maximum PCR positivity",
    x = "Respiratory season starting year",
    y = "Epidemiological week",
    caption = "Source: UK Health Security Agency"
  ) +
  theme_minimal()

peak_timing_plot


# 16. CREATE OUTPUT FOLDERS

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


# 17. SAVE FIGURES

ggsave(
  "outputs/figures/typical_seasonality.png",
  seasonality_plot,
  width = 12,
  height = 9,
  dpi = 300
)

ggsave(
  "outputs/figures/seasonal_heatmap.png",
  heatmap_plot,
  width = 11,
  height = 6,
  dpi = 300
)

ggsave(
  "outputs/figures/seasonal_peak_magnitude.png",
  peak_magnitude_plot,
  width = 12,
  height = 9,
  dpi = 300
)

ggsave(
  "outputs/figures/seasonal_peak_timing.png",
  peak_timing_plot,
  width = 12,
  height = 9,
  dpi = 300
)


# 18. SAVE TABLES

write_csv(
  virus_summary,
  "outputs/tables/virus_summary.csv"
)

write_csv(
  weekly_seasonality,
  "outputs/tables/weekly_seasonality.csv"
)

write_csv(
  seasonal_peaks,
  "outputs/tables/seasonal_peaks.csv"
)

write_csv(
  peak_summary,
  "outputs/tables/peak_summary.csv"
)


# 19. ANALYSIS SUMMARY

cat(
  "\n========================================\n"
)

cat(
  "DESCRIPTIVE SURVEILLANCE ANALYSIS COMPLETE\n"
)

cat(
  "========================================\n"
)

cat(
  "Viruses analysed:",
  n_distinct(respiratory$disease),
  "\n"
)

cat(
  "Weekly observations:",
  nrow(respiratory),
  "\n"
)

cat(
  "Respiratory seasons:",
  n_distinct(respiratory$season),
  "\n"
)

cat(
  "Complete virus-seasons:",
  nrow(seasonal_peaks_complete),
  "\n"
)

cat(
  "========================================\n"
)


# END