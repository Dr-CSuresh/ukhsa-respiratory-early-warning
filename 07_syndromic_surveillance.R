# UKHSA Respiratory Virus Surveillance
# Syndromic surveillance


# 1. PACKAGES

library(httr2)
library(jsonlite)
library(dplyr)
library(tidyr)
library(readr)
library(lubridate)
library(ggplot2)


# 2. CREATE OUTPUT FOLDERS

dir.create(
  "data/raw",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "data/processed",
  recursive = TRUE,
  showWarnings = FALSE
)

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


# 3. UKHSA API SETTINGS

base_url <-
  "https://api.ukhsa-dashboard.data.gov.uk"

theme <-
  "infectious_disease"

sub_theme <-
  "respiratory"

topic <-
  "acute-respiratory-infection"

geography_type <-
  "Nation"

geography <-
  "England"

metric <-
  "acute-respiratory-infection_syndromic_NHS111triagedcalls_countsByDay"


# 4. BUILD API URL

api_url <- paste0(
  base_url,
  "/themes/", theme,
  "/sub_themes/", sub_theme,
  "/topics/", topic,
  "/geography_types/", geography_type,
  "/geographies/", geography,
  "/metrics/", metric,
  "?page_size=365"
)

api_url


# 5. FETCH ALL NHS 111 DATA

message(
  "Downloading NHS 111 acute respiratory infection data..."
)

all_results <- list()

next_url <- api_url

page_number <- 1


while (!is.null(next_url)) {
  
  message(
    "Fetching page ",
    page_number
  )
  
  response <- request(
    next_url
  ) |>
    req_perform()
  
  
  json_data <- resp_body_json(
    response,
    simplifyVector = TRUE
  )
  
  
  page_data <- as_tibble(
    json_data$results
  )
  
  
  all_results[[
    page_number
  ]] <- page_data
  
  
  next_url <-
    json_data[["next"]]
  
  
  page_number <-
    page_number + 1
}


nhs111_raw <- bind_rows(
  all_results
)


message(
  "Download complete ✓"
)


# 6. INSPECT RAW DATA

glimpse(
  nhs111_raw
)

cat(
  "\nRaw rows:",
  nrow(nhs111_raw),
  "\n"
)

cat(
  "Date range:",
  min(as.Date(nhs111_raw$date)),
  "to",
  max(as.Date(nhs111_raw$date)),
  "\n"
)


# 7. CHECK STRATA

cat(
  "\nAge groups:\n"
)

print(
  unique(
    nhs111_raw$age
  )
)


cat(
  "\nSex categories:\n"
)

print(
  unique(
    nhs111_raw$sex
  )
)


cat(
  "\nStrata:\n"
)

print(
  unique(
    nhs111_raw$stratum
  )
)


# 8. SAVE RAW DATA

write_csv(
  nhs111_raw,
  "data/raw/nhs111_acute_respiratory_calls.csv"
)


# 9. CLEAN ALL-AGE SERIES

nhs111 <- nhs111_raw |>
  filter(
    age == "all",
    sex == "all",
    stratum == "default"
  ) |>
  transmute(
    date =
      as.Date(date),
    
    year =
      year(date),
    
    epiweek =
      isoweek(date),
    
    calls =
      as.numeric(metric_value)
  ) |>
  arrange(
    date
  )


glimpse(
  nhs111
)


# 10. DATA QUALITY CHECKS

cat(
  "\nClean rows:",
  nrow(nhs111),
  "\n"
)

cat(
  "Missing call counts:",
  sum(is.na(nhs111$calls)),
  "\n"
)

cat(
  "Duplicate dates:",
  sum(duplicated(nhs111$date)),
  "\n"
)

cat(
  "Minimum calls:",
  min(nhs111$calls, na.rm = TRUE),
  "\n"
)

cat(
  "Maximum calls:",
  max(nhs111$calls, na.rm = TRUE),
  "\n"
)


# 11. CHECK DAILY DATE COMPLETENESS

expected_dates <- tibble(
  date = seq(
    min(nhs111$date),
    max(nhs111$date),
    by = "day"
  )
)


missing_dates <- expected_dates |>
  anti_join(
    nhs111,
    by = "date"
  )


cat(
  "Missing calendar dates:",
  nrow(missing_dates),
  "\n"
)


if (
  nrow(missing_dates) > 0
) {
  
  print(
    missing_dates,
    n = Inf
  )
  
}


# 12. CREATE WEEK START

nhs111 <- nhs111 |>
  mutate(
    week_start =
      floor_date(
        date,
        unit = "week",
        week_start = 1
      )
  )


# 13. AGGREGATE DAILY CALLS TO WEEKLY COUNTS

nhs111_weekly <- nhs111 |>
  group_by(
    week_start
  ) |>
  summarise(
    weekly_calls =
      sum(
        calls,
        na.rm = TRUE
      ),
    
    days_observed =
      n(),
    
    .groups =
      "drop"
  ) |>
  mutate(
    complete_week =
      days_observed == 7
  )


glimpse(
  nhs111_weekly
)


# 14. CHECK WEEKLY COMPLETENESS

weekly_completeness <- nhs111_weekly |>
  count(
    complete_week
  )


print(
  weekly_completeness
)


incomplete_weeks <- nhs111_weekly |>
  filter(
    !complete_week
  )


print(
  incomplete_weeks,
  n = Inf
)


# 15. KEEP COMPLETE WEEKS

nhs111_weekly_complete <- nhs111_weekly |>
  filter(
    complete_week
  )


# 16. IMPORT RESPIRATORY SURVEILLANCE DATA

respiratory <- read_csv(
  "data/processed/respiratory_surveillance_long.csv",
  show_col_types = FALSE
)


# 17. EXTRACT INFLUENZA POSITIVITY

influenza <- respiratory |>
  filter(
    disease == "Influenza"
  ) |>
  transmute(
    week_start =
      as.Date(date),
    
    influenza_positivity =
      positivity
  ) |>
  arrange(
    week_start
  )


glimpse(
  influenza
)


# 18. ALIGN NHS 111 AND INFLUENZA DATA

syndromic_influenza <- influenza |>
  inner_join(
    nhs111_weekly_complete,
    by = "week_start"
  ) |>
  arrange(
    week_start
  )


glimpse(
  syndromic_influenza
)


cat(
  "\nAligned weeks:",
  nrow(syndromic_influenza),
  "\n"
)

cat(
  "Aligned date range:",
  min(syndromic_influenza$week_start),
  "to",
  max(syndromic_influenza$week_start),
  "\n"
)


# 19. CREATE RESPIRATORY SEASON

syndromic_influenza <- syndromic_influenza |>
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


# 20. SUMMARY STATISTICS

syndromic_summary <- syndromic_influenza |>
  summarise(
    weeks =
      n(),
    
    mean_weekly_calls =
      mean(
        weekly_calls,
        na.rm = TRUE
      ),
    
    median_weekly_calls =
      median(
        weekly_calls,
        na.rm = TRUE
      ),
    
    max_weekly_calls =
      max(
        weekly_calls,
        na.rm = TRUE
      ),
    
    mean_influenza_positivity =
      mean(
        influenza_positivity,
        na.rm = TRUE
      ),
    
    max_influenza_positivity =
      max(
        influenza_positivity,
        na.rm = TRUE
      )
  )


print(
  syndromic_summary
)


# 21. PLOT NHS 111 WEEKLY CALLS

nhs111_plot <- ggplot(
  syndromic_influenza,
  aes(
    x = week_start,
    y = weekly_calls
  )
) +
  
  geom_line(
    linewidth = 0.7
  ) +
  
  labs(
    title =
      "NHS 111 Acute Respiratory Infection Calls in England",
    
    subtitle =
      "Weekly total of daily NHS 111 triaged calls",
    
    x =
      NULL,
    
    y =
      "Weekly NHS 111 calls",
    
    caption =
      "Source: UK Health Security Agency"
  ) +
  
  theme_minimal()

nhs111_plot


# 22. STANDARDISE BOTH SERIES

comparison_standardised <- syndromic_influenza |>
  mutate(
    nhs111_z =
      as.numeric(
        scale(
          weekly_calls
        )
      ),
    
    influenza_z =
      as.numeric(
        scale(
          influenza_positivity
        )
      )
  )


# 23. RESHAPE FOR COMPARISON

comparison_long <- comparison_standardised |>
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
    
    names_to =
      "indicator",
    
    values_to =
      "standardised_value"
  ) |>
  mutate(
    indicator =
      recode(
        indicator,
        
        nhs111_z =
          "NHS 111 ARI calls",
        
        influenza_z =
          "Influenza PCR positivity"
      )
  )


# 24. PLOT STANDARDISED SERIES

comparison_plot <- ggplot(
  comparison_long,
  aes(
    x = week_start,
    y = standardised_value,
    linetype = indicator
  )
) +
  
  geom_line(
    linewidth = 0.7
  ) +
  
  labs(
    title =
      "Syndromic and Laboratory Respiratory Surveillance",
    
    subtitle =
      "Standardised NHS 111 acute respiratory calls and influenza PCR positivity",
    
    x =
      NULL,
    
    y =
      "Standardised value (z-score)",
    
    linetype =
      "Indicator",
    
    caption =
      "Source: UK Health Security Agency"
  ) +
  
  theme_minimal()

comparison_plot


# 25. PLOT BY RESPIRATORY SEASON

season_plot <- ggplot(
  comparison_long |>
    left_join(
      syndromic_influenza |>
        select(
          week_start,
          season
        ),
      by = "week_start"
    ),
  
  aes(
    x = week_start,
    y = standardised_value,
    linetype = indicator
  )
) +
  
  geom_line(
    linewidth = 0.7
  ) +
  
  facet_wrap(
    ~ season,
    scales = "free_x",
    ncol = 2
  ) +
  
  labs(
    title =
      "Syndromic and Influenza Activity by Respiratory Season",
    
    subtitle =
      "Standardised NHS 111 ARI calls and influenza PCR positivity",
    
    x =
      NULL,
    
    y =
      "Standardised value (z-score)",
    
    linetype =
      "Indicator",
    
    caption =
      "Source: UK Health Security Agency"
  ) +
  
  theme_minimal()

season_plot


# 26. SAVE PROCESSED DATA

write_csv(
  nhs111_weekly_complete,
  "data/processed/nhs111_ari_weekly.csv"
)

write_csv(
  syndromic_influenza,
  "data/processed/syndromic_influenza_aligned.csv"
)


# 27. SAVE SUMMARY TABLE

write_csv(
  syndromic_summary,
  "outputs/tables/syndromic_summary.csv"
)


# 28. SAVE FIGURES

ggsave(
  "outputs/figures/nhs111_ari_weekly.png",
  nhs111_plot,
  width = 12,
  height = 6,
  dpi = 300
)

ggsave(
  "outputs/figures/syndromic_influenza_comparison.png",
  comparison_plot,
  width = 12,
  height = 7,
  dpi = 300
)

ggsave(
  "outputs/figures/syndromic_influenza_by_season.png",
  season_plot,
  width = 12,
  height = 10,
  dpi = 300
)


# 29. FINAL SUMMARY

cat(
  "\n========================================\n"
)

cat(
  "SYNDROMIC SURVEILLANCE COMPLETE\n"
)

cat(
  "========================================\n"
)

cat(
  "NHS 111 daily observations:",
  nrow(nhs111),
  "\n"
)

cat(
  "Complete NHS 111 weeks:",
  nrow(nhs111_weekly_complete),
  "\n"
)

cat(
  "Weeks aligned with influenza:",
  nrow(syndromic_influenza),
  "\n"
)

cat(
  "Aligned period:",
  min(syndromic_influenza$week_start),
  "to",
  max(syndromic_influenza$week_start),
  "\n"
)

cat(
  "========================================\n"
)


# END