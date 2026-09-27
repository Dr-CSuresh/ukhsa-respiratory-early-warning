# UKHSA Respiratory Virus Surveillance & Early Warning System
# Clean and prepare surveillance data



# 1. PACKAGES


required_packages <- c(
  "dplyr",
  "tidyr",
  "readr",
  "lubridate",
  "ggplot2"
)


new_packages <- required_packages[
  !(required_packages %in% installed.packages()[, "Package"])
]


if (length(new_packages) > 0) {
  install.packages(new_packages)
}


library(dplyr)
library(tidyr)
library(readr)
library(lubridate)
library(ggplot2)



# 2. IMPORT RAW DATA


respiratory_raw <- read_csv(
  "data/raw/ukhsa_respiratory_surveillance.csv",
  show_col_types = FALSE
)


glimpse(respiratory_raw)



# 3. CHECK AVAILABLE AGE GROUPS


respiratory_raw |>
  count(
    disease,
    age
  ) |>
  arrange(
    disease,
    age
  ) |>
  print(n = Inf)



# 4. SELECT MAIN NATIONAL SURVEILLANCE SERIES


# The UKHSA data contain age-specific observations as well
# as an all-age national estimate.
#
# For the primary surveillance analysis we use:
#
#   age = all
#   sex = all
#   stratum = default


respiratory_clean <- respiratory_raw |>
  
  filter(
    age == "all",
    sex == "all",
    stratum == "default"
  ) |>
  
  mutate(
    date = as.Date(date),
    positivity = as.numeric(metric_value)
  ) |>
  
  select(
    date,
    year,
    epiweek,
    disease,
    positivity,
    in_reporting_delay_period
  ) |>
  
  arrange(
    disease,
    date
  )



# 5. INSPECT CLEAN DATA


glimpse(respiratory_clean)


respiratory_clean |>
  count(
    disease,
    sort = TRUE
  ) |>
  print()



# 6. CHECK DATE RANGE FOR EACH VIRUS


date_ranges <- respiratory_clean |>
  
  group_by(disease) |>
  
  summarise(
    
    first_observation =
      min(date, na.rm = TRUE),
    
    last_observation =
      max(date, na.rm = TRUE),
    
    observations =
      n(),
    
    .groups = "drop"
    
  )


print(date_ranges)



# 7. CHECK FOR DUPLICATE VIRUS-WEEKS


duplicates <- respiratory_clean |>
  
  count(
    disease,
    date
  ) |>
  
  filter(
    n > 1
  )


print(duplicates)


cat(
  "\nDuplicate virus-weeks:",
  nrow(duplicates),
  "\n"
)



# 8. CHECK MISSING POSITIVITY VALUES


missing_data <- respiratory_clean |>
  
  group_by(disease) |>
  
  summarise(
    
    observations =
      n(),
    
    missing_positivity =
      sum(is.na(positivity)),
    
    percent_missing =
      round(
        100 * mean(is.na(positivity)),
        2
      ),
    
    .groups = "drop"
    
  )


print(missing_data)



# 9. CHECK POSITIVITY RANGE


positivity_summary <- respiratory_clean |>
  
  group_by(disease) |>
  
  summarise(
    
    minimum =
      min(positivity, na.rm = TRUE),
    
    median =
      median(positivity, na.rm = TRUE),
    
    mean =
      mean(positivity, na.rm = TRUE),
    
    maximum =
      max(positivity, na.rm = TRUE),
    
    .groups = "drop"
    
  )


print(positivity_summary)



# 10. CHECK REPORTING-DELAY FLAGS


reporting_delay_summary <- respiratory_clean |>
  
  group_by(disease) |>
  
  summarise(
    
    observations =
      n(),
    
    reporting_delay =
      sum(
        in_reporting_delay_period == TRUE,
        na.rm = TRUE
      ),
    
    .groups = "drop"
    
  )


print(reporting_delay_summary)



# 11. CREATE WIDE DATASET


# Long format:
#
# date          disease       positivity
# 2025-01-06    Influenza       18.2
# 2025-01-06    RSV              7.4
#
#
# Wide format:
#
# date          Influenza    RSV    hMPV ...
# 2025-01-06       18.2      7.4     ...


respiratory_wide <- respiratory_clean |>
  
  select(
    date,
    disease,
    positivity
  ) |>
  
  pivot_wider(
    names_from = disease,
    values_from = positivity
  ) |>
  
  arrange(date)



# 12. INSPECT WIDE DATASET


glimpse(respiratory_wide)

head(respiratory_wide)



# 13. VISUAL DATA QUALITY CHECK


surveillance_plot <- respiratory_clean |>
  
  ggplot(
    aes(
      x = date,
      y = positivity
    )
  ) +
  
  geom_line() +
  
  facet_wrap(
    ~ disease,
    scales = "free_y",
    ncol = 2
  ) +
  
  labs(
    title = "Respiratory Virus Surveillance in England",
    subtitle = "Weekly PCR positivity by respiratory virus",
    x = NULL,
    y = "PCR positivity (%)",
    caption = "Source: UK Health Security Agency"
  ) +
  
  theme_minimal()


surveillance_plot



# 14. CREATE OUTPUT FOLDER


dir.create(
  "outputs/figures",
  recursive = TRUE,
  showWarnings = FALSE
)



# 15. SAVE SURVEILLANCE FIGURE


ggsave(
  "outputs/figures/respiratory_surveillance_overview.png",
  surveillance_plot,
  width = 12,
  height = 10,
  dpi = 300
)



# 16. SAVE CLEAN LONG DATASET


write_csv(
  respiratory_clean,
  "data/processed/respiratory_surveillance_long.csv"
)



# 17. SAVE CLEAN WIDE DATASET


write_csv(
  respiratory_wide,
  "data/processed/respiratory_surveillance_wide.csv"
)



# 18. PIPELINE SUMMARY


cat("\n")

cat(
  "========================================\n"
)

cat(
  "DATA CLEANING COMPLETE\n"
)

cat(
  "========================================\n"
)


cat(
  "Viruses:",
  n_distinct(respiratory_clean$disease),
  "\n"
)


cat(
  "Clean observations:",
  nrow(respiratory_clean),
  "\n"
)


cat(
  "Earliest observation:",
  as.character(
    min(respiratory_clean$date)
  ),
  "\n"
)


cat(
  "Latest observation:",
  as.character(
    max(respiratory_clean$date)
  ),
  "\n"
)


cat(
  "Duplicate virus-weeks:",
  nrow(duplicates),
  "\n"
)


cat(
  "========================================\n"
)



# END