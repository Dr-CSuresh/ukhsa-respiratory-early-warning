# UKHSA Respiratory Virus Surveillance & Early Warning System
# Script 01: Fetch surveillance data



# 1. PACKAGES


required_packages <- c(
  "httr2",
  "jsonlite",
  "dplyr",
  "purrr",
  "readr"
)


new_packages <- required_packages[
  !(required_packages %in% installed.packages()[, "Package"])
]


if (length(new_packages) > 0) {
  install.packages(new_packages)
}


library(httr2)
library(jsonlite)
library(dplyr)
library(purrr)
library(readr)



# 2. CREATE DATA FOLDERS


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



# 3. UKHSA API BASE URL


base_url <- "https://api.ukhsa-dashboard.data.gov.uk"



# 4. FUNCTION TO DOWNLOAD A UKHSA METRIC


fetch_ukhsa_metric <- function(
    theme,
    sub_theme,
    topic,
    geography_type,
    geography,
    metric
) {
  
  
  # Build API endpoint
  
  endpoint <- paste0(
    base_url,
    "/themes/", theme,
    "/sub_themes/", sub_theme,
    "/topics/", topic,
    "/geography_types/", geography_type,
    "/geographies/", geography,
    "/metrics/", metric,
    "?page_size=365"
  )
  
  
  message(
    "\nDownloading: ",
    topic,
    " | ",
    metric
  )
  
  
  # Create storage for API pages
  
  all_results <- list()
  
  next_url <- endpoint
  
  page_number <- 1
  
  
  # Download every page
  
  while (!is.null(next_url)) {
    
    message(
      "  Page ",
      page_number
    )
    
    
    response <- request(next_url) |>
      req_perform()
    
    
    json_data <- response |>
      resp_body_json(
        simplifyVector = TRUE
      )
    
    
    all_results[[page_number]] <-
      as_tibble(
        json_data$results
      )
    
    
    next_url <-
      json_data[["next"]]
    
    
    page_number <-
      page_number + 1
  }
  
  
  # Combine pages
  
  data <- bind_rows(
    all_results
  )
  
  
  message(
    "  Downloaded ",
    nrow(data),
    " observations."
  )
  
  
  return(data)
}



# 5. COMMON SURVEILLANCE SETTINGS


theme <- "infectious_disease"

sub_theme <- "respiratory"

geography_type <- "Nation"

geography <- "England"



# 6. RESPIRATORY SURVEILLANCE METRICS


surveillance_metrics <- tibble(
  
  disease = c(
    "Influenza",
    "RSV",
    "Adenovirus",
    "hMPV",
    "Parainfluenza",
    "Rhinovirus"
  ),
  
  topic = c(
    "Influenza",
    "RSV",
    "Adenovirus",
    "hMPV",
    "Parainfluenza",
    "Rhinovirus"
  ),
  
  metric = c(
    "influenza_testing_positivityByWeek",
    "RSV_testing_positivityByWeek",
    "adenovirus_testing_positivityByWeek",
    "hMPV_testing_positivityByWeek",
    "parainfluenza_testing_positivityByWeek",
    "rhinovirus_testing_positivityByWeek"
  )
  
)



# 7. DOWNLOAD ALL SURVEILLANCE STREAMS


respiratory_data <- pmap(
  
  surveillance_metrics,
  
  function(
    disease,
    topic,
    metric
  ) {
    
    
    data <- fetch_ukhsa_metric(
      
      theme = theme,
      
      sub_theme = sub_theme,
      
      topic = topic,
      
      geography_type = geography_type,
      
      geography = geography,
      
      metric = metric
      
    )
    
    
    data |>
      mutate(
        disease = disease
      )
    
  }
  
)



# 8. COMBINE ALL RESPIRATORY DATA


respiratory_raw <- bind_rows(
  respiratory_data
)



# 9. INSPECT DATA


glimpse(
  respiratory_raw
)


cat(
  "\n========================================\n"
)

cat(
  "UKHSA DATA DOWNLOAD COMPLETE\n"
)

cat(
  "========================================\n"
)


cat(
  "Total observations:",
  nrow(respiratory_raw),
  "\n"
)


cat(
  "Surveillance streams:",
  n_distinct(respiratory_raw$disease),
  "\n"
)


cat(
  "Earliest observation:",
  min(respiratory_raw$date),
  "\n"
)


cat(
  "Latest observation:",
  max(respiratory_raw$date),
  "\n"
)


cat(
  "========================================\n"
)



# 10. CHECK OBSERVATIONS BY DISEASE


respiratory_raw |>
  
  count(
    disease,
    sort = TRUE
  ) |>
  
  print()



# 11. SAVE RAW DATA


write_csv(
  respiratory_raw,
  "data/raw/ukhsa_respiratory_surveillance.csv"
)



# 12. SAVE DOWNLOAD METADATA


download_metadata <- tibble(
  
  download_time =
    Sys.time(),
  
  number_of_observations =
    nrow(respiratory_raw),
  
  number_of_surveillance_streams =
    n_distinct(respiratory_raw$disease),
  
  earliest_date =
    min(respiratory_raw$date),
  
  latest_date =
    max(respiratory_raw$date)
  
)


write_csv(
  download_metadata,
  "data/raw/download_metadata.csv"
)



# END