# processing results for characteristics for main study for conditions only

# read in results
result_characteristics <- omopgenerics::importSummarisedResult(
  path = here::here("results", "characteristics")
)

results_char_processed <- tableCharacteristics(result_characteristics)$`_data`

# clean column names
cohort_cols <- grepl("\\[header_name\\]", names(results_char_processed))

names(results_char_processed)[cohort_cols] <- sub(
  ".*\\[header_level\\]",
  "",
  names(results_char_processed)[cohort_cols]
)

# tidy up rows and columns for supplement table one
# remove mean values
results_char_processed <- results_char_processed %>%
  dplyr::filter(`Estimate name` != "Mean (SD)")

# remove rows not needed
vars_to_remove <- c(
  "Number records" ,
  "Days in cohort",
  "Future observation",
  "Cohort start date",
  "Cohort end date",
  "Conditions any time prior to index date",
  "Medications 90 days prior to index date",
  "Conditions 365-91 days prior to index date",
  "Medications 365-91 days prior to index date"
)

results_char_processed <- results_char_processed %>%
  dplyr::mutate(`Variable name` = trimws(`Variable name`)) %>%
  dplyr::mutate(`Variable name` = ifelse(`Variable name` == "", NA, `Variable name`)) %>%
  tidyr::fill(`Variable name`, .direction = "down") %>%
  dplyr::filter(!`Variable name` %in% vars_to_remove)

# remove _ between column names
names(results_char_processed) <- gsub(
  "_",
  " ",
  names(results_char_processed)
)

# order the cancers
cohort_order <- c(
  "lip oral cavity",
  "nasopharynx",
  "oropharynx",
  "hypopharynx",
  "larynx",
  "salivary gland",
  "han",
  "esophagus",
  "stomach",
  "colon",
  "rectum",
  "colorectal",
  "anus",
  "liver",
  "pancreas",
  "lung bron trac",
  "bladder",
  "kidney",
  "cervix",
  "corpus uteri",
  "ovary",
  "breast",
  "prostate",
  "thyroid",
  "brain"
)

ordered_cols <- c(
  "Variable name",
  "Variable level",
  "Estimate name",
  as.vector(rbind(
    paste0(cohort_order, " sampled"),
    paste0(cohort_order, " matched")
  ))
)

results_char_processed <- results_char_processed %>%
  dplyr::select(dplyr::any_of(ordered_cols)) %>% 
  filter(`Variable name` != "Medications 365-31 days prior to index date")

# create a pivotted table for main table one
table1 <- results_char_processed %>%
  
  # Convert from wide to long
  pivot_longer(
    cols = -c(`Variable name`, `Variable level`, `Estimate name`),
    names_to = "cohort",
    values_to = "value"
  ) %>%
  
  # Create cancer site and cohort type variables
  mutate(
    cohort_type = case_when(
      grepl("sampled$", cohort) ~ "Cancer cohort",
      grepl("matched$", cohort) ~ "Matched cohort"
    ),
    cancer_site = sub(" (sampled|matched)$", "", cohort)
  ) %>%
  
  # Create characteristic labels
  mutate(
    characteristic = case_when(
      is.na(`Variable level`) | `Variable level` == "" ~ `Variable name`,
      TRUE ~ paste(`Variable name`, `Variable level`, sep = ": ")
    )
  ) %>%
  
  # Keep required columns
  select(
    cancer_site,
    cohort_type,
    characteristic,
    value
  ) %>%
  
  # Convert characteristics to columns
  pivot_wider(
    names_from = characteristic,
    values_from = value
  ) %>%
  
  # Apply ordering
  mutate(
    cancer_site = factor(
      cancer_site,
      levels = cohort_order
    ),
    cohort_type = factor(
      cohort_type,
      levels = c("Cancer cohort", "Matched cohort")
    )
  ) %>%
  
  arrange(cancer_site, cohort_type)

table1_display <- table1 %>%
  mutate(
    cancer_site = ifelse(
      duplicated(cancer_site),
      "",
      as.character(cancer_site)
    )
  )

# remove columns not required
table1_display <- table1_display %>%
  select(-matches("Age Group|Medications|Conditions|observation|Visits")) %>% 
  select(!c("Sex: Male", "Age" ))


# write out results for table 1
write.csv(
  table1_display,
  file = here::here("results", dbName,  "characteristics", "processed", "table1_365.csv"),
  row.names = FALSE
)

# write out results for supplement
write.csv(
  results_char_processed,
  file = here::here("results", dbName, "characteristics", "processed", "supplement_table1_365.csv"),
  row.names = FALSE
)
