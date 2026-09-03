# read in main results
result_lsc <- omopgenerics::bind(results_lsc)

# if reading direct from csv files after analysis
# result_lsc <- omopgenerics::importSummarisedResult(
#   path = here("results", "lsc"))


# get names of cancers
cancer <- names(table(result_lsc$group_level)) %>%
  str_remove("_(matched|sampled)$") %>%
  unique()


# function which extracts each cancer and calculates SMD for main results
lsc_observations_all <- map_dfr(cancer, function(cancer) {
  
  grp_matched <- paste0(cancer, "_matched")
  grp_sampled <- paste0(cancer, "_sampled")
  
  tryCatch({
    lsc_observations <- result_lsc %>%
      omopgenerics::filterSettings(table_name == "observation") %>%
      filter(group_level %in% c(grp_matched, grp_sampled)) %>%
      filter(variable_level == "-365 to -31") %>%
      tableLargeScaleCharacteristics(
        compareBy = "cohort_name",
        smdReference = grp_matched
      )
    
    df <- fromJSON(lsc_observations$x$tag$attribs$data) %>%
      as_tibble() %>%
      rename(
        matched     = !!paste0(grp_matched, " (ref)"),
        sampled     = !!grp_sampled,
        sampled_SMD = !!paste0(grp_sampled, " SMD")
      ) %>%
      mutate(cancer = cancer)
    
    df
  }, error = function(e) {
    warning(sprintf("Skipping %s: %s", cancer, e$message))
    NULL
  })
})


# get results with SMD > 0.1
filtered_observations_results_summary <- lsc_observations_all %>%
  dplyr::filter(
    sampled_SMD > 0.1)


# with wider interval of 365 - 91
# function which extracts each cancer and calculates SMD for main results
lsc_observations_all1 <- map_dfr(cancer, function(cancer) {
  
  grp_matched <- paste0(cancer, "_matched")
  grp_sampled <- paste0(cancer, "_sampled")
  
  tryCatch({
    lsc_observations1 <- result_lsc %>%
      omopgenerics::filterSettings(table_name == "observation") %>%
      filter(group_level %in% c(grp_matched, grp_sampled)) %>%
      filter(variable_level == "-365 to -91") %>%
      tableLargeScaleCharacteristics(
        compareBy = "cohort_name",
        smdReference = grp_matched
      )
    
    df <- fromJSON(lsc_observations1$x$tag$attribs$data) %>%
      as_tibble() %>%
      rename(
        matched     = !!paste0(grp_matched, " (ref)"),
        sampled     = !!grp_sampled,
        sampled_SMD = !!paste0(grp_sampled, " SMD")
      ) %>%
      mutate(cancer = cancer)
    
    df
  }, error = function(e) {
    warning(sprintf("Skipping %s: %s", cancer, e$message))
    NULL
  })
})


# get results with SMD > 0.1
filtered_observations_results_summary1 <- lsc_observations_all1 %>%
  dplyr::filter(
    sampled_SMD > 0.1)


