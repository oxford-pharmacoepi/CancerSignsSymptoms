# code to extract lsc to pull out SMD > 0.1 for review
library(purrr)


# # 365 to 31 days BEFORE - Conditions
# lsc_conditions <- dataFiltered$summarise_large_scale_characteristics %>%
#   omopgenerics::filterSettings(table_name == "condition_occurrence") %>%
#   omopgenerics::filterSettings(analysis == "standard") %>%
#   filter(group_level == "pancreatic_matched" |
#            group_level == "pancreatic_cancer_sampled" ) %>%
#   filter(variable_level == "-365 to -31") %>%
#   tableLargeScaleCharacteristics(compareBy = "cohort_name",
#                                  smdReference = "pancreatic_matched")
# 
# 
# lsc_conditions_df <- jsonlite::fromJSON(lsc_conditions$x$tag$attribs$data)
# # Convert to tibble
# lsc_conditions_df <- as_tibble(lsc_conditions_df)
# 
# #replace NA with 0
# lsc_4plot_conditions_df <- lsc_4plot_conditions_df %>%
#   mutate(across(c(pancreatic_cancer_365_po_sampled, `pancreatic_cancer_365_po_matched (ref)`), ~replace_na(.x, 0)))
# 
# # recalculate SMD for missing ones
# lsc_4plot_conditions_df <- lsc_4plot_conditions_df %>%
#   mutate(across(c(pancreatic_cancer_365_po_sampled, `pancreatic_cancer_365_po_matched (ref)`), ~replace_na(.x, 0))) %>%
#   mutate(
#     `pancreatic_cancer_365_po_sampled SMD` = if_else(
#       is.na(`pancreatic_cancer_365_po_sampled SMD`),
#       if_else(
#         pancreatic_cancer_365_po_sampled == 0 & `pancreatic_cancer_365_po_matched (ref)` == 0,
#         0,
#         round(
#           suppressWarnings(
#             (pancreatic_cancer_365_po_sampled - `pancreatic_cancer_365_po_matched (ref)`) /
#               sqrt(((pancreatic_cancer_365_po_sampled * (100 - pancreatic_cancer_365_po_sampled)) +
#                       (`pancreatic_cancer_365_po_matched (ref)` * (100 - `pancreatic_cancer_365_po_matched (ref)`))) / 2)
#           ),
#           4
#         )
#       ),
#       `pancreatic_cancer_365_po_sampled SMD`
#     )
#   )
# 
# 




library(dplyr)
library(stringr)
library(purrr)
library(jsonlite)
library(glue)
library(tibble)
library(omopgenerics)

# --- User config: put your 30 cancer identifiers here (use the short identifier you used in group_level)
cancers <- dataFiltered$summarise_large_scale_characteristics$group_level %>%
  unique() %>%
  # keep only those WITHOUT "matched" or "sampled"
  .[!str_detect(., regex("matched|sampled", ignore_case = TRUE))]

# Where to save outputs
out_dir <- "lsc_by_cancer_outputs"
if (!dir.exists(out_dir)) dir.create(out_dir)

process_one_cancer <- function(cancer) {
  message("Processing: ", cancer)
  # build expected group_level variants (adjust if your naming differs)
  expected_groups <- c(glue("{cancer}_matched"),
                       glue("{cancer}_sampled")) # a few alternatives just in case
  
  # run the pipeline (adapt the source dataframe name if different)
  lsc_conditions <- dataFiltered$summarise_large_scale_characteristics %>%
    omopgenerics::filterSettings(table_name == "condition_occurrence") %>%
    omopgenerics::filterSettings(analysis == "standard") %>%
    filter(group_level %in% expected_groups) %>%
    filter(variable_level == "-365 to -31") %>%
    tableLargeScaleCharacteristics(compareBy = "cohort_name",
                                   smdReference = glue("{cancer}_matched"))
  
  # parse JSON -> tibble
  # protect against empty result
  if (is.null(lsc_conditions$x$tag$attribs$data) || lsc_conditions$x$tag$attribs$data == "") {
    warning(glue("No LSC JSON data for {cancer}. Returning NULL."))
    return(NULL)
  }
  
  df <- jsonlite::fromJSON(lsc_conditions$x$tag$attribs$data) %>% as_tibble()
  
  # find sampled and matched percentage columns dynamically
  # Sampled: contains "sampled" AND cancer
  sampled_col <- names(df)[str_detect(names(df), regex(paste0(cancer, ".*sampled|sampled.*", cancer), ignore_case = TRUE))]
  # Matched: contains "matched" AND cancer or contains "(ref)"
  matched_col <- names(df)[str_detect(names(df), regex(paste0(cancer, ".*matched|matched.*", cancer), ignore_case = TRUE)) |
                             str_detect(names(df), regex("\\(ref\\)", ignore_case = TRUE))]
  
  # If not found, fallback to any column containing 'sampled' or 'matched'
  if (length(sampled_col) == 0) sampled_col <- names(df)[str_detect(names(df), regex("sampled", ignore_case = TRUE))]
  if (length(matched_col) == 0) matched_col <- names(df)[str_detect(names(df), regex("matched|\\(ref\\)", ignore_case = TRUE))]
  
  # if ambiguous (multiple matches), try to pick the best that contains the cancer name; otherwise take first
  if (length(sampled_col) > 1) {
    sampled_col <- sampled_col[str_detect(sampled_col, regex(cancer, ignore_case = TRUE))]
    if (length(sampled_col) == 0) sampled_col <- names(df)[str_detect(names(df), regex("sampled", ignore_case = TRUE))][1]
  }
  if (length(matched_col) > 1) {
    matched_col <- matched_col[str_detect(matched_col, regex(cancer, ignore_case = TRUE))]
    if (length(matched_col) == 0) matched_col <- names(df)[str_detect(names(df), regex("matched|\\(ref\\)", ignore_case = TRUE))][1]
  }
  
  if (length(sampled_col) == 0 || length(matched_col) == 0) {
    warning(glue("Could not confidently identify sampled/matched columns for '{cancer}'.\nIdentified sampled_col: {toString(sampled_col)}; matched_col: {toString(matched_col)}"))
    # return df so you can inspect
    return(df)
  }
  
  sampled_col <- sampled_col[1]
  matched_col <- matched_col[1]
  
  message("Found columns: sampled = ", sampled_col, " ; matched = ", matched_col)
  
  # identify SMD column name for sampled side (if exists) or create new name
  sampled_smd_col <- names(df)[str_detect(names(df), regex(paste0("^", sampled_col, ".*SMD$"), ignore_case = TRUE))]
  if (length(sampled_smd_col) == 0) {
    sampled_smd_col <- paste0(sampled_col, " SMD")
  } else sampled_smd_col <- sampled_smd_col[1]
  
  # Make sure sampled and matched columns are numeric (they may be strings)
  df <- df %>%
    mutate(across(all_of(c(sampled_col, matched_col)), ~ as.numeric(.x)))
  
  # replace NA with 0 in the two columns
  df <- df %>%
    mutate(
      "{sampled_col}" := replace_na(.data[[sampled_col]], 0),
      "{matched_col}" := replace_na(.data[[matched_col]], 0)
    )
  
  # Recalculate sampled SMD where missing
  # Calculation matches your original: (p1 - p2) / sqrt(((p1*(100-p1))+(p2*(100-p2)))/2)
  df <- df %>%
    mutate(
      !!sampled_smd_col := if_else(
        is.na(.data[[sampled_smd_col]]),
        if_else(
          .data[[sampled_col]] == 0 & .data[[matched_col]] == 0,
          0,
          round(
            suppressWarnings(
              (.data[[sampled_col]] - .data[[matched_col]]) /
                sqrt(((.data[[sampled_col]] * (100 - .data[[sampled_col]])) +
                        (.data[[matched_col]] * (100 - .data[[matched_col]]))) / 2)
            ),
            4
          )
        ),
        .data[[sampled_smd_col]]
      )
    )
  
  # write csv for inspection
  out_file <- file.path(out_dir, glue("{cancer}_lsc_conditions_-365_to_-31.csv"))
  write.csv(df, out_file, row.names = FALSE)
  message("Wrote ", out_file)
  
  return(df)
}

# run for all cancers and collect results in a named list
results_list <- map(cancers, safely(process_one_cancer)) %>%
  set_names(cancers)

# Unpack results and collect successful outputs
successful <- map(results_list, "result")
errors <- map(results_list, "error")

# Print any errors
imap(errors, ~ if (!is.null(.x)) message(glue("Error processing { .y }: { .x$message }")) )

# Optionally bind all non-null results into one big table with a column indicating cancer
all_lsc <- imap_dfr(successful, function(df, cancer) {
  if (is.null(df)) return(NULL)
  df %>% mutate(cancer = cancer)
})

# Save combined file
write.csv(all_lsc, file.path(out_dir, "all_cancers_lsc_conditions_-365_to_-31.csv"), row.names = FALSE)

# need to make this into long format
