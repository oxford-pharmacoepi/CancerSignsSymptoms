# processing results for lsc for main study and sensitivity for conditions only

result_lsc <- omopgenerics::bind(results_lsc)

# if reading direct from csv files after analysis
# result_lsc <- omopgenerics::importSummarisedResult(
#   path = here("results", "lsc"))


# get names of cancers
cancer <- names(table(result_lsc$group_level)) %>%
  str_remove("_(matched|sampled)$") %>%
  unique()


# function which extracts each cancer and calculates SMD for main results
lsc_drugs_all <- map_dfr(cancer, function(cancer) {
  
  grp_matched <- paste0(cancer, "_matched")
  grp_sampled <- paste0(cancer, "_sampled")
  
  tryCatch({
    lsc_drugs <- result_lsc %>%
      omopgenerics::filterSettings(table_name == "drug_exposure") %>%
      filter(group_level %in% c(grp_matched, grp_sampled)) %>%
      filter(variable_level == "-365 to -31") %>%
      tableLargeScaleCharacteristics(
        compareBy = "cohort_name",
        smdReference = grp_matched
      )
    
    df <- fromJSON(lsc_drugs$x$tag$attribs$data) %>%
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
filtered_drugs_results_summary <- lsc_drugs_all %>%
  dplyr::filter(
    sampled_SMD > 0.1)


# get ordering for cancers
preferred_order <- c(
  # CNS
  "brain",
  # Head & neck
  "han",
  "nasopharynx",
  "oropharynx",
  "hypopharynx",
  "larynx",
  "lip_oral_cavity",
  "salivary_gland",
  # Thyroid
  "thyroid",
  # Thoracic
  "lung_bron_trac",
  # Upper GI
  "esophagus",
  "stomach",
  "pancreas",
  "liver",
  # Lower GI
  "colorectal",
  "colon",
  "rectum",
  "anus",
  # Female-specific
  "breast",
  "cervix",
  "corpus_uteri",
  "ovary",
  # Genitourinary
  "prostate",
  "bladder",
  "kidney"
)

# refining cancer names
pretty_order <- c(
  "Brain",
  "Head and Neck",
  "Nasopharynx",
  "Oropharynx",
  "Hypopharynx",
  "Larynx",
  "Lip/Oral Cavity",
  "Salivary Gland",
  "Thyroid",
  "Lung",
  "Esophagus",
  "Stomach",
  "Pancreas",
  "Liver",
  "Colorectal",
  "Colon",
  "Rectum",
  "Anus",
  "Breast",
  "Cervix",
  "Corpus Uteri",
  "Ovary",
  "Prostate",
  "Bladder",
  "Kidney"
)





# start from here
# run code for extracting ancestors if not already in the results file
# this code extracts the ancestors for each concept SMD > 0.1 if the file is not already in results folder
if (file.exists(here("concept_mapping", "concept_mapping_4_classification_drugs.csv"))) {
  
  concept_mapping_4_classification <-
    read.csv(here("concept_mapping", "concept_mapping_4_classification_drugs.csv"))
  
} else {
  
  # getting ancestor ids for all significant concepts to help group significant concepts
  # get concept ids
  concept_ids_drugs <- filtered_drugs_results_summary %>%
    distinct(concept_id) %>%
    pull(concept_id)
  
  # get ancestors for concepts
  immediate_ancestors_drugs <- cdm$concept_ancestor %>%
    filter(
      descendant_concept_id %in% concept_ids,
      min_levels_of_separation > 0
    ) %>%
    inner_join(
      cdm$concept %>%
        select(concept_id, ancestor_concept_name = concept_name),
      by = c("ancestor_concept_id" = "concept_id")
    ) %>%
    inner_join(
      cdm$concept %>%
        select(concept_id, descendant_concept_name = concept_name),
      by = c("descendant_concept_id" = "concept_id")
    ) %>%
    collect() %>%
    group_by(descendant_concept_id, descendant_concept_name) %>%
    summarise(
      ancestor_concept_ids   = paste(sort(unique(ancestor_concept_id)), collapse = " | "),
      ancestor_concept_names = paste(sort(unique(ancestor_concept_name)), collapse = " | "),
      .groups = "drop"
    )
  
  # join ancestor summary back onto the filtered condition results
  concept_mapping_drugs <- filtered_drugs_results_summary %>%
    mutate(concept_id = as.numeric(concept_id)) %>%
    distinct(concept_id, concept_name) %>%
    left_join(
      immediate_ancestors,
      by = c("concept_id" = "descendant_concept_id")
    ) %>%
    select(
      concept_id,
      concept_name,
      ancestor_concept_ids,
      ancestor_concept_names
    )
  
  
  # read the file out containing significant concepts and their ancestors
  readr::write_csv(
    concept_mapping,
    here("concept_mapping", "concept_mapping_4_classification_drugs.csv"
    ))
  
}

# read in the file with concepts grouped into categories you will need this for plotting

# drug concepts were grouped into broader clinically meaningful therapeutic categories informed by ATC classification
# and intended clinical use, with related subclasses combined to improve interpretability

# Group	Includes
# Analgesics	Paracetamol, opioids, NSAIDs, combination analgesics
# Antibiotics	All systemic antibiotics, urinary antibiotics
# Antifungals	Topical and systemic antifungals
# Upper GI therapies	PPIs, H2 antagonists, antacids, alginates
# Bowel therapies	Laxatives, antidiarrhoeals, bowel preparation
# Anorectal preparations	Proctofoam HC, Proctosedyl, Rectogesic, Anusol, diltiazem cream, etc.
# Dermatological preparations	Fucibet, Timodine, Sudocrem, emollients
# Respiratory drugs	Bronchodilators, inhalers, COPD/asthma drugs
# Cardiovascular drugs	ACE inhibitors, beta blockers, statins, anticoagulants, antiplatelets, diuretics
# Psychotropic drugs	Antidepressants, anxiolytics, hypnotics, antipsychotics
# Immunomodulators	Corticosteroids, DMARDs, biologics, azathioprine
# Hormonal therapies	Thyroid drugs, HRT, endocrine therapies
# Nutritional supplements	Iron, folate, B12, vitamin D, calcium
# Other	Anything genuinely uncommon

# Gastrointestinal
# Pain & Inflammation
# Anti-infectives
# Cardiovascular
# Respiratory
# Genitourinary
# Neurological
# Immunomodulators
# Dermatological
# Endocrine & Metabolic
# Vitamins & Supplements
# Other (only if truly unavoidable)

if (file.exists(here("concept_mapping", "concept_mapping_4_classification_drugs_with_domains.csv"))) {
  
  concept_mapping_4_classification_with_domains <-
    read.csv(here("concept_mapping", "concept_mapping_4_classification_drugs_with_domains.csv"))
  
} else {
  
  message(
    "You do not have the mapped concept file with groupings: concept_mapping_4_classification_drugs_with_domains.csv"
  )
  
}


# bind domain groupings to results for visulisation
filtered_drugs_results_summary <- filtered_drugs_results_summary %>%
  left_join(
    concept_mapping_4_classification_drugs_with_domains %>%
      mutate(concept_id = as.character(concept_id)) %>%
      select(concept_id, assigned_domain),
    by = "concept_id"
  )


# remove combined phenotypes
filtered_drugs_results_summary2 <- filtered_drugs_results_summary %>%
  dplyr::filter(
    ! cancer %in% c("colorectal", "han")
  )

# figure 1 (all cancer sites and conditions/presentation grouped into groups)
filtered_drugs_results_summary3 <- filtered_drugs_results_summary2 %>%
  mutate(
    assigned_domain = case_when(
      assigned_domain %in% c(
        "Genitourinary >> Urological",
        "Genitourinary >> Male specific"
      ) ~ "Genitourinary",
      
      assigned_domain %in% c(
        "Breast symptoms",
        "Female specific"
      ) ~ "Female specific/breast",
      
      TRUE ~ assigned_domain
    )
  )


domain_cols <- c(
  "Central nervous system" = "#6A3D9A",
  "Lump/mass/lymph node" = "#CAB2D6",
  "Respiratory" = "#1F78B4",
  "Upper abdominal" = "#2A9D8F",    
  "Lower abdominal" = "#B2DF8A",      
  "Female specific/breast" = "#E7298A",
  "Genitourinary" = "#e41a1c",
  "Skin Lesion" = "#D95F02",
  "Ulceration" = "#FB9A99",
  "Musculoskeletal" = "#ff7f00",
  "Non-specific" = "#ffd92f",
  "Other clinical findings" = "#7F7F7F"
)

domain_order <- c(
  "Central nervous system",
  "Lump/mass/lymph node",
  "Respiratory",
  "Upper abdominal",
  "Lower abdominal",
  "Female specific/breast",
  "Genitourinary",
  "Skin Lesion",
  "Ulceration",
  "Musculoskeletal",
  "Non-specific",
  "Other clinical findings"
)

#background for each facet
domain_bg_cols <- c(
  "Central nervous system" = "#EEE8F7",
  "Lump/mass/lymph node" = "#F5F0FA",
  "Respiratory" = "#EAF3FB",
  "Upper abdominal" = "#E8F6F4",      
  "Lower abdominal" = "#F2F8EA",      
  "Female specific/breast" = "#FCEAF3",
  "Genitourinary" = "#FCECEC",
  "Skin Lesion" = "#FDF0E8",
  "Ulceration" = "#FDF3F2",
  "Musculoskeletal" = "#FFF2E6",
  "Non-specific" = "#FFFBE6",
  "Other clinical findings" = "#F2F2F2"
)

domain_cols <- domain_cols[domain_order]


plot_df <- filtered_drugs_results_summary3 %>%
  select(concept_name, assigned_domain, cancer, sampled_SMD) %>%
  mutate(
    sampled_SMD = as.numeric(sampled_SMD),
    assigned_domain = factor(assigned_domain, levels = domain_order),
    cancer = factor(
      cancer,
      levels = preferred_order_figure,
      labels = pretty_cancer_names
    )
  ) %>%
  group_by(assigned_domain, concept_name, cancer) %>%
  summarise(
    sampled_SMD = max(sampled_SMD, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(assigned_domain) %>%
  tidyr::complete(concept_name, cancer) %>%
  ungroup() %>%
  arrange(assigned_domain, concept_name) %>%
  mutate(
    concept_name = factor(concept_name, levels = unique(concept_name))
  )

n_cancer <- nlevels(plot_df$cancer)

# Full panel background rectangles
bg_df <- distinct(plot_df, assigned_domain) %>%
  mutate(
    fill = unname(domain_bg_cols[as.character(assigned_domain)]),
    xmin = 0.5,
    xmax = n_cancer + 0.5
  )

p <- ggplot(
  plot_df,
  aes(
    x = cancer,
    y = concept_name,
    fill = sampled_SMD
  )
) +
  geom_rect(
    data = bg_df,
    inherit.aes = FALSE,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = -Inf,
      ymax = Inf,
      fill = fill
    ),
    colour = NA
  ) +
  scale_fill_identity(guide = "none") +
  ggnewscale::new_scale_fill() +
  
  geom_tile(
    data = subset(plot_df, is.na(sampled_SMD)),
    fill = NA,
    colour = NA
  ) +
  geom_tile(
    data = subset(plot_df, !is.na(sampled_SMD)),
    aes(fill = sampled_SMD),
    colour = "black",
    linewidth = 0.1
  ) +
  
  facet_grid2(
    rows = vars(assigned_domain),
    scales = "free_y",
    space = "free_y",
    switch = "y",
    strip = strip_themed(
      background_y = elem_list_rect(
        fill = domain_cols,
        colour = NA
      ),
      text_y = elem_list_text(
        colour = "black",
        size = 10,
        face = "bold"
      )
    )
  ) +
  
  scale_fill_gradientn(
    name = "SMD",
    colours = c("#1a9850", "#fdae61", "#b2182b"),
    limits = c(0.1, 0.4),
    breaks = c(0.1, 0.2, 0.3, 0.4),
    oob = scales::squish,
    na.value = "grey95"
  ) +
  scale_x_discrete(expand = expansion(add = c(0.4, 0.4))) +
  scale_y_discrete(expand = expansion(add = 0.4)) +
  labs(
    x = "Cancer",
    y = "Domain"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    panel.spacing.y = unit(0.06, "lines"),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, size = 12),
    strip.switch.pad.grid = unit(0, "pt")
  )

p

# save plot
ggsave(
  here("results", "lsc", "processed", "drug_domain_heatmap.png"),
  plot = p,
  width = 12,
  height = 14,
  units = "in",
  dpi = 600,
  device = "png",
  type = "cairo",
  bg = "white"
)



# sensitivity 365-91 window
# function which extracts each cancer and calculates SMD for main results
lsc_drugs_all1 <- map_dfr(cancer, function(cancer) {
  
  grp_matched <- paste0(cancer, "_matched")
  grp_sampled <- paste0(cancer, "_sampled")
  
  tryCatch({
    lsc_drugs1 <- result_lsc %>%
      omopgenerics::filterSettings(table_name == "drug_exposure") %>%
      filter(group_level %in% c(grp_matched, grp_sampled)) %>%
      filter(variable_level == "-365 to -91") %>%
      tableLargeScaleCharacteristics(
        compareBy = "cohort_name",
        smdReference = grp_matched
      )
    
    df <- fromJSON(lsc_drugs1$x$tag$attribs$data) %>%
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
filtered_drugs_results_summary1 <- lsc_drugs_all1 %>%
  dplyr::filter(
    sampled_SMD > 0.1)
