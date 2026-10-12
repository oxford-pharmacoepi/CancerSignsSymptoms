# processing results for lsc for main study and sensitivity for conditions only

# read in main results
result_lsc <- omopgenerics::bind(results_lsc)

# if reading direct from csv files after analysis
result_lsc <- omopgenerics::importSummarisedResult(
  path = here("results", dbName, "lsc"))


# get names of cancers
cancer <- names(table(result_lsc$group_level)) %>%
  str_remove("_(matched|sampled)$") %>%
  unique()


# function which extracts each cancer and calculates SMD for main results
lsc_conditions_all <- map_dfr(cancer, function(cancer) {
  
  grp_matched <- paste0(cancer, "_matched")
  grp_sampled <- paste0(cancer, "_sampled")
  
  tryCatch({
    lsc_conditions <- result_lsc %>%
      omopgenerics::filterSettings(table_name == "condition_occurrence") %>%
      filter(group_level %in% c(grp_matched, grp_sampled)) %>%
      filter(variable_level == "-365 to -31") %>%
      tableLargeScaleCharacteristics(
        compareBy = "cohort_name",
        smdReference = grp_matched
      )
    
    df <- fromJSON(lsc_conditions$x$tag$attribs$data) %>%
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
filtered_condition_results_summary <- lsc_conditions_all %>%
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


# run code for extracting ancestors if not already in the results file
# this code extracts the ancestors for each concept SMD > 0.1 if the file is not already in results folder
if (file.exists(here("concept_mapping", "concept_mapping_4_classification_conditions.csv"))) {
  
  concept_mapping_4_classification <-
    read.csv(here("concept_mapping", "concept_mapping_4_classification_conditions.csv"))
  
} else {
  
  # getting ancestor ids for all significant concepts to help group significant concepts
  # get concept ids
  concept_ids <- filtered_condition_results_summary %>%
    distinct(concept_id) %>%
    pull(concept_id)
  
  # get ancestors for concepts
  immediate_ancestors <- cdm$concept_ancestor %>%
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
  concept_mapping <- filtered_condition_results_summary %>%
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
    here("concept_mapping", "concept_mapping_4_classification_conditions.csv"
    ))
  
}

# read in the file with concepts grouped into categories you will need this for plotting
if (file.exists(here("concept_mapping", "concept_mapping_4_classification_conditions_with_domains.csv"))) {
  
  concept_mapping_4_classification_conditions_with_domains <-
    read.csv(here("concept_mapping", "concept_mapping_4_classification_conditions_with_domains.csv"))
  
} else {
  
  message(
    "You do not have the mapped concept file with groupings: concept_mapping_4_classification_conditions_with_domains.csv"
  )
  
}


# bind domain groupings to results for visulisation
filtered_condition_results_summary <- filtered_condition_results_summary %>%
  left_join(
    concept_mapping_4_classification_conditions_with_domains %>%
      mutate(concept_id = as.character(concept_id)) %>%
      select(concept_id, assigned_domain),
    by = "concept_id"
  )


# remove combined phenotypes
filtered_condition_results_summary2 <- filtered_condition_results_summary %>%
  dplyr::filter(
    ! cancer %in% c("colorectal", "han")
  )

# figure 1 (all cancer sites and conditions/presentation grouped into groups)
filtered_condition_results_summary3 <- filtered_condition_results_summary2 %>%
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


# assign red, orange and yellow 
# RED > 0.4, 
# ORANGE 0.2-0.4
# YELLOW: 0.1-0.2
filtered_condition_results_summary3 <- filtered_condition_results_summary3 %>%
  mutate(
    flag = case_when(
      sampled_SMD > 0.4 ~ "RED",
      sampled_SMD > 0.2 & sampled_SMD <= 0.4 ~ "ORANGE",
      sampled_SMD > 0.1 & sampled_SMD <= 0.2 ~ "YELLOW",
      TRUE ~ NA_character_   # or "GREEN" if desired
    )
  )


# write the filtered result out
readr::write_csv(
  filtered_condition_results_summary3,
  here::here("results", dbName, "lsc", "processed", "condition_smds_by_cancer.csv")
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

# get ordering for cancers
preferred_order_figure <- c(
  # CNS
  "brain",
  # Head & neck
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
pretty_cancer_names <- c(
  "Brain",
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

plot_df <- filtered_condition_results_summary3 %>%
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
    name = "ASMD",
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

# save plot
ggsave(
  here("results", dbName, "lsc", "processed", "FIGURE1.png"),
  plot = p,
  width = 12,
  height = 14,
  units = "in",
  dpi = 600,
  device = "png",
  type = "cairo",
  bg = "white"
)


# other figures for specific cancers
# this will create figures with concepts grouped per domain but no domains are shown

# example for lower GI cancers
plot_df_lower_gi <- filtered_condition_results_summary3 %>%
  filter(cancer %in% c("colon", "rectum", "anus")) %>%
  mutate(
    cancer = factor(
      cancer,
      levels = c("colon", "rectum", "anus"),
      labels = c("Colon", "Rectum", "Anus")
    ),
    assigned_domain = factor(assigned_domain, levels = domain_order)
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
    concept_name = factor(
      concept_name,
      levels = rev(unique(concept_name))
    )
  )

p_lower_gi <- ggplot(
  plot_df_lower_gi,
  aes(
    x = cancer,
    y = concept_name,
    fill = sampled_SMD
  )
) +
  geom_tile(
    colour = "grey70",
    linewidth = 0.25
  ) +
  scale_fill_gradientn(
    name = "ASMD",
    colours = c("#1a9850", "#fdae61", "#b2182b"),
    limits = c(0.1, 0.4),
    breaks = c(0.1, 0.2, 0.3, 0.4),
    oob = scales::squish,
    na.value = "grey95"
  ) +
  scale_x_discrete(position = "top") +
  labs(
    x = "Cancer",
    y = "Clinical Concept"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 11, face = "bold"),
    axis.text.y = element_text(size = 9),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.position = "right",
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.switch.pad.grid = unit(0, "pt"),
    panel.spacing.y = unit(0.06, "lines")
  )

ggsave(
  here("results", dbName, "lsc", "processed", "lower_gi_cancers_heatmap_with_concepts.png"),
  plot = p_lower_gi,
  width = 6,
  height = 6,
  units = "in",
  dpi = 600,
  device = "png",
  type = "cairo",
  bg = "white"
)

# upper gi cancers -----------------------------
plot_df_upper_gi <- filtered_condition_results_summary3 %>%
  filter(cancer %in% c("esophagus", "stomach", "pancreas", "liver")) %>%
  mutate(
    cancer = factor(
      cancer,
      levels = c("esophagus", "stomach", "pancreas", "liver"),
      labels = c("Esophagus", "Stomach", "Pancreas", "Liver")
    ),
    assigned_domain = factor(assigned_domain, levels = domain_order)
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
    concept_name = factor(
      concept_name,
      levels = rev(unique(concept_name))
    )
  )

p_upper_gi <- ggplot(
  plot_df_upper_gi,
  aes(
    x = cancer,
    y = concept_name,
    fill = sampled_SMD
  )
) +
  geom_tile(
    colour = "grey70",
    linewidth = 0.25
  ) +
  scale_fill_gradientn(
    name = "ASMD",
    colours = c("#1a9850", "#fdae61", "#b2182b"),
    limits = c(0.1, 0.4),
    breaks = c(0.1, 0.2, 0.3, 0.4),
    oob = scales::squish,
    na.value = "grey95"
  ) +
  scale_x_discrete(position = "top") +
  labs(
    x = "Cancer",
    y = "Clinical Concept"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 11, face = "bold"),
    axis.text.y = element_text(size = 9),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.position = "right",
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.switch.pad.grid = unit(0, "pt"),
    panel.spacing.y = unit(0.06, "lines")
  )

ggsave(
  here("results", dbName, "lsc", "processed", "upper_gi_cancers_heatmap_with_concepts.png"),
  plot = p_upper_gi,
  width = 8,
  height = 9,
  units = "in",
  dpi = 600,
  device = "png",
  type = "cairo",
  bg = "white"
)


# han cancers -----------------------------------
plot_df_han <- filtered_condition_results_summary3 %>%
  filter(cancer %in% c( "hypopharynx", 
                        "larynx", 
                        "oropharynx" ,
                        "salivary_gland" ,
                        "lip_oral_cavity",
                        "nasopharynx")) %>%
  mutate(
    cancer = factor(
      cancer,
      levels = c("nasopharynx",
                 "oropharynx",
                 "hypopharynx",
                 "larynx",
                 "lip_oral_cavity",
                 "salivary_gland"),
      labels = c("Nasopharynx",
                 "Oropharynx",
                 "Hypopharynx",
                 "Larynx",
                 "Lip/oral cavity",
                 "Salivary gland")
    ),
    assigned_domain = factor(assigned_domain, levels = domain_order)
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
    concept_name = factor(
      concept_name,
      levels = rev(unique(concept_name))
    )
  )

p_han <- ggplot(
  plot_df_han,
  aes(
    x = cancer,
    y = concept_name,
    fill = sampled_SMD
  )
) +
  geom_tile(
    colour = "grey70",
    linewidth = 0.25
  ) +
  scale_fill_gradientn(
    name = "ASMD",
    colours = c("#1a9850", "#fdae61", "#b2182b"),
    limits = c(0.1, 0.4),
    breaks = c(0.1, 0.2, 0.3, 0.4),
    oob = scales::squish,
    na.value = "grey95"
  ) +
  scale_x_discrete(position = "top") +
  labs(
    x = "Cancer",
    y = "Clinical Concept"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 11, face = "bold"),
    axis.text.y = element_text(size = 9),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.position = "right",
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.switch.pad.grid = unit(0, "pt"),
    panel.spacing.y = unit(0.06, "lines")
  )

ggsave(
  here("results", dbName, "lsc", "processed", "han_cancers_heatmap_with_concepts.png"),
  plot = p_han,
  width = 12,
  height = 14,
  units = "in",
  dpi = 600,
  device = "png",
  type = "cairo",
  bg = "white"
)

# lung cancer  -----------------------------
plot_df_lung <- filtered_condition_results_summary3 %>%
  filter(cancer %in% c( "lung_bron_trac")) %>%
  mutate(
    cancer = factor(
      cancer,
      levels = c("lung_bron_trac"),
      labels = c("Lung/Bronchus/Trachea")
    ),
    assigned_domain = factor(assigned_domain, levels = domain_order)
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
    concept_name = factor(
      concept_name,
      levels = rev(unique(concept_name))
    )
  )

p_lung <- ggplot(
  plot_df_lung,
  aes(
    x = cancer,
    y = concept_name,
    fill = sampled_SMD
  )
) +
  geom_tile(
    colour = "grey70",
    linewidth = 0.25
  ) +
  scale_fill_gradientn(
    name = "ASMD",
    colours = c("#1a9850", "#fdae61", "#b2182b"),
    limits = c(0.1, 0.4),
    breaks = c(0.1, 0.2, 0.3, 0.4),
    oob = scales::squish,
    na.value = "grey95"
  ) +
  scale_x_discrete(position = "top") +
  labs(
    x = "Cancer",
    y = "Clinical Concept"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 11, face = "bold"),
    axis.text.y = element_text(size = 9),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.position = "right",
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.switch.pad.grid = unit(0, "pt"),
    panel.spacing.y = unit(0.06, "lines")
  )

ggsave(
  here("results", dbName, "lsc", "processed", "lung_cancers_heatmap_with_concepts.png"),
  plot = p_lung,
  width = 6,
  height = 6,
  units = "in",
  dpi = 600,
  device = "png",
  type = "cairo",
  bg = "white"
)


# female specific/breast ------
plot_df_female <- filtered_condition_results_summary3 %>%
  filter(cancer %in% c( "breast", 
                        "cervix", 
                        "corpus_uteri" ,
                        "ovary")) %>%
  mutate(
    cancer = factor(
      cancer,
      levels = c("cervix",
                 "corpus_uteri" ,
                 "ovary",
                 "breast"),
      labels = c("Cervix",
                 "Corpus Uteri",
                 "Ovary",
                 "Breast")
    ),
    assigned_domain = factor(assigned_domain, levels = domain_order)
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
    concept_name = factor(
      concept_name,
      levels = rev(unique(concept_name))
    )
  )

p_female <- ggplot(
  plot_df_female,
  aes(
    x = cancer,
    y = concept_name,
    fill = sampled_SMD
  )
) +
  geom_tile(
    colour = "grey70",
    linewidth = 0.25
  ) +
  scale_fill_gradientn(
    name = "ASMD",
    colours = c("#1a9850", "#fdae61", "#b2182b"),
    limits = c(0.1, 0.4),
    breaks = c(0.1, 0.2, 0.3, 0.4),
    oob = scales::squish,
    na.value = "grey95"
  ) +
  scale_x_discrete(position = "top") +
  labs(
    x = "Cancer",
    y = "Clinical Concept"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 11, face = "bold"),
    axis.text.y = element_text(size = 9),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.position = "right",
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.switch.pad.grid = unit(0, "pt"),
    panel.spacing.y = unit(0.06, "lines")
  )

ggsave(
  here("results", dbName, "lsc", "processed", "female_and_breast_cancers_heatmap_with_concepts.png"),
  plot = p_female,
  width = 9,
  height = 8,
  units = "in",
  dpi = 600,
  device = "png",
  type = "cairo",
  bg = "white"
)


# genitourinary ------------------------
plot_df_genitourinary <- filtered_condition_results_summary3 %>%
  filter(cancer %in% c( "bladder", 
                        "kidney", 
                        "prostate" )) %>%
  mutate(
    cancer = factor(
      cancer,
      levels = c( "kidney", 
                  "bladder",
                  "prostate"),
      labels = c("Kidney", 
                 "Bladder",
                 "Prostate")
    ),
    assigned_domain = factor(assigned_domain, levels = domain_order)
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
    concept_name = factor(
      concept_name,
      levels = rev(unique(concept_name))
    )
  )

p_genitourinary <- ggplot(
  plot_df_genitourinary,
  aes(
    x = cancer,
    y = concept_name,
    fill = sampled_SMD
  )
) +
  geom_tile(
    colour = "grey70",
    linewidth = 0.25
  ) +
  scale_fill_gradientn(
    name = "ASMD",
    colours = c("#1a9850", "#fdae61", "#b2182b"),
    limits = c(0.1, 0.4),
    breaks = c(0.1, 0.2, 0.3, 0.4),
    oob = scales::squish,
    na.value = "grey95"
  ) +
  scale_x_discrete(position = "top") +
  labs(
    x = "Cancer",
    y = "Clinical Concept"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 11, face = "bold"),
    axis.text.y = element_text(size = 9),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.position = "right",
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.switch.pad.grid = unit(0, "pt"),
    panel.spacing.y = unit(0.06, "lines")
  )

ggsave(
  here("results", dbName, "lsc", "processed", "genitourinary_cancers_heatmap_with_concepts.png"),
  plot = p_genitourinary,
  width = 6,
  height = 6,
  units = "in",
  dpi = 600,
  device = "png",
  type = "cairo",
  bg = "white"
)

# other ----------------------------------
plot_df_other <- filtered_condition_results_summary3 %>%
  filter(cancer %in% c( "brain", 
                        "thyroid" )) %>%
  mutate(
    cancer = factor(
      cancer,
      levels = c( "brain", 
                  "thyroid"),
      labels = c("Brain", 
                 "Thyroid")
    ),
    assigned_domain = factor(assigned_domain, levels = domain_order)
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
    concept_name = factor(
      concept_name,
      levels = rev(unique(concept_name))
    )
  )

p_other <- ggplot(
  plot_df_other,
  aes(
    x = cancer,
    y = concept_name,
    fill = sampled_SMD
  )
) +
  geom_tile(
    colour = "grey70",
    linewidth = 0.25
  ) +
  scale_fill_gradientn(
    name = "ASMD",
    colours = c("#1a9850", "#fdae61", "#b2182b"),
    limits = c(0.1, 0.4),
    breaks = c(0.1, 0.2, 0.3, 0.4),
    oob = scales::squish,
    na.value = "grey95"
  ) +
  scale_x_discrete(position = "top") +
  labs(
    x = "Cancer",
    y = "Clinical Concept"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 11, face = "bold"),
    axis.text.y = element_text(size = 9),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.position = "right",
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.switch.pad.grid = unit(0, "pt"),
    panel.spacing.y = unit(0.06, "lines")
  )

ggsave(
  here("results", dbName, "lsc", "processed", "other_cancers_heatmap_with_concepts.png"),
  plot = p_other,
  width = 6,
  height = 6,
  units = "in",
  dpi = 600,
  device = "png",
  type = "cairo",
  bg = "white"
)


# sensitivity 365-91 window #######
lsc_conditions_all1 <- map_dfr(cancer, function(cancer) {
  
  grp_matched <- paste0(cancer, "_matched")
  grp_sampled <- paste0(cancer, "_sampled")
  
  tryCatch({
    lsc_conditions1 <- result_lsc %>%
      omopgenerics::filterSettings(table_name == "condition_occurrence") %>%
      filter(group_level %in% c(grp_matched, grp_sampled)) %>%
      filter(variable_level == "-365 to -91") %>%
      tableLargeScaleCharacteristics(
        compareBy = "cohort_name",
        smdReference = grp_matched
      )
    
    df <- fromJSON(lsc_conditions1$x$tag$attribs$data) %>%
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
filtered_condition_results_summary_sens <- lsc_conditions_all1 %>%
  dplyr::filter(
    sampled_SMD > 0.1)


# run code for extracting ancestors if not already in the results file
# this code extracts the ancestors for each concept SMD > 0.1 if the file is not already in results folder
if (file.exists(here("concept_mapping", "concept_mapping_4_classification_conditions_sens.csv"))) {
  
  concept_mapping_4_classification_sens <-
    read.csv(here("concept_mapping", "concept_mapping_4_classification_conditions_sens.csv"))
  
} else {
  
  # getting ancestor ids for all significant concepts to help group significant concepts
  # get concept ids
  concept_ids <- filtered_condition_results_summary_sens %>%
    distinct(concept_id) %>%
    pull(concept_id)
  
  # get ancestors for concepts
  immediate_ancestors <- cdm$concept_ancestor %>%
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
  concept_mapping_sens <- filtered_condition_results_summary_sens %>%
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
    concept_mapping_sens,
    here("concept_mapping", "concept_mapping_4_classification_conditions_sens.csv"
    ))
  
}

# read in the file with concepts grouped into categories you will need this for plotting
if (file.exists(here("concept_mapping", "concept_mapping_4_classification_conditions_with_domains_sens.csv"))) {
  
  concept_mapping_4_classification_conditions_with_domains_sens <-
    read.csv(here("concept_mapping", "concept_mapping_4_classification_conditions_with_domains_sens.csv"))
  
} else {
  
  message(
    "You do not have the mapped concept file with groupings: concept_mapping_4_classification_conditions_with_domains.csv"
  )
  
}


# bind domain groupings to results for visulisation
filtered_condition_results_summary_sens <- filtered_condition_results_summary_sens %>%
  left_join(
    concept_mapping_4_classification_conditions_with_domains_sens %>%
      mutate(concept_id = as.character(concept_id)) %>%
      select(concept_id, assigned_domain),
    by = "concept_id"
  )


# remove combined phenotypes
filtered_condition_results_summary2_sens <- filtered_condition_results_summary_sens %>%
  dplyr::filter(
    ! cancer %in% c("colorectal", "han")
  )

# figure 1 (all cancer sites and conditions/presentation grouped into groups)
filtered_condition_results_summary3_sens <- filtered_condition_results_summary2_sens %>%
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


# assign red, orange and yellow 
# RED > 0.4, 
# ORANGE 0.2-0.4
# YELLOW: 0.1-0.2
filtered_condition_results_summary3_sens <- filtered_condition_results_summary3_sens %>%
  mutate(
    flag = case_when(
      sampled_SMD > 0.4 ~ "RED",
      sampled_SMD > 0.2 & sampled_SMD <= 0.4 ~ "ORANGE",
      sampled_SMD > 0.1 & sampled_SMD <= 0.2 ~ "YELLOW",
      TRUE ~ NA_character_   
    )
  )

# write the filtered result out
readr::write_csv(
  filtered_condition_results_summary3_sens,
  here::here("results", dbName, "lsc", "processed", "condition_smds_by_cancer_sens.csv")
)


# figures
plot_df1 <- filtered_condition_results_summary3_sens %>%
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

n_cancer <- nlevels(plot_df1$cancer)

# Full panel background rectangles
bg_df1 <- distinct(plot_df1, assigned_domain) %>%
  mutate(
    fill = unname(domain_bg_cols[as.character(assigned_domain)]),
    xmin = 0.5,
    xmax = n_cancer + 0.5
  )

p1 <- ggplot(
  plot_df1,
  aes(
    x = cancer,
    y = concept_name,
    fill = sampled_SMD
  )
) +
  geom_rect(
    data = bg_df1,
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
    data = subset(plot_df1, is.na(sampled_SMD)),
    fill = NA,
    colour = NA
  ) +
  geom_tile(
    data = subset(plot_df1, !is.na(sampled_SMD)),
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
    name = "ASMD",
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

# save plot
ggsave(
  here("results", dbName, "lsc", "processed", "symptom_domain_heatmap_sens.png"),
  plot = p1,
  width = 12,
  height = 14,
  units = "in",
  dpi = 600,
  device = "png",
  type = "cairo",
  bg = "white"
)


# other figures for specific cancers
# this will create figures with concepts grouped per domain but no domains are shown

# example for lower GI cancers
plot_df_lower_gi1 <- filtered_condition_results_summary3_sens %>%
  filter(cancer %in% c("colon", "rectum", "anus")) %>%
  mutate(
    cancer = factor(
      cancer,
      levels = c("colon", "rectum", "anus"),
      labels = c("Colon", "Rectum", "Anus")
    ),
    assigned_domain = factor(assigned_domain, levels = domain_order)
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
    concept_name = factor(
      concept_name,
      levels = rev(unique(concept_name))
    )
  )

p_lower_gi1 <- ggplot(
  plot_df_lower_gi1,
  aes(
    x = cancer,
    y = concept_name,
    fill = sampled_SMD
  )
) +
  geom_tile(
    colour = "grey70",
    linewidth = 0.25
  ) +
  scale_fill_gradientn(
    name = "ASMD",
    colours = c("#1a9850", "#fdae61", "#b2182b"),
    limits = c(0.1, 0.4),
    breaks = c(0.1, 0.2, 0.3, 0.4),
    oob = scales::squish,
    na.value = "grey95"
  ) +
  scale_x_discrete(position = "top") +
  labs(
    x = "Cancer",
    y = "Clinical Concept"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 11, face = "bold"),
    axis.text.y = element_text(size = 9),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.position = "right",
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.switch.pad.grid = unit(0, "pt"),
    panel.spacing.y = unit(0.06, "lines")
  )

ggsave(
  here("results", dbName, "lsc", "processed", "lower_gi_cancers_heatmap_with_concepts_sens.png"),
  plot = p_lower_gi1,
  width = 6,
  height = 6,
  units = "in",
  dpi = 600,
  device = "png",
  type = "cairo",
  bg = "white"
)

# upper gi cancers -----------------------------
plot_df_upper_gi1 <- filtered_condition_results_summary3_sens %>%
  filter(cancer %in% c("esophagus", "stomach", "pancreas", "liver")) %>%
  mutate(
    cancer = factor(
      cancer,
      levels = c("esophagus", "stomach", "pancreas", "liver"),
      labels = c("Esophagus", "Stomach", "Pancreas", "Liver")
    ),
    assigned_domain = factor(assigned_domain, levels = domain_order)
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
    concept_name = factor(
      concept_name,
      levels = rev(unique(concept_name))
    )
  )

p_upper_gi1 <- ggplot(
  plot_df_upper_gi1,
  aes(
    x = cancer,
    y = concept_name,
    fill = sampled_SMD
  )
) +
  geom_tile(
    colour = "grey70",
    linewidth = 0.25
  ) +
  scale_fill_gradientn(
    name = "ASMD",
    colours = c("#1a9850", "#fdae61", "#b2182b"),
    limits = c(0.1, 0.4),
    breaks = c(0.1, 0.2, 0.3, 0.4),
    oob = scales::squish,
    na.value = "grey95"
  ) +
  scale_x_discrete(position = "top") +
  labs(
    x = "Cancer",
    y = "Clinical Concept"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 11, face = "bold"),
    axis.text.y = element_text(size = 9),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.position = "right",
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.switch.pad.grid = unit(0, "pt"),
    panel.spacing.y = unit(0.06, "lines")
  )

ggsave(
  here("results", dbName, "lsc", "processed", "upper_gi_cancers_heatmap_with_concepts_sens.png"),
  plot = p_upper_gi1,
  width = 8,
  height = 9,
  units = "in",
  dpi = 600,
  device = "png",
  type = "cairo",
  bg = "white"
)


# han cancers -----------------------------------
plot_df_han1 <- filtered_condition_results_summary3_sens %>%
  filter(cancer %in% c( "hypopharynx", 
                        "larynx", 
                        "oropharynx" ,
                        "salivary_gland" ,
                        "lip_oral_cavity",
                        "nasopharynx")) %>%
  mutate(
    cancer = factor(
      cancer,
      levels = c("nasopharynx",
                 "oropharynx",
                 "hypopharynx",
                 "larynx",
                 "lip_oral_cavity",
                 "salivary_gland"),
      labels = c("Nasopharynx",
                 "Oropharynx",
                 "Hypopharynx",
                 "Larynx",
                 "Lip/oral cavity",
                 "Salivary gland")
    ),
    assigned_domain = factor(assigned_domain, levels = domain_order)
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
    concept_name = factor(
      concept_name,
      levels = rev(unique(concept_name))
    )
  )

p_han1 <- ggplot(
  plot_df_han1,
  aes(
    x = cancer,
    y = concept_name,
    fill = sampled_SMD
  )
) +
  geom_tile(
    colour = "grey70",
    linewidth = 0.25
  ) +
  scale_fill_gradientn(
    name = "ASMD",
    colours = c("#1a9850", "#fdae61", "#b2182b"),
    limits = c(0.1, 0.4),
    breaks = c(0.1, 0.2, 0.3, 0.4),
    oob = scales::squish,
    na.value = "grey95"
  ) +
  scale_x_discrete(position = "top") +
  labs(
    x = "Cancer",
    y = "Clinical Concept"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 11, face = "bold"),
    axis.text.y = element_text(size = 9),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.position = "right",
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.switch.pad.grid = unit(0, "pt"),
    panel.spacing.y = unit(0.06, "lines")
  )

ggsave(
  here("results", dbName, "lsc", "processed", "han_cancers_heatmap_with_concepts_sens.png"),
  plot = p_han1,
  width = 12,
  height = 14,
  units = "in",
  dpi = 600,
  device = "png",
  type = "cairo",
  bg = "white"
)

# lung cancer  -----------------------------
plot_df_lung1 <- filtered_condition_results_summary3_sens %>%
  filter(cancer %in% c( "lung_bron_trac")) %>%
  mutate(
    cancer = factor(
      cancer,
      levels = c("lung_bron_trac"),
      labels = c("Lung/Bronchus/Trachea")
    ),
    assigned_domain = factor(assigned_domain, levels = domain_order)
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
    concept_name = factor(
      concept_name,
      levels = rev(unique(concept_name))
    )
  )

p_lung1 <- ggplot(
  plot_df_lung1,
  aes(
    x = cancer,
    y = concept_name,
    fill = sampled_SMD
  )
) +
  geom_tile(
    colour = "grey70",
    linewidth = 0.25
  ) +
  scale_fill_gradientn(
    name = "ASMD",
    colours = c("#1a9850", "#fdae61", "#b2182b"),
    limits = c(0.1, 0.4),
    breaks = c(0.1, 0.2, 0.3, 0.4),
    oob = scales::squish,
    na.value = "grey95"
  ) +
  scale_x_discrete(position = "top") +
  labs(
    x = "Cancer",
    y = "Clinical Concept"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 11, face = "bold"),
    axis.text.y = element_text(size = 9),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.position = "right",
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.switch.pad.grid = unit(0, "pt"),
    panel.spacing.y = unit(0.06, "lines")
  )

ggsave(
  here("results", dbName, "lsc", "processed", "lung_cancers_heatmap_with_concepts_sens.png"),
  plot = p_lung1,
  width = 6,
  height = 6,
  units = "in",
  dpi = 600,
  device = "png",
  type = "cairo",
  bg = "white"
)


# female specific/breast ------
plot_df_female1 <- filtered_condition_results_summary3_sens %>%
  filter(cancer %in% c( "breast", 
                        "cervix", 
                        "corpus_uteri" ,
                        "ovary")) %>%
  mutate(
    cancer = factor(
      cancer,
      levels = c("cervix",
                 "corpus_uteri" ,
                 "ovary",
                 "breast"),
      labels = c("Cervix",
                 "Corpus Uteri",
                 "Ovary",
                 "Breast")
    ),
    assigned_domain = factor(assigned_domain, levels = domain_order)
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
    concept_name = factor(
      concept_name,
      levels = rev(unique(concept_name))
    )
  )

p_female1 <- ggplot(
  plot_df_female1,
  aes(
    x = cancer,
    y = concept_name,
    fill = sampled_SMD
  )
) +
  geom_tile(
    colour = "grey70",
    linewidth = 0.25
  ) +
  scale_fill_gradientn(
    name = "ASMD",
    colours = c("#1a9850", "#fdae61", "#b2182b"),
    limits = c(0.1, 0.4),
    breaks = c(0.1, 0.2, 0.3, 0.4),
    oob = scales::squish,
    na.value = "grey95"
  ) +
  scale_x_discrete(position = "top") +
  labs(
    x = "Cancer",
    y = "Clinical Concept"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 11, face = "bold"),
    axis.text.y = element_text(size = 9),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.position = "right",
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.switch.pad.grid = unit(0, "pt"),
    panel.spacing.y = unit(0.06, "lines")
  )

p_female1

ggsave(
  here("results", dbName, "lsc", "processed", "female_and_breast_cancers_heatmap_with_concepts_sens.png"),
  plot = p_female1,
  width = 9,
  height = 8,
  units = "in",
  dpi = 600,
  device = "png",
  type = "cairo",
  bg = "white"
)


# genitourinary ------------------------
plot_df_genitourinary1 <- filtered_condition_results_summary3_sens %>%
  filter(cancer %in% c( "bladder", 
                        "kidney", 
                        "prostate" )) %>%
  mutate(
    cancer = factor(
      cancer,
      levels = c( "kidney", 
                  "bladder",
                  "prostate"),
      labels = c("Kidney", 
                 "Bladder",
                 "Prostate")
    ),
    assigned_domain = factor(assigned_domain, levels = domain_order)
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
    concept_name = factor(
      concept_name,
      levels = rev(unique(concept_name))
    )
  )

p_genitourinary1 <- ggplot(
  plot_df_genitourinary1,
  aes(
    x = cancer,
    y = concept_name,
    fill = sampled_SMD
  )
) +
  geom_tile(
    colour = "grey70",
    linewidth = 0.25
  ) +
  scale_fill_gradientn(
    name = "ASMD",
    colours = c("#1a9850", "#fdae61", "#b2182b"),
    limits = c(0.1, 0.4),
    breaks = c(0.1, 0.2, 0.3, 0.4),
    oob = scales::squish,
    na.value = "grey95"
  ) +
  scale_x_discrete(position = "top") +
  labs(
    x = "Cancer",
    y = "Clinical Concept"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 11, face = "bold"),
    axis.text.y = element_text(size = 9),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.position = "right",
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.switch.pad.grid = unit(0, "pt"),
    panel.spacing.y = unit(0.06, "lines")
  )

ggsave(
  here("results", dbName, "lsc", "processed", "genitourinary_cancers_heatmap_with_concepts_sens.png"),
  plot = p_genitourinary1,
  width = 6,
  height = 6,
  units = "in",
  dpi = 600,
  device = "png",
  type = "cairo",
  bg = "white"
)

# other ----------------------------------
plot_df_other1 <- filtered_condition_results_summary3_sens %>%
  filter(cancer %in% c( "brain", 
                        "thyroid" )) %>%
  mutate(
    cancer = factor(
      cancer,
      levels = c( "brain", 
                  "thyroid"),
      labels = c("Brain", 
                 "Thyroid")
    ),
    assigned_domain = factor(assigned_domain, levels = domain_order)
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
    concept_name = factor(
      concept_name,
      levels = rev(unique(concept_name))
    )
  )

p_other1 <- ggplot(
  plot_df_other1,
  aes(
    x = cancer,
    y = concept_name,
    fill = sampled_SMD
  )
) +
  geom_tile(
    colour = "grey70",
    linewidth = 0.25
  ) +
  scale_fill_gradientn(
    name = "ASMD",
    colours = c("#1a9850", "#fdae61", "#b2182b"),
    limits = c(0.1, 0.4),
    breaks = c(0.1, 0.2, 0.3, 0.4),
    oob = scales::squish,
    na.value = "grey95"
  ) +
  scale_x_discrete(position = "top") +
  labs(
    x = "Cancer",
    y = "Clinical Concept"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 11, face = "bold"),
    axis.text.y = element_text(size = 9),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.position = "right",
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.switch.pad.grid = unit(0, "pt"),
    panel.spacing.y = unit(0.06, "lines")
  )

ggsave(
  here("results", dbName, "lsc", "processed", "other_cancers_heatmap_with_concepts_sens.png"),
  plot = p_other1,
  width = 6,
  height = 6,
  units = "in",
  dpi = 600,
  device = "png",
  type = "cairo",
  bg = "white"
)
