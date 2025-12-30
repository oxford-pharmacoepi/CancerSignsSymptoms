set.seed(1234)

# read in and create a codelist and only use those in cdm
cancer_codelists <- read_csv(paste0(here("codelists"), "/prelim_codelists_reviewed.csv")) %>%
  filter(keep == "y") %>%
  group_by(cancer) %>%
  summarise(concept_id = list(concept_id)) %>%
  deframe() %>%
  newCodelist() 

#han combined
han_codelist <- read_csv(paste0(here("codelists"), "/prelim_codelists_reviewed.csv")) %>%
  filter(cancer %in% c("hypopharynx",
                       "larynx",
                       "lip",
                       "nasopharynx",
                       "oral cavity",
                       "oropharynx" ,
                       "salivary gland",
                       "head and neck") ) %>% 
  filter(keep == "y")

#create han codelist using codelist by combining all concepts together
# # create a combination for han cancer
han <- list(han = c(han_codelist$concept_id)) %>%
  omopgenerics::newCodelist() 

# create a codelist for lip and oral cavity combined
lip_oral_cavity_codelist <- read_csv(paste0(here("codelists"), "/prelim_codelists_reviewed.csv")) %>%
  filter(cancer %in% c("lip",
                       "oral cavity" ) ) %>% 
  filter(keep == "y")

lip_oral_cavity <- list(lip_oral_cavity = c(lip_oral_cavity_codelist$concept_id)) %>%
  omopgenerics::newCodelist() 

# create a combined lung, trachea and bronchus codelist
lung_trachea_bronchus_codelist <- read_csv(paste0(here("codelists"), "/prelim_codelists_reviewed.csv")) %>%
  filter(cancer %in% c("lung",
                       "bronchus" ,
                       "trachea"  ,
                       "respiratory tract" ) ) %>% 
  filter(keep == "y")

lung_bron_trac <- list(lung_bron_trac = c(lung_trachea_bronchus_codelist$concept_id)) %>%
  omopgenerics::newCodelist()


# create a combined colorectal cancer codelist
colorectal_codelist <- read_csv(paste0(here("codelists"), "/prelim_codelists_reviewed.csv")) %>%
  filter(cancer %in% c("colon",
                       "rectum",
                       "colorectal") ) %>% 
  filter(keep == "y")

# create colorectal combined codelist
colorectal <- list(colorectal = c(colorectal_codelist$concept_id)) %>%
  omopgenerics::newCodelist()

# create a combined brain cancer codelist
brain_codelist <- read_csv(paste0(here("codelists"), "/prelim_codelists_reviewed.csv")) %>%
  filter(cancer %in% c("brain",
                       "glioma") ) %>% 
  filter(keep == "y")

# create brain combined codelist
brain <- list(brain = c(brain_codelist$concept_id)) %>%
  omopgenerics::newCodelist()

# remove codelists from cancer_codelists that are no longer required (individual components so lip, oral cavity, lung, bronc etc)
cancer_codelists <- cancer_codelists[!names(cancer_codelists) %in% c("lung",
                                                                     "bronchus",
                                                                     "trachea",
                                                                     "respiratory tract",
                                                                     "colorectal",
                                                                     "head and neck",
                                                                     "lip",
                                                                     "brain",
                                                                     "glioma",
                                                                     "oral cavity")] %>% 
  newCodelist()

# # combine all together into one codelist
cancer_codelists <- bind(
  cancer_codelists,
  han,
  lip_oral_cavity,
  lung_bron_trac,
  colorectal,
  brain
) %>%
  newCodelist()

# any cancer (to make sure cases are primary cases)
any_cancer_codelists <- CodelistGenerator::codesFromConceptSet(
  path = here::here("codelists", "exclusions"),
  cdm = cdm
)

# create cohort of any cancer (excl skin cancer)
cdm$any_cancer <- CohortConstructor::conceptCohort(
  cdm = cdm,
  conceptSet = any_cancer_codelists,
  name = "any_cancer",
)

# create the cancer cohorts ------
cdm[["pss_365_po"]] <- conceptCohort(cdm,
                                     conceptSet = cancer_codelists,
                                     exit = "event_start_date",
                                     name = "pss_365_po") %>% 
  requireIsFirstEntry() %>%
  requirePriorObservation(minPriorObservation = 365) %>% 
  requireInDateRange(indexDate = "cohort_start_date",
                     dateRange = as.Date(c(study_start_date, study_end_date))) %>% 
  requireAge(ageRange = list(c(18, 150)) ) 

# exclude those with other cancers the year before cancer of interest
cdm[["pss_365_po_clean"]] <- cdm[["pss_365_po"]] %>%
  requireCohortIntersect(
    intersections = 0,
    name = "pss_365_po_clean",
    targetCohortTable = "any_cancer",
    window = c(-365, -1)       
  )

# run phenotypeR
result <- phenotypeDiagnostics(
  cdm[["pss_365_po_clean"]],
  diagnostics = c("databaseDiagnostics",
                  "codelistDiagnostics",
                  "cohortDiagnostics"),
  survival = FALSE,
  matchedSample = NULL,
  populationSample = NULL,
  populationDateRange = as.Date(c(NA, NA))
)

# # add in baseline characteristics
# codelistConditions <- CodelistGenerator::codesFromConceptSet(here("codelists", "Conditions"), cdm)
#
# cdm <- CDMConnector::generateConceptCohortSet(cdm = cdm,
#                                               conceptSet = codelistConditions,
#                                               name = "conditions",
#                                               overwrite = TRUE)
#
#
#
# # medications -----
# codelistMedications <- CodelistGenerator::codesFromConceptSet(here("codelists", "Medications"), cdm)
#
# cdm <- DrugUtilisation::generateDrugUtilisationCohortSet(cdm = cdm,
#                                                          conceptSet = codelistMedications,
#                                                          name = "medications")
#
# # get matched cohort for table one
# cdm$matched_pss <- matchCohorts(
#   cohort = cdm$pss_365_po_clean,
#   name = "matched_pss")
#
# # get nsaid cohort and add in inclusion/exclusion ---
# cdm$pss_characteristics <- cdm$matched_pss %>%
#   PatientProfiles::addDemographics()
#
# cdm$pss_characteristics <- cdm$pss_characteristics %>%
#   PatientProfiles::addDemographics(
#     ageGroup = list(
#       "age_group" =
#         list(
#           "18 to 49" = c(18, 49),
#           "50 to 59" = c(50, 59),
#           "60 to 69" = c(60, 69),
#           "70 to 79" = c(70, 79),
#           "80 to 89" = c(80, 89),
#           "90+" = c(90, 150)
#         )
#     ))
#
#
#
# # characterise comorbs and medications
# suppressWarnings(
#
#   summarycharacteristics <- cdm$pss_characteristics %>%
#     CohortCharacteristics::summariseCharacteristics(
#       ageGroup = list( "18 to 49" = c(18, 49),
#                        "50 to 59" = c(50, 59),
#                        "60 to 69" = c(60, 69),
#                        "70 to 79" = c(70, 79),
#                        "80 to 89" = c(80, 89),
#                        "90+" = c(90, 150)),
#       cohortIntersectFlag = list(
#         "Conditions any time prior to index date" = list(
#           targetCohortTable = "conditions",
#           window = c(-Inf, -1)
#         ),
#
#         "Conditions 365 days prior to index date" = list(
#           targetCohortTable = "conditions",
#           window = c(-365, -1)
#         ),
#
#
#         "Medications 365 days prior to index date" = list(
#           targetCohortTable = "medications",
#           window = c(-365, -1)
#         ),
#
#         "Medications 90 days prior to index date" = list(
#           targetCohortTable = "medications",
#           window = c(-90, -1)
#         )
#
#
#       )
#
#     )
# )
#
# omopgenerics::exportSummarisedResult(summarycharacteristics,
#                                      minCellCount = 5,
#                                      path = here("results"),
#                                      fileName = paste0(cdmName(cdm),
#                                                        "_summary_characteristics.csv")
# )

# export the phenotypeR results
exportSummarisedResult(result, path = here::here("results"), minCellCount = minCellCount, fileName = paste0(dbName, "_result.csv"))
