# create results folders
dir.create(here::here("results", dbName, "characteristics", "processed"),
           recursive = TRUE,
           showWarnings = FALSE)

dir.create(here::here("results", dbName, "lsc", "processed"),
           recursive = TRUE,
           showWarnings = FALSE)


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

# anyone with a prior cancer (apart from skin cancer) is excluded
# cdm[["pss_365_po_clean_inf"]] <- cdm[["pss_365_po"]] %>%
#   requireCohortIntersect(
#     intersections = 0,
#     name = "pss_365_po_clean_inf",
#     targetCohortTable = "any_cancer",
#     window = c(-Inf, -1)       
#   )

# baseline characteristics -------

# conditions -------
codelistConditions <- CodelistGenerator::codesFromConceptSet(here("codelists", "Conditions"), cdm)

cdm <- CDMConnector::generateConceptCohortSet(cdm = cdm,
                                              conceptSet = codelistConditions,
                                              name = "conditions",
                                              overwrite = TRUE)



# medications -----
codelistMedications <- CodelistGenerator::codesFromConceptSet(here("codelists", "Medications"), cdm)

cdm <- DrugUtilisation::generateDrugUtilisationCohortSet(cdm = cdm,
                                                         conceptSet = codelistMedications,
                                                         name = "medications")




# updated match cohort function adapted from CohortConstructor
source(here::here("match_cohort_po.R"))

#set up for lsc
lscWindows <- list(c(-Inf, -366), c(-365, -31), c(-365, -91))
lscTableEvents <-c("condition_occurrence", 
                   "procedure_occurrence" ,
                   "observation")
lscTableEpisodes <- c("drug_exposure")
lscMminimumFrequency <- 0


# capture results for processing
results_lsc <- list()

# cdm <- cdmFromCon(
#   con = db,
#   cdmSchema = cdmSchema,
#   writeSchema = writeSchema,
#   writePrefix = writePrefix,
#   achillesSchema = writeSchema,
#   cohortTables = c("pss_365_po" ,
#                    "conditions",
#                    "medications",
#                    #"pss_365_po_clean_inf",
#                    "any_cancer",
#                    "pss_365_po_clean"),
#   cdmName = dbName
# )


# create a copy of the cohort and filter it per cancer then do the matching and create a loop of analysis
# you could do it all in one step however if you have a slow server it is better to do per cancer
for (i in 1:length(names(cancer_codelists))) {
  
#create cohort name (remove and spaces with underscores)
cancer_name <- paste0(
  "pss_365_",
  gsub(" ", "_", names(cancer_codelists)[[i]])
)

message("Starting: ", cancer_name)

# create a copy of the cancer cohort
cdm[[cancer_name]] <- copyCohorts(cdm[["pss_365_po_clean"]], cohortId = i, name = cancer_name)

# exclude males for females based cancer and females from male based cancers
cdm[[cancer_name]] <- cdm[[cancer_name]] %>%
  addSex() %>%
  filter(
    !(sex == "Male" & cancer_name %in% c("pss_365_ovary", "pss_365_cervix", "pss_365_corpus_uteri")),
    !(sex == "Female" & cancer_name == "pss_365_prostate")
  )

# create cohort name for matched cohorts
cancer_name_match <- paste0(
  "pss_365_matched_",
  gsub(" ", "_", names(cancer_codelists)[[i]])
)

set.seed(digest::digest2int(cancer_name_match))

# carry out matching
cdm[[cancer_name_match]] <- matchCohortsPo(
  cohort = cdm[[cancer_name]],
  matchSex = TRUE,
  matchYearOfBirth = TRUE,
  ratio = 1,
  keepOriginalCohorts = FALSE,
  name = cancer_name_match
)

# baseline characteristics
cdm[[cancer_name_match]] <- cdm[[cancer_name_match]] %>%
  PatientProfiles::addDemographics(
    ageGroup = list(
      "age_group" =
        list(
          "18 to 49" = c(18, 49),
          "50 to 59" = c(50, 59),
          "60 to 69" = c(60, 69),
          "70 to 79" = c(70, 79),
          "80 to 89" = c(80, 89),
          "90+" = c(90, 150)
        )
    ))



# characterise comorbs and medications
suppressWarnings(

  summarycharacteristics <- cdm[[cancer_name_match]] %>%
    CohortCharacteristics::summariseCharacteristics(
      ageGroup = list( "18 to 49" = c(18, 49),
                       "50 to 59" = c(50, 59),
                       "60 to 69" = c(60, 69),
                       "70 to 79" = c(70, 79),
                       "80 to 89" = c(80, 89),
                       "90+" = c(90, 150)),
      cohortIntersectFlag = list(
        "Conditions any time prior to index date" = list(
          targetCohortTable = "conditions",
          window = c(-Inf, -1)
        ),

        "Conditions 365-31 days prior to index date" = list(
          targetCohortTable = "conditions",
          window = c(-365, -31)
        ),
        
        "Conditions 365-91 days prior to index date" = list(
          targetCohortTable = "conditions",
          window = c(-365, -91)
        ),
        
        "Any Cancer 365-31 days prior to index date" = list(
          targetCohortTable = "any_cancer",
          window = c(-365, -31)
        ),


        "Medications 365-31 days prior to index date" = list(
          targetCohortTable = "medications",
          window = c(-365, -31)
        ),
        
        "Medications 365-91 days prior to index date" = list(
          targetCohortTable = "medications",
          window = c(-365, -91)
        ),

        "Medications 90 days prior to index date" = list(
          targetCohortTable = "medications",
          window = c(-90, -1)
        )


      ) ,

      tableIntersectCount = list(
        "Visits in 365-31 days prior index date" = list(
          tableName = "visit_occurrence",
          window = c(-365, -31)
        )
      ) 
      
    )
)


exportSummarisedResult(summarycharacteristics,
                       path = here::here("results", dbName, "characteristics"),
                       minCellCount = minCellCount,
                       fileName = paste0(dbName, "_", cancer_name, "_baseline_characteristics.csv"))


# lsc
lsc1 <- CohortCharacteristics::summariseLargeScaleCharacteristics(
  cohort = cdm[[cancer_name_match]],
  window = lscWindows,
  eventInWindow = lscTableEvents,
  episodeInWindow = lscTableEpisodes,
  minimumFrequency = lscMminimumFrequency,
  includeSource = c(FALSE),
  excludedCodes = NULL
)

exportSummarisedResult(lsc1, path = here::here("results", dbName, "lsc"), 
                       minCellCount = 0, 
                       fileName = paste0(dbName, "_", cancer_name, "_lsc.csv"))

# save result for processing
results_lsc[[cancer_name]] <- lsc1

rm(lsc1)

message("Finished: ", cancer_name)


}


# script for processing lsc to calculate SMD and creating figures for manuscripts
source(here::here("result_processing_lsc_conditions.R"))
source(here::here("result_processing_lsc_drugs.R"))
source(here::here("result_processing_lsc_observation.R"))

# script for processing characteristics for paper
source(here::here("result_processing_characteristics.R"))




# run phenotypeR
# result <- phenotypeDiagnostics(
#   cdm[["pss_365_po_clean_inf"]],
#   diagnostics = c("databaseDiagnostics",
#                   "codelistDiagnostics",
#                   "cohortDiagnostics"),
#   survival = FALSE,
#   matchedSample = NULL,
#   populationSample = NULL,
#   populationDateRange = as.Date(c(NA, NA))
# )

# export the phenotypeR results
#exportSummarisedResult(result, path = here::here("results"), minCellCount = minCellCount, fileName = paste0(dbName, "_result.csv"))
