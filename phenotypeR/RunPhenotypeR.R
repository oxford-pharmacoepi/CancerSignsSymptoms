# Get cancer cohorts

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

# create basic cohort
cdm[["pss_phenotyper"]] <- conceptCohort(cdm,
                                     conceptSet = cancer_codelists,
                                     exit = "event_start_date",
                                     name = "pss_phenotyper")



# run phenotypeR
result <- phenotypeDiagnostics(cdm$pss_phenotyper)

# export the phenotypeR results
exportSummarisedResult(result, path = here::here("phenotypeR"), minCellCount = minCellCount, fileName = paste0(dbName, "_result.csv"))



