# call the packages you need for the study
library(CDMConnector)
library(DBI)
library(dplyr)
library(here)
library(stringr)
library(CodelistGenerator)
library(odbc)
library(RPostgres)
library(readr)

# database metadata and connection details
# The name/ acronym for the database
dbName <- "..."

server_dbi <- "..."
user       <- "..."
password   <- "..."
port       <- "..."
host       <- "..."


# database connection details
db <- dbConnect(RPostgres::Postgres(),
                dbname = server_dbi,
                host   = host,
                user   = user,
                password = password)


# The name of the schema that contains the OMOP CDM with patient-level data
cdmSchema <- "public"

# A prefix for all permanent tables in the database
writePrefix <- "css_codelists"

# The name of the schema where results tables will be created
writeSchema <- "results"


# Create cdm object ----
cdm <- cdmFromCon(
  con = db,
  cdmSchema = cdmSchema,
  writeSchema = writeSchema,
  writePrefix = writePrefix,
  achillesSchema = writeSchema,
  cdmName = dbName
)


# create a vector of cancers
x <- c( "bladder",
        "kidney",
        "liver",
        "colon",
        "rectum",
        "anus",
        "brain",
        "breast",
        "cervix",
        "corpus uteri",
        "hypopharynx",
        "larynx",
        "lip",
        "oral cavity",
        "lung",
        "nasopharynx",
        "oropharynx",
        "ovary",
        "pancreas",
        "prostate",
        "salivary gland",
        "stomach",
        "thyroid",
        "esophagus")

codes <- list()

for (i in 1:length(x)) {

keywords <- c(
  paste("malignant neoplasm of", x[i]),
  paste("malignant tumor of", x[i])
)

codes[[i]] <- getCandidateCodes(
    cdm = cdm,
    keywords = keywords,
    exclude = c("benign",
                "skin",
                "recurrence",
                "lymphoma",
                "secondary") ,
    domains = "Condition",
    standardConcept = "Standard",
    searchInSynonyms = FALSE,
    searchNonStandard = FALSE,
    includeDescendants = TRUE,
    includeAncestor = FALSE
  ) %>%
  mutate(cancer = x[i])

}
# bind codes together
codes_together <- bind_rows(codes)

#export codelists
write_csv(codes_together, paste0(here("codelists"), "/prelim_codelists.csv"))
