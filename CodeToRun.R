# Characterisation of signs and symptoms of cancer

# install renv and restore the environment
# install.packages("renv")
# if its your first time using the project need to activate it (uncomment and run line below)
# renv::activate()
renv::restore()

# call the packages you need for the study
library(CDMConnector)
library(DBI)
library(log4r)
library(dplyr)
library(here)
library(OmopSketch)
library(omopgenerics)
library(CohortConstructor)
library(PhenotypeR)
library(CohortCharacteristics)
library(PatientProfiles)
library(visOmopResults)
library(stringr)
library(CodelistGenerator)
library(DrugUtilisation)
library(odbc)
library(RPostgres)
library(readr)
library(tibble)

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
cdmSchema <- "..."

# A prefix for all permanent tables in the database
writePrefix <- "..."

# The name of the schema where results tables will be created
writeSchema <- "..."

# minimum counts that can be displayed according to data governance
minCellCount <- 5

# start date for study period
study_start_date <- as.Date("2005-01-01")

# end date for the study period
study_end_date <- as.Date("2024-12-31")

# Create cdm object ----
cdm <- cdmFromCon(
  con = db,
  cdmSchema = cdmSchema,
  writeSchema = writeSchema,
  writePrefix = writePrefix,
  achillesSchema = writeSchema,
  cdmName = dbName
)

# run study
source(here("RunStudy.R"))
