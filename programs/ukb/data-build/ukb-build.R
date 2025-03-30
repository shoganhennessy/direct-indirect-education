#!/usr/bin/R
## Senan Hogan-Hennessy, 11 March 2025.
## Script to ingest raw UKB data, output familiar tabular format.
print(Sys.time())
set.seed(47)
## Packages:
# functions for data manipulation and visualisation
library(tidyverse)
# Functions for fast data manipulation (adapting to data.table code from Kweon+)
library(data.table)

# Define folder paths (1) where data are (2) input data (3) intermed files.
data.folder <- file.path("..", "..", "..", "data", "ukb-restricted")
input.folder <- file.path(data.folder, "input")
output.folder <- file.path(data.folder, "intermediate-files")


################################################################################
## Load raw data

# Load the phenotype data (raw format from UKB RAP servers).
ukb_pheno.data <- input.folder %>%
    file.path("phenotype-extract.csv") %>%
    read_csv()

# Load the Ed PGI (Okbay+ 2022) data, raw format from UKB RAP servers.
ukb_edpgi.data <- input.folder %>%
    file.path("ed-pgi-score.sscore") %>%
    read_tsv()

#TODO:
# Load the imputed Ed PGI (Young+ 2022), raw format from UKB RAP servers.
ukb_imputed.data <- input.folder %>%
    file.path("ed-pgi-score.sscore") %>%
    read_tsv()

# Load the occupation-coded income data (thanks to Kweon Koellinger+ 2025).
load(file.path(input.folder, "earnings-imputed", "data_input.Rdata"))



# CPI-U from (WHERE?)
#! Placeholder.


################################################################################
## Clean Pheno + PGI files.

# Clean the phenotype data.
cleaned_pheno.data <- ukb_pheno.data %>%
    # Restrict data types
    transmute(
        eid                    = as.integer(participant.eid),
        sex_male               = as.integer(participant.p31),
        birthyear              = as.integer(participant.p34),
        birthmonth             = as.integer(participant.p52),
        datevisit              = date(participant.p53_i0),
        # Income categories,
        householdincome_cat    = as.integer(participant.p738_i0),
        recruitedage           = as.integer(participant.p21022),
        datelastcontact        = date(participant.p20143),
        edquals                = as.character(participant.p6138_i0),
        agefinishededuc        = as.integer(participant.p845_i0),
        ethnicity              = as.integer(participant.p21000_i0),
        genetic_race           = as.integer(participant.p22006),
        jobcode                = as.integer(participant.p20277_i0),
        numjobs                = as.integer(participant.p22599),
        numjobgaps             = as.integer(participant.p22661),
        jobcode_online         = as.integer(participant.p22601_a0),
        jobcode_onlineSOC2000  = as.integer(participant.p22617_a0),
        jobyearstarted_online  = as.integer(participant.p22602_a0),
        jobhoursworked_online  = as.integer(participant.p22604_a0),
        jobhoursperweek_online = as.numeric(participant.p22605_a0),
        yearendededuc_online   = as.integer(participant.p22501),
        datedeath              = date(participant.p40000_i0),
        PCA                    = as.integer(participant.p22020)) %>%
    # Clean the resulting columns
    mutate(
        # householdincome_cat,    p738_i0,    https://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=738
        householdincome_string = 
            ifelse(householdincome_cat == 1, "< 18k",
            ifelse(householdincome_cat == 2, "18 - 40k",
            ifelse(householdincome_cat == 3, "31 - 60k",
            ifelse(householdincome_cat == 4, "52 - 100k",
            ifelse(householdincome_cat == 5, "100k < ", "Unknown"))))),
        # Missing values to missing.
        agefinishededuc    = ifelse(agefinishededuc < 0, NA, agefinishededuc),
        # Replace missing values with bianry values.
        genetic_euroancest = ifelse(
            is.na(genetic_race), 0, ifelse(genetic_race == 1, 1, 0)),
        inPCA = ifelse(is.na(PCA), 0, ifelse(PCA == 1, 1, 0)))

# Restrict Ed PGI to identifiers, and the score
cleaned_edpgi.data <- ukb_edpgi.data %>%
    transmute(
        eid = IID,
        edpgi_raw = SCORE1_AVG)

# Restrict imputed variables identifiers, and the score
cleaned_imputed.data <- ukb_imputed.data %>%
    transmute(
        #TODO: adjust to output of PRS impute
        eid = IID,
        edpgi_raw = SCORE1_AVG)


################################################################################
## Calculate OCC coded wages, thanks to Kweon Koellinger provided code.

#TODO: code here that adapts the data.table code they provided to
#TODO: the format I have here.
ukb_income.data <- ukb_pheno.data %>% data.table()

################################################################################
## Merge UKB data files.

# Merge the Phenotype data with Ed PGI.
cleaned_pheno.data <- cleaned_pheno.data %>%
    left_join(cleaned_edpgi.data, by = "eid")

# Merge the Phenotype data with parental imputed Ed PGI.
cleaned_pheno.data <- cleaned_pheno.data %>%
    left_join(cleaned_edpgi.data, by = "eid")



################################################################################
## Clean the Ed PGI data fields.



# Standardise the Ed PGI.
cleaned_pheno.data <- cleaned_pheno.data %>%
    # PGI reliable for EUro ancestry, and reliable PCA data.
    mutate(edpgi_norm = edpgi_raw * inPCA * genetic_euroancest) %>%
    # Normalise to have 0 mean, 1 SD.
    mutate(edpgi_norm = (edpgi_norm - mean(edpgi_norm, na.rm = TRUE)
        ) / sd(edpgi_norm, na.rm = TRUE))

# Ensure outliers in PGI do not persist.
hist(cleaned_pheno.data$edpgi_raw)
hist(cleaned_pheno.data$edpgi_norm)


################################################################################
## Save workable data file.

cleaned_pheno.data %>%
    write_csv(file.path(output.folder, "ukb-intermed-pheno.csv"))
